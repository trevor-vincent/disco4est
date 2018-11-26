#include <d4est_ghost.h>
#include <d4est_ghost_data_ext.h>
#include <d4est_element_data.h>
#include <string.h>

//#define D4EST_GHOST_DATA_EXT_DEBUG

static void
d4est_ghost_data_ext_compute_ids_and_strides
(
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_ext_t* d4est_ghost_data_ext
)
{
  int stride_global_receive = 0;
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ghost_elem_data = &d4est_ghost->ghost_elements[gid];
    ghost_elem_data->id = gid;

    int stride_local_receive = 0;
    for (int i = 0; i < d4est_ghost_data_ext->num_vecs; i++){
      /* d4est_field_type_t type = d4est_ghost_data_ext->transfer_types[i]; */
      int field_elem_size = d4est_ghost_data_ext->field_size_of_ghost_fcn
                            (
                             d4est_ghost,
                             gid,
                             i,
                             d4est_ghost_data_ext->user_ctx
                            );
      
      d4est_ghost_data_ext->receive_strides[gid][i] = stride_global_receive + stride_local_receive;     
      stride_local_receive += field_elem_size;
    }    
    stride_global_receive += stride_local_receive;
  }
}


d4est_ghost_data_ext_t*
d4est_ghost_data_ext_init
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 int num_vecs,
 d4est_ghost_data_ext_fcn_t field_size_of_ghost_fcn,
 d4est_ghost_data_ext_fcn_t field_size_of_mirror_fcn,
 d4est_ghost_data_ext_fcn_t field_stride_of_mirror_fcn,
 void* user_ctx
)
{
#ifdef D4EST_GHOST_DATA_EXT_DEBUG
  printf("Starting d4est_ghost_data_ext_init\n");
#endif
 
  d4est_ghost_data_ext_t* d4est_ghost_data_ext = D4EST_ALLOC(d4est_ghost_data_ext_t, 1); 

  d4est_ghost_data_ext->field_size_of_ghost_fcn = field_size_of_ghost_fcn;
  d4est_ghost_data_ext->field_size_of_mirror_fcn = field_size_of_mirror_fcn;
  d4est_ghost_data_ext->field_stride_of_mirror_fcn = field_stride_of_mirror_fcn;
  d4est_ghost_data_ext->user_ctx = user_ctx;
  
  d4est_ghost_data_ext->num_ghosts = d4est_ghost->ghost->ghosts.elem_count;
  d4est_ghost_data_ext->num_vecs = num_vecs;

  d4est_ghost_data_ext->receive_strides = D4EST_ALLOC(int*, d4est_ghost_data_ext->num_ghosts);
  for (int i = 0; i < d4est_ghost_data_ext->num_ghosts; i++){
    d4est_ghost_data_ext->receive_strides[i] = D4EST_ALLOC(int, d4est_ghost_data_ext->num_vecs);
  }
  

  d4est_ghost_data_ext->ghost_data_sizes = D4EST_ALLOC_ZERO(int, p4est->mpisize);
  d4est_ghost_data_ext->receive_size = 0;

  for (int tn = 0; tn < num_vecs; tn++){
    for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
      d4est_element_data_t* ghost_elem_data = &d4est_ghost->ghost_elements[gid];
      int ghost_rank = d4est_ghost->ghost_elements[gid].mpirank;
      int size = d4est_ghost_data_ext->field_size_of_ghost_fcn
                            (
                             d4est_ghost,
                             gid,
                             tn,
                             d4est_ghost_data_ext->user_ctx
                            );
      d4est_ghost_data_ext->ghost_data_sizes[ghost_rank] += size;
      d4est_ghost_data_ext->receive_size += size;
    }
  }

  d4est_ghost_data_ext->receive_data = D4EST_ALLOC(char, d4est_ghost_data_ext->receive_size);
  d4est_ghost_data_ext_compute_ids_and_strides
    (
     d4est_ghost,
     d4est_ghost_data_ext
    );
  
  return d4est_ghost_data_ext;
}

void
d4est_ghost_data_ext_destroy
(
 d4est_ghost_data_ext_t* d4est_ghost_data_ext
)
{
  if (d4est_ghost_data_ext != NULL){
    /* D4EST_FREE(d4est_ghost_data_ext->transfer_types); */
    /* D4EST_FREE(d4est_ghost_data_ext->transfer_strides); */
    D4EST_FREE(d4est_ghost_data_ext->ghost_data_sizes);
  
    for (int i = 0; i < d4est_ghost_data_ext->num_ghosts; i++){
      D4EST_FREE(d4est_ghost_data_ext->receive_strides[i]);
    }
    D4EST_FREE(d4est_ghost_data_ext->receive_strides);
    D4EST_FREE(d4est_ghost_data_ext->receive_data);

    d4est_ghost_data_ext->receive_size = -1;
    d4est_ghost_data_ext->num_ghosts = -1;
    d4est_ghost_data_ext->num_vecs = -1;
  
    D4EST_FREE(d4est_ghost_data_ext);
    d4est_ghost_data_ext = NULL;
  }
}





void
d4est_ghost_data_ext_exchange
(
 p4est_t * p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_ext_t* d4est_ghost_data_ext,
 char** transfer_vecs
)
{
#ifdef D4EST_GHOST_DATA_EXT_DEBUG
  printf("Starting d4est_ghost_data_ext_exchange\n");
#endif
  
  const int           num_procs = p4est->mpisize;
  int                 mpiret;
  int                 q;
  char               *mem, **sbuf;
  size_t              zz;
  sc_array_t          requests, sbuffers;
  p4est_locidx_t      ng_excl, ng_incl, ng, theg;
  p4est_locidx_t      mirr;
  sc_MPI_Request     *r;

  sc_array_init (&requests, sizeof (sc_MPI_Request));
  sc_array_init (&sbuffers, sizeof (char *));

  /* receive data from other processors */
  int data_size_offset = 0;
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = d4est_ghost->ghost->proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    D4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      r = (sc_MPI_Request *) sc_array_push (&requests);

      /* this is wrong, check if ghost_data_sizes is sizeof(double) in it, otherwise its wrong */
      mpiret = sc_MPI_Irecv ((char *) (d4est_ghost_data_ext->receive_data) + data_size_offset,
                             d4est_ghost_data_ext->ghost_data_sizes[q], sc_MPI_BYTE, q,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      ng_excl = ng_incl;
      data_size_offset += d4est_ghost_data_ext->ghost_data_sizes[q];
    }
  }
  D4EST_ASSERT (ng_excl == (p4est_locidx_t) d4est_ghost->ghost->ghosts.elem_count);

  /* send data to other processors */
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = d4est_ghost->ghost->mirror_proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    D4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      /* calculate total data size to send */
      int total_memsize = 0;
      for (theg = 0; theg < ng; ++theg){
        mirr = d4est_ghost->ghost->mirror_proc_mirrors[ng_excl + theg];
        D4EST_ASSERT (0 <= mirr && (size_t) mirr < d4est_ghost->ghost->mirrors.elem_count);
        for (int tn = 0; tn < d4est_ghost_data_ext->num_vecs; tn++){
          total_memsize +=  d4est_ghost_data_ext->field_size_of_mirror_fcn
                            (
                             d4est_ghost,
                             mirr,
                             tn,
                             d4est_ghost_data_ext->user_ctx
                            );
          /* printf("total_memsize = %d on proc %d\n", total_memsize, p4est->mpirank); */
        }
      }
      /* every peer populates its own send buffer */
      sbuf = (char **) sc_array_push (&sbuffers);
      mem = *sbuf = D4EST_ALLOC (char, total_memsize);

      for (theg = 0; theg < ng; ++theg) {
        mirr = d4est_ghost->ghost->mirror_proc_mirrors[ng_excl + theg];
        D4EST_ASSERT (0 <= mirr && (size_t) mirr < d4est_ghost->ghost->mirrors.elem_count);
        for (int tn = 0; tn < d4est_ghost_data_ext->num_vecs; tn++){
          int mem_size = d4est_ghost_data_ext->field_size_of_mirror_fcn
                            (
                             d4est_ghost,
                             mirr,
                             tn,
                             d4est_ghost_data_ext->user_ctx
                            );

          char* field = transfer_vecs[tn];
          D4EST_ASSERT(field != NULL);
          int stride = d4est_ghost_data_ext->field_stride_of_mirror_fcn
                            (
                             d4est_ghost,
                             mirr,
                             tn,
                             d4est_ghost_data_ext->user_ctx
                            );
          
          memcpy (mem, &field[stride], mem_size);
          mem += mem_size;
        }
      }
      r = (sc_MPI_Request *) sc_array_push (&requests);
      mpiret = sc_MPI_Isend (*sbuf, total_memsize, sc_MPI_BYTE, q,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      ng_excl = ng_incl;
    }
  }

  /* wait and clean up */
  mpiret = sc_MPI_Waitall (requests.elem_count, (sc_MPI_Request *)
                           requests.array, sc_MPI_STATUSES_IGNORE);
  /* SC_CHECK_MPI (mpiret); */
  sc_array_reset (&requests);
  for (zz = 0; zz < sbuffers.elem_count; ++zz) {
    sbuf = (char **) sc_array_index (&sbuffers, zz);
    D4EST_FREE (*sbuf);
  }
  sc_array_reset (&sbuffers);
  
}


void* d4est_ghost_data_ext_get_field_on_element
(
 d4est_element_data_t* ed,
 int ghost_name_id,
 d4est_ghost_data_ext_t* dgd
){
  D4EST_ASSERT(dgd != NULL);
  D4EST_ASSERT(ghost_name_id >= 0 && ghost_name_id < dgd->num_vecs);
  int stride = dgd->receive_strides[ed->id][ghost_name_id];
  return (void*)&dgd->receive_data[stride];
}
