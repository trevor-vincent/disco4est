#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_element_data.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_transfer_ghost_data.h>
#include <time.h>

static int field_size_of_ghost_fcn
(
 d4est_ghost_t* d4est_ghost,
 int gid,
 int tn,
 void* user_ctx
)
{
  d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx;
  d4est_solver_schwarz_subdomain_metadata_t* sub_data =
    (d4est_solver_schwarz_subdomain_metadata_t*)
    d4est_ghost_data_ext_get_field_on_element
    (
     &d4est_ghost->ghost_elements[gid],
     0,
     schwarz_metadata->subdomain_ghostdata
    );

  int point_size = sub_data->nodal_size;
  return sizeof(double)*point_size;
}

static int field_size_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirror,
 int tn,
 void* user_ctx
)
{
  d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx;
  d4est_element_data_t* ed = &d4est_ghost->mirror_elements[mirror];
  int point_size = schwarz_metadata->subdomain_metadata[ed->id].nodal_size;
  return sizeof(double)*point_size;
}

static int field_stride_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx;
  d4est_element_data_t* ed = &d4est_ghost->mirror_elements[mirr];
  int point_stride = schwarz_metadata->subdomain_metadata[ed->id].nodal_stride;
  return sizeof(double)*point_stride;
}


void
d4est_solver_schwarz_transfer_ghost_data_and_add_corrections
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_ghost_data_ext_t** d4est_ghost_data_ext,
 double* u,
 double* du_over_subdomains
)
{
  zlog_category_t* c_default = zlog_get_category("d4est_schwarz_solver");
  clock_t begin = clock();
  
  if (p4est->mpisize > 1){
    if (*d4est_ghost_data_ext == NULL){
      *d4est_ghost_data_ext =
        d4est_ghost_data_ext_init
        (
         p4est,
         d4est_ghost,
         1,
         field_size_of_ghost_fcn,
         field_size_of_mirror_fcn,
         field_stride_of_mirror_fcn,
         schwarz_metadata
        );
    }
  
    d4est_ghost_data_ext_exchange
      (
       p4est,
       d4est_ghost,
       *d4est_ghost_data_ext,
       (char**)&du_over_subdomains
      );
  }
  
  /* add local correction du to u */
  for (int i = 0; i < schwarz_metadata->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[i];
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[j];

      if (schwarz_ed->mpirank == p4est->mpirank){
        d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr
                                        (
                                         p4est,
                                         schwarz_ed->tree,
                                         schwarz_ed->tree_quadid
                                        );
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), mesh_ed->deg);
        int local_nodal_stride = mesh_ed->nodal_stride;
        int sub_nodal_stride = sub_data->nodal_stride + schwarz_ed->nodal_stride;
        for (int k = 0; k < volume_nodes; k++){
          u[local_nodal_stride + k] += du_over_subdomains[sub_nodal_stride + k];
          /* printf("du_over_subdomains[sub_nodal_stride + k] = %.15f\n", */
                 /* du_over_subdomains[sub_nodal_stride + k]); */
        }
      }
    }
  }


  int num_ghost_subdomains = schwarz_metadata->subdomain_ghostdata->num_ghosts;

  if (p4est->mpisize > 1){
    /* add ghost correction du to u */
    for (int gid = 0; gid < num_ghost_subdomains; gid++){
      /* get ghost subdomain data */
      d4est_solver_schwarz_subdomain_metadata_t* sub_data =
        (d4est_solver_schwarz_subdomain_metadata_t*)
        d4est_ghost_data_ext_get_field_on_element
        (
         &d4est_ghost->ghost_elements[gid],
         0,
         schwarz_metadata->subdomain_ghostdata
        );
    
      for (int el = 0; el < sub_data->num_elements; el++){
        d4est_solver_schwarz_element_metadata_t* schwarz_ed
          = &sub_data->element_metadata[el];

        /* if ghost subdomain has an element which is local to this 
       * processor we will add the contribution */
        if (schwarz_ed->mpirank == p4est->mpirank){

          /* get ghost subdomain correction data */
          double* du_data =
            (double*)
            d4est_ghost_data_ext_get_field_on_element
            (
             &d4est_ghost->ghost_elements[gid],
             0,
             *d4est_ghost_data_ext
            );

          int schwarz_stride = schwarz_ed->nodal_stride;
          d4est_element_data_t* mesh_ed = d4est_element_data_get_ptr
                                          (p4est,
                                           schwarz_ed->tree,
                                           schwarz_ed->tree_quadid
                                          );
          int mesh_stride = mesh_ed->nodal_stride;
          int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), mesh_ed->deg);
          for (int i = 0; i < volume_nodes; i++){
            u[mesh_stride + i] += du_data[schwarz_stride + i];
          }
        }
      
      }

    }
  }

  clock_t end = clock();
  double time_spent = (double)(end-begin)/CLOCKS_PER_SEC;

  zlog_info(c_default, "Gathering schwarz contributions in %f seconds", time_spent);
}


