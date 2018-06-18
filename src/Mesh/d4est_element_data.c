#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_xyz_functions.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_field.h>
#include <sc_reduce.h>

/* static void */
/* d4est_element_data_copy_from_vec_to_storage_callback */
/* ( */
/*  p4est_iter_volume_info_t * info, */
/*  void *user_data */
/* ) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   d4est_element_data_t* d4est_element_data = (d4est_element_data_t *) q->p.user_data; */
/*   double* u = (double*) info->p4est->user_pointer; */
/*   int* stride = (int*) user_data; */
/*   int dim = (P4EST_DIM); */
/*   int deg = d4est_element_data->deg; */
/*   int volume_nodes = d4est_lgl_get_nodes(dim,deg); */
  
/*   d4est_util_copy_1st_to_2nd */
/*     ( */
/*      &u[*stride], */
/*      &(d4est_element_data->u_elem)[0], */
/*      volume_nodes */
/*     ); */

/*   *stride += volume_nodes; */
/* } */

/* void */
/* d4est_element_data_copy_from_vec_to_storage( */
/*                                              p4est_t* p4est, */
/*                                              double* vec */
/* ) */
/* { */
/*   void* tmp = p4est->user_pointer; */
/*   p4est->user_pointer = vec; */
/*   int stride = 0; */
  
/*   p4est_iterate(p4est, */
/*                 NULL, */
/*                 (void*)&stride, */
/*                 d4est_element_data_copy_from_vec_to_storage_callback, */
/*                 NULL, */
/* #if (P4EST_DIM)==3 */
/*                 NULL,        */
/* #endif                 */
/*                 NULL); */

/*   p4est->user_pointer = tmp; */
/* } */



/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
/* static void */
/* d4est_element_data_copy_from_storage_to_vec_callback */
/* ( */
/*  p4est_iter_volume_info_t * info, */
/*  void *user_data */
/* ) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   d4est_element_data_t* d4est_element_data = (d4est_element_data_t *) q->p.user_data; */
/*   double* u = (double*) info->p4est->user_pointer; */
/*   int* stride = (int*) user_data; */
  
/*   int dim = (P4EST_DIM); */
/*   int deg = d4est_element_data->deg; */
/*   int volume_nodes = d4est_lgl_get_nodes(dim,deg); */
  
/*   d4est_util_copy_1st_to_2nd */
/*     ( */
/*      &(d4est_element_data->u_elem)[0], */
/*      &u[*stride], */
/*      volume_nodes */
/*     ); */

/*   *stride += volume_nodes; */
/* } */

/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
/* void */
/* d4est_element_data_copy_from_storage_to_vec */
/* ( */
/*  p4est_t* p4est, */
/*  double* vec */
/* ) */
/* { */
/*   void* tmp = p4est->user_pointer; */
/*   p4est->user_pointer = vec; */
/*   int stride = 0; */
  
/*   p4est_iterate(p4est, */
/*                 NULL, */
/*                 (void*)&stride, */
/*                 d4est_element_data_copy_from_storage_to_vec_callback, */
/*                 NULL, */
/* #if (P4EST_DIM)==3 */
/*                 NULL,        */
/* #endif                 */
/*                 NULL); */

/*   p4est->user_pointer = tmp; */
/* } */


/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
static void
d4est_element_data_store_element_scalar_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) q->p.user_data;
  double* vertex_array_vec = (double*)user_data;
  double (*get_local_scalar_fcn)(d4est_element_data_t*) = (double (*)(d4est_element_data_t*))info->p4est->user_pointer;  

  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid; 
  p4est_tree_t       *tree;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (info->p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = (P4EST_CHILDREN) * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in */

  double est = get_local_scalar_fcn(elem_data);
  /* int quadid = elem_data->id; */

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    vertex_array_vec[arrayoffset + i] = est;
  }
}


/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
void
d4est_element_data_store_element_scalar_in_vertex_array
(
 p4est_t* p4est,
 double* vertex_array,
 double (*get_local_scalar_fcn)(d4est_element_data_t*)
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = (void*)get_local_scalar_fcn;
  
  p4est_iterate(p4est,
                NULL,
                (void*)vertex_array,
                d4est_element_data_store_element_scalar_in_vertex_array_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);

  p4est->user_pointer = tmp;
}


/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
static void
d4est_element_data_store_nodal_vec_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) q->p.user_data;

  double* nodal_vec = (double*)user_data;
  double* vertex_vec = (double*) info->p4est->user_pointer;
  /* int id = elem_data->id; */
  int stride = elem_data->nodal_stride;

  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid; 
  p4est_tree_t       *tree;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (info->p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = (P4EST_CHILDREN) * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in */
  
  int c;
  for (c = 0; c < (P4EST_CHILDREN); c++){
    int corner_to_node = d4est_reference_corner_to_node((P4EST_DIM), elem_data->deg, c);
    vertex_vec[arrayoffset + c] = nodal_vec[stride + corner_to_node];
  }
}


/**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
void
d4est_element_data_store_nodal_vec_in_vertex_array
(
 p4est_t* p4est,
 double* nodal_vec,
 double* corner_vec
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = corner_vec;
  
  p4est_iterate(p4est,
                NULL,
                nodal_vec,
                d4est_element_data_store_nodal_vec_in_vertex_array_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
                NULL);

  p4est->user_pointer = tmp;
}

 /**
* WILL BE DEPRECATED EVENTUALLY
* 
*/
void
d4est_element_data_reorient_f_p_elements_to_f_m_order
(
 d4est_element_data_t** e_p,
 int face_dim,
 int f_m,
 int f_p,
 int o,
 int faces_p,
 d4est_element_data_t* e_p_oriented [(P4EST_HALF)]
)
{

  if (faces_p == 1){
    e_p_oriented[0] = e_p[0];
    return;
  }
  for (int i = 0; i < faces_p; i++){
    int inew = d4est_reference_reorient_face_order(face_dim, f_m, f_p, o, i);
    e_p_oriented[i] = e_p[inew];
  }
}


int
d4est_element_data_get_size_of_field
(
 d4est_element_data_t* ed,
 d4est_field_type_t type
)
{
  D4EST_FIELD_CHECK_TYPE(type);
  if (type == VOLUME_NODAL){
    return d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
  }
  else if (type == VOLUME){
    return 1;
  }
  else if (type == FACE){
    return (P4EST_FACES);
  }
  else {
    D4EST_ABORT("not a supported type");
  }
}

int d4est_element_data_get_stride_for_field
(
 d4est_element_data_t* ed,
 d4est_field_type_t type
)
{
  D4EST_FIELD_CHECK_TYPE(type);
  if (type == VOLUME_NODAL){
    return ed->nodal_stride;
  }
  else if (type == VOLUME){
    return ed->id;
  }
  else if (type == FACE){
    return ed->id*(P4EST_FACES);
  }
  else {
    D4EST_ABORT("not a supported type");
  }
}
