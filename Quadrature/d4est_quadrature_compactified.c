#include <d4est_operators.h>
#include <d4est_element_data.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_linalg.h>
#include <curved_compute_flux.h>
#include <util.h>

/* store only weights and abscissas and compute rst on the spot */

typedef struct {

  /* Used for volume quantity storage */
  int* volume_strides_1d;
  int* volume_strides_2d;
  int* volume_strides_3d;
  int* volume_deg_lobatto;
  int* volume_deg_quad;
  
  double* volume_interp;
  double* volume_interp_trans;
  double* volume_rst [(P4EST_DIM)];
  double* volume_weights;
  double* volume_abscissas;

  /* Used for trace quantity storage */
  double* mortar_rst [(P4EST_DIM)-1];
  double* mortar_interp;
  double* mortar_interp_trans;
  double* mortar_weights;
  double* mortar_abscissas;

  int* mortar_strides_1d;
  int* mortar_strides_2d;

  int total_mortar_sides;
  int total_mortar_quad_nodes_1d;
  int total_mortar_quad_nodes_2d;

} d4est_quadrature_compactified_storage_t;


static int
check_if_direction_is_compactified
(
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int rst_direction,
 int compactified_direction
){
  if (object_type == QUAD_VOLUME){
    return ((rst_direction == compactified_direction) && (integrand_type == QUAD_JAC_TIMES_POLY_INTEGRAND));
  }
  else if (object_type == QUAD_MORTAR){
    d4est_geometry_face_info_t face_info;
    d4est_geometry_get_face_info(((d4est_quadrature_mortar_t*)object)->face, &face_info);
    if (rst_direction == 0){
      return (face_info.a == compactified_direction) && (integrand_type == QUAD_JAC_TIMES_POLY_INTEGRAND);
    }
    else if (rst_direction == 1){
      return (face_info.b == compactified_direction) && (integrand_type == QUAD_JAC_TIMES_POLY_INTEGRAND);
    }
    else {
      mpi_abort("[D4EST_ERROR]: rst_direction < (P4EST_DIM)-1 when type == ELEMENT_FACE");
      return -1;
    }
  }
  else {
    mpi_abort("[D4EST_ERROR]: type == ELEMENT_FACE or type == ELEMENT_VOLUME");
    return -1;
  }
}

static void
d4est_quadrature_compactified_compute_storage_for_volume
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  int stride_1d_quad = 0;
  int stride_2d_quad = 0;
  int stride_3d_quad = 0;
  int k = 0;

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        
        storage->volume_strides_1d[k] = stride_1d_quad;
        storage->volume_strides_2d[k] = stride_2d_quad;
        storage->volume_strides_3d[k] = stride_3d_quad;
        storage->volume_deg_lobatto[k] = ed->deg;
        storage->volume_deg_quad[k] = ed->deg_quad;

        int nodes_quad_1d = ed->deg_quad + 1;
        
        double* eye = P4EST_ALLOC_ZERO(double, nodes_quad_1d);
        for (int i = 0; i < nodes_quad_1d; i++) eye[i] = 1.;

        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           &storage->volume_abscissas[stride_1d_quad],
           &storage->volume_weights[stride_1d_quad],
           ed->tree,
           ed->q,
           ed->dq,
           ed->deg_quad
          );

        double* rst_volume [(P4EST_DIM)];
        for (int i = 0; i < (P4EST_DIM); i++){
          rst_volume[i] = &storage->volume_rst[i][stride_3d_quad];
        }


        for (int i = 0; i < (P4EST_DIM); i++){
          d4est_quadrature_compactified_compute_rst_volume
            (
             &storage->volume_abscissas[stride_1d_quad],
             ed->deg_quad,
             rst_volume[i],
             i
            );
        }
        
        
        d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
          (
           d4est_ops,
           &storage->volume_abscissas[stride_1d_quad],
           &storage->volume_interp[stride_2d_quad],
           &storage->volume_interp_trans[stride_2d_quad],
           ed->deg,
           ed->deg_quad
          );
       
        int nodes_1d = ed->deg_quad + 1;
        stride_1d_quad += nodes_1d;
        stride_2d_quad += nodes_1d*nodes_1d;
        stride_3d_quad += nodes_1d*nodes_1d*nodes_1d;

        P4EST_FREE(eye);
        k++;
      }
    }    
}

static double*
d4est_quadrature_compactified_get_weights
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int degree,
 int rst_direction
){

  int compactified_direction = 2;
  int is_it_compactified = check_if_direction_is_compactified(object,
                                                              object_type,
                                                              integrand_type,
                                                              rst_direction,
                                                              compactified_direction);
  mpi_assert((rst_direction < (P4EST_DIM) && rst_direction >= 0) || (rst_direction < (P4EST_DIM) - 1 && rst_direction >= 0));

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (!is_it_compactified){
    return d4est_operators_fetch_GL_weights_1d(d4est_ops, degree);
  }
  else{
    if (object_type == QUAD_VOLUME){
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_1d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      mpi_assert(degree == deg_quad_stored);      
      return &storage->volume_weights[stride];
    }
    else if (object_type == QUAD_MORTAR){
      int mortar_side_id = ((d4est_quadrature_mortar_t*)object)->mortar_side_id;
      int mortar_subface_id = ((d4est_quadrature_mortar_t*)object)->mortar_subface_id;
      int stride = storage->mortar_strides_1d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
      return &storage->mortar_weights[stride];
    }
    else {
      mpi_abort("[D4EST_ERROR]: object_type must be QUAD_VOLUME OR QUAD_MORTAR");
    }
    
  }
}

static double*
d4est_quadrature_compactified_get_rst
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int degree,
 int rst_direction
){

  int compactified_direction = 2;
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (object_type == QUAD_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);

    if (!is_it_compactified){
      return d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), degree, rst_direction);
    }
    else{
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_3d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      mpi_assert(degree == deg_quad_stored);      
      return &(storage->volume_rst[rst_direction][stride]);
    }
  }
  else if (object_type == QUAD_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM)-1, degree, rst_direction);      
    }
    else {
      mpi_assert(rst_direction == 1 && compactified_direction == 2); /* sanity check for now */
      int mortar_side_id = ((d4est_quadrature_mortar_t*)object)->mortar_side_id;
      int mortar_subface_id = ((d4est_quadrature_mortar_t*)object)->mortar_subface_id;
#if (P4EST_DIM)==3
      int stride = storage->mortar_strides_2d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
#else
      int stride = storage->mortar_strides_1d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
#endif
      return &(storage->mortar_rst[rst_direction][stride]);
    }
  } 
 else {
   mpi_abort("[D4EST_ERROR]: ELEMENT_TYPE not supported");
 }
}


static double*
d4est_quadrature_compactified_get_interp
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int deg_lobatto,
 int deg_quad,
 int rst_direction
){

  int compactified_direction = 2;
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;

  if (object_type == QUAD_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);
    if (!is_it_compactified){
      return d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else{
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_2d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      int deg_lobatto_stored = storage->volume_deg_lobatto[local_id];
      mpi_assert(deg_quad == deg_quad_stored);      
      mpi_assert(deg_lobatto == deg_lobatto_stored);      
      return &(storage->volume_interp[stride]);
    }
  }
  else if (object_type == QUAD_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else {
      mpi_assert(rst_direction == 1 && compactified_direction == 2); /* sanity check for now */
      int mortar_side_id = ((d4est_quadrature_mortar_t*)object)->mortar_side_id;
      int mortar_subface_id = ((d4est_quadrature_mortar_t*)object)->mortar_subface_id;
      int stride = storage->mortar_strides_2d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
      return &(storage->mortar_interp[stride]);
    }
  } 
 else {
   mpi_abort("[D4EST_ERROR]: ELEMENT_TYPE not supported");
 }
  

}

static double*
d4est_quadrature_compactified_get_interp_trans
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int deg_lobatto,
 int deg_quad,
 int rst_direction
 
){
  int compactified_direction = 2;
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (object_type == QUAD_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);
    if (!is_it_compactified){
      return d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else{
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_2d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      int deg_lobatto_stored = storage->volume_deg_lobatto[local_id];
      mpi_assert(deg_quad == deg_quad_stored);      
      mpi_assert(deg_lobatto == deg_lobatto_stored);      
      return &(storage->volume_interp_trans[stride]);
    }
  }
  else if (object_type == QUAD_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else {
      mpi_assert(rst_direction == 1 && compactified_direction == 2); /* sanity check for now */
      int mortar_side_id = ((d4est_quadrature_mortar_t*)object)->mortar_side_id;
      int mortar_subface_id = ((d4est_quadrature_mortar_t*)object)->mortar_subface_id;
      int stride = storage->mortar_strides_2d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
      return &(storage->mortar_interp_trans[stride]);
    }
  } 
 else {
   mpi_abort("[D4EST_ERROR]: ELEMENT_TYPE not supported");
 }
}


void
d4est_quadrature_compactified_compute_abscissas_and_weights
(
 d4est_geometry_t* d4est_geom,
 double* abscissas,
 double* weights,
 p4est_topidx_t tree,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int degree
)
{
  d4est_geometry_cubed_sphere_outer_shell_block_get_custom_quadrature
    (
     d4est_geom,
     tree,
     q,
     dq,
     degree,
     abscissas,
     weights,
     0 /* do not test moments */
    );
}

void
d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
(
 d4est_operators_t* d4est_ops,
 double* abscissas,
 double* interp,
 double* interp_trans,
 int deg,
 int deg_quad
){
  d4est_operators_build_custom_GLL_interp_1d
    (
     d4est_ops,
     interp,
     deg,
     deg_quad,
     abscissas
    );

  d4est_linalg_mat_transpose_nonsqr
    (
     interp,
     interp_trans,
     deg_quad + 1,
     deg + 1
    );  
}

void
d4est_quadrature_compactified_compute_rst_volume
(
 double* abscissas,
 int deg_quad,
 double* rst_volume,
 int rst_direction
)
{
  int nodes_quad_1d = deg_quad + 1;
  double* eye = P4EST_ALLOC_ZERO(double, nodes_quad_1d);
  for (int i = 0; i < nodes_quad_1d; i++) eye[i] = 1.;


  if (rst_direction == 2){
  d4est_linalg_kron_AoBoC
    (
     abscissas,
     eye,
     eye,
     rst_volume,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1
    );
  }
  else if (rst_direction == 1){
  d4est_linalg_kron_AoBoC
    (
     eye,
     abscissas,
     eye,
     rst_volume,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1
    );
  }
  else if (rst_direction == 0){
  d4est_linalg_kron_AoBoC
    (
     eye,
     eye,
     abscissas,
     rst_volume,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1,
     nodes_quad_1d, 1
    );
  }
  else {
    mpi_abort("[D4EST_ERROR]: rst_direction != 0,1 or 2");
  }

  P4EST_FREE(eye);
}

void
d4est_quadrature_compactified_compute_weights_and_abscissas
(
 d4est_geometry_t* d4est_geom,
 double* abscissas,
 double* weights,
 p4est_topidx_t tree,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int degree
)
{
  d4est_geometry_cubed_sphere_outer_shell_block_get_custom_quadrature
    (
     d4est_geom,
     tree,
     q,
     dq,
     degree,
     abscissas,
     weights,
     0 /* do not test moments */
    );
}


void
d4est_quadrature_compactified_compute_rst_face
(
 double* abscissas,
 int deg_quad,
 double* rst_face [(P4EST_DIM)-1]
)
{
  int nodes_quad_1d = deg_quad + 1;
  double* eye = P4EST_ALLOC_ZERO(double, nodes_quad_1d);
  for (int i = 0; i < nodes_quad_1d; i++) eye[i] = 1.;
  d4est_linalg_kron_AoB(eye,
                  abscissas,
                  rst_face[0],
                  nodes_quad_1d,
                  1,
                  nodes_quad_1d,
                  1
                 );
        
  d4est_linalg_kron_AoB(abscissas,
                  eye,
                  rst_face[1],
                  nodes_quad_1d,
                  1,
                  nodes_quad_1d,
                  1
                 );

  P4EST_FREE(eye);
}       




static d4est_quadrature_compactified_storage_t*
d4est_quadrature_compactified_storage_init()
{
  d4est_quadrature_compactified_storage_t* storage = P4EST_ALLOC(d4est_quadrature_compactified_storage_t, 1);
  
  storage->volume_strides_1d = NULL;
  storage->volume_strides_2d = NULL;
  storage->volume_strides_3d = NULL;
  storage->volume_interp = NULL;
  storage->volume_interp_trans = NULL;
  storage->volume_weights = NULL;
  storage->volume_abscissas = NULL;
  storage->volume_deg_lobatto = NULL;
  storage->volume_deg_quad = NULL;
  for (int i = 0; i < (P4EST_DIM); i++){
    storage->volume_rst[i] = NULL;
  }

  storage->mortar_abscissas = NULL;
  storage->mortar_strides_1d = NULL;
  storage->mortar_strides_2d = NULL;
  storage->mortar_weights = NULL;
  storage->mortar_interp = NULL;
  storage->mortar_interp_trans = NULL;
  for (int i = 0; i < (P4EST_DIM)-1; i++){
    storage->mortar_rst[i] = NULL;
  }
  return storage;
}

void
d4est_quadrature_compactified_storage_destroy
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  P4EST_FREE(storage->volume_strides_1d);
  P4EST_FREE(storage->volume_strides_2d);
  P4EST_FREE(storage->volume_strides_3d);
  P4EST_FREE(storage->volume_deg_lobatto);
  P4EST_FREE(storage->volume_deg_quad);
  P4EST_FREE(storage->volume_interp);
  P4EST_FREE(storage->volume_interp_trans);
  P4EST_FREE(storage->volume_weights);
  P4EST_FREE(storage->volume_abscissas);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(storage->volume_rst[i]);
  }

  P4EST_FREE(storage->mortar_abscissas);
  P4EST_FREE(storage->mortar_strides_1d);
  P4EST_FREE(storage->mortar_strides_2d);
  P4EST_FREE(storage->mortar_weights);
  P4EST_FREE(storage->mortar_interp);
  P4EST_FREE(storage->mortar_interp_trans);
  for (int i = 0; i < (P4EST_DIM)-1; i++){
    P4EST_FREE(storage->mortar_rst[i]);
  }
  P4EST_FREE(storage);
}


void
d4est_quadrature_compactified_new
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 const char* input_file
)
{
  mpi_assert(d4est_geom->geom_type == GEOM_CUBED_SPHERE_OUTER_SHELL);
  
  d4est_quad->get_weights = d4est_quadrature_compactified_get_weights;
  d4est_quad->get_rst = d4est_quadrature_compactified_get_rst;
  d4est_quad->get_interp = d4est_quadrature_compactified_get_interp;
  d4est_quad->get_interp_trans = d4est_quadrature_compactified_get_interp_trans;

  d4est_quad->user_destroy = d4est_quadrature_compactified_storage_destroy;
  d4est_quad->user_reinit = d4est_quadrature_compactified_setup_storage;

  d4est_quad->user = d4est_quadrature_compactified_storage_init();
}





static void
d4est_quadrature_compactified_compute_mortar_strides_and_sizes_boundary
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  /* WE WILL NEED TO SKIP CERTAIN TREES IN THE FUTURE WHEN WE GO TO
   TO CUBED SPHERES*/

  
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;

  storage->mortar_strides_1d[mortar_side_id_m*(P4EST_HALF)] = storage->total_mortar_quad_nodes_1d;
  storage->mortar_strides_2d[mortar_side_id_m*(P4EST_HALF)] = storage->total_mortar_quad_nodes_2d;
  
  storage->total_mortar_quad_nodes_1d += d4est_operators_get_nodes(1, e_m->deg_quad);  
  storage->total_mortar_quad_nodes_2d += d4est_operators_get_nodes(2, e_m->deg_quad);
  storage->total_mortar_sides += 1;
}

static void
d4est_quadrature_compactified_store_data_for_boundary
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  /* WE WILL NEED TO SKIP CERTAIN TREES IN THE FUTURE WHEN WE GO TO
   TO CUBED SPHERES*/
  
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  int mortar_stride_2d = storage->mortar_strides_2d[mortar_side_id_m*(P4EST_HALF)];
  int mortar_stride_1d = storage->mortar_strides_1d[mortar_side_id_m*(P4EST_HALF)];
  int mortar_stride_face = ((P4EST_DIM)==3) ? mortar_stride_2d : mortar_stride_1d;
  
  d4est_quadrature_compactified_compute_abscissas_and_weights
    (
     d4est_geom,
     &storage->mortar_abscissas[mortar_stride_1d],
     &storage->mortar_weights[mortar_stride_1d],
     e_m->tree,
     e_m->q,
     e_m->dq,
     e_m->deg_quad
    );

  double* rst_interface [(P4EST_DIM)-1];
  for (int i = 0; i < (P4EST_DIM)-1; i++){
    rst_interface[i] = &storage->mortar_rst[i][mortar_stride_face];
  }
  
  d4est_quadrature_compactified_compute_rst_face
    (
     &storage->mortar_abscissas[mortar_stride_1d],
     e_m->deg_quad,
     rst_interface
    );

  d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
    (
     d4est_ops,
     &storage->mortar_abscissas[mortar_stride_1d],
     &storage->mortar_interp[mortar_stride_2d],
     &storage->mortar_interp_trans[mortar_stride_2d],
     e_m->deg,
     e_m->deg_quad
    );
}

static void
d4est_quadrature_compactified_compute_mortar_strides_and_sizes_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
      
  for (int i = 0; i < faces_m; i++){
    for (int j = 0; j < faces_p; j++){
      int deg_mortar_quad = util_max_int(e_m[i]->deg_quad,e_p_oriented[j]->deg_quad);
      storage->mortar_strides_1d[mortar_side_id_m*(P4EST_HALF) + i + j] = storage->total_mortar_quad_nodes_1d;
      storage->mortar_strides_2d[mortar_side_id_m*(P4EST_HALF) + i + j] = storage->total_mortar_quad_nodes_2d;
      storage->total_mortar_quad_nodes_2d += d4est_operators_get_nodes(2, deg_mortar_quad);
      storage->total_mortar_quad_nodes_1d += d4est_operators_get_nodes(1, deg_mortar_quad);
    }
  }

  storage->total_mortar_sides += 1;  
}


static void
d4est_quadrature_compactified_store_data_for_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  int deg_mortar_quad [(P4EST_HALF)];
  int num_faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
    
  for (int i = 0; i < faces_m; i++){
    for (int j = 0; j < faces_p; j++){
      deg_mortar_quad[i+j] = util_max_int(e_m[i]->deg_quad,e_p_oriented[j]->deg_quad);
    }
  }

  p4est_qcoord_t mortar_q0 [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t mortar_dq;

  d4est_geometry_compute_qcoords_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q, /* qcoord of first element of side */
     e_m[0]->dq, /* qcoord vector spanning first element of side */
     faces_m, 
     num_faces_mortar,
     f_m,
     mortar_q0,
     &mortar_dq
    );

  p4est_qcoord_t q [(P4EST_DIM)];
  
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){

    int mortar_stride_1d = storage->mortar_strides_1d[(P4EST_HALF)*mortar_side_id_m + face_mortar];
    int mortar_stride_2d = storage->mortar_strides_2d[(P4EST_HALF)*mortar_side_id_m + face_mortar];
    int mortar_stride_face = ((P4EST_DIM)==3) ? mortar_stride_2d : mortar_stride_1d;

    for (int d = 0; d < (P4EST_DIM); d++){
      q[d] = mortar_q0[face_mortar][d];
    }
    
  d4est_quadrature_compactified_compute_abscissas_and_weights
    (
     d4est_geom,
     &storage->mortar_abscissas[mortar_stride_1d],
     &storage->mortar_weights[mortar_stride_1d],
     e_m[0]->tree,
     q,
     mortar_dq,
     deg_mortar_quad[face_mortar]
    );

  double* rst_interface [(P4EST_DIM)-1];
  for (int i = 0; i < (P4EST_DIM)-1; i++){
    rst_interface[i] = &storage->mortar_rst[i][mortar_stride_face];
  }
  
  d4est_quadrature_compactified_compute_rst_face
    (
     &storage->mortar_abscissas[mortar_stride_1d],
     deg_mortar_quad[face_mortar],
     rst_interface
    );

  d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
    (
     d4est_ops,
     &storage->mortar_abscissas[mortar_stride_1d],
     &storage->mortar_interp[mortar_stride_2d],
     &storage->mortar_interp_trans[mortar_stride_2d],
     deg_mortar_quad[face_mortar],
     deg_mortar_quad[face_mortar]
    );

  }  
}

static void
d4est_quadrature_compactified_compute_storage_for_mortar_space
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  storage->mortar_strides_1d = P4EST_REALLOC(storage->mortar_strides_1d, int, (P4EST_HALF)*(P4EST_FACES)*p4est->local_num_quadrants);
  storage->mortar_strides_2d = P4EST_REALLOC(storage->mortar_strides_2d, int, (P4EST_HALF)*(P4EST_FACES)*p4est->local_num_quadrants);
  
  storage->total_mortar_quad_nodes_1d = 0;
  storage->total_mortar_quad_nodes_2d = 0;
  storage->total_mortar_sides = 0;
  
  curved_flux_fcn_ptrs_t compute_strides_and_sizes;
  compute_strides_and_sizes.flux_interface_fcn = d4est_quadrature_compactified_compute_mortar_strides_and_sizes_interface;
  compute_strides_and_sizes.flux_boundary_fcn = d4est_quadrature_compactified_compute_mortar_strides_and_sizes_boundary;
  compute_strides_and_sizes.bndry_fcn = NULL;
  compute_strides_and_sizes.params = NULL;

  curved_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &compute_strides_and_sizes,
     EXCHANGE_GHOST_DATA
    );

  storage->mortar_interp = P4EST_REALLOC(storage->mortar_interp, double, storage->total_mortar_quad_nodes_2d);
  storage->mortar_interp_trans = P4EST_REALLOC(storage->mortar_interp_trans, double, storage->total_mortar_quad_nodes_2d);
  storage->mortar_weights = P4EST_REALLOC(storage->mortar_weights, double, storage->total_mortar_quad_nodes_1d);
  storage->mortar_abscissas = P4EST_REALLOC(storage->mortar_abscissas, double, storage->total_mortar_quad_nodes_1d);
  
  for (int i = 0; i < (P4EST_DIM)-1; i++){
#if (P4EST_DIM)==3
    storage->mortar_rst[i] = P4EST_REALLOC(storage->mortar_rst[i], double, storage->total_mortar_quad_nodes_2d);
#else
    storage->mortar_rst[i] = P4EST_REALLOC(storage->mortar_rst[i], double, storage->total_mortar_quad_nodes_1d);
#endif
  }
  
  
  curved_flux_fcn_ptrs_t compute_storage;
  compute_storage.flux_interface_fcn = d4est_quadrature_compactified_store_data_for_interface;
  compute_storage.flux_boundary_fcn = d4est_quadrature_compactified_store_data_for_boundary;
  compute_storage.bndry_fcn = NULL;
  compute_storage.params = NULL;  

  curved_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &compute_storage,
     DO_NOT_EXCHANGE_GHOST_DATA
    ); 
}


void
d4est_quadrature_compactified_setup_storage
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  int nodes_1d_quad = 0;
  int nodes_2d_quad = 0;
  int nodes_3d_quad = 0;
  int k = 0;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        k++;
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;

        int nodes_1d = ed->deg_quad + 1;
        nodes_1d_quad += nodes_1d;
        nodes_2d_quad += nodes_1d*nodes_1d;
        nodes_3d_quad += nodes_1d*nodes_1d*nodes_1d;
      }
    }    

  storage->volume_strides_1d = P4EST_REALLOC(storage->volume_strides_1d, int, p4est->local_num_quadrants);
  storage->volume_strides_2d = P4EST_REALLOC(storage->volume_strides_2d, int, p4est->local_num_quadrants);
  storage->volume_strides_3d = P4EST_REALLOC(storage->volume_strides_3d, int, p4est->local_num_quadrants);
  storage->volume_deg_quad = P4EST_REALLOC(storage->volume_deg_quad, int, p4est->local_num_quadrants);
  storage->volume_deg_lobatto = P4EST_REALLOC(storage->volume_deg_lobatto, int, p4est->local_num_quadrants);

  storage->volume_interp = P4EST_REALLOC(storage->volume_interp, double, nodes_2d_quad);
  storage->volume_interp_trans = P4EST_REALLOC(storage->volume_interp_trans, double, nodes_2d_quad);
  storage->volume_weights = P4EST_REALLOC(storage->volume_weights, double, nodes_1d_quad);
  storage->volume_abscissas = P4EST_REALLOC(storage->volume_abscissas, double, nodes_1d_quad);

  for (int i = 0; i < (P4EST_DIM); i++){
    storage->volume_rst[i] = P4EST_REALLOC(storage->volume_rst[i], double, nodes_3d_quad);
  }
  
  d4est_quadrature_compactified_compute_storage_for_volume(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_quadrature_compactified_compute_storage_for_mortar_space(p4est, ghost, ghost_data, d4est_ops, d4est_geom, d4est_quad);
  
}
