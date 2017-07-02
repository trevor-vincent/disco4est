#include <d4est_operators.h>
#include <d4est_element_data.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_geometry_disk.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_mortars.h>
#include <d4est_linalg.h>
#include <arbquad.h>
#include <util.h>

/* store only weights and abscissas and compute rst on the spot */

typedef struct {

  long double c1;
  long double c2;

} d4est_quadrature_compactified_params_t;

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
get_compactified_direction
(
 d4est_geometry_t* d4est_geom,
 void* object,
 d4est_quadrature_object_type_t object_type
){

  int tree;
  if (object_type == QUAD_OBJECT_MORTAR){
    tree = ((d4est_quadrature_mortar_t*)object)->tree;
  }
  else if (object_type == QUAD_OBJECT_VOLUME){
    tree = ((d4est_quadrature_volume_t*)object)->tree;
  }
  else {
    printf("[D4EST_ERROR]: Not a supported object_type");
  }
  
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_OUTER_SHELL){
    mpi_assert(tree == 0);
    return 2;
  }
  else if (d4est_geom->geom_type == GEOM_DISK_OUTER_WEDGE){
    mpi_assert(tree == 0);
    return 1;
  }
  else {
    printf("[D4EST_ERROR]: compactified quadrature does not support this type");
    return -1;
  }
}

static int
check_if_direction_is_compactified
(
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int rst_direction,
 int compactified_direction
){
  if (object_type == QUAD_OBJECT_VOLUME){
    return ((rst_direction == compactified_direction) && (integrand_type == QUAD_INTEGRAND_UNKNOWN));
  }
  else if (object_type == QUAD_OBJECT_MORTAR){
    d4est_geometry_face_info_t face_info;
    d4est_geometry_get_face_info(((d4est_quadrature_mortar_t*)object)->face, &face_info);
    if (rst_direction == 0){
      return (face_info.a == compactified_direction) && (integrand_type == QUAD_INTEGRAND_UNKNOWN);
    }
    else if (rst_direction == 1){
      return (face_info.b == compactified_direction) && (integrand_type == QUAD_INTEGRAND_UNKNOWN);
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
           d4est_quad,
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

  int compactified_direction = get_compactified_direction(d4est_geom, object, object_type);
  int is_it_compactified = check_if_direction_is_compactified(object,
                                                              object_type,
                                                              integrand_type,
                                                              rst_direction,
                                                              compactified_direction);
  mpi_assert((rst_direction < (P4EST_DIM) && rst_direction >= 0) || (rst_direction < (P4EST_DIM) - 1 && rst_direction >= 0));

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (!is_it_compactified){
    return d4est_operators_fetch_gauss_weights_1d(d4est_ops, degree);
  }
  else{
    if (object_type == QUAD_OBJECT_VOLUME){
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_1d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      mpi_assert(degree == deg_quad_stored);      
      return &storage->volume_weights[stride];
    }
    else if (object_type == QUAD_OBJECT_MORTAR){
      int mortar_side_id = ((d4est_quadrature_mortar_t*)object)->mortar_side_id;
      int mortar_subface_id = ((d4est_quadrature_mortar_t*)object)->mortar_subface_id;
      int stride = storage->mortar_strides_1d[(P4EST_HALF)*mortar_side_id + mortar_subface_id];
      return &storage->mortar_weights[stride];
    }
    else {
      mpi_abort("[D4EST_ERROR]: object_type must be QUAD_OBJECT_VOLUME OR QUAD_OBJECT_MORTAR");
    }
    
  }
}

int
d4est_quadrature_compactified_check_type
(
 d4est_quadrature_t* d4est_quad
)
{
    if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG1){
      return 1;
    }
    else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG2){
      return 1;
    }
    else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG3){
      return 1;
    }
    else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG4){
      return 1;
    }
    else {
      return 0;
    }
}
 
  /* Should be set to static after debugging is done */
double*
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

  int compactified_direction = get_compactified_direction(d4est_geom, object, object_type);
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);
  
  /* if(!check_if_using_compactified_quadrature(d4est_quad)){ */
  /*   return NULL; */
  /* } */
    
  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (object_type == QUAD_OBJECT_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);

    if (!is_it_compactified){
      return d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), degree, rst_direction);
    }
    else{
      int local_id = ((d4est_quadrature_volume_t*)object)->element_id;
      int stride = storage->volume_strides_3d[local_id];
      int deg_quad_stored = storage->volume_deg_quad[local_id];
      mpi_assert(degree == deg_quad_stored);      
      return &(storage->volume_rst[rst_direction][stride]);
    }
  }
  else if (object_type == QUAD_OBJECT_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, degree, rst_direction);      
    }
    else {
      /* mpi_assert(rst_direction == 1 && compactified_direction == 2); /\* sanity check for now *\/ */
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

  int compactified_direction = get_compactified_direction(d4est_geom, object, object_type);
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;

  if (object_type == QUAD_OBJECT_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);
    if (!is_it_compactified){
      return d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg_lobatto, deg_quad);
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
  else if (object_type == QUAD_OBJECT_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else {
      /* mpi_assert(rst_direction == 1 && compactified_direction == 2); /\* sanity check for now *\/ */
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
  int compactified_direction = get_compactified_direction(d4est_geom, object, object_type);
  int is_it_compactified = check_if_direction_is_compactified(object, object_type, integrand_type, rst_direction, compactified_direction);

  d4est_quadrature_compactified_storage_t* storage = d4est_quad->user;
  
  if (object_type == QUAD_OBJECT_VOLUME){
    mpi_assert(rst_direction < (P4EST_DIM) && rst_direction >= 0);
    if (!is_it_compactified){
      return d4est_operators_fetch_lobatto_to_gauss_interp_trans_1d(d4est_ops, deg_lobatto, deg_quad);
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
  else if (object_type == QUAD_OBJECT_MORTAR){

    if (!is_it_compactified){
      return d4est_operators_fetch_lobatto_to_gauss_interp_trans_1d(d4est_ops, deg_lobatto, deg_quad);
    }
    else {
      /* mpi_assert(rst_direction == 1 && compactified_direction == 2); /\* sanity check for now *\/ */
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
d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
(
 d4est_operators_t* d4est_ops,
 double* abscissas,
 double* interp,
 double* interp_trans,
 int deg,
 int deg_quad
){
  d4est_operators_build_custom_lobatto_interp_1d
    (
     d4est_ops,
     interp,
     deg,
     deg_quad,
     abscissas
    );

  if (interp_trans != NULL){
    d4est_linalg_mat_transpose_nonsqr
      (
       interp,
       interp_trans,
       deg_quad + 1,
       deg + 1
      );
  }
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
d4est_quadrature_compactified_compute_rst_face
(
 double* abscissas,
 int deg_quad,
 double* rst_face [(P4EST_DIM)-1]
)
{
  int nodes_quad_1d = deg_quad + 1;
#if(P4EST_DIM)==3
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
#else
  d4est_linalg_copy_1st_to_2nd(abscissas, rst_face[0], nodes_quad_1d);
#endif
}       




static d4est_quadrature_compactified_storage_t*
d4est_quadrature_compactified_storage_init()
{
  d4est_quadrature_compactified_storage_t* storage = P4EST_ALLOC(d4est_quadrature_compactified_storage_t, 1);

  /* if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_OUTER_SHELL){ */
  /*   int num_trees = 1; */
  /*   storage->compactified_directions = P4EST_ALLOC(int,num_trees); */
  /*   storage->compactified_directions[0] = 2; */
  /* } */
  /* else if (d4est_geom->geom_type == GEOM_DISK_OUTER_WEDGE){ */
  /*   int num_trees = 1; */
  /*   storage->compactified_directions = P4EST_ALLOC(int,num_trees); */
  /*   storage->compactified_directions[0] = 1;     */
  /* } */
  /* else { */
  /*   printf("[D4EST_ERROR]: compactified quadrature does not support this type"); */
  /* } */
  
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
 const char* input_file,
 const char* input_section
)
{
  mpi_assert(d4est_geom->geom_type == GEOM_CUBED_SPHERE_OUTER_SHELL
             ||
             d4est_geom->geom_type == GEOM_DISK_OUTER_WEDGE
            );

  
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
 d4est_grid_fcn_t bndry_fcn,
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
  
  storage->total_mortar_quad_nodes_1d += d4est_lgl_get_nodes(1, e_m->deg_quad);  
  storage->total_mortar_quad_nodes_2d += d4est_lgl_get_nodes(2, e_m->deg_quad);
  storage->total_mortar_sides += 1;
}

static void
d4est_quadrature_compactified_store_data_for_boundary
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_grid_fcn_t bndry_fcn,
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
     d4est_quad,
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

  double* tmp = &storage->mortar_rst[0][mortar_stride_face];
  printf("storaage for f_m = %d\n", f_m);
  DEBUG_PRINT_ARR_DBL(rst_interface[0], d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg_quad));
  DEBUG_PRINT_ARR_DBL(tmp, d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg_quad));
  
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
      storage->total_mortar_quad_nodes_2d += d4est_lgl_get_nodes(2, deg_mortar_quad);
      storage->total_mortar_quad_nodes_1d += d4est_lgl_get_nodes(1, deg_mortar_quad);
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

  d4est_mortars_compute_qcoords_on_mortar
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
       d4est_quad,
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
  
  d4est_mortar_fcn_ptrs_t compute_strides_and_sizes;
  compute_strides_and_sizes.flux_interface_fcn = d4est_quadrature_compactified_compute_mortar_strides_and_sizes_interface;
  compute_strides_and_sizes.flux_boundary_fcn = d4est_quadrature_compactified_compute_mortar_strides_and_sizes_boundary;
  compute_strides_and_sizes.bndry_fcn = NULL;
  compute_strides_and_sizes.params = NULL;

  d4est_mortar_compute_flux_on_local_elements
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
  
  
  d4est_mortar_fcn_ptrs_t compute_storage;
  compute_storage.flux_interface_fcn = d4est_quadrature_compactified_store_data_for_interface;
  compute_storage.flux_boundary_fcn = d4est_quadrature_compactified_store_data_for_boundary;
  compute_storage.bndry_fcn = NULL;
  compute_storage.params = NULL;  

  d4est_mortar_compute_flux_on_local_elements
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




void
d4est_quadrature_compactified_c1tpc2_neg4_aa_and_bb
(
 int n,
 long double* aa,
 long double* bb,
 void* user
)
{
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  
  if (n == 1){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    bb[0] = 0;
  }
  else if (n == 2){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
  }
  else if (n == 3){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*powl(c2,3)*(powl(c1,2) - 2*powl(c2,2)) + (27*powl(c1,4)*powl(c2,2) - 39*powl(c1,2)*powl(c2,4) - 2*powl(c2,6))*atanhl(c2/c1) + (-27*powl(c1,5)*c2 + 24*powl(c1,3)*powl(c2,3) + 7*c1*powl(c2,5))*powl(atanhl(c2/c1),2) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,4)*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2)));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
    bb[2] = -(((c1 - c2)*(c1 + c2)*(3*powl(c1,2) + powl(c2,2))*(3*powl(c1,2)*powl(c2,2) - 4*powl(c2,4) + (-6*powl(c1,3)*c2 + 6*c1*powl(c2,3))*atanhl(c2/c1) + (3*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) - powl(c2,4))*powl(atanhl(c2/c1),2)))/powl(c2,8));
  }
  else if (n == 4){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*powl(c2,3)*(powl(c1,2) - 2*powl(c2,2)) + (27*powl(c1,4)*powl(c2,2) - 39*powl(c1,2)*powl(c2,4) - 2*powl(c2,6))*atanhl(c2/c1) + (-27*powl(c1,5)*c2 + 24*powl(c1,3)*powl(c2,3) + 7*c1*powl(c2,5))*powl(atanhl(c2/c1),2) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,4)*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2)));
    aa[3] = -(((c1 - c2)*(c1 + c2)*(-9*(3*powl(c1,5)*powl(c2,5) + 16*powl(c1,3)*powl(c2,7) - 16*c1*powl(c2,9)) + powl(c2,4)*(135*powl(c1,6) + 351*powl(c1,4)*powl(c2,2) - 576*powl(c1,2)*powl(c2,4) + 112*powl(c2,6))*atanhl(c2/c1) - 2*(135*powl(c1,7)*powl(c2,3) - 189*powl(c1,3)*powl(c2,7) + 56*c1*powl(c2,9))*powl(atanhl(c2/c1),2) + 2*powl(c2,2)*(135*powl(c1,8) - 315*powl(c1,6)*powl(c2,2) + 369*powl(c1,4)*powl(c2,4) - 245*powl(c1,2)*powl(c2,6) + 56*powl(c2,8))*powl(atanhl(c2/c1),3) + (-135*powl(c1,9)*c2 + 540*powl(c1,7)*powl(c2,3) - 882*powl(c1,5)*powl(c2,5) + 652*powl(c1,3)*powl(c2,7) - 175*c1*powl(c2,9))*powl(atanhl(c2/c1),4) + 27*powl(powl(c1,2) - powl(c2,2),4)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),5) - 36*c1*c2*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),6)))/((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
    bb[2] = -(((c1 - c2)*(c1 + c2)*(3*powl(c1,2) + powl(c2,2))*(3*powl(c1,2)*powl(c2,2) - 4*powl(c2,4) + (-6*powl(c1,3)*c2 + 6*c1*powl(c2,3))*atanhl(c2/c1) + (3*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) - powl(c2,4))*powl(atanhl(c2/c1),2)))/powl(c2,8));
    bb[3] = (-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))/(3.*powl(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2),2));
  }

  else if (n == 5){

    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));

    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));

    aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*powl(c2,3)*(powl(c1,2) - 2*powl(c2,2)) + (27*powl(c1,4)*powl(c2,2) - 39*powl(c1,2)*powl(c2,4) - 2*powl(c2,6))*atanhl(c2/c1) + (-27*powl(c1,5)*c2 + 24*powl(c1,3)*powl(c2,3) + 7*c1*powl(c2,5))*powl(atanhl(c2/c1),2) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,4)*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2)));

    aa[3] = -(((c1 - c2)*(c1 + c2)*(-9*(3*powl(c1,5)*powl(c2,5) + 16*powl(c1,3)*powl(c2,7) - 16*c1*powl(c2,9)) + powl(c2,4)*(135*powl(c1,6) + 351*powl(c1,4)*powl(c2,2) - 576*powl(c1,2)*powl(c2,4) + 112*powl(c2,6))*atanhl(c2/c1) - 2*(135*powl(c1,7)*powl(c2,3) - 189*powl(c1,3)*powl(c2,7) + 56*c1*powl(c2,9))*powl(atanhl(c2/c1),2) + 2*powl(c2,2)*(135*powl(c1,8) - 315*powl(c1,6)*powl(c2,2) + 369*powl(c1,4)*powl(c2,4) - 245*powl(c1,2)*powl(c2,6) + 56*powl(c2,8))*powl(atanhl(c2/c1),3) + (-135*powl(c1,9)*c2 + 540*powl(c1,7)*powl(c2,3) - 882*powl(c1,5)*powl(c2,5) + 652*powl(c1,3)*powl(c2,7) - 175*c1*powl(c2,9))*powl(atanhl(c2/c1),4) + 27*powl(powl(c1,2) - powl(c2,2),4)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),5) - 36*c1*c2*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),6)))/((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))));

    aa[4] = -((81*powl(c2,2)*powl(-powl(c1,2) + powl(c2,2),5)*(45315*powl(c1,10) - 55485*powl(c1,8)*powl(c2,2) + 23902*powl(c1,6)*powl(c2,4) - 426*powl(c1,4)*powl(c2,6) + 111*powl(c1,2)*powl(c2,8) + 23*powl(c2,10))*powl(atanhl(c2/c1),7) + 306180*powl(c1,6)*powl(powl(c1,2) - powl(c2,2),8)*powl(atanhl(c2/c1),9) + 486*c1*powl(powl(c1,2) - powl(c2,2),8)*powl(atanhl(c2/c1),8)*(-10*powl(c1,2)*powl(c2,3) - 6*powl(c2,5) + 315*powl(c1,5)*logl(c1 - c2) - 315*powl(c1,5)*logl(c1 + c2)) + 9*powl(c2,3)*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),5)*(1029105*powl(c1,14)*c2 - 2534220*powl(c1,12)*powl(c2,3) + 2920905*powl(c1,10)*powl(c2,5) - 1812408*powl(c1,8)*powl(c2,7) + 526595*powl(c1,6)*powl(c2,9) - 22092*powl(c1,4)*powl(c2,11) - 1869*powl(c1,2)*powl(c2,13) + 1504*powl(c2,15) + 90720*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*(3*powl(c1,2) - powl(c2,2))*logl(c1 - c2) - 90720*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*(3*powl(c1,2) - powl(c2,2))*logl(c1 + c2)) + 3*c1*powl(c2,8)*(-1215*powl(c1,12)*c2 - 2835*powl(c1,10)*powl(c2,3) + 35775*powl(c1,8)*powl(c2,5) - 94221*powl(c1,6)*powl(c2,7) + 111792*powl(c1,4)*powl(c2,9) - 62608*powl(c1,2)*powl(c2,11) + 13312*powl(c2,13) + 630*powl(c1,5)*powl(-27*powl(c1,4) + 12*powl(c1,2)*powl(c2,2) + 16*powl(c2,4),2)*logl(c1 - c2) - 630*powl(c1,5)*powl(-27*powl(c1,4) + 12*powl(c1,2)*powl(c2,2) + 16*powl(c2,4),2)*logl(c1 + c2)) - 27*c1*powl(c2,2)*powl(powl(c1,2) - powl(c2,2),5)*powl(atanhl(c2/c1),6)*(-180495*powl(c1,8)*c2 + 58500*powl(c1,6)*powl(c2,3) + 906*powl(c1,4)*powl(c2,5) + 1772*powl(c1,2)*powl(c2,7) - 1083*powl(c2,9) + 2520*(27*powl(c1,9) - 33*powl(c1,7)*powl(c2,2) + 14*powl(c1,5)*powl(c2,4))*logl(c1 - c2) - 2520*(27*powl(c1,9) - 33*powl(c1,7)*powl(c2,2) + 14*powl(c1,5)*powl(c2,4))*logl(c1 + c2)) + c1*powl(c2,6)*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)*(-14773185*powl(c1,12)*c2 + 11484180*powl(c1,10)*powl(c2,3) + 7136262*powl(c1,8)*powl(c2,5) - 4352724*powl(c1,6)*powl(c2,7) + 1293147*powl(c1,4)*powl(c2,9) - 472704*powl(c1,2)*powl(c2,11) + 47888*powl(c2,13) + 7560*powl(c1,5)*(2025*powl(c1,8) - 3375*powl(c1,6)*powl(c2,2) + 1350*powl(c1,4)*powl(c2,4) + 216*powl(c1,2)*powl(c2,6) - 224*powl(c2,8))*logl(c1 - c2) - 7560*powl(c1,5)*(2025*powl(c1,8) - 3375*powl(c1,6)*powl(c2,2) + 1350*powl(c1,4)*powl(c2,4) + 216*powl(c1,2)*powl(c2,6) - 224*powl(c2,8))*logl(c1 + c2)) + powl(c2,7)*atanhl(c2/c1)*(2781135*powl(c1,14)*c2 - 2442150*powl(c1,12)*powl(c2,3) - 3187188*powl(c1,10)*powl(c2,5) + 2795094*powl(c1,8)*powl(c2,7) - 735831*powl(c1,6)*powl(c2,9) + 1074816*powl(c1,4)*powl(c2,11) - 307696*powl(c1,2)*powl(c2,13) + 25600*powl(c2,15) - 90720*powl(c1,7)*(81*powl(c1,8) - 144*powl(c1,6)*powl(c2,2) + 27*powl(c1,4)*powl(c2,4) + 52*powl(c1,2)*powl(c2,6) - 16*powl(c2,8))*logl(c1 - c2) + 90720*powl(c1,7)*(81*powl(c1,8) - 144*powl(c1,6)*powl(c2,2) + 27*powl(c1,4)*powl(c2,4) + 52*powl(c1,2)*powl(c2,6) - 16*powl(c2,8))*logl(c1 + c2)) - powl(c2,5)*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),3)*(-30745575*powl(c1,14)*c2 + 51359265*powl(c1,12)*powl(c2,3) - 20519406*powl(c1,10)*powl(c2,5) - 3545478*powl(c1,8)*powl(c2,7) + 3445749*powl(c1,6)*powl(c2,9) + 299109*powl(c1,4)*powl(c2,11) - 205136*powl(c1,2)*powl(c2,13) + 32432*powl(c2,15) + 181440*powl(c1,7)*(81*powl(c1,8) - 207*powl(c1,6)*powl(c2,2) + 201*powl(c1,4)*powl(c2,4) - 89*powl(c1,2)*powl(c2,6) + 14*powl(c2,8))*logl(c1 - c2) - 181440*powl(c1,7)*(81*powl(c1,8) - 207*powl(c1,6)*powl(c2,2) + 201*powl(c1,4)*powl(c2,4) - 89*powl(c1,2)*powl(c2,6) + 14*powl(c2,8))*logl(c1 + c2)) + 3*c1*powl(c2,4)*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4)*(-9840285*powl(c1,12)*c2 + 15380685*powl(c1,10)*powl(c2,3) - 9305010*powl(c1,8)*powl(c2,5) + 2008170*powl(c1,6)*powl(c2,7) - 304065*powl(c1,4)*powl(c2,9) + 151897*powl(c1,2)*powl(c2,11) - 26880*powl(c2,13) + 1260*powl(c1,5)*(1215*powl(c1,8) - 2970*powl(c1,6)*powl(c2,2) + 3375*powl(c1,4)*powl(c2,4) - 2028*powl(c1,2)*powl(c2,6) + 536*powl(c2,8))*logl(c1 - c2) - 1260*powl(c1,5)*(1215*powl(c1,8) - 2970*powl(c1,6)*powl(c2,2) + 3375*powl(c1,4)*powl(c2,4) - 2028*powl(c1,2)*powl(c2,6) + 536*powl(c2,8))*logl(c1 + c2)))/(powl(c2,2)*(-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))*(powl(c2,4)*(405*powl(c1,8) - 2970*powl(c1,6)*powl(c2,2) + 5157*powl(c1,4)*powl(c2,4) - 3612*powl(c1,2)*powl(c2,6) + 1024*powl(c2,8)) - 12*c1*powl(c2,3)*(135*powl(c1,8) - 810*powl(c1,6)*powl(c2,2) + 1431*powl(c1,4)*powl(c2,4) - 1018*powl(c1,2)*powl(c2,6) + 262*powl(c2,8))*atanhl(c2/c1) + 18*powl(c2,2)*(135*powl(c1,10) - 630*powl(c1,8)*powl(c2,2) + 1008*powl(c1,6)*powl(c2,4) - 640*powl(c1,4)*powl(c2,6) + 93*powl(c1,2)*powl(c2,8) + 34*powl(c2,10))*powl(atanhl(c2/c1),2) - 540*c1*c2*powl(powl(c1,2) - powl(c2,2),4)*(3*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),3) + 81*powl(powl(c1,2) - powl(c2,2),4)*(5*powl(c1,4) + 10*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),4))));

    bb[0] = 0;

    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);

    bb[2] = -(((c1 - c2)*(c1 + c2)*(3*powl(c1,2) + powl(c2,2))*(3*powl(c1,2)*powl(c2,2) - 4*powl(c2,4) + (-6*powl(c1,3)*c2 + 6*c1*powl(c2,3))*atanhl(c2/c1) + (3*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) - powl(c2,4))*powl(atanhl(c2/c1),2)))/powl(c2,8));

    bb[3] = (-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))/(3.*powl(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2),2));

    bb[4] = ((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(powl(c2,4)*(405*powl(c1,8) - 2970*powl(c1,6)*powl(c2,2) + 5157*powl(c1,4)*powl(c2,4) - 3612*powl(c1,2)*powl(c2,6) + 1024*powl(c2,8)) - 12*c1*powl(c2,3)*(135*powl(c1,8) - 810*powl(c1,6)*powl(c2,2) + 1431*powl(c1,4)*powl(c2,4) - 1018*powl(c1,2)*powl(c2,6) + 262*powl(c2,8))*atanhl(c2/c1) + 18*powl(c2,2)*(135*powl(c1,10) - 630*powl(c1,8)*powl(c2,2) + 1008*powl(c1,6)*powl(c2,4) - 640*powl(c1,4)*powl(c2,6) + 93*powl(c1,2)*powl(c2,8) + 34*powl(c2,10))*powl(atanhl(c2/c1),2) - 540*c1*c2*powl(powl(c1,2) - powl(c2,2),4)*(3*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),3) + 81*powl(powl(c1,2) - powl(c2,2),4)*(5*powl(c1,4) + 10*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),4)))/(15.*powl(-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4),2));

  }
  else {
    mpi_abort("[D4EST_ERROR]: Do not support n >= 5 yet\n");
  }

}

long double
d4est_quadrature_compactified_c1tpc2_neg4_weight_fcn(long double x, void* user){
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  return 1.0l/powl(c1 + c2*x,4);
}

long double
d4est_quadrature_compactified_c1tpc2_neg3_weight_fcn(long double x, void* user){
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  return 1.0l/powl(c1 + c2*x,3);
}

long double
d4est_quadrature_compactified_c1tpc2_neg2_weight_fcn
(
 long double x,
 void* user
){
  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  return 1.0l/powl(c1 + c2*x,2);
}

long double
d4est_quadrature_compactified_c1tpc2_neg1_weight_fcn
(
 long double x,
 void* user
){
  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  return 1.0l/(c1 + c2*x);
}


long double
d4est_quadrature_compactified_c1tpc2_neg3_moment_fcn(int n, void* user)
{
  
  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;

  if (n == 0) return (2*c1)/(powl(c1 - c2,2)*powl(c1 + c2,2));

  else if (n == 1) return (-2*c2)/(powl(c1 - c2,2)*powl(c1 + c2,2));

  else if (n == 2) return (2*((-(powl(c1,3)*c2) + 2*c1*powl(c2,3))/powl(powl(c1,2) - powl(c2,2),2) + atanhl(c2/c1)))/powl(c2,3);

  else if (n == 3) return (2*(3*powl(c1,4)*c2 - 5*powl(c1,2)*powl(c2,3) + powl(c2,5) - 3*c1*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(powl(c2,4)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 4) return (2*c1*(-6*powl(c1,4)*c2 + 10*powl(c1,2)*powl(c2,3) - 3*powl(c2,5) + 6*c1*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(powl(c2,5)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 5) return (2*(30*powl(c1,6)*c2 - 50*powl(c1,4)*powl(c2,3) + 16*powl(c1,2)*powl(c2,5) + powl(c2,7) - 30*powl(c1,3)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(3.*powl(c2,6)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 6) return (2*c1*(-((c2*(15*powl(c1,6) - 25*powl(c1,4)*powl(c2,2) + 8*powl(c1,2)*powl(c2,4) + powl(c2,6)))/powl(powl(c1,2) - powl(c2,2),2)) + 15*powl(c1,3)*atanhl(c2/c1)))/powl(c2,7);

  else if (n == 7) return (2*(105*powl(c1,8)*c2 - 175*powl(c1,6)*powl(c2,3) + 56*powl(c1,4)*powl(c2,5) + 8*powl(c1,2)*powl(c2,7) + powl(c2,9) - 105*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(5.*powl(c2,8)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 8) return (-2*(420*powl(c1,9) - 700*powl(c1,7)*powl(c2,2) + 224*powl(c1,5)*powl(c2,4) + 32*powl(c1,3)*powl(c2,6) + 9*c1*powl(c2,8)))/(15.*powl(c2,8)*powl(powl(c1,2) - powl(c2,2),2)) + (56*powl(c1,6)*atanhl(c2/c1))/powl(c2,9);

  else if (n == 9) return (2*(1260*powl(c1,10)*c2 - 2100*powl(c1,8)*powl(c2,3) + 672*powl(c1,6)*powl(c2,5) + 96*powl(c1,4)*powl(c2,7) + 32*powl(c1,2)*powl(c2,9) + 5*powl(c2,11) + 630*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*logl(c1 - c2) - 630*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*logl(c1 + c2)))/(35.*powl(c1 - c2,2)*powl(c2,10)*powl(c1 + c2,2));

  else if (n == 10) return (2*c1*(-((c2*(315*powl(c1,10) - 525*powl(c1,8)*powl(c2,2) + 168*powl(c1,6)*powl(c2,4) + 24*powl(c1,4)*powl(c2,6) + 8*powl(c1,2)*powl(c2,8) + 3*powl(c2,10)))/powl(powl(c1,2) - powl(c2,2),2)) + 315*powl(c1,7)*atanhl(c2/c1)))/(7.*powl(c2,11));

  else if (n == 11) return (2*(3465*powl(c1,12)*c2 - 5775*powl(c1,10)*powl(c2,3) + 1848*powl(c1,8)*powl(c2,5) + 264*powl(c1,6)*powl(c2,7) + 88*powl(c1,4)*powl(c2,9) + 40*powl(c1,2)*powl(c2,11) + 7*powl(c2,13) - 3465*powl(c1,9)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(63.*powl(c2,12)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 12) return (2*c1*(-((c2*(6930*powl(c1,12) - 11550*powl(c1,10)*powl(c2,2) + 3696*powl(c1,8)*powl(c2,4) + 528*powl(c1,6)*powl(c2,6) + 176*powl(c1,4)*powl(c2,8) + 80*powl(c1,2)*powl(c2,10) + 35*powl(c2,12)))/powl(powl(c1,2) - powl(c2,2),2)) + 6930*powl(c1,9)*atanhl(c2/c1)))/(105.*powl(c2,13));

  else if (n == 13) return (2*(90090*powl(c1,14)*c2 - 150150*powl(c1,12)*powl(c2,3) + 48048*powl(c1,10)*powl(c2,5) + 6864*powl(c1,8)*powl(c2,7) + 2288*powl(c1,6)*powl(c2,9) + 1040*powl(c1,4)*powl(c2,11) + 560*powl(c1,2)*powl(c2,13) + 105*powl(c2,15) - 90090*powl(c1,11)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(1155.*powl(c2,14)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 14) return (-2*(45045*powl(c1,15) - 75075*powl(c1,13)*powl(c2,2) + 24024*powl(c1,11)*powl(c2,4) + 3432*powl(c1,9)*powl(c2,6) + 1144*powl(c1,7)*powl(c2,8) + 520*powl(c1,5)*powl(c2,10) + 280*powl(c1,3)*powl(c2,12) + 135*c1*powl(c2,14)))/(495.*powl(c2,14)*powl(powl(c1,2) - powl(c2,2),2)) + (182*powl(c1,12)*atanhl(c2/c1))/powl(c2,15);

  else if (n == 15) return (2*((45045*powl(c1,16)*c2 - 75075*powl(c1,14)*powl(c2,3) + 24024*powl(c1,12)*powl(c2,5) + 3432*powl(c1,10)*powl(c2,7) + 1144*powl(c1,8)*powl(c2,9) + 520*powl(c1,6)*powl(c2,11) + 280*powl(c1,4)*powl(c2,13) + 168*powl(c1,2)*powl(c2,15) + 33*powl(c2,17))/powl(powl(c1,2) - powl(c2,2),2) - 45045*powl(c1,13)*atanhl(c2/c1)))/(429.*powl(c2,16));

  else if (n == 16) return (-2*(360360*powl(c1,17) - 600600*powl(c1,15)*powl(c2,2) + 192192*powl(c1,13)*powl(c2,4) + 27456*powl(c1,11)*powl(c2,6) + 9152*powl(c1,9)*powl(c2,8) + 4160*powl(c1,7)*powl(c2,10) + 2240*powl(c1,5)*powl(c2,12) + 1344*powl(c1,3)*powl(c2,14) + 693*c1*powl(c2,16)))/(3003.*powl(c2,16)*powl(powl(c1,2) - powl(c2,2),2)) + (240*powl(c1,14)*atanhl(c2/c1))/powl(c2,17);

  else if (n == 17) return (2*(6126120*powl(c1,18)*c2 - 10210200*powl(c1,16)*powl(c2,3) + 3267264*powl(c1,14)*powl(c2,5) + 466752*powl(c1,12)*powl(c2,7) + 155584*powl(c1,10)*powl(c2,9) + 70720*powl(c1,8)*powl(c2,11) + 38080*powl(c1,6)*powl(c2,13) + 22848*powl(c1,4)*powl(c2,15) + 14784*powl(c1,2)*powl(c2,17) + 3003*powl(c2,19) - 6126120*powl(c1,15)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(45045.*powl(c2,18)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 18) return (-2*(765765*powl(c1,19) - 1276275*powl(c1,17)*powl(c2,2) + 408408*powl(c1,15)*powl(c2,4) + 58344*powl(c1,13)*powl(c2,6) + 19448*powl(c1,11)*powl(c2,8) + 8840*powl(c1,9)*powl(c2,10) + 4760*powl(c1,7)*powl(c2,12) + 2856*powl(c1,5)*powl(c2,14) + 1848*powl(c1,3)*powl(c2,16) + 1001*c1*powl(c2,18)))/(5005.*powl(c2,18)*powl(powl(c1,2) - powl(c2,2),2)) + (306*powl(c1,16)*atanhl(c2/c1))/powl(c2,19);

  else if (n == 19) return (2*(14549535*powl(c1,20)*c2 - 24249225*powl(c1,18)*powl(c2,3) + 7759752*powl(c1,16)*powl(c2,5) + 1108536*powl(c1,14)*powl(c2,7) + 369512*powl(c1,12)*powl(c2,9) + 167960*powl(c1,10)*powl(c2,11) + 90440*powl(c1,8)*powl(c2,13) + 54264*powl(c1,6)*powl(c2,15) + 35112*powl(c1,4)*powl(c2,17) + 24024*powl(c1,2)*powl(c2,19) + 5005*powl(c2,21) - 14549535*powl(c1,17)*powl(powl(c1,2) - powl(c2,2),2)*atanhl(c2/c1)))/(85085.*powl(c2,20)*powl(powl(c1,2) - powl(c2,2),2));

  else if (n == 20) return (-2*(29099070*powl(c1,21) - 48498450*powl(c1,19)*powl(c2,2) + 15519504*powl(c1,17)*powl(c2,4) + 2217072*powl(c1,15)*powl(c2,6) + 739024*powl(c1,13)*powl(c2,8) + 335920*powl(c1,11)*powl(c2,10) + 180880*powl(c1,9)*powl(c2,12) + 108528*powl(c1,7)*powl(c2,14) + 70224*powl(c1,5)*powl(c2,16) + 48048*powl(c1,3)*powl(c2,18) + 27027*c1*powl(c2,20)))/(153153.*powl(c2,20)*powl(powl(c1,2) - powl(c2,2),2)) + (380*powl(c1,18)*atanhl(c2/c1))/powl(c2,21);

  else {
    mpi_abort("n must be <= 20");
    return NAN;
  }
    
  
}

long double
d4est_quadrature_compactified_c1tpc2_neg2_moment_fcn(int n, void* user)
{
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  

  if (n == 0) return 2/(powl(c1,2) - powl(c2,2));

  else if (n == 1) return (2*(-((c1*c2)/(powl(c1,2) - powl(c2,2))) + atanhl(c2/c1)))/powl(c2,2);

  else if (n == 2) return (2*((2*powl(c1,2)*c2 - powl(c2,3))/(powl(c1,2) - powl(c2,2)) - 2*c1*atanhl(c2/c1)))/powl(c2,3);

  else if (n == 3) return (2*c1*((-3*powl(c1,2)*c2 + 2*powl(c2,3))/(powl(c1,2) - powl(c2,2)) + 3*c1*atanhl(c2/c1)))/powl(c2,4);

  else if (n == 4) return (2*(-((-12*powl(c1,4)*c2 + 8*powl(c1,2)*powl(c2,3) + powl(c2,5))/(powl(c1,2) - powl(c2,2))) - 12*powl(c1,3)*atanhl(c2/c1)))/(3.*powl(c2,5));

  else if (n == 5) return (2*c1*((-15*powl(c1,4)*c2 + 10*powl(c1,2)*powl(c2,3) + 2*powl(c2,5))/(powl(c1,2) - powl(c2,2)) + 15*powl(c1,3)*atanhl(c2/c1)))/(3.*powl(c2,6));

  else if (n == 6) return (2*(-((-30*powl(c1,6)*c2 + 20*powl(c1,4)*powl(c2,3) + 4*powl(c1,2)*powl(c2,5) + powl(c2,7))/(powl(c1,2) - powl(c2,2))) - 30*powl(c1,5)*atanhl(c2/c1)))/(5.*powl(c2,7));

  else if (n == 7) return (2*(-105*powl(c1,7)*c2 + 70*powl(c1,5)*powl(c2,3) + 14*powl(c1,3)*powl(c2,5) + 6*c1*powl(c2,7) + 105*(powl(c1,8) - powl(c1,6)*powl(c2,2))*atanhl(c2/c1)))/(15.*(c1 - c2)*powl(c2,8)*(c1 + c2));

  else if (n == 8) return (2*(840*powl(c1,8)*c2 - 560*powl(c1,6)*powl(c2,3) - 112*powl(c1,4)*powl(c2,5) - 48*powl(c1,2)*powl(c2,7) - 15*powl(c2,9) + 840*powl(c1,7)*(-powl(c1,2) + powl(c2,2))*atanhl(c2/c1)))/(105.*(c1 - c2)*powl(c2,9)*(c1 + c2));

  else if (n == 9) return (-630*powl(c1,9)*c2 + 420*powl(c1,7)*powl(c2,3) + 84*powl(c1,5)*powl(c2,5) + 36*powl(c1,3)*powl(c2,7) + 20*c1*powl(c2,9) + 630*powl(c1,8)*(c1 - c2)*(c1 + c2)*atanhl(c2/c1))/(35.*(c1 - c2)*powl(c2,10)*(c1 + c2));

  else if (n == 10) return (2*(630*powl(c1,10)*c2 - 420*powl(c1,8)*powl(c2,3) - 84*powl(c1,6)*powl(c2,5) - 36*powl(c1,4)*powl(c2,7) - 20*powl(c1,2)*powl(c2,9) - 7*powl(c2,11) + 315*powl(c1,9)*(powl(c1,2) - powl(c2,2))*logl(c1 - c2) - 315*powl(c1,9)*(powl(c1,2) - powl(c2,2))*logl(c1 + c2)))/(63.*(c1 - c2)*powl(c2,11)*(c1 + c2));

  else if (n == 11) return (2*c1*(-3465*powl(c1,10)*c2 + 2310*powl(c1,8)*powl(c2,3) + 462*powl(c1,6)*powl(c2,5) + 198*powl(c1,4)*powl(c2,7) + 110*powl(c1,2)*powl(c2,9) + 70*powl(c2,11) + 3465*powl(c1,9)*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1)))/(315.*(c1 - c2)*powl(c2,12)*(c1 + c2));

  else if (n == 12) return (24*powl(c1,10))/powl(c2,12) + (8*powl(c1,8))/powl(c2,10) + (24*powl(c1,6))/(5.*powl(c2,8)) + (24*powl(c1,4))/(7.*powl(c2,6)) + (8*powl(c1,2))/(3.*powl(c2,4)) + 24/(11.*powl(c2,2)) + 2/(powl(c1,2) - powl(c2,2)) - (24*powl(c1,11)*atanhl(c2/c1))/powl(c2,13);

  else if (n == 13) return (2*c1*(-45045*powl(c1,12) + 30030*powl(c1,10)*powl(c2,2) + 6006*powl(c1,8)*powl(c2,4) + 2574*powl(c1,6)*powl(c2,6) + 1430*powl(c1,4)*powl(c2,8) + 910*powl(c1,2)*powl(c2,10) + 630*powl(c2,12)))/(3465.*(c1 - c2)*powl(c2,13)*(c1 + c2)) + (26*powl(c1,12)*atanhl(c2/c1))/powl(c2,14);

  else if (n == 14) return (28*powl(c1,12))/powl(c2,14) + (28*powl(c1,10))/(3.*powl(c2,12)) + (28*powl(c1,8))/(5.*powl(c2,10)) + (4*powl(c1,6))/powl(c2,8) + (28*powl(c1,4))/(9.*powl(c2,6)) + (28*powl(c1,2))/(11.*powl(c2,4)) + 28/(13.*powl(c2,2)) + 2/(powl(c1,2) - powl(c2,2)) - (28*powl(c1,13)*atanhl(c2/c1))/powl(c2,15);

  else if (n == 15) return (2*c1*(-45045*powl(c1,14) + 30030*powl(c1,12)*powl(c2,2) + 6006*powl(c1,10)*powl(c2,4) + 2574*powl(c1,8)*powl(c2,6) + 1430*powl(c1,6)*powl(c2,8) + 910*powl(c1,4)*powl(c2,10) + 630*powl(c1,2)*powl(c2,12) + 462*powl(c2,14)))/(3003.*(c1 - c2)*powl(c2,15)*(c1 + c2)) + (30*powl(c1,14)*atanhl(c2/c1))/powl(c2,16);

  else if (n == 16) return (32*powl(c1,14))/powl(c2,16) + (32*powl(c1,12))/(3.*powl(c2,14)) + (32*powl(c1,10))/(5.*powl(c2,12)) + (32*powl(c1,8))/(7.*powl(c2,10)) + (32*powl(c1,6))/(9.*powl(c2,8)) + (32*powl(c1,4))/(11.*powl(c2,6)) + (32*powl(c1,2))/(13.*powl(c2,4)) + 32/(15.*powl(c2,2)) + 2/(powl(c1,2) - powl(c2,2)) - (32*powl(c1,15)*atanhl(c2/c1))/powl(c2,17);

  else if (n == 17) return (2*c1*(-765765*powl(c1,16) + 510510*powl(c1,14)*powl(c2,2) + 102102*powl(c1,12)*powl(c2,4) + 43758*powl(c1,10)*powl(c2,6) + 24310*powl(c1,8)*powl(c2,8) + 15470*powl(c1,6)*powl(c2,10) + 10710*powl(c1,4)*powl(c2,12) + 7854*powl(c1,2)*powl(c2,14) + 6006*powl(c2,16)))/(45045.*(c1 - c2)*powl(c2,17)*(c1 + c2)) + (17*powl(c1,16)*logl((c1 + c2)/(c1 - c2)))/powl(c2,18);

  else if (n == 18) return (36*powl(c1,16))/powl(c2,18) + (12*powl(c1,14))/powl(c2,16) + (36*powl(c1,12))/(5.*powl(c2,14)) + (36*powl(c1,10))/(7.*powl(c2,12)) + (4*powl(c1,8))/powl(c2,10) + (36*powl(c1,6))/(11.*powl(c2,8)) + (36*powl(c1,4))/(13.*powl(c2,6)) + (12*powl(c1,2))/(5.*powl(c2,4)) + 36/(17.*powl(c2,2)) + 2/(powl(c1,2) - powl(c2,2)) - (36*powl(c1,17)*atanhl(c2/c1))/powl(c2,19);

  else if (n == 19) return (-38*powl(c1,17))/powl(c2,19) - (38*powl(c1,15))/(3.*powl(c2,17)) - (38*powl(c1,13))/(5.*powl(c2,15)) - (38*powl(c1,11))/(7.*powl(c2,13)) - (38*powl(c1,9))/(9.*powl(c2,11)) - (38*powl(c1,7))/(11.*powl(c2,9)) - (38*powl(c1,5))/(13.*powl(c2,7)) - (38*powl(c1,3))/(15.*powl(c2,5)) - (38*c1)/(17.*powl(c2,3)) + 1/(powl(c1,2) + c1*c2) + (-2/c2 + 1/(-c1 + c2))/c1 + (19*powl(c1,18)*logl((c1 + c2)/(c1 - c2)))/powl(c2,20);

  else if (n == 20) return (40*powl(c1,18))/powl(c2,20) + (40*powl(c1,16))/(3.*powl(c2,18)) + (8*powl(c1,14))/powl(c2,16) + (40*powl(c1,12))/(7.*powl(c2,14)) + (40*powl(c1,10))/(9.*powl(c2,12)) + (40*powl(c1,8))/(11.*powl(c2,10)) + (40*powl(c1,6))/(13.*powl(c2,8)) + (8*powl(c1,4))/(3.*powl(c2,6)) + (40*powl(c1,2))/(17.*powl(c2,4)) + 40/(19.*powl(c2,2)) + 2/(powl(c1,2) - powl(c2,2)) - (40*powl(c1,19)*atanhl(c2/c1))/powl(c2,21);
  else {
    mpi_abort("n must be <= 20");
    return NAN;
  }
}

long double
d4est_quadrature_compactified_c1tpc2_neg4_moment_fcn(int n, void* user)
{
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  
  if (!(c1 + c2 > 0.0l && c1 - c2 > 0.0l && n <= 20)){
    printf(" c1 = %Le\n", c1);
    printf(" c2 = %Le\n", c2);
    /* printf(" R2 = %Le\n", R2); */
    /* printf(" R1 = %Le\n", R1); */
    /* printf(" cmax = %Le\n", cmax); */
    /* printf(" cmin == %Le\n", cmin); */
    /* printf(" n == %d\n", n); */
    mpi_abort("[D4EST_ERROR]: condition c1 + c2 > 0 && c1 - c2 > 0 && n <= 20 failed\n");
  }
  
  if (n == 0) return (2*(3*powl(c1,2) + powl(c2,2)))/(3.*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 1) return (-8*c1*c2)/(3.*powl(c1 - c2,3)*powl(c1 + c2,3));
  else if (n == 2) return (2*(powl(c1,2) + 3*powl(c2,2)))/(3.*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 3) return (2*((-3*powl(c1,5)*c2 + 8*powl(c1,3)*powl(c2,3) - 9*c1*powl(c2,5))/powl(powl(c1,2) - powl(c2,2),3) + 3*atanhl(c2/c1)))/(3.*powl(c2,4));
  else if (n == 4) return (-24*powl(c1,6)*c2 + 64*powl(c1,4)*powl(c2,3) - 54*powl(c1,2)*powl(c2,5) + 6*powl(c2,7) + 24*c1*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1))/(3.*powl(c2,5)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 5) return (2*c1*(-30*powl(c1,6)*c2 + 80*powl(c1,4)*powl(c2,3) - 66*powl(c1,2)*powl(c2,5) + 12*powl(c2,7) + 30*c1*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,6)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 6) return (2*(-60*powl(c1,8)*c2 + 160*powl(c1,6)*powl(c2,3) - 132*powl(c1,4)*powl(c2,5) + 27*powl(c1,2)*powl(c2,7) + powl(c2,9) + 60*powl(powl(c1,3) - c1*powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,7)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 7) return (c1*(-210*powl(c1,8)*c2 + 560*powl(c1,6)*powl(c2,3) - 462*powl(c1,4)*powl(c2,5) + 96*powl(c1,2)*powl(c2,7) + 8*powl(c2,9) + 210*powl(powl(c1,3) - c1*powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,8)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 8) return (840*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 2*(-840*powl(c1,10)*c2 + 2240*powl(c1,8)*powl(c2,3) - 1848*powl(c1,6)*powl(c2,5) + 384*powl(c1,4)*powl(c2,7) + 41*powl(c1,2)*powl(c2,9) + 3*powl(c2,11) + 420*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15.*powl(c1 - c2,3)*powl(c2,9)*powl(c1 + c2,3));
  else if (n == 9) return (-4*c1*(-2*c2*(-315*powl(c1,10) + 840*powl(c1,8)*powl(c2,2) - 693*powl(c1,6)*powl(c2,4) + 144*powl(c1,4)*powl(c2,6) + 16*powl(c1,2)*powl(c2,8) + 3*powl(c2,10)) + 315*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 315*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15.*powl(c1 - c2,3)*powl(c2,10)*powl(c1 + c2,3));
  else if (n == 10) return (2*(-2520*powl(c1,12)*c2 + 6720*powl(c1,10)*powl(c2,3) - 5544*powl(c1,8)*powl(c2,5) + 1152*powl(c1,6)*powl(c2,7) + 128*powl(c1,4)*powl(c2,9) + 33*powl(c1,2)*powl(c2,11) + 3*powl(c2,13) + 2520*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(21.*powl(c2,11)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 11) return (-2*c1*(-3465*powl(c1,12)*c2 + 9240*powl(c1,10)*powl(c2,3) - 7623*powl(c1,8)*powl(c2,5) + 1584*powl(c1,6)*powl(c2,7) + 176*powl(c1,4)*powl(c2,9) + 48*powl(c1,2)*powl(c2,11) + 12*powl(c2,13) + 3465*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(21.*powl(c2,12)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 12) return (2*(-((-13860*powl(c1,14)*c2 + 36960*powl(c1,12)*powl(c2,3) - 30492*powl(c1,10)*powl(c2,5) + 6336*powl(c1,8)*powl(c2,7) + 704*powl(c1,6)*powl(c2,9) + 192*powl(c1,4)*powl(c2,11) + 69*powl(c1,2)*powl(c2,13) + 7*powl(c2,15))/powl(powl(c1,2) - powl(c2,2),3)) - 13860*powl(c1,9)*atanhl(c2/c1)))/(63.*powl(c2,13));
  else if (n == 13) return (4*c1*((-45045*powl(c1,14)*c2 + 120120*powl(c1,12)*powl(c2,3) - 99099*powl(c1,10)*powl(c2,5) + 20592*powl(c1,8)*powl(c2,7) + 2288*powl(c1,6)*powl(c2,9) + 624*powl(c1,4)*powl(c2,11) + 240*powl(c1,2)*powl(c2,13) + 70*powl(c2,15))/powl(powl(c1,2) - powl(c2,2),3) + 45045*powl(c1,9)*atanhl(c2/c1)))/(315.*powl(c2,14));
  else if (n == 14) return (2*(-180180*powl(c1,16) + 480480*powl(c1,14)*powl(c2,2) - 396396*powl(c1,12)*powl(c2,4) + 82368*powl(c1,10)*powl(c2,6) + 9152*powl(c1,8)*powl(c2,8) + 2496*powl(c1,6)*powl(c2,10) + 960*powl(c1,4)*powl(c2,12) + 415*powl(c1,2)*powl(c2,14) + 45*powl(c2,16)))/(495.*powl(c2,14)*powl(-powl(c1,2) + powl(c2,2),3)) - (728*powl(c1,11)*atanhl(c2/c1))/powl(c2,15);
  else if (n == 15) return (-2*c1*(-45045*powl(c1,16)*c2 + 120120*powl(c1,14)*powl(c2,3) - 99099*powl(c1,12)*powl(c2,5) + 20592*powl(c1,10)*powl(c2,7) + 2288*powl(c1,8)*powl(c2,9) + 624*powl(c1,6)*powl(c2,11) + 240*powl(c1,4)*powl(c2,13) + 112*powl(c1,2)*powl(c2,15) + 36*powl(c2,17) + 45045*powl(c1,11)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(99.*powl(c2,16)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 16) return (2*(720720*powl(c1,18)*c2 - 1921920*powl(c1,16)*powl(c2,3) + 1585584*powl(c1,14)*powl(c2,5) - 329472*powl(c1,12)*powl(c2,7) - 36608*powl(c1,10)*powl(c2,9) - 9984*powl(c1,8)*powl(c2,11) - 3840*powl(c1,6)*powl(c2,13) - 1792*powl(c1,4)*powl(c2,15) - 873*powl(c1,2)*powl(c2,17) - 99*powl(c2,19) + 360360*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 360360*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(1287.*powl(c1 - c2,3)*powl(c2,17)*powl(c1 + c2,3));
  else if (n == 17) return (8*c1*(-1531530*powl(c1,18)*c2 + 4084080*powl(c1,16)*powl(c2,3) - 3369366*powl(c1,14)*powl(c2,5) + 700128*powl(c1,12)*powl(c2,7) + 77792*powl(c1,10)*powl(c2,9) + 21216*powl(c1,8)*powl(c2,11) + 8160*powl(c1,6)*powl(c2,13) + 3808*powl(c1,4)*powl(c2,15) + 2016*powl(c1,2)*powl(c2,17) + 693*powl(c2,19) + 1531530*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(9009.*powl(c2,18)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 18) return (2*(12252240*powl(c1,20)*c2 - 32672640*powl(c1,18)*powl(c2,3) + 26954928*powl(c1,16)*powl(c2,5) - 5601024*powl(c1,14)*powl(c2,7) - 622336*powl(c1,12)*powl(c2,9) - 169728*powl(c1,10)*powl(c2,11) - 65280*powl(c1,8)*powl(c2,13) - 30464*powl(c1,6)*powl(c2,15) - 16128*powl(c1,4)*powl(c2,17) - 8547*powl(c1,2)*powl(c2,19) - 1001*powl(c2,21) + 6126120*powl(c1,15)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 6126120*powl(c1,15)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15015.*powl(c1 - c2,3)*powl(c2,19)*powl(c1 + c2,3));
  else if (n == 19) return (2*c1*c2*(-14549535*powl(c1,20) + 38798760*powl(c1,18)*powl(c2,2) - 32008977*powl(c1,16)*powl(c2,4) + 6651216*powl(c1,14)*powl(c2,6) + 739024*powl(c1,12)*powl(c2,8) + 201552*powl(c1,10)*powl(c2,10) + 77520*powl(c1,8)*powl(c2,12) + 36176*powl(c1,6)*powl(c2,14) + 19152*powl(c1,4)*powl(c2,16) + 11088*powl(c1,2)*powl(c2,18) + 4004*powl(c2,20)) - 14549535*powl(c1,16)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) + 14549535*powl(c1,16)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2))/(15015.*powl(c1 - c2,3)*powl(c2,20)*powl(c1 + c2,3));
  else if (n == 20) return (2*(-58198140*powl(c1,22) + 155195040*powl(c1,20)*powl(c2,2) - 128035908*powl(c1,18)*powl(c2,4) + 26604864*powl(c1,16)*powl(c2,6) + 2956096*powl(c1,14)*powl(c2,8) + 806208*powl(c1,12)*powl(c2,10) + 310080*powl(c1,10)*powl(c2,12) + 144704*powl(c1,8)*powl(c2,14) + 76608*powl(c1,6)*powl(c2,16) + 44352*powl(c1,4)*powl(c2,18) + 25025*powl(c1,2)*powl(c2,20) + 3003*powl(c2,22)))/(51051.*powl(c2,20)*powl(-powl(c1,2) + powl(c2,2),3)) - (2280*powl(c1,17)*atanhl(c2/c1))/powl(c2,21);
  else {
    mpi_abort("n must be <= 20");
    return NAN;
  }
}

long double
d4est_quadrature_compactified_c1tpc2_neg1_moment_fcn(int n, void* user)
{
  /* long double cmin == ((d4est_quadrature_compactified_params_t*)user)->cmin; */
  /* long double cmax = ((d4est_quadrature_compactified_params_t*)user)->cmax; */
  /* long double R1 = ((d4est_quadrature_compactified_params_t*)user)->R1; */
  /* long double R2 = ((d4est_quadrature_compactified_params_t*)user)->R2; */
  /* long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1); */
  /* long double c2 = -((R2-R1)*(cmax-cmin)); */

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  
if (n == 0) return logl((c1 + c2)/(c1 - c2))/c2;

else if (n == 1) return (2*(c2 - c1*atanhl(c2/c1)))/powl(c2,2);

else if (n == 2) return (2*c1*(-c2 + c1*atanhl(c2/c1)))/powl(c2,3);

else if (n == 3) return (2*(3*powl(c1,2)*c2 + powl(c2,3) - 3*powl(c1,3)*atanhl(c2/c1)))/(3.*powl(c2,4));

else if (n == 4) return (-2*(3*powl(c1,3)*c2 + c1*powl(c2,3) - 3*powl(c1,4)*atanhl(c2/c1)))/(3.*powl(c2,5));

else if (n == 5) return (2*(15*powl(c1,4)*c2 + 5*powl(c1,2)*powl(c2,3) + 3*powl(c2,5) - 15*powl(c1,5)*atanhl(c2/c1)))/(15.*powl(c2,6));

else if (n == 6) return (-2*(15*powl(c1,5)*c2 + 5*powl(c1,3)*powl(c2,3) + 3*c1*powl(c2,5) - 15*powl(c1,6)*atanhl(c2/c1)))/(15.*powl(c2,7));

else if (n == 7) return (2*powl(c1,6))/powl(c2,7) + (2*powl(c1,4))/(3.*powl(c2,5)) + (2*powl(c1,2))/(5.*powl(c2,3)) + 2/(7.*c2) - (2*powl(c1,7)*atanhl(c2/c1))/powl(c2,8);

else if (n == 8) return -(210*powl(c1,7)*c2 + 70*powl(c1,5)*powl(c2,3) + 42*powl(c1,3)*powl(c2,5) + 30*c1*powl(c2,7) - 210*powl(c1,8)*atanhl(c2/c1))/(105.*powl(c2,9));

else if (n == 9) return (2*powl(c1,8))/powl(c2,9) + (2*powl(c1,6))/(3.*powl(c2,7)) + (2*powl(c1,4))/(5.*powl(c2,5)) + (2*powl(c1,2))/(7.*powl(c2,3)) + 2/(9.*c2) - (2*powl(c1,9)*atanhl(c2/c1))/powl(c2,10);

else if (n == 10) return -(c1*(630*powl(c1,8)*c2 + 210*powl(c1,6)*powl(c2,3) + 126*powl(c1,4)*powl(c2,5) + 90*powl(c1,2)*powl(c2,7) + 70*powl(c2,9) - 630*powl(c1,9)*atanhl(c2/c1)))/(315.*powl(c2,11));

else if (n == 11) return (2*powl(c1,10))/powl(c2,11) + (2*powl(c1,8))/(3.*powl(c2,9)) + (2*powl(c1,6))/(5.*powl(c2,7)) + (2*powl(c1,4))/(7.*powl(c2,5)) + (2*powl(c1,2))/(9.*powl(c2,3)) + 2/(11.*c2) - (2*powl(c1,11)*atanhl(c2/c1))/powl(c2,12);

else if (n == 12) return (-2*(3465*powl(c1,11)*c2 + 1155*powl(c1,9)*powl(c2,3) + 693*powl(c1,7)*powl(c2,5) + 495*powl(c1,5)*powl(c2,7) + 385*powl(c1,3)*powl(c2,9) + 315*c1*powl(c2,11) - 3465*powl(c1,12)*atanhl(c2/c1)))/(3465.*powl(c2,13));

else if (n == 13) return (2*powl(c1,12))/powl(c2,13) + (2*powl(c1,10))/(3.*powl(c2,11)) + (2*powl(c1,8))/(5.*powl(c2,9)) + (2*powl(c1,6))/(7.*powl(c2,7)) + (2*powl(c1,4))/(9.*powl(c2,5)) + (2*powl(c1,2))/(11.*powl(c2,3)) + 2/(13.*c2) - (2*powl(c1,13)*atanhl(c2/c1))/powl(c2,14);

else if (n == 14) return (-2*powl(c1,13))/powl(c2,14) - (2*powl(c1,11))/(3.*powl(c2,12)) - (2*powl(c1,9))/(5.*powl(c2,10)) - (2*powl(c1,7))/(7.*powl(c2,8)) - (2*powl(c1,5))/(9.*powl(c2,6)) - (2*powl(c1,3))/(11.*powl(c2,4)) - (2*c1)/(13.*powl(c2,2)) + (2*powl(c1,14)*atanhl(c2/c1))/powl(c2,15);

else if (n == 15) return (2*powl(c1,14))/powl(c2,15) + (2*powl(c1,12))/(3.*powl(c2,13)) + (2*powl(c1,10))/(5.*powl(c2,11)) + (2*powl(c1,8))/(7.*powl(c2,9)) + (2*powl(c1,6))/(9.*powl(c2,7)) + (2*powl(c1,4))/(11.*powl(c2,5)) + (2*powl(c1,2))/(13.*powl(c2,3)) + 2/(15.*c2) - (2*powl(c1,15)*atanhl(c2/c1))/powl(c2,16);

else if (n == 16) return (-2*powl(c1,15))/powl(c2,16) - (2*powl(c1,13))/(3.*powl(c2,14)) - (2*powl(c1,11))/(5.*powl(c2,12)) - (2*powl(c1,9))/(7.*powl(c2,10)) - (2*powl(c1,7))/(9.*powl(c2,8)) - (2*powl(c1,5))/(11.*powl(c2,6)) - (2*powl(c1,3))/(13.*powl(c2,4)) - (2*c1)/(15.*powl(c2,2)) + (2*powl(c1,16)*atanhl(c2/c1))/powl(c2,17);

else if (n == 17) return (2*powl(c1,16))/powl(c2,17) + (2*powl(c1,14))/(3.*powl(c2,15)) + (2*powl(c1,12))/(5.*powl(c2,13)) + (2*powl(c1,10))/(7.*powl(c2,11)) + (2*powl(c1,8))/(9.*powl(c2,9)) + (2*powl(c1,6))/(11.*powl(c2,7)) + (2*powl(c1,4))/(13.*powl(c2,5)) + (2*powl(c1,2))/(15.*powl(c2,3)) + 2/(17.*c2) - (2*powl(c1,17)*atanhl(c2/c1))/powl(c2,18);

else if (n == 18) return (-2*powl(c1,17))/powl(c2,18) - (2*powl(c1,15))/(3.*powl(c2,16)) - (2*powl(c1,13))/(5.*powl(c2,14)) - (2*powl(c1,11))/(7.*powl(c2,12)) - (2*powl(c1,9))/(9.*powl(c2,10)) - (2*powl(c1,7))/(11.*powl(c2,8)) - (2*powl(c1,5))/(13.*powl(c2,6)) - (2*powl(c1,3))/(15.*powl(c2,4)) - (2*c1)/(17.*powl(c2,2)) + (2*powl(c1,18)*atanhl(c2/c1))/powl(c2,19);

else if (n == 19) return (2*powl(c1,18))/powl(c2,19) + (2*powl(c1,16))/(3.*powl(c2,17)) + (2*powl(c1,14))/(5.*powl(c2,15)) + (2*powl(c1,12))/(7.*powl(c2,13)) + (2*powl(c1,10))/(9.*powl(c2,11)) + (2*powl(c1,8))/(11.*powl(c2,9)) + (2*powl(c1,6))/(13.*powl(c2,7)) + (2*powl(c1,4))/(15.*powl(c2,5)) + (2*powl(c1,2))/(17.*powl(c2,3)) + 2/(19.*c2) - (2*powl(c1,19)*atanhl(c2/c1))/powl(c2,20);

else if (n == 20) return (-2*powl(c1,19))/powl(c2,20) - (2*powl(c1,17))/(3.*powl(c2,18)) - (2*powl(c1,15))/(5.*powl(c2,16)) - (2*powl(c1,13))/(7.*powl(c2,14)) - (2*powl(c1,11))/(9.*powl(c2,12)) - (2*powl(c1,9))/(11.*powl(c2,10)) - (2*powl(c1,7))/(13.*powl(c2,8)) - (2*powl(c1,5))/(15.*powl(c2,6)) - (2*powl(c1,3))/(17.*powl(c2,4)) - (2*c1)/(19.*powl(c2,2)) + (2*powl(c1,20)*atanhl(c2/c1))/powl(c2,21);

 else {
   mpi_abort("[D4EST_ERROR]: Do not support n>=20 moments for neg1 compactified quadrature points");
   return NAN;
 }

 
}

void
d4est_quadrature_compactified_c1tpc2_neg1_aa_and_bb
(
 int n,
 long double* aa,
 long double* bb,
 void* user
){

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  

if (n == 1){

aa[0] = (2*(c2 - c1*atanhl(c2/c1)))/(c2*logl((c1 + c2)/(c1 - c2)));

bb[0] = 0;

}
 
else if (n == 2){

aa[0] = (2*(c2 - c1*atanhl(c2/c1)))/(c2*logl((c1 + c2)/(c1 - c2)));

aa[1] = (3*powl(c1,2)*c2 + powl(c2,3) - 3*powl(c1,3)*atanhl(c2/c1) + (12*powl(c2 - c1*atanhl(c2/c1),3))/powl(logl((c1 + c2)/(c1 - c2)),2) + (12*c1*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2)))/(3.*c2*(c1*(-c2 + c1*atanhl(c2/c1)) - (2*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2))));

bb[0] = 0;

bb[1] = (-2*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(powl(c2,2)*powl(logl((c1 + c2)/(c1 - c2)),2));

}

else if (n == 3){

aa[0] = (2*(c2 - c1*atanhl(c2/c1)))/(c2*logl((c1 + c2)/(c1 - c2)));

aa[1] = (3*powl(c1,2)*c2 + powl(c2,3) - 3*powl(c1,3)*atanhl(c2/c1) + (12*powl(c2 - c1*atanhl(c2/c1),3))/powl(logl((c1 + c2)/(c1 - c2)),2) + (12*c1*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2)))/(3.*c2*(c1*(-c2 + c1*atanhl(c2/c1)) - (2*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2))));

aa[2] = (c2*(108*powl(c1,4)*powl(atanhl(c2/c1),4) - 12*powl(c1,3)*powl(atanhl(c2/c1),3)*(31*c2 + 9*c1*logl((c1 + c2)/(c1 - c2))) + 9*powl(c1,2)*powl(atanhl(c2/c1),2)*(52*powl(c2,2) + 36*c1*c2*logl((c1 + c2)/(c1 - c2)) + 3*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) - 3*c1*c2*atanhl(c2/c1)*(84*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + 23*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(48*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + (42*powl(c1,2) + 5*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(15.*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))));

bb[0] = 0;

bb[1] = (-2*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(powl(c2,2)*powl(logl((c1 + c2)/(c1 - c2)),2));

bb[2] = -(c2*logl((c1 + c2)/(c1 - c2))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))))/(9.*powl(c2 - c1*atanhl(c2/c1),2)*powl(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)),2));

}

else if (n == 4){

aa[0] = (2*(c2 - c1*atanhl(c2/c1)))/(c2*logl((c1 + c2)/(c1 - c2)));

aa[1] = (3*powl(c1,2)*c2 + powl(c2,3) - 3*powl(c1,3)*atanhl(c2/c1) + (12*powl(c2 - c1*atanhl(c2/c1),3))/powl(logl((c1 + c2)/(c1 - c2)),2) + (12*c1*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2)))/(3.*c2*(c1*(-c2 + c1*atanhl(c2/c1)) - (2*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2))));

aa[2] = (c2*(108*powl(c1,4)*powl(atanhl(c2/c1),4) - 12*powl(c1,3)*powl(atanhl(c2/c1),3)*(31*c2 + 9*c1*logl((c1 + c2)/(c1 - c2))) + 9*powl(c1,2)*powl(atanhl(c2/c1),2)*(52*powl(c2,2) + 36*c1*c2*logl((c1 + c2)/(c1 - c2)) + 3*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) - 3*c1*c2*atanhl(c2/c1)*(84*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + 23*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(48*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + (42*powl(c1,2) + 5*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(15.*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))));

aa[3] = (c2*(-324*(85*powl(c1,6) - 63*powl(c1,4)*powl(c2,2))*powl(atanhl(c2/c1),4) + 108*powl(c1,3)*powl(atanhl(c2/c1),3)*(670*powl(c1,2)*c2 - 441*powl(c2,3) + 3*(85*powl(c1,3) - 63*c1*powl(c2,2))*logl((c1 + c2)/(c1 - c2))) + c1*c2*atanhl(c2/c1)*(32*(555*powl(c1,2)*powl(c2,2) - 434*powl(c2,4)) + 36*(1245*powl(c1,3)*c2 - 736*c1*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) + 27*(510*powl(c1,4) - 313*powl(c1,2)*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(atanhl(c2/c1),2)*(-62340*powl(c1,4)*powl(c2,2) + 39312*powl(c1,2)*powl(c2,4) - 108*(590*powl(c1,5)*c2 - 377*powl(c1,3)*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) - 81*(85*powl(c1,6) - 63*powl(c1,4)*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(-240*powl(c1,2)*powl(c2,2) + 1792*powl(c2,4) - 144*(60*powl(c1,3)*c2 - 43*c1*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) - 9*(765*powl(c1,4) - 372*powl(c1,2)*powl(c2,2) - 20*powl(c2,4))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(35.*(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1))*(-8*c2 + 18*c1*atanhl(c2/c1) - 9*c1*logl((c1 + c2)/(c1 - c2)))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))));

bb[0] = 0;

bb[1] = (-2*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(powl(c2,2)*powl(logl((c1 + c2)/(c1 - c2)),2));

bb[2] = -(c2*logl((c1 + c2)/(c1 - c2))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))))/(9.*powl(c2 - c1*atanhl(c2/c1),2)*powl(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)),2));

bb[3] = ((c2 - c1*atanhl(c2/c1))*(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1))*(-8*c2 + 18*c1*atanhl(c2/c1) - 9*c1*logl((c1 + c2)/(c1 - c2)))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(25.*powl(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2))),2));

}

else if (n == 5){

aa[0] = (2*(c2 - c1*atanhl(c2/c1)))/(c2*logl((c1 + c2)/(c1 - c2)));

aa[1] = (3*powl(c1,2)*c2 + powl(c2,3) - 3*powl(c1,3)*atanhl(c2/c1) + (12*powl(c2 - c1*atanhl(c2/c1),3))/powl(logl((c1 + c2)/(c1 - c2)),2) + (12*c1*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2)))/(3.*c2*(c1*(-c2 + c1*atanhl(c2/c1)) - (2*powl(c2 - c1*atanhl(c2/c1),2))/logl((c1 + c2)/(c1 - c2))));

aa[2] = (c2*(108*powl(c1,4)*powl(atanhl(c2/c1),4) - 12*powl(c1,3)*powl(atanhl(c2/c1),3)*(31*c2 + 9*c1*logl((c1 + c2)/(c1 - c2))) + 9*powl(c1,2)*powl(atanhl(c2/c1),2)*(52*powl(c2,2) + 36*c1*c2*logl((c1 + c2)/(c1 - c2)) + 3*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) - 3*c1*c2*atanhl(c2/c1)*(84*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + 23*powl(c1,2)*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(48*powl(c2,2) + 108*c1*c2*logl((c1 + c2)/(c1 - c2)) + (42*powl(c1,2) + 5*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(15.*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))));

aa[3] = (c2*(-324*(85*powl(c1,6) - 63*powl(c1,4)*powl(c2,2))*powl(atanhl(c2/c1),4) + 108*powl(c1,3)*powl(atanhl(c2/c1),3)*(670*powl(c1,2)*c2 - 441*powl(c2,3) + 3*(85*powl(c1,3) - 63*c1*powl(c2,2))*logl((c1 + c2)/(c1 - c2))) + c1*c2*atanhl(c2/c1)*(32*(555*powl(c1,2)*powl(c2,2) - 434*powl(c2,4)) + 36*(1245*powl(c1,3)*c2 - 736*c1*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) + 27*(510*powl(c1,4) - 313*powl(c1,2)*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(atanhl(c2/c1),2)*(-62340*powl(c1,4)*powl(c2,2) + 39312*powl(c1,2)*powl(c2,4) - 108*(590*powl(c1,5)*c2 - 377*powl(c1,3)*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) - 81*(85*powl(c1,6) - 63*powl(c1,4)*powl(c2,2))*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(-240*powl(c1,2)*powl(c2,2) + 1792*powl(c2,4) - 144*(60*powl(c1,3)*c2 - 43*c1*powl(c2,3))*logl((c1 + c2)/(c1 - c2)) - 9*(765*powl(c1,4) - 372*powl(c1,2)*powl(c2,2) - 20*powl(c2,4))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(35.*(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1))*(-8*c2 + 18*c1*atanhl(c2/c1) - 9*c1*logl((c1 + c2)/(c1 - c2)))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))));

aa[4] = (c2*(14580*(623*powl(c1,8) - 870*powl(c1,6)*powl(c2,2) + 315*powl(c1,4)*powl(c2,4))*powl(atanhl(c2/c1),4) - 324*powl(c1,3)*powl(atanhl(c2/c1),3)*(68390*powl(c1,4)*c2 - 77850*powl(c1,2)*powl(c2,3) + 22932*powl(c2,5) + 45*(623*powl(c1,5) - 870*powl(c1,3)*powl(c2,2) + 315*c1*powl(c2,4))*logl((c1 + c2)/(c1 - c2))) + 9*powl(c1,2)*powl(atanhl(c2/c1),2)*(4*powl(c2,2)*(473515*powl(c1,4) - 438300*powl(c1,2)*powl(c2,2) + 123984*powl(c2,4)) + 72*(31115*powl(c1,5)*c2 - 34365*powl(c1,3)*powl(c2,3) + 9666*c1*powl(c2,5))*logl((c1 + c2)/(c1 - c2)) + 405*(623*powl(c1,6) - 870*powl(c1,4)*powl(c2,2) + 315*powl(c1,2)*powl(c2,4))*powl(logl((c1 + c2)/(c1 - c2)),2)) - 6*c1*c2*atanhl(c2/c1)*(128*(5145*powl(c1,4)*powl(c2,2) - 4225*powl(c1,2)*powl(c2,4) + 1533*powl(c2,6)) + 18*(121065*powl(c1,5)*c2 - 103780*powl(c1,3)*powl(c2,3) + 27792*c1*powl(c2,5))*logl((c1 + c2)/(c1 - c2)) + 81*(9345*powl(c1,6) - 9935*powl(c1,4)*powl(c2,2) + 2622*powl(c1,2)*powl(c2,4))*powl(logl((c1 + c2)/(c1 - c2)),2)) + powl(c2,2)*(-64*(315*powl(c1,4)*powl(c2,2) + 75*powl(c1,2)*powl(c2,4) - 1792*powl(c2,6)) + 576*(3465*powl(c1,5)*c2 - 2820*powl(c1,3)*powl(c2,3) + 943*c1*powl(c2,5))*logl((c1 + c2)/(c1 - c2)) + 81*(28035*powl(c1,6) - 20460*powl(c1,4)*powl(c2,2) + 4672*powl(c1,2)*powl(c2,4) + 144*powl(c2,6))*powl(logl((c1 + c2)/(c1 - c2)),2))))/(63.*(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1))*(8*c2 - 18*c1*atanhl(c2/c1) + 9*c1*logl((c1 + c2)/(c1 - c2)))*(-270*(7*powl(c1,5) - 5*powl(c1,3)*powl(c2,2))*powl(atanhl(c2/c1),2) + 15*powl(c1,2)*atanhl(c2/c1)*(182*powl(c1,2)*c2 - 96*powl(c2,3) + 9*(7*powl(c1,3) - 5*c1*powl(c2,2))*logl((c1 + c2)/(c1 - c2))) + c2*(-840*powl(c1,3)*c2 + 440*c1*powl(c2,3) + (-945*powl(c1,4) + 360*powl(c1,2)*powl(c2,2) + 36*powl(c2,4))*logl((c1 + c2)/(c1 - c2)))));

bb[0] = 0;

bb[1] = (-2*(c2 - c1*atanhl(c2/c1))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(powl(c2,2)*powl(logl((c1 + c2)/(c1 - c2)),2));

bb[2] = -(c2*logl((c1 + c2)/(c1 - c2))*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2)))))/(9.*powl(c2 - c1*atanhl(c2/c1),2)*powl(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2)),2));

bb[3] = ((c2 - c1*atanhl(c2/c1))*(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1))*(-8*c2 + 18*c1*atanhl(c2/c1) - 9*c1*logl((c1 + c2)/(c1 - c2)))*(2*c2 - 2*c1*atanhl(c2/c1) + c1*logl((c1 + c2)/(c1 - c2))))/(25.*powl(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2))),2));

bb[4] = (36*(6*powl(c1,3)*powl(atanhl(c2/c1),2) - 3*powl(c1,2)*atanhl(c2/c1)*(4*c2 + c1*logl((c1 + c2)/(c1 - c2))) + c2*(6*c1*c2 + (3*powl(c1,2) + powl(c2,2))*logl((c1 + c2)/(c1 - c2))))*(270*(7*powl(c1,5) - 5*powl(c1,3)*powl(c2,2))*powl(atanhl(c2/c1),2) - 15*powl(c1,2)*atanhl(c2/c1)*(182*powl(c1,2)*c2 - 96*powl(c2,3) + 9*(7*powl(c1,3) - 5*c1*powl(c2,2))*logl((c1 + c2)/(c1 - c2))) + c2*(840*powl(c1,3)*c2 - 440*c1*powl(c2,3) + 9*(105*powl(c1,4) - 40*powl(c1,2)*powl(c2,2) - 4*powl(c2,4))*logl((c1 + c2)/(c1 - c2)))))/(49.*powl(-15*powl(c1,2)*c2 + 4*powl(c2,3) + 3*(5*powl(c1,3) - 3*c1*powl(c2,2))*atanhl(c2/c1),2)*powl(8*c2 - 18*c1*atanhl(c2/c1) + 9*c1*logl((c1 + c2)/(c1 - c2)),2));

}

  else {
    mpi_abort("[D4EST_ERROR]: Do not support n >= 5 yet\n");
  }
 
}


void
d4est_quadrature_compactified_c1tpc2_neg2_aa_and_bb
(
 int n,
 long double* aa,
 long double* bb,
 void* user
){

  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  
  if (n == 1){
    aa[0] = (-(c1*c2) + (powl(c1,2) - powl(c2,2))*atanhl(c2/c1))/powl(c2,2);
    bb[0] = 0;
  }

  if (n == 2){
    aa[0] = (-(c1*c2) + (powl(c1,2) - powl(c2,2))*atanhl(c2/c1))/powl(c2,2);
    aa[1] = ((powl(c1,2) - powl(c2,2))*atanhl(c2/c1)*(-2*powl(c2,2) + c1*c2*atanhl(c2/c1) + (powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)))/(powl(c2,4) + (-(powl(c1,2)*powl(c2,2)) + powl(c2,4))*powl(atanhl(c2/c1),2));
    bb[0] = 0;
    bb[1] = -(((-powl(c1,2) + powl(c2,2))*(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2)))/powl(c2,4));
  }

  if (n == 3){

    aa[0] = (-(c1*c2) + (powl(c1,2) - powl(c2,2))*atanhl(c2/c1))/powl(c2,2);

    aa[1] = ((powl(c1,2) - powl(c2,2))*atanhl(c2/c1)*(-2*powl(c2,2) + c1*c2*atanhl(c2/c1) + (powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)))/(powl(c2,4) + (-(powl(c1,2)*powl(c2,2)) + powl(c2,4))*powl(atanhl(c2/c1),2));

    aa[2] = -(((powl(c1,2) - powl(c2,2))*(-3*c1*powl(c2,3) + (9*powl(c1,2)*powl(c2,2) - 5*powl(c2,4))*atanhl(c2/c1) + (-9*powl(c1,3)*c2 + 8*c1*powl(c2,3))*powl(atanhl(c2/c1),2) + (3*powl(c1,4) - powl(c1,2)*powl(c2,2) - 2*powl(c2,4))*powl(atanhl(c2/c1),3) + (-2*powl(c1,3)*c2 + 2*c1*powl(c2,3))*powl(atanhl(c2/c1),4)))/((powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2))*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))));

    bb[0] = 0;

    bb[1] = -(((-powl(c1,2) + powl(c2,2))*(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2)))/powl(c2,4));

    bb[2] = (-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))/(3.*powl(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2),2));

  }

  if (n == 4){

    aa[0] = (-(c1*c2) + (powl(c1,2) - powl(c2,2))*atanhl(c2/c1))/powl(c2,2);

    aa[1] = ((powl(c1,2) - powl(c2,2))*atanhl(c2/c1)*(-2*powl(c2,2) + c1*c2*atanhl(c2/c1) + (powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)))/(powl(c2,4) + (-(powl(c1,2)*powl(c2,2)) + powl(c2,4))*powl(atanhl(c2/c1),2));

    aa[2] = -(((powl(c1,2) - powl(c2,2))*(-3*c1*powl(c2,3) + (9*powl(c1,2)*powl(c2,2) - 5*powl(c2,4))*atanhl(c2/c1) + (-9*powl(c1,3)*c2 + 8*c1*powl(c2,3))*powl(atanhl(c2/c1),2) + (3*powl(c1,4) - powl(c1,2)*powl(c2,2) - 2*powl(c2,4))*powl(atanhl(c2/c1),3) + (-2*powl(c1,3)*c2 + 2*c1*powl(c2,3))*powl(atanhl(c2/c1),4)))/((powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2))*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))));

    aa[3] = -(((c1 - c2)*c2*(c1 + c2)*(-117*powl(c1,3)*powl(c2,4) + 96*c1*powl(c2,6) + (333*powl(c1,4)*powl(c2,3) - 345*powl(c1,2)*powl(c2,5) + 32*powl(c2,7))*atanhl(c2/c1) + (-297*powl(c1,5)*powl(c2,2) + 324*powl(c1,3)*powl(c2,4) - 31*c1*powl(c2,6))*powl(atanhl(c2/c1),2) + (63*powl(c1,6)*c2 + 3*powl(c1,4)*powl(c2,3) - 79*powl(c1,2)*powl(c2,5) + 13*powl(c2,7))*powl(atanhl(c2/c1),3) + 6*(3*powl(c1,7) - 13*powl(c1,5)*powl(c2,2) + 13*powl(c1,3)*powl(c2,4) - 3*c1*powl(c2,6))*powl(atanhl(c2/c1),4)))/((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2))));

    bb[0] = 0;

    bb[1] = -(((-powl(c1,2) + powl(c2,2))*(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2)))/powl(c2,4));

    bb[2] = (-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))/(3.*powl(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2),2));

    bb[3] = ((powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2))*(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2)))/(15.*powl(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2),2));

  }

  if (n == 5){

    aa[0] = (-(c1*c2) + (powl(c1,2) - powl(c2,2))*atanhl(c2/c1))/powl(c2,2);

    aa[1] = ((powl(c1,2) - powl(c2,2))*atanhl(c2/c1)*(-2*powl(c2,2) + c1*c2*atanhl(c2/c1) + (powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)))/(powl(c2,4) + (-(powl(c1,2)*powl(c2,2)) + powl(c2,4))*powl(atanhl(c2/c1),2));

    aa[2] = -(((powl(c1,2) - powl(c2,2))*(-3*c1*powl(c2,3) + (9*powl(c1,2)*powl(c2,2) - 5*powl(c2,4))*atanhl(c2/c1) + (-9*powl(c1,3)*c2 + 8*c1*powl(c2,3))*powl(atanhl(c2/c1),2) + (3*powl(c1,4) - powl(c1,2)*powl(c2,2) - 2*powl(c2,4))*powl(atanhl(c2/c1),3) + (-2*powl(c1,3)*c2 + 2*c1*powl(c2,3))*powl(atanhl(c2/c1),4)))/((powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2))*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))));

    aa[3] = -(((c1 - c2)*c2*(c1 + c2)*(-117*powl(c1,3)*powl(c2,4) + 96*c1*powl(c2,6) + (333*powl(c1,4)*powl(c2,3) - 345*powl(c1,2)*powl(c2,5) + 32*powl(c2,7))*atanhl(c2/c1) + (-297*powl(c1,5)*powl(c2,2) + 324*powl(c1,3)*powl(c2,4) - 31*c1*powl(c2,6))*powl(atanhl(c2/c1),2) + (63*powl(c1,6)*c2 + 3*powl(c1,4)*powl(c2,3) - 79*powl(c1,2)*powl(c2,5) + 13*powl(c2,7))*powl(atanhl(c2/c1),3) + 6*(3*powl(c1,7) - 13*powl(c1,5)*powl(c2,2) + 13*powl(c1,3)*powl(c2,4) - 3*c1*powl(c2,6))*powl(atanhl(c2/c1),4)))/((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2))));

    aa[4] = (2*(170100*powl(c1,8)*powl(powl(c1,2) - powl(c2,2),2)*(75*powl(c1,8) - 80*powl(c1,6)*powl(c2,2) + 50*powl(c1,4)*powl(c2,4) - 16*powl(c1,2)*powl(c2,6) + 3*powl(c2,8))*powl(atanhl(c2/c1),5) + 18*c1*powl(c2,2)*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)*(-2835000*powl(c1,14)*c2 + 3024000*powl(c1,12)*powl(c2,3) + 94500*powl(c1,10)*powl(c2,5) - 149325*powl(c1,8)*powl(c2,7) + 13360*powl(c1,6)*powl(c2,9) - 2885*powl(c1,4)*powl(c2,11) + 1016*powl(c1,2)*powl(c2,13) - 62*powl(c2,15) + 1050*powl(c1,7)*(2025*powl(c1,8) - 2835*powl(c1,6)*powl(c2,2) + 900*powl(c1,4)*powl(c2,4) - 93*powl(c1,2)*powl(c2,6) - 25*powl(c2,8))*logl(c1 - c2) - 1050*powl(c1,7)*(2025*powl(c1,8) - 2835*powl(c1,6)*powl(c2,2) + 900*powl(c1,4)*powl(c2,4) - 93*powl(c1,2)*powl(c2,6) - 25*powl(c2,8))*logl(c1 + c2)) + 27*c1*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4)*(-1890000*powl(c1,14)*c2 + 1386000*powl(c1,12)*powl(c2,3) - 567000*powl(c1,10)*powl(c2,5) + 83975*powl(c1,8)*powl(c2,7) - 4260*powl(c1,6)*powl(c2,9) + 270*powl(c1,4)*powl(c2,11) - 180*powl(c1,2)*powl(c2,13) + 27*powl(c2,15) + 3150*powl(c1,7)*(75*powl(c1,8) - 80*powl(c1,6)*powl(c2,2) + 50*powl(c1,4)*powl(c2,4) - 16*powl(c1,2)*powl(c2,6) + 3*powl(c2,8))*logl(c1 - c2) - 3150*powl(c1,7)*(75*powl(c1,8) - 80*powl(c1,6)*powl(c2,2) + 50*powl(c1,4)*powl(c2,4) - 16*powl(c1,2)*powl(c2,6) + 3*powl(c2,8))*logl(c1 + c2)) + c1*powl(c2,4)*(powl(c2,7)*(-675*powl(c1,8) + 14580*powl(c1,6)*powl(c2,2) - 29925*powl(c1,4)*powl(c2,4) + 19316*powl(c1,2)*powl(c2,6) - 3296*powl(c2,8)) + 3150*powl(c1,7)*(2025*powl(c1,8) - 3510*powl(c1,6)*powl(c2,2) + 675*powl(c1,4)*powl(c2,4) + 786*powl(c1,2)*powl(c2,6) + 32*powl(c2,8))*logl(c1 - c2) - 3150*powl(c1,7)*(2025*powl(c1,8) - 3510*powl(c1,6)*powl(c2,2) + 675*powl(c1,4)*powl(c2,4) + 786*powl(c1,2)*powl(c2,6) + 32*powl(c2,8))*logl(c1 + c2)) - 54*c2*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),3)*(-1417500*powl(c1,16)*c2 + 1984500*powl(c1,14)*powl(c2,3) - 630000*powl(c1,12)*powl(c2,5) + 65050*powl(c1,10)*powl(c2,7) + 17705*powl(c1,8)*powl(c2,9) - 85*powl(c1,6)*powl(c2,11) - 151*powl(c1,4)*powl(c2,13) + 87*powl(c1,2)*powl(c2,15) - 6*powl(c2,17) + 1050*powl(c1,7)*(450*powl(c1,10) - 780*powl(c1,8)*powl(c2,2) + 465*powl(c1,6)*powl(c2,4) - 155*powl(c1,4)*powl(c2,6) + 21*powl(c1,2)*powl(c2,8) - powl(c2,10))*logl(c1 - c2) - 1050*powl(c1,7)*(450*powl(c1,10) - 780*powl(c1,8)*powl(c2,2) + 465*powl(c1,6)*powl(c2,4) - 155*powl(c1,4)*powl(c2,6) + 21*powl(c1,2)*powl(c2,8) - powl(c2,10))*logl(c1 + c2)) - 2*powl(c2,3)*atanhl(c2/c1)*(-(c2*(6378750*powl(c1,16) - 11056500*powl(c1,14)*powl(c2,2) + 2126250*powl(c1,12)*powl(c2,4) + 2477250*powl(c1,10)*powl(c2,6) + 79065*powl(c1,8)*powl(c2,8) + 47610*powl(c1,6)*powl(c2,10) - 36019*powl(c1,4)*powl(c2,12) + 9194*powl(c1,2)*powl(c2,14) - 400*powl(c2,16))) + 3150*powl(c1,7)*(4050*powl(c1,10) - 8370*powl(c1,8)*powl(c2,2) + 4185*powl(c1,6)*powl(c2,4) + 348*powl(c1,4)*powl(c2,6) - 229*powl(c1,2)*powl(c2,8) + 16*powl(c2,10))*logl(c1 - c2) - 3150*powl(c1,7)*(4050*powl(c1,10) - 8370*powl(c1,8)*powl(c2,2) + 4185*powl(c1,6)*powl(c2,4) + 348*powl(c1,4)*powl(c2,6) - 229*powl(c1,2)*powl(c2,8) + 16*powl(c2,10))*logl(c1 + c2))))/(c2*(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2))*(-28350*powl(c1,7)*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),3) + powl(c2,2)*(-1575*powl(c1,6)*powl(c2,5) + 2010*powl(c1,4)*powl(c2,7) - 655*powl(c1,2)*powl(c2,9) + 256*powl(c2,11) - 1575*powl(c1,7)*(45*powl(c1,4) - 33*powl(c1,2)*powl(c2,2) - 16*powl(c2,4))*logl(c1 - c2) + 1575*powl(c1,7)*(45*powl(c1,4) - 33*powl(c1,2)*powl(c2,2) - 16*powl(c2,4))*logl(c1 + c2)) + 30*c1*c2*atanhl(c2/c1)*(-4725*powl(c1,10)*c2 + 3465*powl(c1,8)*powl(c2,3) + 1785*powl(c1,6)*powl(c2,5) - 169*powl(c1,4)*powl(c2,7) + 79*powl(c1,2)*powl(c2,9) - 15*powl(c2,11) + 315*powl(c1,7)*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 - c2) - 315*powl(c1,7)*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 + c2)) - 9*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)*(-31500*powl(c1,10)*c2 + 2100*powl(c1,8)*powl(c2,3) + 175*powl(c1,6)*powl(c2,5) - 165*powl(c1,4)*powl(c2,7) + 45*powl(c1,2)*powl(c2,9) + 9*powl(c2,11) + 1575*powl(c1,7)*(5*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 - c2) - 1575*powl(c1,7)*(5*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 + c2))));

    bb[0] = 0;

    bb[1] = -(((-powl(c1,2) + powl(c2,2))*(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2)))/powl(c2,4));

    bb[2] = (-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))/(3.*powl(powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2),2));

    bb[3] = ((powl(c2,2) + (-powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),2))*(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2)))/(15.*powl(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2),2));

    bb[4] = (4*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-28350*powl(c1,7)*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),3) + powl(c2,2)*(-1575*powl(c1,6)*powl(c2,5) + 2010*powl(c1,4)*powl(c2,7) - 655*powl(c1,2)*powl(c2,9) + 256*powl(c2,11) - 1575*powl(c1,7)*(45*powl(c1,4) - 33*powl(c1,2)*powl(c2,2) - 16*powl(c2,4))*logl(c1 - c2) + 1575*powl(c1,7)*(45*powl(c1,4) - 33*powl(c1,2)*powl(c2,2) - 16*powl(c2,4))*logl(c1 + c2)) + 30*c1*c2*atanhl(c2/c1)*(-4725*powl(c1,10)*c2 + 3465*powl(c1,8)*powl(c2,3) + 1785*powl(c1,6)*powl(c2,5) - 169*powl(c1,4)*powl(c2,7) + 79*powl(c1,2)*powl(c2,9) - 15*powl(c2,11) + 315*powl(c1,7)*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 - c2) - 315*powl(c1,7)*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 + c2)) - 9*(powl(c1,2) - powl(c2,2))*powl(atanhl(c2/c1),2)*(-31500*powl(c1,10)*c2 + 2100*powl(c1,8)*powl(c2,3) + 175*powl(c1,6)*powl(c2,5) - 165*powl(c1,4)*powl(c2,7) + 45*powl(c1,2)*powl(c2,9) + 9*powl(c2,11) + 1575*powl(c1,7)*(5*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 - c2) - 1575*powl(c1,7)*(5*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*logl(c1 + c2))))/(35.*powl(c2,5)*powl(-45*powl(c1,4)*powl(c2,2) + 33*powl(c1,2)*powl(c2,4) + 16*powl(c2,6) + 6*c1*c2*(15*powl(c1,4) - 16*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) - 9*(5*powl(c1,6) - 7*powl(c1,4)*powl(c2,2) + 3*powl(c1,2)*powl(c2,4) - powl(c2,6))*powl(atanhl(c2/c1),2),2));

  }
}


void
d4est_quadrature_compactified_c1tpc2_neg3_aa_and_bb
(
 int n,
 long double* aa,
 long double* bb,
 void* user
){
  long double c1 = ((d4est_quadrature_compactified_params_t*)user)->c1;
  long double c2 = ((d4est_quadrature_compactified_params_t*)user)->c2;
  

  if (n == 1){

    aa[0] = -(c2/c1);

    bb[0] = 0;

  }

  else if (n == 2){

    aa[0] = -(c2/c1);

    aa[1] = (3*powl(c1,2)*c2 - powl(c2,3) + (-3*powl(c1,3) + 2*c1*powl(c2,2))*atanhl(c2/c1))/(c1*c2*(-c2 + c1*atanhl(c2/c1)));

    bb[0] = 0;

    bb[1] = (powl(c1 - c2,2)*powl(c1 + c2,2)*(-c2 + c1*atanhl(c2/c1)))/(powl(c1,2)*powl(c2,3));

  }

  else if (n == 3){

    aa[0] = -(c2/c1);

    aa[1] = (3*powl(c1,2)*c2 - powl(c2,3) + (-3*powl(c1,3) + 2*c1*powl(c2,2))*atanhl(c2/c1))/(c1*c2*(-c2 + c1*atanhl(c2/c1)));

    aa[2] = (3*powl(c1,2)*powl(c2,3) + 4*powl(c2,5) + (-9*powl(c1,3)*powl(c2,2) + 4*c1*powl(c2,4))*atanhl(c2/c1) + (9*powl(c1,4)*c2 - 17*powl(c1,2)*powl(c2,3) + 9*powl(c2,5))*powl(atanhl(c2/c1),2) - 3*c1*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4))/(3.*(c2 - c1*atanhl(c2/c1))*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)));

    bb[0] = 0;

    bb[1] = (powl(c1 - c2,2)*powl(c1 + c2,2)*(-c2 + c1*atanhl(c2/c1)))/(powl(c1,2)*powl(c2,3));

    bb[2] = -((c1*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,3)*powl(c2 - c1*atanhl(c2/c1),2)));

  }

  else if (n == 4){

    aa[0] = -(c2/c1);

    aa[1] = (3*powl(c1,2)*c2 - powl(c2,3) + (-3*powl(c1,3) + 2*c1*powl(c2,2))*atanhl(c2/c1))/(c1*c2*(-c2 + c1*atanhl(c2/c1)));

    aa[2] = (3*powl(c1,2)*powl(c2,3) + 4*powl(c2,5) + (-9*powl(c1,3)*powl(c2,2) + 4*c1*powl(c2,4))*atanhl(c2/c1) + (9*powl(c1,4)*c2 - 17*powl(c1,2)*powl(c2,3) + 9*powl(c2,5))*powl(atanhl(c2/c1),2) - 3*c1*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4))/(3.*(c2 - c1*atanhl(c2/c1))*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)));

    aa[3] = (powl(c2,5)*(135*powl(c1,6) - 432*powl(c1,2)*powl(c2,4) + 320*powl(c2,6)) + (-675*powl(c1,7)*powl(c2,4) + 810*powl(c1,5)*powl(c2,6) + 396*powl(c1,3)*powl(c2,8) - 544*c1*powl(c2,10))*atanhl(c2/c1) + powl(c2,3)*(1350*powl(c1,8) - 3105*powl(c1,6)*powl(c2,2) + 2673*powl(c1,4)*powl(c2,4) - 1444*powl(c1,2)*powl(c2,6) + 528*powl(c2,8))*powl(atanhl(c2/c1),2) - 54*c1*powl(c2,2)*powl(powl(c1,2) - powl(c2,2),2)*(25*powl(c1,4) - 30*powl(c1,2)*powl(c2,2) + 18*powl(c2,4))*powl(atanhl(c2/c1),3) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(225*powl(c1,6) - 360*powl(c1,4)*powl(c2,2) + 111*powl(c1,2)*powl(c2,4) + 76*powl(c2,6))*powl(atanhl(c2/c1),4) - 135*c1*powl(powl(c1,2) - powl(c2,2),4)*(powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),5) + 27*c2*powl(powl(c1,2) - powl(c2,2),4)*(5*powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),6))/(15.*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3))*(-9*powl(c1,4)*powl(c2,3) + 30*powl(c1,2)*powl(c2,5) - 16*powl(c2,7) + (27*powl(c1,5)*powl(c2,2) - 69*powl(c1,3)*powl(c2,4) + 40*c1*powl(c2,6))*atanhl(c2/c1) - 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),2) + 9*c1*powl(powl(c1,2) - powl(c2,2),2)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),3)));

    bb[0] = 0;

    bb[1] = (powl(c1 - c2,2)*powl(c1 + c2,2)*(-c2 + c1*atanhl(c2/c1)))/(powl(c1,2)*powl(c2,3));

    bb[2] = -((c1*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,3)*powl(c2 - c1*atanhl(c2/c1),2)));

    bb[3] = -((c2 - c1*atanhl(c2/c1))*(9*powl(c1,4)*powl(c2,3) - 30*powl(c1,2)*powl(c2,5) + 16*powl(c2,7) + (-27*powl(c1,5)*powl(c2,2) + 69*powl(c1,3)*powl(c2,4) - 40*c1*powl(c2,6))*atanhl(c2/c1) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),2) - 9*c1*powl(powl(c1,2) - powl(c2,2),2)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),3)))/(9.*powl(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3),2));

  }

  else if (n == 5){

    aa[0] = -(c2/c1);

    aa[1] = (3*powl(c1,2)*c2 - powl(c2,3) + (-3*powl(c1,3) + 2*c1*powl(c2,2))*atanhl(c2/c1))/(c1*c2*(-c2 + c1*atanhl(c2/c1)));

    aa[2] = (3*powl(c1,2)*powl(c2,3) + 4*powl(c2,5) + (-9*powl(c1,3)*powl(c2,2) + 4*c1*powl(c2,4))*atanhl(c2/c1) + (9*powl(c1,4)*c2 - 17*powl(c1,2)*powl(c2,3) + 9*powl(c2,5))*powl(atanhl(c2/c1),2) - 3*c1*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4))/(3.*(c2 - c1*atanhl(c2/c1))*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)));

    aa[3] = (powl(c2,5)*(135*powl(c1,6) - 432*powl(c1,2)*powl(c2,4) + 320*powl(c2,6)) + (-675*powl(c1,7)*powl(c2,4) + 810*powl(c1,5)*powl(c2,6) + 396*powl(c1,3)*powl(c2,8) - 544*c1*powl(c2,10))*atanhl(c2/c1) + powl(c2,3)*(1350*powl(c1,8) - 3105*powl(c1,6)*powl(c2,2) + 2673*powl(c1,4)*powl(c2,4) - 1444*powl(c1,2)*powl(c2,6) + 528*powl(c2,8))*powl(atanhl(c2/c1),2) - 54*c1*powl(c2,2)*powl(powl(c1,2) - powl(c2,2),2)*(25*powl(c1,4) - 30*powl(c1,2)*powl(c2,2) + 18*powl(c2,4))*powl(atanhl(c2/c1),3) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(225*powl(c1,6) - 360*powl(c1,4)*powl(c2,2) + 111*powl(c1,2)*powl(c2,4) + 76*powl(c2,6))*powl(atanhl(c2/c1),4) - 135*c1*powl(powl(c1,2) - powl(c2,2),4)*(powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),5) + 27*c2*powl(powl(c1,2) - powl(c2,2),4)*(5*powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),6))/(15.*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3))*(-9*powl(c1,4)*powl(c2,3) + 30*powl(c1,2)*powl(c2,5) - 16*powl(c2,7) + (27*powl(c1,5)*powl(c2,2) - 69*powl(c1,3)*powl(c2,4) + 40*c1*powl(c2,6))*atanhl(c2/c1) - 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),2) + 9*c1*powl(powl(c1,2) - powl(c2,2),2)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),3)));

    aa[4] = ((-3*powl(c1,2)*powl(c2,2) + 8*powl(c2,4) + (6*powl(c1,3)*c2 - 8*c1*powl(c2,3))*atanhl(c2/c1) - 3*(powl(c1,4) - powl(c2,4))*powl(atanhl(c2/c1),2))*(22963500*powl(c1,9)*powl(powl(c1,2) - powl(c2,2),4)*powl(powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),7) + 243*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),6)*(-567000*powl(c1,12)*c2 - 693000*powl(c1,10)*powl(c2,3) - 126175*powl(c1,8)*powl(c2,5) + 900*powl(c1,6)*powl(c2,7) - 230*powl(c1,4)*powl(c2,9) + 60*powl(c1,2)*powl(c2,11) + 21*powl(c2,13) + 47250*powl(c1,9)*powl(powl(c1,2) + powl(c2,2),2)*logl(c1 - c2) - 47250*powl(c1,9)*powl(powl(c1,2) + powl(c2,2),2)*logl(c1 + c2)) - 1620*c1*c2*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),5)*(-212625*powl(c1,14)*c2 + 330750*powl(c1,12)*powl(c2,3) + 70980*powl(c1,10)*powl(c2,5) - 169920*powl(c1,8)*powl(c2,7) - 7028*powl(c1,6)*powl(c2,9) + 522*powl(c1,4)*powl(c2,11) - 105*powl(c1,2)*powl(c2,13) + 26*powl(c2,15) + 4725*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,4) + 11*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*logl(c1 - c2) - 4725*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,4) + 11*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*logl(c1 + c2)) + powl(c2,6)*(powl(c2,5)*(382725*powl(c1,10) - 1283850*powl(c1,8)*powl(c2,2) + 1638765*powl(c1,6)*powl(c2,4) - 878580*powl(c1,4)*powl(c2,6) + 59712*powl(c1,2)*powl(c2,8) + 81920*powl(c2,10)) + 141750*powl(c1,7)*powl(9*powl(c1,4) - 30*powl(c1,2)*powl(c2,2) + 16*powl(c2,4),2)*logl(c1 - c2) - 141750*powl(c1,7)*powl(9*powl(c1,4) - 30*powl(c1,2)*powl(c2,2) + 16*powl(c2,4),2)*logl(c1 + c2)) - 36*c1*powl(c2,3)*powl(atanhl(c2/c1),3)*(-9568125*powl(c1,16)*c2 + 48903750*powl(c1,14)*powl(c2,3) - 86018625*powl(c1,12)*powl(c2,5) + 61436250*powl(c1,10)*powl(c2,7) - 13008900*powl(c1,8)*powl(c2,9) - 1865860*powl(c1,6)*powl(c2,11) + 98488*powl(c1,4)*powl(c2,13) - 13046*powl(c1,2)*powl(c2,15) + 4568*powl(c2,17) + 47250*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*(135*powl(c1,6) - 315*powl(c1,4)*powl(c2,2) + 90*powl(c1,2)*powl(c2,4) + 64*powl(c2,6))*logl(c1 - c2) - 47250*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),2)*(135*powl(c1,6) - 315*powl(c1,4)*powl(c2,2) + 90*powl(c1,2)*powl(c2,4) + 64*powl(c2,6))*logl(c1 + c2)) + 27*powl(c2,2)*powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),4)*(-17010000*powl(c1,14)*c2 + 39690000*powl(c1,12)*powl(c2,3) - 11284875*powl(c1,10)*powl(c2,5) - 8152500*powl(c1,8)*powl(c2,7) + 14350*powl(c1,6)*powl(c2,9) + 19280*powl(c1,4)*powl(c2,11) - 1059*powl(c1,2)*powl(c2,13) + 1556*powl(c2,15) + 47250*powl(c1,7)*(135*powl(c1,8) - 210*powl(c1,6)*powl(c2,2) - 45*powl(c1,4)*powl(c2,4) + 108*powl(c1,2)*powl(c2,6) + 4*powl(c2,8))*logl(c1 - c2) - 47250*powl(c1,7)*(135*powl(c1,8) - 210*powl(c1,6)*powl(c2,2) - 45*powl(c1,4)*powl(c2,4) + 108*powl(c1,2)*powl(c2,6) + 4*powl(c2,8))*logl(c1 + c2)) - 4*c1*powl(c2,5)*atanhl(c2/c1)*(-5740875*powl(c1,14)*c2 + 38272500*powl(c1,12)*powl(c2,3) - 83731725*powl(c1,10)*powl(c2,5) + 66355200*powl(c1,8)*powl(c2,7) - 15835635*powl(c1,6)*powl(c2,9) - 1400580*powl(c1,4)*powl(c2,11) + 276144*powl(c1,2)*powl(c2,13) + 33184*powl(c2,15) + 70875*powl(c1,7)*(243*powl(c1,8) - 1431*powl(c1,6)*powl(c2,2) + 2862*powl(c1,4)*powl(c2,4) - 2304*powl(c1,2)*powl(c2,6) + 640*powl(c2,8))*logl(c1 - c2) - 70875*powl(c1,7)*(243*powl(c1,8) - 1431*powl(c1,6)*powl(c2,2) + 2862*powl(c1,4)*powl(c2,4) - 2304*powl(c1,2)*powl(c2,6) + 640*powl(c2,8))*logl(c1 + c2)) + 3*powl(c2,4)*powl(atanhl(c2/c1),2)*(-45927000*powl(c1,16)*c2 + 270459000*powl(c1,14)*powl(c2,3) - 539713125*powl(c1,12)*powl(c2,5) + 430879500*powl(c1,10)*powl(c2,7) - 114349050*powl(c1,8)*powl(c2,9) - 4347780*powl(c1,6)*powl(c2,11) + 1179279*powl(c1,4)*powl(c2,13) - 106072*powl(c1,2)*powl(c2,15) + 35264*powl(c2,17) + 47250*powl(c1,7)*(1215*powl(c1,10) - 6210*powl(c1,8)*powl(c2,2) + 10935*powl(c1,6)*powl(c2,4) - 7848*powl(c1,4)*powl(c2,6) + 1720*powl(c1,2)*powl(c2,8) + 192*powl(c2,10))*logl(c1 - c2) - 47250*powl(c1,7)*(1215*powl(c1,10) - 6210*powl(c1,8)*powl(c2,2) + 10935*powl(c1,6)*powl(c2,4) - 7848*powl(c1,4)*powl(c2,6) + 1720*powl(c1,2)*powl(c2,8) + 192*powl(c2,10))*logl(c1 + c2))))/(35.*powl(c2,4)*(-3*powl(c1,2)*c2 + 2*powl(c2,3) + 3*(powl(c1,3) - c1*powl(c2,2))*atanhl(c2/c1))*powl(3*powl(c1,2)*powl(c2,2) - 8*powl(c2,4) + (-6*powl(c1,3)*c2 + 8*c1*powl(c2,3))*atanhl(c2/c1) + 3*(powl(c1,4) - powl(c2,4))*powl(atanhl(c2/c1),2),2)*(675*powl(c1,7)*powl(c2,3) - 810*powl(c1,5)*powl(c2,5) - 765*powl(c1,3)*powl(c2,7) + 832*c1*powl(c2,9) + 3*powl(c2,2)*(-675*powl(c1,8) + 1035*powl(c1,6)*powl(c2,2) + 135*powl(c1,4)*powl(c2,4) - 551*powl(c1,2)*powl(c2,6) + 64*powl(c2,8))*atanhl(c2/c1) + 45*c1*c2*powl(powl(c1,2) - powl(c2,2),2)*(45*powl(c1,4) + 6*powl(c1,2)*powl(c2,2) + 5*powl(c2,4))*powl(atanhl(c2/c1),2) - 27*powl(powl(c1,2) - powl(c2,2),2)*(25*powl(c1,6) - 5*powl(c1,4)*powl(c2,2) + 15*powl(c1,2)*powl(c2,4) - 3*powl(c2,6))*powl(atanhl(c2/c1),3)));

    bb[0] = 0;

    bb[1] = (powl(c1 - c2,2)*powl(c1 + c2,2)*(-c2 + c1*atanhl(c2/c1)))/(powl(c1,2)*powl(c2,3));

    bb[2] = -((c1*(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,3)*powl(c2 - c1*atanhl(c2/c1),2)));

    bb[3] = -((c2 - c1*atanhl(c2/c1))*(9*powl(c1,4)*powl(c2,3) - 30*powl(c1,2)*powl(c2,5) + 16*powl(c2,7) + (-27*powl(c1,5)*powl(c2,2) + 69*powl(c1,3)*powl(c2,4) - 40*c1*powl(c2,6))*atanhl(c2/c1) + 3*c2*powl(powl(c1,2) - powl(c2,2),2)*(9*powl(c1,2) + 2*powl(c2,2))*powl(atanhl(c2/c1),2) - 9*c1*powl(powl(c1,2) - powl(c2,2),2)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),3)))/(9.*powl(2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3),2));

    bb[4] = ((2*c1*powl(c2,3) + (-3*powl(c1,2)*powl(c2,2) + 2*powl(c2,4))*atanhl(c2/c1) + powl(powl(c1,2) - powl(c2,2),2)*powl(atanhl(c2/c1),3))*(-675*powl(c1,7)*powl(c2,3) + 810*powl(c1,5)*powl(c2,5) + 765*powl(c1,3)*powl(c2,7) - 832*c1*powl(c2,9) - 3*powl(c2,2)*(-675*powl(c1,8) + 1035*powl(c1,6)*powl(c2,2) + 135*powl(c1,4)*powl(c2,4) - 551*powl(c1,2)*powl(c2,6) + 64*powl(c2,8))*atanhl(c2/c1) - 45*c1*c2*powl(powl(c1,2) - powl(c2,2),2)*(45*powl(c1,4) + 6*powl(c1,2)*powl(c2,2) + 5*powl(c2,4))*powl(atanhl(c2/c1),2) + 27*powl(powl(c1,2) - powl(c2,2),2)*(25*powl(c1,6) - 5*powl(c1,4)*powl(c2,2) + 15*powl(c1,2)*powl(c2,4) - 3*powl(c2,6))*powl(atanhl(c2/c1),3)))/(25.*powl(-3*powl(c1,2)*c2 + 2*powl(c2,3) + 3*(powl(c1,3) - c1*powl(c2,2))*atanhl(c2/c1),2)*powl(3*powl(c1,2)*powl(c2,2) - 8*powl(c2,4) + (-6*powl(c1,3)*c2 + 8*c1*powl(c2,3))*atanhl(c2/c1) + 3*(powl(c1,4) - powl(c2,4))*powl(atanhl(c2/c1),2),2));

  }


  else {
    mpi_abort("[D4EST_ERROR]: Do not support n >= 5 yet\n");
  }
  
}

typedef struct {

  long double weight;
  long double abscissa;

} weight_and_abscissa_pair_t;


static int
sort_weight_and_abscissa_pairs_callback
(
 const void* p,
 const void* q
)
{
  const weight_and_abscissa_pair_t x = *(const weight_and_abscissa_pair_t *)p;
  const weight_and_abscissa_pair_t y = *(const weight_and_abscissa_pair_t *)q;

  /* Avoid return x - y, which can cause undefined behaviour
       because of signed integer overflow. */
  if (x.abscissa == y.abscissa)
    return 0;
  else if (x.abscissa < y.abscissa)
    return -1;
  else
    return 1;
}

void
d4est_quadrature_compactified_compute_abscissas_and_weights
(
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* abscissas,
 double* weights,
 p4est_topidx_t tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int degree
)
{
  long double* custom_weights = P4EST_ALLOC(long double, degree + 1);
  long double* custom_abscissas = P4EST_ALLOC(long double, degree + 1);
  long double c1;
  long double c2;
  
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_OUTER_SHELL){
    /* transform element corners to [-1,1]^2 x [1,2] topological space */

    double amin = q0[0];
    double amax = q0[0] + dq;
    double bmin = q0[1];
    double bmax = q0[1] + dq;
    double cmin = q0[2];
    double cmax = q0[2] + dq;

    /* transform element corners to [0,1]^3 topological space */
    amin /= (double)P4EST_ROOT_LEN;
    amax /= (double)P4EST_ROOT_LEN;
    bmin /= (double)P4EST_ROOT_LEN;
    bmax /= (double)P4EST_ROOT_LEN;
    cmin /= (double)P4EST_ROOT_LEN;
    cmax /= (double)P4EST_ROOT_LEN;
    
    amin = 2.*amin - 1.;
    amax = 2.*amax - 1.;
    bmin = 2.*bmin - 1.;
    bmax = 2.*bmax - 1.;
    cmin = cmin + 1.;
    cmax = cmax + 1.;
    
    int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
    mpi_assert(compactify_outer_shell == 1);
    double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
    double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;
    c1 = ((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1);
    c2 = ((R2-R1)*(cmax-cmin));    
  }

  else if (d4est_geom->geom_type == GEOM_DISK_OUTER_WEDGE){
    /* topological coordinates of element corners */
    double amin = q0[0];
    double amax = q0[0] + dq;
    double bmin = q0[1];
    double bmax = q0[1] + dq;

    /* transform element corners to [0,1]^3 topological space */
    amin /= (double)P4EST_ROOT_LEN;
    amax /= (double)P4EST_ROOT_LEN;
    bmin /= (double)P4EST_ROOT_LEN;
    bmax /= (double)P4EST_ROOT_LEN;

    amin = 2.*amin - 1.;
    amax = 2.*amax - 1.;
    bmin = bmin + 1.;
    bmax = bmax + 1.;

    printf("amin,amax,bmin,bmax = %.15f,%.15f,%.15f,%.15f\n", amin, amax, bmin, bmax);
    
    int compactify_outer_wedge = ((d4est_geometry_disk_attr_t*)(d4est_geom->user))->compactify_outer_wedge;
    mpi_assert(compactify_outer_wedge == 1);
    double R1 = ((d4est_geometry_disk_attr_t*)(d4est_geom->user))->R1;
    double R2 = ((d4est_geometry_disk_attr_t*)(d4est_geom->user))->R2;
    c1 = ((R2-R1)*(bmax+bmin) - 4.0l*R2 + 2.0l*R1);
    c2 = ((R2-R1)*(bmax-bmin));
  }
  else {
    mpi_abort("[D4EST_ERROR]: d4est_quadrature_compactified does not support this type at the moment\n");
  }
  
  d4est_quadrature_compactified_params_t params;
  params.c1 = c1;
  params.c2 = c2;
  if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG1){
    arbquad_get_abscissas_and_weights_use_aa_and_bb(degree+1,
                                                    custom_weights,
                                                    custom_abscissas,
                                                    d4est_quadrature_compactified_c1tpc2_neg1_moment_fcn,
                                                    d4est_quadrature_compactified_c1tpc2_neg1_aa_and_bb,
                                                    &params,
                                                    DIVIDE_WEIGHTS_BY_WEIGHT_FCN,
                                                    d4est_quadrature_compactified_c1tpc2_neg1_weight_fcn
                                                   );
  }
  else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG2){
    arbquad_get_abscissas_and_weights_use_aa_and_bb(degree+1,
                                                    custom_weights,
                                                    custom_abscissas,
                                                    d4est_quadrature_compactified_c1tpc2_neg2_moment_fcn,
                                                    d4est_quadrature_compactified_c1tpc2_neg2_aa_and_bb,
                                                    &params,
                                                    DIVIDE_WEIGHTS_BY_WEIGHT_FCN,
                                                    d4est_quadrature_compactified_c1tpc2_neg2_weight_fcn
                                                   );
  }
  else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG3){
    arbquad_get_abscissas_and_weights_use_aa_and_bb(degree+1,
                                                    custom_weights,
                                                    custom_abscissas,
                                                    d4est_quadrature_compactified_c1tpc2_neg3_moment_fcn,
                                                    d4est_quadrature_compactified_c1tpc2_neg3_aa_and_bb,
                                                    &params,
                                                    DIVIDE_WEIGHTS_BY_WEIGHT_FCN,
                                                    d4est_quadrature_compactified_c1tpc2_neg3_weight_fcn
                                                   );
  }
  else if (d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG4){
    arbquad_get_abscissas_and_weights_use_aa_and_bb(degree+1,
                                                    custom_weights,
                                                    custom_abscissas,
                                                    d4est_quadrature_compactified_c1tpc2_neg4_moment_fcn,
                                                    d4est_quadrature_compactified_c1tpc2_neg4_aa_and_bb,
                                                    &params,
                                                    DIVIDE_WEIGHTS_BY_WEIGHT_FCN,
                                                    d4est_quadrature_compactified_c1tpc2_neg4_weight_fcn
                                                   );
  }
  else {
    mpi_abort("[D4EST_ERROR]: Not a supported type");
  }
    
  weight_and_abscissa_pair_t* pairs = P4EST_ALLOC(weight_and_abscissa_pair_t, degree + 1);
    
  for (int i = 0; i < degree + 1; i++){
    pairs[i].weight = custom_weights[i];
    pairs[i].abscissa = custom_abscissas[i];
  }

  qsort(pairs, degree + 1, sizeof(weight_and_abscissa_pair_t), sort_weight_and_abscissa_pairs_callback);
  
  for (int i = 0; i < degree + 1; i++){
    weights[i] = pairs[i].weight;
    abscissas[i] = pairs[i].abscissa;
  }
  
  P4EST_FREE(custom_weights);
  P4EST_FREE(custom_abscissas);
  P4EST_FREE(pairs);
  
}

