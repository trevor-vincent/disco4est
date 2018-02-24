#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_geometry.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_util.h>
#include <ini.h>
#include <sc_reduce.h>
#include <d4est_linalg.h>
#include <d4est_ghost_data.h>
static
int d4est_mesh_initial_extents_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_mesh_initial_extents_t* pconfig = (d4est_mesh_initial_extents_t*) user;
  if (d4est_util_match_couple(section,"initial_mesh",name,"min_quadrants")) {
    D4EST_ASSERT(pconfig->min_quadrants == -2);
    pconfig->min_quadrants = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name, "min_level")) {
    D4EST_ASSERT(pconfig->min_level == -2);
    pconfig->min_level = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"fill_uniform")) {
    D4EST_ASSERT(pconfig->fill_uniform == -2);
    pconfig->fill_uniform = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"load_from_checkpoint")) {
    D4EST_ASSERT(pconfig->load_from_checkpoint == 0);
    pconfig->load_from_checkpoint = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"checkpoint_prefix")) {
    D4EST_ASSERT(pconfig->checkpoint_prefix == NULL);
    asprintf(&pconfig->checkpoint_prefix,"%s", value);
  }
  else {
    for (int i = 0; i < pconfig->number_of_regions; i++){
      char* deg_name;
      char* deg_quad_name;
      asprintf(&deg_name,"region%d_deg", i);
      asprintf(&deg_quad_name,"region%d_deg_quad_inc", i);
      int hit = 0;
      if (d4est_util_match_couple(section,"initial_mesh",name,deg_name)) {
        D4EST_ASSERT(pconfig->deg[i] == -1);
        pconfig->deg[i] = atoi(value);
        hit++;
      }
      if (d4est_util_match_couple(section,"initial_mesh",name,deg_quad_name)) {
        D4EST_ASSERT(pconfig->deg_quad_inc[i] == -1);
        pconfig->deg_quad_inc[i] = atoi(value);
        hit++;
      }
      free(deg_name);
      free(deg_quad_name);
      if (hit)
        return 1;
    }
    return 0;
  }
  return 1;
}


void
d4est_mesh_set_initial_extents
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  d4est_mesh_initial_extents_t* input = user_ctx;
  int region = elem_data->region;

  if (input->load_from_checkpoint == 0)
    elem_data->deg = input->deg[region];
  else{
    D4EST_ASSERT(input->checkpoint_deg_array != NULL);
    elem_data->deg = input->checkpoint_deg_array[elem_data->id];
  }
  
  elem_data->deg_quad = input->deg_quad_inc[region] +
                            elem_data->deg;
  /* printf("[BEFORE_AMR]: deg, deg_quad, deg_quad_inc, region = %d, %d, %d, %d\n", elem_data->deg, elem_data->deg_quad, input->deg_quad_inc[region], region); */
}


void
d4est_mesh_set_quadratures_after_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  d4est_mesh_initial_extents_t* input = user_ctx;
  int region = elem_data->region;
  elem_data->deg_quad = input->deg_quad_inc[region] +
                            elem_data->deg;

  /* printf("[AFTER_AMR]: deg, deg_quad, deg_quad_inc, region = %d, %d, %d, %d\n", elem_data->deg, elem_data->deg_quad, input->deg_quad_inc[region], region); */
               
}

void
d4est_mesh_initial_extents_destroy
(
 d4est_mesh_initial_extents_t* initial_extents
){
  P4EST_FREE(initial_extents->deg);
  P4EST_FREE(initial_extents->deg_quad_inc);
  free(initial_extents->checkpoint_prefix);
  P4EST_FREE(initial_extents);
}

d4est_mesh_initial_extents_t*
d4est_mesh_initial_extents_parse
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
  d4est_mesh_initial_extents_t* initial_extents = P4EST_ALLOC(d4est_mesh_initial_extents_t, 1);
  initial_extents->min_quadrants = -2;
  initial_extents->min_level = -2;
  initial_extents->fill_uniform = -2;
  initial_extents->number_of_regions = d4est_geom->get_number_of_regions(d4est_geom);
  initial_extents->deg = P4EST_ALLOC(int, initial_extents->number_of_regions);
  initial_extents->deg_quad_inc = P4EST_ALLOC(int, initial_extents->number_of_regions);
  initial_extents->load_from_checkpoint = 0;
  initial_extents->checkpoint_prefix = NULL;
  
 for (int i = 0; i < initial_extents->number_of_regions; i++){
    initial_extents->deg[i] = -1;
    initial_extents->deg_quad_inc[i] = -1;
  }
 
  if (ini_parse(input_file, d4est_mesh_initial_extents_handler, initial_extents) < 0) {
    printf("[D4EST_ERROR]: pXest input_file = %s\n", input_file);
    D4EST_ABORT("[D4EST_ERROR]: Can't load pXest input file");
  }


  for (int i = 0; i < initial_extents->number_of_regions; i++){
    D4EST_CHECK_INPUT("initial_mesh", initial_extents->deg_quad_inc[i], -1);
  }

  if (initial_extents->load_from_checkpoint == 0){

    D4EST_CHECK_INPUT("initial_mesh", initial_extents->min_quadrants, -2);
    D4EST_CHECK_INPUT("initial_mesh", initial_extents->min_level, -2);
    D4EST_CHECK_INPUT("initial_mesh", initial_extents->fill_uniform, -2);
    
    for (int i = 0; i < initial_extents->number_of_regions; i++){
      D4EST_CHECK_INPUT("initial_mesh", initial_extents->deg[i], -1);
    }
    
  }
  else {
    D4EST_CHECK_INPUT("initial_mesh", initial_extents->checkpoint_prefix, NULL);
  }
  
  if(
     proc_size != 1
     &&
     d4est_util_int_pow_int((P4EST_CHILDREN), initial_extents->min_level) < proc_size
     &&
     initial_extents->min_quadrants < proc_size
  ){
    if (proc_rank == 0){
      printf("[D4EST_ERROR]: proc_size = %d\n", proc_size);
      printf("[D4EST_ERROR]: min_level = %d\n", initial_extents->min_level);
      D4EST_ABORT("[D4EST_ERROR]: Starting p4est with elements < processes\n");
    }
  }

  return initial_extents;
}


void
d4est_mesh_compute_mortar_quadrature_quantities_boundary_callback
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params
)
{
  int total_mortar_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg_quad);
  int scalar_stride = e_m->mortar_quad_scalar_stride[f_m];
  int vector_stride = e_m->mortar_quad_vector_stride[f_m];
  int matrix_stride = e_m->mortar_quad_matrix_stride[f_m];
  int vector_boundary_stride = e_m->boundary_quad_vector_stride[f_m];

  double* n_m_mortar_quad [(P4EST_DIM)];
  double* xyz_m_mortar_quad [(P4EST_DIM)];
  double* xyz_m_mortar_lobatto [(P4EST_DIM)];
  double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  double* sj_m_mortar_quad = &d4est_factors->sj_m_mortar_quad[scalar_stride];
  
  for (int i = 0; i < (P4EST_DIM); i++){
    n_m_mortar_quad[i] = &d4est_factors->n_m_mortar_quad[vector_stride + i*total_mortar_nodes_quad];
    xyz_m_mortar_quad[i] = &d4est_factors->xyz_m_mortar_quad[vector_boundary_stride + i*total_mortar_nodes_quad];
    xyz_m_mortar_lobatto[i] = &d4est_factors->xyz_m_mortar_lobatto[vector_boundary_stride + i*total_mortar_nodes_quad];
    for (int j = 0; j < (P4EST_DIM); j++){
      drst_dxyz_m_mortar_quad[i][j] = &d4est_factors->drst_dxyz_m_mortar_quad[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
      drst_dxyz_p_mortar_quad_porder[i][j] = &d4est_factors->drst_dxyz_p_mortar_quad_porder[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
    }
  }

  for (int d = 0; d < (P4EST_DIM); d++){

    d4est_operators_apply_slicer(d4est_ops,
                                 e_m->xyz[d],
                                 (P4EST_DIM),
                                 f_m,
                                 e_m->deg,
                                 xyz_m_mortar_lobatto[d]);

  }
 
  
  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_m->tree,
     e_m->q,
     e_m->dq,
     mortar_side_id_m,
     1,
     1,
     &e_m->deg,
     &e_m->deg_quad,
     f_m,
     xyz_m_mortar_quad,
     drst_dxyz_m_mortar_quad,
     sj_m_mortar_quad,
     n_m_mortar_quad,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
}
  
void
d4est_mesh_compute_mortar_quadrature_quantities_interface_callback
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
 d4est_mesh_data_t* d4est_factors,
 void* params
){

  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int deg_m_quad [(P4EST_HALF)];
  int deg_p_quad [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_quad [(P4EST_HALF)];
  int mortar_nodes_quad [(P4EST_HALF)];
  int mortar_nodes_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    deg_m_quad[i] = e_m[i]->deg_quad;
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_m_quad[i]);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_quad[i] = e_p_oriented[i]->deg_quad;
    deg_p_lobatto_porder[i] = e_p[i]->deg;

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_p_quad[i]);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }

  /* calculate degs and nodes of the mortar faces */
  int total_mortar_nodes_quad = 0;
  int total_mortar_nodes_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = d4est_util_max_int(deg_m_quad[i],
                                          deg_p_quad[j]);
      deg_mortar_lobatto[i+j] = d4est_util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );
      mortar_nodes_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );
      mortar_nodes_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );
      total_mortar_nodes_quad += mortar_nodes_quad[i+j];
      total_mortar_nodes_lobatto += mortar_nodes_lobatto[i+j];

    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int deg_mortar_lobatto_porder [(P4EST_HALF)];
  int mortar_nodes_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    deg_mortar_lobatto_porder[inew] = deg_mortar_lobatto[i];
    mortar_nodes_quad_porder[inew] = mortar_nodes_quad[i];
  }
  
  int scalar_stride = e_m[0]->mortar_quad_scalar_stride[f_m];
  int vector_stride = e_m[0]->mortar_quad_vector_stride[f_m];
  int matrix_stride = e_m[0]->mortar_quad_matrix_stride[f_m];
  
  double* n_m_mortar_quad [(P4EST_DIM)];
  double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  double* sj_m_mortar_quad = &d4est_factors->sj_m_mortar_quad[scalar_stride];
  
  for (int i = 0; i < (P4EST_DIM); i++){
    n_m_mortar_quad[i] = &d4est_factors->n_m_mortar_quad[vector_stride + i*total_mortar_nodes_quad];
    for (int j = 0; j < (P4EST_DIM); j++){
      drst_dxyz_m_mortar_quad[i][j] = &d4est_factors->drst_dxyz_m_mortar_quad[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
      drst_dxyz_p_mortar_quad_porder[i][j] = &d4est_factors->drst_dxyz_p_mortar_quad_porder[matrix_stride + (i + j*(P4EST_DIM))*total_mortar_nodes_quad];
    }
  }

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     mortar_side_id_m,
     faces_m,
     faces_mortar,
     &deg_mortar_lobatto[0],
     &deg_mortar_quad[0],
     f_m,
     NULL,
     drst_dxyz_m_mortar_quad,
     sj_m_mortar_quad,
     n_m_mortar_quad,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     mortar_side_id_p,
     faces_p,
     faces_mortar,
     &deg_mortar_lobatto_porder[0],
     &deg_mortar_quad_porder[0],
     f_p,
     NULL,
     drst_dxyz_p_mortar_quad_porder,
     NULL,
     NULL,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
}

void
d4est_mesh_compute_mortar_quadrature_quantities
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{

  /* printf("[D4EST_INFO]: Computing quadrature quantities\n"); */
  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_mesh_compute_mortar_quadrature_quantities_interface_callback;
  /* flux_fcns.flux_interface_fcn = NULL; */
  flux_fcns.flux_boundary_fcn = d4est_mesh_compute_mortar_quadrature_quantities_boundary_callback;
  flux_fcns.user_ctx = (void*)d4est_factors;
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns,
     EXCHANGE_GHOST_DATA /* might not be needed because it's done for size computations */
    );

}


static void
d4est_mesh_compute_mortar_quadrature_sizes_boundary_callback
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params
)
{
  d4est_local_sizes_t* sizes = params;

  e_m->boundary_quad_vector_stride[f_m] = (P4EST_DIM)*sizes->local_boundary_nodes_quad;
  e_m->mortar_quad_scalar_stride[f_m] = sizes->local_mortar_nodes_quad;
  e_m->mortar_quad_vector_stride[f_m] = sizes->local_mortar_nodes_quad*(P4EST_DIM);
  e_m->mortar_quad_matrix_stride[f_m] = sizes->local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  
  sizes->local_mortar_nodes_quad += d4est_lgl_get_nodes(
                                                        (P4EST_DIM)-1,
                                                        e_m->deg_quad);

  sizes->local_boundary_nodes_quad += d4est_lgl_get_nodes(
                                                        (P4EST_DIM)-1,
                                                        e_m->deg_quad);
}

static void
d4est_mesh_compute_mortar_quadrature_sizes_interface_callback
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
 d4est_mesh_data_t* d4est_factors,
 void* params
)
{
  d4est_local_sizes_t* sizes = params;
  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of the mortar faces */
  int total_mortar_nodes_quad = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      int deg_mortar_quad = d4est_util_max_int
                            (
                             e_m[i]->deg_quad,
                             e_p_oriented[j]->deg_quad
                            );
      total_mortar_nodes_quad += d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad );
    }

  e_m[0]->mortar_quad_scalar_stride[f_m] = sizes->local_mortar_nodes_quad;
  e_m[0]->mortar_quad_vector_stride[f_m] = sizes->local_mortar_nodes_quad*(P4EST_DIM);
  e_m[0]->mortar_quad_matrix_stride[f_m] = sizes->local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);

  for (int f = 1; f < faces_m; f++){
    e_m[f]->mortar_quad_scalar_stride[f_m] = -1;
    e_m[f]->mortar_quad_vector_stride[f_m] = -1;
    e_m[f]->mortar_quad_matrix_stride[f_m] = -1;
  }
  
  sizes->local_mortar_nodes_quad += total_mortar_nodes_quad;
}



void
d4est_mesh_compute_mortar_quadrature_sizes
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_local_sizes_t* local_sizes
){

  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_mesh_compute_mortar_quadrature_sizes_interface_callback;
  flux_fcns.flux_boundary_fcn = d4est_mesh_compute_mortar_quadrature_sizes_boundary_callback;
  flux_fcns.user_ctx = (void*)local_sizes;
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns,
     EXCHANGE_GHOST_DATA
    );
}

d4est_mesh_data_t*
d4est_mesh_geometry_storage_init()
{
  d4est_mesh_data_t* d4est_factors = P4EST_ALLOC(d4est_mesh_data_t, 1);
  d4est_factors->J_quad = NULL;
  d4est_factors->xyz = NULL;
  d4est_factors->xyz_quad = NULL;
  d4est_factors->xyz_rst_quad = NULL;
  d4est_factors->rst_xyz_quad = NULL;

  d4est_factors->drst_dxyz_m_mortar_quad = NULL;
  d4est_factors->drst_dxyz_p_mortar_quad_porder = NULL;
  d4est_factors->sj_m_mortar_quad = NULL;
  d4est_factors->n_m_mortar_quad = NULL;
  d4est_factors->xyz_m_mortar_quad = NULL;
  d4est_factors->xyz_m_mortar_lobatto = NULL;
  
  return d4est_factors;
}

static void
d4est_mesh_geometry_storage_realloc
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_local_sizes_t local_sizes
)
{
  d4est_factors->local_sizes = local_sizes;

  if (p4est->mpirank == 0)
    printf("[D4EST_INFO]: Reallocing storage for geometry data\n");
  int local_nodes = local_sizes.local_nodes;
  int local_nodes_quad = local_sizes.local_nodes_quad;

  int vector_nodes = local_nodes*(P4EST_DIM);
  d4est_factors->xyz = P4EST_REALLOC(d4est_factors->xyz,double,vector_nodes);
  d4est_factors->xyz_quad = P4EST_REALLOC(d4est_factors->xyz_quad, double, (P4EST_DIM)*local_nodes_quad);

  int matrix_nodes_quad = local_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  d4est_factors->J_quad = P4EST_REALLOC(d4est_factors->J_quad,double,local_nodes_quad);
  d4est_factors->xyz_rst_quad = P4EST_REALLOC(d4est_factors->xyz_rst_quad,double,matrix_nodes_quad);
  d4est_factors->rst_xyz_quad = P4EST_REALLOC(d4est_factors->rst_xyz_quad,double,matrix_nodes_quad);

  int local_matrix_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);

  int local_vector_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM);
  int local_vector_boundary_nodes_quad = local_sizes.local_boundary_nodes_quad*(P4EST_DIM);
  
  d4est_factors->drst_dxyz_m_mortar_quad = P4EST_REALLOC(d4est_factors->drst_dxyz_m_mortar_quad,
                                                             double,
                                                             local_matrix_mortar_nodes_quad);

  d4est_factors->drst_dxyz_p_mortar_quad_porder = P4EST_REALLOC(d4est_factors->drst_dxyz_p_mortar_quad_porder,
                                                                    double,
                                                                    local_matrix_mortar_nodes_quad);


  d4est_factors->sj_m_mortar_quad = P4EST_REALLOC(d4est_factors->sj_m_mortar_quad,
                                                      double,
                                                      local_sizes.local_mortar_nodes_quad);


  d4est_factors->n_m_mortar_quad = P4EST_REALLOC(d4est_factors->n_m_mortar_quad,
                                                     double,
                                                     local_vector_mortar_nodes_quad);
  
  d4est_factors->xyz_m_mortar_quad = P4EST_REALLOC(d4est_factors->xyz_m_mortar_quad,
                                                   double,
                                                   local_vector_boundary_nodes_quad);

  d4est_factors->xyz_m_mortar_lobatto = P4EST_REALLOC(d4est_factors->xyz_m_mortar_lobatto,
                                                      double,
                                                      local_vector_boundary_nodes_quad);
  
}


void
d4est_mesh_geometry_storage_printout
(
 d4est_mesh_data_t* d4est_factors
)
{
  /* int local_nodes = d4est_factors->local_sizes.local_nodes; */
  /* int local_nodes_quad = d4est_factors->local_sizes.local_nodes_quad; */
  /* int vector_nodes = local_nodes*(P4EST_DIM);  */
  /* int matrix_nodes_quad = local_nodes_quad*(P4EST_DIM)*(P4EST_DIM); */
  int local_matrix_mortar_nodes_quad = d4est_factors->local_sizes.local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_vector_mortar_nodes_quad = d4est_factors->local_sizes.local_mortar_nodes_quad*(P4EST_DIM);
  int local_vector_boundary_nodes_quad = d4est_factors->local_sizes.local_boundary_nodes_quad*(P4EST_DIM);

  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->drst_dxyz_p_mortar_quad_porder, local_matrix_mortar_nodes_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->drst_dxyz_m_mortar_quad, local_matrix_mortar_nodes_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->xyz_m_mortar_lobatto, local_vector_boundary_nodes_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->xyz_m_mortar_quad, local_vector_boundary_nodes_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->n_m_mortar_quad, local_vector_mortar_nodes_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(d4est_factors->sj_m_mortar_quad, d4est_factors->local_sizes.local_mortar_nodes_quad); */
}


void
d4est_mesh_geometry_storage_destroy
(
 d4est_mesh_data_t* d4est_factors
)
{
  if (d4est_factors != NULL){
    P4EST_FREE(d4est_factors->J_quad);
    P4EST_FREE(d4est_factors->xyz);
    P4EST_FREE(d4est_factors->xyz_quad);
    P4EST_FREE(d4est_factors->xyz_rst_quad);
    P4EST_FREE(d4est_factors->rst_xyz_quad);
    P4EST_FREE(d4est_factors->drst_dxyz_m_mortar_quad);
    P4EST_FREE(d4est_factors->drst_dxyz_p_mortar_quad_porder);
    P4EST_FREE(d4est_factors->n_m_mortar_quad);
    P4EST_FREE(d4est_factors->xyz_m_mortar_quad);
    P4EST_FREE(d4est_factors->xyz_m_mortar_lobatto);
    P4EST_FREE(d4est_factors->sj_m_mortar_quad);
    P4EST_FREE(d4est_factors);
  }
}

static void
d4est_mesh_init_element_size_parameters
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t,1);
  d4est_quadrature_lobatto_new(d4est_quad, NULL, NULL);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree; tt++)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);

        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto = d4est_quadrature_get_rst_points
                             (
                              d4est_ops,
                              d4est_quad,
                              d4est_geom,
                              NULL,
                              QUAD_OBJECT_VOLUME,
                              QUAD_INTEGRAND_UNKNOWN,
                              ed->deg
                             );

        int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        int face_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM)-1, ed->deg);
  
        double* xyz_lobatto [(P4EST_DIM)];
        D4EST_ALLOC_DIM_VEC(xyz_lobatto, volume_nodes_lobatto);
        double* dxyz_drst_lobatto [(P4EST_DIM)][(P4EST_DIM)];
        D4EST_ALLOC_DBYD_MAT(dxyz_drst_lobatto, volume_nodes_lobatto);
        double* j_lobatto = P4EST_ALLOC(double, volume_nodes_lobatto);

        double* xyz_on_f_lobatto [(P4EST_DIM)];
        D4EST_ALLOC_DIM_VEC(xyz_on_f_lobatto, face_nodes_lobatto);
        double* sj_on_f_lobatto = P4EST_ALLOC(double, face_nodes_lobatto);
        double* j_div_sj_on_f_lobatto = P4EST_ALLOC(double, face_nodes_lobatto);
  
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           ed->tree,
           ed->deg,
           ed->q,
           ed->dq,
           xyz_lobatto
          );

        double diam_volume = -1;
        for (int i = 0; i < volume_nodes_lobatto; i++){
          for (int j = i + 1; j < volume_nodes_lobatto; j++){
            double diam_temp = 0.;
            for (int d = 0; d < (P4EST_DIM); d++){
              diam_temp += pow(fabs(xyz_lobatto[d][i] - xyz_lobatto[d][j]), 2);
            }
            diam_temp = sqrt(diam_temp);
            if (diam_temp > diam_volume)
              diam_volume = diam_temp;
          }
        }
        ed->diam_volume = diam_volume;
  
        d4est_geometry_compute_dxyz_drst
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           ed->tree,
           ed->q,
           ed->dq,
           ed->deg,
           dxyz_drst_lobatto
          );

    
        d4est_geometry_compute_jacobian
          (
           dxyz_drst_lobatto,
           j_lobatto,
           volume_nodes_lobatto
          );
  
        double* ones = P4EST_ALLOC(double, volume_nodes_lobatto);
        for (int i = 0; i < volume_nodes_lobatto; i++)
          ones[i] = 1.;

        ed->volume = d4est_quadrature_innerproduct
                     (
                      d4est_ops,
                      d4est_geom,
                      d4est_quad,
                      NULL,
                      QUAD_OBJECT_VOLUME,
                      QUAD_INTEGRAND_UNKNOWN,
                      ones,
                      NULL,
                      j_lobatto,
                      ed->deg
                     );


        for (int f = 0; f < (P4EST_FACES); f++){

          d4est_mortars_compute_geometric_data_on_mortar
            (
             d4est_ops,
             d4est_geom,
             d4est_quad,
             QUAD_INTEGRAND_UNKNOWN,
             ed->tree,
             ed->q,
             ed->dq,
             -1,
             1,
             1,
             &ed->deg,
             &ed->deg,
             f,
             xyz_on_f_lobatto,
             NULL,
             sj_on_f_lobatto,
             NULL,
             NULL,
             j_div_sj_on_f_lobatto,
             COMPUTE_NORMAL_USING_JACOBIAN
            );

          ed->area[f] = d4est_quadrature_innerproduct
                        (
                         d4est_ops,
                         d4est_geom,
                         d4est_quad,
                         NULL,
                         QUAD_OBJECT_MORTAR,
                         QUAD_INTEGRAND_UNKNOWN,
                         ones,
                         NULL,
                         sj_on_f_lobatto,
                         ed->deg
                        );


          double diam_face = -1;
          for (int i = 0; i < face_nodes_lobatto; i++){
            for (int j = i + 1; j < face_nodes_lobatto; j++){
              double diam_temp = 0.;
              for (int d = 0; d < (P4EST_DIM); d++){
                diam_temp += pow(fabs(xyz_on_f_lobatto[d][i] - xyz_on_f_lobatto[d][j]), 2);
              }
              diam_temp = sqrt(diam_temp);
              if (diam_temp > diam_face)
                diam_face = diam_temp;
            }
          }

          ed->diam_face[f] = diam_face;

          ed->j_div_sj_min[f] = d4est_util_min_dbl_array(j_div_sj_on_f_lobatto, face_nodes_lobatto);
 
        }
     
        P4EST_FREE(ones);

        D4EST_FREE_DIM_VEC(xyz_lobatto);
        D4EST_FREE_DBYD_MAT(dxyz_drst_lobatto);
        D4EST_FREE_DIM_VEC(xyz_on_f_lobatto);
        D4EST_FREE(sj_on_f_lobatto);
        D4EST_FREE(j_div_sj_on_f_lobatto);
        D4EST_FREE(j_lobatto);
      }
    }


  P4EST_FREE(d4est_quad);
}

void
d4est_mesh_get_array_of_degrees
(
 p4est_t* p4est,
 void* deg_array,
 d4est_builtin_t type
)
{
  int stride = 0;
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
        if (type == D4EST_INT)
          ((int*)deg_array)[stride] = ed->deg;
        else if (type == D4EST_DOUBLE)
          ((double*)deg_array)[stride] = (double)ed->deg;
        else
          D4EST_ABORT("Not a supported builtin-type");
        stride++;
      }
    }
}

void
d4est_mesh_get_array_of_quadrature_degrees
(
 p4est_t* p4est,
 void* deg_array,
 d4est_builtin_t type
)
{
  int stride = 0;
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
        if (type == D4EST_INT)
          ((int*)deg_array)[stride] = ed->deg_quad;
        else if (type == D4EST_DOUBLE)
          ((double*)deg_array)[stride] = (double)ed->deg_quad;
        else
          D4EST_ABORT("Not a supported builtin-type");
        stride++;
      }
    }
}


void
d4est_mesh_get_array_of_estimators
(
 p4est_t* p4est,
 double* eta2_array
)
{
  int stride = 0;
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
        eta2_array[stride] = ed->local_estimator;
        stride++;
      }
    }
}

void
d4est_mesh_compute_jacobian_on_lgl_grid
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double* jacobian_lgl
)
{
  int nodal_stride = 0;
  double* xyz_rst [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(xyz_rst, d4est_lgl_get_nodes((P4EST_DIM),(MAX_DEGREE)));

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg);


        d4est_quadrature_volume_t mesh_object = {.dq = elem_data->dq,
                                                 .tree = elem_data->tree,
                                                 .q[0] = elem_data->q[0],
                                                 .q[1] = elem_data->q[1],
#if (P4EST_DIM)==3
                                                 .q[2] = elem_data->q[2],
#endif
                                                 .element_id = elem_data->id
                                                };
      
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg,
                                                                0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL, NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg, 2);
#endif

        d4est_geometry_compute_dxyz_drst
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           elem_data->tree,
           elem_data->q,
           elem_data->dq,
           elem_data->deg,
           xyz_rst
          );


        d4est_geometry_compute_jacobian
          (
           xyz_rst,
           &jacobian_lgl[nodal_stride],
           volume_nodes
          );

        
        nodal_stride += volume_nodes;
      }
    }
                       
  D4EST_FREE_DBYD_MAT(xyz_rst);
}

int d4est_mesh_get_local_matrix_nodes(p4est_t* p4est){

  int local_matrix_nodes = 0;
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        local_matrix_nodes += volume_nodes*volume_nodes;
      }
    }
  return local_matrix_nodes;
}

int d4est_mesh_get_local_quad_nodes(p4est_t* p4est){

  int local_quad_nodes = 0;
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
        local_quad_nodes += volume_nodes;
      }
    }
  return local_quad_nodes;
}


void
d4est_mesh_print_number_of_elements_per_tree
(
 p4est_t* p4est
)
{
  int num_trees = p4est->connectivity->num_trees;
  int* elements_per_tree_local = P4EST_ALLOC_ZERO(int, num_trees);
  int* elements_per_tree_global = P4EST_ALLOC_ZERO(int, num_trees);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      elements_per_tree_local[tt] = Q;
    }

    sc_reduce
      (
       &elements_per_tree_local[0],
       &elements_per_tree_global[0],
       num_trees,
       sc_MPI_INT,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    if (p4est->mpirank == 0){
      for (int i = 0; i < num_trees; i++){
        printf(" Tree %d: Number of Elements = %d\n", i, elements_per_tree_global[i]);
      }
    }
  
  P4EST_FREE(elements_per_tree_local);
  P4EST_FREE(elements_per_tree_global);
}

d4est_element_data_t*
d4est_mesh_get_element_data
(
 p4est_t* p4est,
 int local_element_id
)
{
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
        if (ed->id == local_element_id)
          return ed;
      }
    }
  return NULL;
}

int
d4est_mesh_global_node_to_local_node
(
 p4est_t* p4est,
 int global_node
)
{
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        if(global_node >= ed->nodal_stride && global_node < (ed->nodal_stride + volume_nodes))
          return global_node - ed->nodal_stride;
      }
    }
  return -1;
}

/* only for serial use */
int
d4est_mesh_debug_find_node
(
 p4est_t* p4est,
 int node
)
{
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        if(node >= ed->nodal_stride && node < (ed->nodal_stride + volume_nodes))
          return ed->id;
      }
    }
  return -1;
}

void
d4est_mesh_print_element_data_debug
(
 p4est_t* p4est
)
{
/* #ifndef D4EST_DEBUG */
/*   D4EST_ABORT("compile with the debug flag if you want to print curved element data"); */
/* #endif */

  
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

        printf("** Element %d **\n", ed->id);
        printf("deg, deg_quad = %d, %d\n", ed->deg, ed->deg_quad);
        printf("tree, tree_quadid = %d, %d\n", ed->tree, ed->tree_quadid);

        
        int volume_nodes = d4est_lgl_get_nodes( (P4EST_DIM), ed->deg );
       
#ifndef NDEBUG
#if (P4EST_DIM)==2
        printf("q = %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->dq);
        DEBUG_PRINT_2ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes);
        /* DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes); */
#elif (P4EST_DIM)==3
        printf("q = %d, %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->q[2], ed->dq);
        DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], ed->xyz[2], volume_nodes);
#else
        D4EST_ABORT("DIM = 2 or 3");
#endif
#endif
/* #else */
        /* D4EST_ABORT("DEBUG flag must be set"); */
/* #endif */
      }
    }
}

double
d4est_mesh_compute_l2_norm_sqr
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* nodal_vec,
 int local_nodes,
 norm_storage_option_t store_local,
 int (*skip_element_fcn)(d4est_element_data_t*),
 double* l2_array
)
{
  int k = 0;
  double* Mvec = P4EST_ALLOC(double, local_nodes);
  double l2_norm_sqr = 0.;
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        int skip_element = (skip_element_fcn != NULL) ? skip_element_fcn(ed) : 0;
        
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif

        d4est_quadrature_apply_mass_matrix
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &nodal_vec[ed->nodal_stride],
           ed->deg,
           ed->J_quad,
           ed->deg_quad,
           &Mvec[ed->nodal_stride]
          );
      
        double norm2
          = d4est_linalg_vec_dot(&nodal_vec[ed->nodal_stride],
                                 &Mvec[ed->nodal_stride],
                                 volume_nodes);
        
        if (store_local == STORE_LOCALLY){
          ed->local_estimator = norm2;
        }

        if (!skip_element)
          l2_norm_sqr += norm2;

        if(l2_array != NULL)
          l2_array[k] = norm2;
        
        k++;
      }
    }
  P4EST_FREE(Mvec);
  return l2_norm_sqr;
}


d4est_local_sizes_t
d4est_mesh_init_element_data
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void(*user_fcn)(d4est_element_data_t*, void*),
 void* user_ctx
)
{
  /* sizes */
  int local_nodes = 0;
  int local_sqr_nodes = 0;
  int local_nodes_quad = 0;
  /* int local_sqr_nodes_invM = 0; */

  /* strides */
  int sqr_nodal_stride = 0;
  int sqr_mortar_stride = 0;
  int nodal_stride = 0;
  int quad_stride = 0;
  int id_stride = 0;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);
        elem_data->tree = tt;
        elem_data->tree_quadid = q;
        elem_data->dq = P4EST_QUADRANT_LEN(quad->level);
        elem_data->q[0] = quad->x;
        elem_data->q[1] = quad->y;
#if (P4EST_DIM)==3
        elem_data->q[2] = quad->z;
#endif


        
        elem_data->id = id_stride;
        elem_data->sqr_nodal_stride = sqr_nodal_stride;
        /* elem_data->sqr_mortar_stride = sqr_mortar_stride; */
        elem_data->nodal_stride = nodal_stride;
        elem_data->quad_stride = quad_stride;

        elem_data->region = d4est_geom->get_region(d4est_geom, elem_data->q, elem_data->dq, elem_data->tree);

        
        /* user_fcn should set degree,
           or the degree will be assumed to be set */
        if (user_fcn != NULL){
          /* problem_set_degrees_donald_trump(elem_data, user_ctx); */
          user_fcn(elem_data, user_ctx);
        }

        /* printf("[UPDATE]: deg, deg_quad, region = %d, %d, %d\n", elem_data->deg, elem_data->deg_quad, elem_data->region); */

        /* printf("elem_data->deg = %d\n", elem_data->deg); */
        /* printf("elem_data->deg_quad = %d\n", elem_data->deg_quad); */
        
        /* elem_data->deg = 2; */
        /* elem_data->deg_quad = 2; */

        D4EST_ASSERT(elem_data->deg > 0
                   &&
                   elem_data->deg_quad > 0
                   &&
                   elem_data->deg < MAX_DEGREE
                  );
        
        int nodes = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg);
        int nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg_quad);
        int face_nodes = d4est_lgl_get_nodes((P4EST_DIM)-1, elem_data->deg);
        local_nodes += nodes;
        local_sqr_nodes += nodes*nodes;
        local_nodes_quad += d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg_quad);
        
        sqr_nodal_stride += nodes*nodes;
        sqr_mortar_stride += face_nodes*face_nodes*(P4EST_FACES);
        nodal_stride += nodes;
        quad_stride += nodes_quad;
        id_stride += 1;
      }
    }

  d4est_local_sizes_t local_sizes;
  local_sizes.local_nodes = local_nodes;
  local_sizes.local_sqr_nodes = local_sqr_nodes;
  local_sizes.local_nodes_quad = local_nodes_quad;
  local_sizes.local_mortar_nodes_quad = 0.;
  local_sizes.local_boundary_nodes_quad = 0.;

  if (ghost != NULL)
    d4est_mesh_compute_mortar_quadrature_sizes
      (
       p4est,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &local_sizes
      );
  
  /* local_sizes.local_sqr_nodes_invM = local_sqr_nodes_invM; */
  return local_sizes;
}

static void
d4est_mesh_geometry_storage_initialize_data
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
 
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
      
        d4est_rst_t rst_points_quad;
        rst_points_quad = d4est_quadrature_get_rst_points
                          (
                           d4est_ops,
                           d4est_quad,
                           d4est_geom,
                           &mesh_object,
                           QUAD_OBJECT_VOLUME,
                           QUAD_INTEGRAND_UNKNOWN,
                           ed->deg_quad
                          );


        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 2);
#endif
              
        int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           tt,
           ed->deg,
           ed->q,
           ed->dq,
           ed->xyz
          );

        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_quad,
           tt,
           ed->deg_quad,
           ed->q,
           ed->dq,
           ed->xyz_quad
          );

        
        if (d4est_geom->DX_compute_method == GEOM_COMPUTE_NUMERICAL){
            double* tmp = P4EST_ALLOC(double, volume_nodes);
            for (int d = 0; d < (P4EST_DIM); d++){
              for (int d1 = 0; d1 < (P4EST_DIM); d1++){
                d4est_operators_apply_dij(d4est_ops, &ed->xyz[d][0], (P4EST_DIM), ed->deg, d1, tmp);
                d4est_quadrature_interpolate(d4est_ops, d4est_quad, d4est_geom, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, tmp, ed->deg, &ed->xyz_rst_quad[d][d1][0], ed->deg_quad);
                /* double *dxyz_print = &ed->xyz_rst_quad[d][d1][0]; */
                /* DEBUG_PRINT_ARR_DBL(dxyz_print, volume_nodes_quad); */
              }
            }
            P4EST_FREE(tmp);
        }
        else if (d4est_geom->DX_compute_method == GEOM_COMPUTE_ANALYTIC){
          d4est_geometry_compute_dxyz_drst_analytic
            (
             d4est_ops,
             d4est_geom,
             rst_points_quad,
             ed->tree,
             ed->q,
             ed->dq,
             ed->deg_quad,
             ed->xyz_rst_quad
            );
        }
        else {
          D4EST_ABORT("Not a supported compute method for DX");
        }

        /* for (int d = 0; d < (P4EST_DIM); d++){ */
        /*   for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
        /*     double *dxyz_print = &ed->xyz_rst_quad[d][d1][0]; */
        /*     DEBUG_PRINT_ARR_DBL(dxyz_print, volume_nodes_quad); */
        /*   } */
        /* } */
        
    
        d4est_geometry_compute_jacobian
          (
           ed->xyz_rst_quad,
           ed->J_quad,
           volume_nodes_quad
          );

        d4est_geometry_compute_drst_dxyz
          (
           ed->xyz_rst_quad,
           ed->J_quad,
           ed->rst_xyz_quad,
           volume_nodes_quad
          );
      }
    }
}



void
d4est_mesh_geometry_storage_initialize_aliases
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_local_sizes_t local_sizes
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);

        elem_data->J_quad = &d4est_factors->J_quad[elem_data->quad_stride];
        for (int i = 0; i < (P4EST_DIM); i++){
          elem_data->xyz[i] = &d4est_factors->xyz[i*local_sizes.local_nodes + elem_data->nodal_stride];
          elem_data->xyz_quad[i] = &d4est_factors->xyz_quad[i*local_sizes.local_nodes_quad + elem_data->quad_stride];
          for (int j = 0; j < (P4EST_DIM); j++){
            elem_data->xyz_rst_quad[i][j] = &d4est_factors->xyz_rst_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride];
            elem_data->rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride];
          }
        }

      }
    }
}

int
d4est_mesh_update
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_quadrature_data_init_option_t quad_init_option,
 d4est_mesh_geometry_data_init_option_t geom_init_option,
 d4est_mesh_geometry_aliases_init_option_t alias_init_option,
 void(*user_fcn)(d4est_element_data_t*, void*),
 void* user_ctx
)
{
  d4est_local_sizes_t local_sizes =
    d4est_mesh_init_element_data(p4est,
                                 ghost, ghost_data,
                                 d4est_ops,
                                 d4est_geom,
                                 d4est_quad,
                                 d4est_factors,
                                 user_fcn,//problem_set_degrees_donald_trump,
                                 user_ctx);
  
  if (quad_init_option == INITIALIZE_QUADRATURE_DATA)
    {
      d4est_quadrature_reinit(
                              p4est,
                              ghost,
                              ghost_data,
                              d4est_ops,
                              d4est_geom,
                              d4est_quad
      );
    }

    
  
  if (geom_init_option == INITIALIZE_GEOMETRY_DATA)
    {
      d4est_mesh_geometry_storage_realloc(
                                          p4est,
                                          d4est_factors,
                                          local_sizes
      );
      d4est_mesh_geometry_storage_initialize_aliases(p4est, d4est_factors, local_sizes);
      d4est_mesh_geometry_storage_initialize_data(
                                                  p4est,
                                                  d4est_ops,
                                                  d4est_geom,
                                                  d4est_quad,
                                                  d4est_factors
      );
      d4est_mesh_compute_mortar_quadrature_quantities
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         d4est_factors
        );


      
    }

  if (alias_init_option == INITIALIZE_GEOMETRY_ALIASES &&
      geom_init_option != INITIALIZE_GEOMETRY_DATA /* if this is false, then this will be a waste of time */
     ){
    d4est_mesh_geometry_storage_initialize_aliases(p4est, d4est_factors, local_sizes);
  }


  if (alias_init_option == INITIALIZE_GEOMETRY_ALIASES ||
      geom_init_option == INITIALIZE_GEOMETRY_DATA /* if this is false, then this will be a waste of time */
     ){

    d4est_mesh_init_element_size_parameters
      (
       p4est,
       d4est_ops,
       d4est_geom
      );

  }
  
  return local_sizes.local_nodes;
}

void
d4est_mesh_init_field
(
 p4est_t* p4est,
 double* node_vec,
 d4est_xyz_fcn_t init_fcn,
 d4est_operators_t* d4est_ops, // TODO: unused, remove?
 d4est_geometry_t* d4est_geom, // TODO: unused, remove?
 d4est_mesh_init_field_option_t option,
 void* user
)
{
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

        if (option == INIT_FIELD_ON_LOBATTO){
          int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
          for (int i = 0; i < volume_nodes; i++){
            node_vec[ed->nodal_stride + i] = init_fcn(ed->xyz[0][i],
                                                      ed->xyz[1][i],
#if (P4EST_DIM)==3
                                                      ed->xyz[2][i],
#endif
                                                      user
                                                     );
          }
        }
        else if (option == INIT_FIELD_ON_QUAD){
          int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
          for (int i = 0; i < volume_nodes_quad; i++){
            node_vec[ed->quad_stride + i] = init_fcn(ed->xyz_quad[0][i],
                                                      ed->xyz_quad[1][i],
#if (P4EST_DIM)==3
                                                      ed->xyz_quad[2][i],
#endif
                                                      user
                                                     );
          }
        }
        else {
          D4EST_ABORT("Not a supported init option");
        }
      }
    }
  
}

void
d4est_mesh_compute_point_error
(
 double* v1,
 double* v2,
 double* error,
 int local_nodes
)
{
  for (int i = 0; i < local_nodes; i++){
    error[i] = fabs(v1[i] - v2[i]);
  }
}
 

void
d4est_mesh_init_field_ext
(
 p4est_t* p4est,
 double* node_vec,
 d4est_xyz_fcn_ext_t xyz_fcn,
 void* user,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{

  double* xyz_temp [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    xyz_temp[d] = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), (MAX_DEGREE)));
  }
  

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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);


          d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
      
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg,
                                                                0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL, NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 2);
#endif

        
        for (int i = 0; i < volume_nodes; i++){
          node_vec[ed->nodal_stride + i] = xyz_fcn(xyz_temp[0][i],
                                                    xyz_temp[1][i],
#if (P4EST_DIM)==3
                                                    xyz_temp[2][i],
#endif
                                                    user,
                                               d4est_geom,
                                               ed
                                              );

        }
      }
    }


  for (int d = 0; d < (P4EST_DIM); d++) {
    P4EST_FREE(xyz_temp[d]);
  }
  
}

void
d4est_mesh_get_local_nodes_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) q->p.user_data;
  int* local_nodes = (int*) user_data;
  *local_nodes += d4est_lgl_get_nodes( (P4EST_DIM),
                                    elem_data->deg);
}

int d4est_mesh_get_ghost_nodes
(
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data
)
{
  int ghost_nodes = 0;
  for (int gid = 0; gid < ghost->ghosts.elem_count; gid++){
    ghost_nodes += d4est_lgl_get_nodes((P4EST_DIM),ghost_data[gid].deg);
  }
  return ghost_nodes;
}

int d4est_mesh_get_local_nodes(p4est_t* p4est)
{
  int local_nodes = 0;
  
  p4est_iterate(p4est,
                NULL,
                (void *) (&local_nodes),
                d4est_mesh_get_local_nodes_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
                NULL);
  
  return local_nodes;
}


double
d4est_mesh_compare_two_fields
(
 p4est_t* p4est,
 double* field1,
 double* field2,
 const char* msg,
 d4est_mesh_boundary_option_t boundary_option,
 d4est_mesh_print_option_t print_option,
 double eps
)
{
  printf("%s -> \n", msg);
  int fail = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        for (int i = 0; i < volume_nodes; i++){
          int on_bndry = d4est_geometry_does_element_touch_boundary(p4est, quad, tt);

          if ((boundary_option == DISCARD_BOUNDARY && !on_bndry) ||
              (boundary_option == DISCARD_INTERIOR && on_bndry) ||
              (boundary_option == DISCARD_NOTHING))
            {
              int faili = fabs(field1[ed->nodal_stride + i] - field2[ed->nodal_stride + i]) > eps;
              fail += faili;

              if (print_option == PRINT || (print_option == PRINT_ON_ERROR && faili))
                {
                  printf("Element Info: tree %d, id %d, deg %d, on_bndry %d x %f y %f z %f field 1 %.25f field 2 %.25f error %.25f error_>_eps %d\n",
                         ed->tree,
                         ed->id,
                         ed->deg,
                         on_bndry,
                         ed->xyz[0][i],
                         ed->xyz[1][i],
                         ((P4EST_DIM)==3) ? ed->xyz[2][i] : 0.,
                         field1[ed->nodal_stride + i],
                         field2[ed->nodal_stride + i],
                         fabs(field1[ed->nodal_stride + i] - field2[ed->nodal_stride + i]),
                         faili
                        );
                }
            }
        }
      }
    }

  return (fail == 0);
}
 

double
d4est_mesh_surface_integral
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 int (*is_it_on_surface)(d4est_element_data_t*, int, void*),
 double (*compute_face_integral)(d4est_element_data_t*,
                                 int, double*,
                                 double* [(P4EST_DIM)],
                                 double* [(P4EST_DIM)],
                                 double* [(P4EST_DIM)][(P4EST_DIM)],
                                 void*),
 void* user
)
{

  double integral = 0.;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        int face_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, ed->deg_quad);


        double* normal_quad [(P4EST_DIM)];
        double* xyz_quad [(P4EST_DIM)];
        double* sj_quad = P4EST_ALLOC(double, face_nodes_quad);
        double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];

        D4EST_ALLOC_DIM_VEC(normal_quad, face_nodes_quad);
        D4EST_ALLOC_DIM_VEC(xyz_quad, face_nodes_quad);
        D4EST_ALLOC_DBYD_MAT(drst_dxyz_quad, face_nodes_quad);
        
        for (int f = 0; f < (P4EST_FACES); f++){
          int on_surface = is_it_on_surface(ed, f, user);

          
          d4est_mortars_compute_geometric_data_on_mortar
            (
             d4est_ops,
             d4est_geom,
             d4est_quad,
             QUAD_INTEGRAND_UNKNOWN,
             ed->tree,
             ed->q,
             ed->dq,
             ed->id*(P4EST_FACES) + f,
             1,
             1,
             &ed->deg,
             &ed->deg_quad,
             f,
             xyz_quad,
             drst_dxyz_quad,
             sj_quad,
             normal_quad,
             NULL,
             NULL,
             COMPUTE_NORMAL_USING_JACOBIAN
            );
          
          if (on_surface){
            integral += compute_face_integral(ed, f, sj_quad, normal_quad, xyz_quad, drst_dxyz_quad, user);
          }
        }

        D4EST_FREE_DIM_VEC(xyz_quad);
        D4EST_FREE_DBYD_MAT(drst_dxyz_quad);
        D4EST_FREE_DIM_VEC(normal_quad);
        P4EST_FREE(sj_quad);
      }
    }
  return integral;
}



double
d4est_mesh_volume_integral
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 int (*is_it_in_volume)(d4est_element_data_t*, void*),
 double (*compute_volume_integral)(d4est_element_data_t*, void*),
 void* user
)
{
  double integral = 0.;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        int in_volume = is_it_in_volume(ed, user);
        if (in_volume){
          integral += compute_volume_integral(ed, user);
        }
      }
    }
  return integral;
}

d4est_mesh_interpolate_data_t
d4est_mesh_interpolate_at_tree_coord
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double abc [(P4EST_DIM)],
 int tree_id,
 double* f,
 int print
){

  d4est_mesh_interpolate_data_t data;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      if (tt == tree_id){
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;

      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        int check = 0;
        for (int d = 0; d < (P4EST_DIM); d++){
          double amin = ed->q[d];
          double amax = ed->q[d] + ed->dq;
          amin /= (double)P4EST_ROOT_LEN;
          amax /= (double)P4EST_ROOT_LEN;
          check += (abc[d] <= amax) && (abc[d] >= amin);
        }
        if (check == (P4EST_DIM)){
          double xyz[(P4EST_DIM)];
          double rst [(P4EST_DIM)];

          d4est_geom->X(d4est_geom, tree_id, ed->q, ed->dq, abc, COORDS_TREE_UNITCUBE, xyz);
          
          for (int d = 0; d < (P4EST_DIM); d++){
            double amin = ed->q[d];
            double amax = ed->q[d] + ed->dq;
            amin /= (double)P4EST_ROOT_LEN;
            amax /= (double)P4EST_ROOT_LEN;
            rst[d] = 2*(abc[d] - amin)/(amax - amin) - 1;
            data.rst[d] = rst[d];
            data.xyz[d] = xyz[d];
            data.q[d] = ed->q[d];
            data.abc[d] = abc[d];
          }
          
          data.dq = ed->dq;
          data.id = ed->id;
          data.nodal_stride = ed->nodal_stride;
          data.f_at_xyz = d4est_operators_interpolate(d4est_ops, &rst[0], &f[ed->nodal_stride], (P4EST_DIM), ed->deg);

          if (print){
            printf("abc, xyz, f = %f, %f, %f, %f,%f,%f, %.15f\n", abc[0], abc[1], (P4EST_DIM)==3 ? abc[2] : 0, xyz[0], xyz[1], (P4EST_DIM)==3 ? xyz[2] : 0, data.f_at_xyz);
          }
          data.err = 0;
          return data;
        }
      }
      }
    }
  data.err = 1;
  return data;
}


