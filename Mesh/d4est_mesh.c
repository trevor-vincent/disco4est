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
  d4est__local_sizes_t local_sizes =
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



