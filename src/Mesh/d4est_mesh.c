#define _GNU_SOURCE
#include <stdio.h>
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

static void
d4est_mesh_data_get_string_from_face_h_type
(
 d4est_mesh_face_h_t h_calc,
 char* string
)
{
  if (h_calc == FACE_H_EQ_J_DIV_SJ_QUAD){
    strcpy(string, "FACE_H_EQ_J_DIV_SJ_QUAD");
  }
  else if (h_calc == FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO){
    strcpy(string,"FACE_H_EQ_J_DIV_SJ_MIN_lOBATTO");
  }
  else if (h_calc == FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO){
    strcpy(string,"FACE_H_EQ_J_DIV_SJ_MEAN_lOBATTO");
  }
  else if (h_calc == FACE_H_EQ_VOLUME_DIV_AREA){
    strcpy(string,"FACE_H_EQ_VOLUME_DIV_AREA");
  }
  else if (h_calc == FACE_H_EQ_VOLUME_DIV_AREA){
    strcpy(string,"FACE_H_EQ_FACE_DIAM");
  }
  else if (h_calc == FACE_H_EQ_TREE_H){
    strcpy(string,"FACE_H_EQ_TREE_H");
  }
  else {
    strcpy(string,"NOT_SET");
  }
}


static void
d4est_mesh_data_get_string_from_vol_h_type
(
 d4est_mesh_volume_h_t h_calc,
 char* string
)
{
  if (h_calc == VOL_H_EQ_DIAM){
    strcpy(string, "VOL_H_EQ_DIAM");
  }
  else if (h_calc == VOL_H_EQ_CUBE_APPROX){
    strcpy(string,"VOL_H_EQ_CUBE_APPROX");
  }
  else {
    strcpy(string,"NOT_SET");
  }
}




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
  else if (d4est_util_match_couple(section,"initial_mesh",name,"load_newton_checkpoint")) {
    D4EST_ASSERT(pconfig->load_newton_checkpoint == 0);
    pconfig->load_newton_checkpoint = atoi(value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"load_krylov_checkpoint")) {
    D4EST_ASSERT(pconfig->load_krylov_checkpoint == 0);
    pconfig->load_krylov_checkpoint = atoi(value);
  }  
  else if (d4est_util_match_couple(section,"initial_mesh",name,"checkpoint_prefix")) {
    D4EST_ASSERT(pconfig->checkpoint_prefix == NULL);
    asprintf(&pconfig->checkpoint_prefix,"%s", value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"newton_checkpoint_prefix")) {
    D4EST_ASSERT(pconfig->newton_checkpoint_prefix == NULL);
    asprintf(&pconfig->newton_checkpoint_prefix,"%s", value);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"krylov_checkpoint_prefix")) {
    D4EST_ASSERT(pconfig->krylov_checkpoint_prefix == NULL);
    asprintf(&pconfig->krylov_checkpoint_prefix,"%s", value);
  }  
  else if (d4est_util_match_couple(section,"mesh_parameters",name,"max_degree")) {
    D4EST_ASSERT(pconfig->max_degree == -1);
    pconfig->max_degree = atoi(value);
    D4EST_ASSERT(atoi(value) > 0);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"checkpoint_number")) {
    D4EST_ASSERT(pconfig->checkpoint_number == -1);
    pconfig->checkpoint_number = atoi(value);
    D4EST_ASSERT(atoi(value) > 0);
  }
  else if (d4est_util_match_couple(section,"initial_mesh",name,"initial_checkpoint_number")) {
    D4EST_ASSERT(pconfig->initial_checkpoint_number == -1);
    pconfig->initial_checkpoint_number = atoi(value);
    D4EST_ASSERT(atoi(value) > 0);
  }
 else if (d4est_util_match_couple(section,"initial_mesh",name,"checkpoint_type")) {  
    D4EST_ASSERT(pconfig->checkpoint_type == D4EST_CHKPT_NOT_SET);
    if(d4est_util_match(value, "D4EST_CHKPT_HISTORY_H5")){
      pconfig->checkpoint_type = D4EST_CHKPT_HISTORY_H5;
    }
    else if(d4est_util_match(value, "D4EST_CHKPT_P4EST_H5")){
      pconfig->checkpoint_type = D4EST_CHKPT_P4EST_H5;
    }
    else {
      printf("checkpoint_type = %s\n", value);
      D4EST_ABORT("checkpoint_type is not set to a supported value\n");
    }
  }
 else if (d4est_util_match_couple(section,"mesh_parameters",name,"face_h_type")) {
    D4EST_ASSERT(pconfig->face_h_type == FACE_H_EQ_NOT_SET);
    if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_QUAD")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_QUAD;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_VOLUME_DIV_AREA")){
      pconfig->face_h_type = FACE_H_EQ_VOLUME_DIV_AREA;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_FACE_DIAM")){
      pconfig->face_h_type = FACE_H_EQ_FACE_DIAM;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_TREE_H")){
      pconfig->face_h_type = FACE_H_EQ_TREE_H;
    }
    else {
      printf("face_h_type = %s\n", value);
      D4EST_ABORT("face_h_type is not set to a supported value\n");
    }
  }
 else if (d4est_util_match_couple(section,"mesh_parameters",name,"volume_h_type")) {
    D4EST_ASSERT(pconfig->volume_h_type == VOL_H_EQ_NOT_SET);
    if(d4est_util_match(value, "VOL_H_EQ_DIAM")){
      pconfig->volume_h_type = VOL_H_EQ_DIAM;
    }
    else if(d4est_util_match(value, "VOL_H_EQ_CUBE_APPROX")){
      pconfig->volume_h_type = VOL_H_EQ_CUBE_APPROX;
    }
    else {
      printf("volume_h_type = %s\n", value);
      D4EST_ABORT("volume_h_type is not set to a supported value\n");
    }
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
d4est_mesh_set_initial_extents_using_input_file_regions
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  d4est_mesh_initial_extents_t* input = user_ctx;
  int region = elem_data->region;
  elem_data->deg = input->deg[region];
  elem_data->deg_quad = input->deg_quad_inc[region] +
                            elem_data->deg;
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
  if (input->load_from_checkpoint == 0){
    elem_data->deg = input->deg[region];
  }
  else{
    if (input->checkpoint_type == D4EST_CHKPT_P4EST_H5){
      D4EST_ASSERT(input->checkpoint_deg_array != NULL);
      elem_data->deg = input->checkpoint_deg_array[elem_data->id];
    }
    else if (input->checkpoint_type == D4EST_CHKPT_HISTORY_H5){
      /* elem_data->deg already set */
      /* do nothing */
    }
    else {
      D4EST_ABORT("set_initial_extents ran into unsupported else branch");
    }
  }  
  elem_data->deg_quad = input->deg_quad_inc[region] +
                            elem_data->deg;
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
}

void
d4est_mesh_initial_extents_destroy
(
 d4est_mesh_initial_extents_t* initial_extents
){
  P4EST_FREE(initial_extents->deg);
  P4EST_FREE(initial_extents->deg_quad_inc);
  free(initial_extents->checkpoint_prefix);
  free(initial_extents->newton_checkpoint_prefix);
  free(initial_extents->krylov_checkpoint_prefix);
  P4EST_FREE(initial_extents);
}




d4est_mesh_initial_extents_t*
d4est_mesh_initial_extents_parse
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  zlog_category_t* c_default = zlog_get_category("d4est_mesh_initial_extents");
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
  d4est_mesh_initial_extents_t* initial_extents = P4EST_ALLOC(d4est_mesh_initial_extents_t, 1);
  initial_extents->min_quadrants = -2;
  initial_extents->min_level = -2;
  initial_extents->max_degree = -1;
  initial_extents->fill_uniform = -2;
  initial_extents->number_of_regions = d4est_geom->get_number_of_regions(d4est_geom);
  initial_extents->deg = P4EST_ALLOC(int, initial_extents->number_of_regions);
  initial_extents->deg_quad_inc = P4EST_ALLOC(int, initial_extents->number_of_regions);
  initial_extents->load_from_checkpoint = 0;
  initial_extents->load_newton_checkpoint = 0;
  initial_extents->load_krylov_checkpoint = 0;
  initial_extents->checkpoint_prefix = NULL;
  initial_extents->newton_checkpoint_prefix = NULL;
  initial_extents->krylov_checkpoint_prefix = NULL;
  initial_extents->checkpoint_number = -1;
  initial_extents->initial_checkpoint_number = -1;
  initial_extents->checkpoint_type = D4EST_CHKPT_NOT_SET;
  initial_extents->volume_h_type = VOL_H_EQ_NOT_SET;
  initial_extents->face_h_type = FACE_H_EQ_NOT_SET;
  
 for (int i = 0; i < initial_extents->number_of_regions; i++){
    initial_extents->deg[i] = -1;
    initial_extents->deg_quad_inc[i] = -1;
  }
 
  if (ini_parse(input_file, d4est_mesh_initial_extents_handler, initial_extents) < 0) {
    zlog_error(c_default,"[D4EST_ERROR]: pXest input_file = %s\n", input_file);
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
    D4EST_CHECK_INPUT("initial_mesh", initial_extents->checkpoint_type, D4EST_CHKPT_NOT_SET);
    if (initial_extents->checkpoint_type == D4EST_CHKPT_HISTORY_H5){
      D4EST_CHECK_INPUT("initial_mesh", initial_extents->checkpoint_number, -1);
      D4EST_CHECK_INPUT("initial_mesh", initial_extents->initial_checkpoint_number, -1);
    }
  }

  D4EST_CHECK_INPUT("mesh_parameters", initial_extents->volume_h_type, VOL_H_EQ_NOT_SET);
  D4EST_CHECK_INPUT("mesh_parameters", initial_extents->max_degree, -1);
  D4EST_CHECK_INPUT("mesh_parameters", initial_extents->face_h_type, FACE_H_EQ_NOT_SET);


  char face_h_string [50];
  char vol_h_string [50];
  d4est_mesh_data_get_string_from_face_h_type( initial_extents->face_h_type,&face_h_string[0]);
  d4est_mesh_data_get_string_from_vol_h_type(initial_extents->volume_h_type, &vol_h_string[0]);
  if(proc_rank == 0){
    zlog_debug(c_default,"face_h = %s", face_h_string);
    zlog_debug(c_default,"vol_h = %s", vol_h_string);
    zlog_debug(c_default,"max_degree = %d", initial_extents->max_degree);
  }
  int trees = d4est_geom->p4est_conn->num_trees;
  
  if(
     proc_size != 1
     &&
     trees*d4est_util_int_pow_int((P4EST_CHILDREN), initial_extents->min_level) < proc_size
     &&
     initial_extents->min_quadrants < proc_size
  ){
    if (proc_rank == 0){
      zlog_error(c_default,"proc_size = %d", proc_size);
      zlog_error(c_default,"min_level = %d", initial_extents->min_level);
      zlog_error(c_default,"elements = %d",  trees*d4est_util_int_pow_int((P4EST_CHILDREN), initial_extents->min_level));
      zlog_error(c_default,"Starting p4est with elements < processes");
      D4EST_ABORT("");
    }
  }

  /* printf("%s: d4est_laplacian_flux_sipg_params_h_calc = %s\n", printf_prefix, h_eq); */
  /* d4est_laplacian_flux_sipg_params_get_string_from_h_calc (input->sipg_flux_h,h_eq); */
  
  return initial_extents;
}


void
d4est_mesh_compute_mortar_quadrature_quantities_boundary_callback
(
 p4est_t* p4est,
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

  d4est_factors->element_touches_boundary[e_m->id] = 1;
  double* h_quad = &d4est_factors->hm_mortar_quad[scalar_stride];
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


  double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, total_mortar_nodes_quad);
  
  d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                         (
                                          d4est_factors,
                                          e_m
                                         );


  int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), e_m->deg);
  int face_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg);
  
  double* node_on_face_vol = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* node_on_face = P4EST_ALLOC(double, face_nodes_lobatto);
  int* node_touches_boundary = &d4est_factors->node_touches_boundary[e_m->nodal_stride];
  
  d4est_util_fill_array(node_on_face, 1., face_nodes_lobatto);

  d4est_operators_apply_lift(d4est_ops,
                             node_on_face,
                             (P4EST_DIM),
                             e_m->deg,
                             f_m,
                             node_on_face_vol);

  for (int i = 0; i < volume_nodes_lobatto; i++){
    if (node_on_face_vol[i] > 0.){
      node_touches_boundary[i] = 1;
    }
  }
  
  P4EST_FREE(node_on_face);
  P4EST_FREE(node_on_face_vol);
  
  for (int d = 0; d < (P4EST_DIM); d++){

    d4est_operators_apply_slicer(d4est_ops,
                                 md_on_e.xyz[d],
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
     j_div_sj_m_mortar_quad,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_mesh_calculate_mortar_h
    (
     p4est,
     &e_m,
     f_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     *(d4est_mesh_face_h_t*)params,
     j_div_sj_m_mortar_quad,
     h_quad,
     1,
     1,
     &total_mortar_nodes_quad,
     d4est_mesh_get_size_parameters(d4est_factors)
    );

  P4EST_FREE(j_div_sj_m_mortar_quad);
  
}

void
d4est_mesh_calculate_mortar_h
(
 p4est_t* p4est,
 d4est_element_data_t** elems_side,
 int face_side,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_face_h_t face_h_type,
 double* j_div_sj_quad_mortar,
 double* h_quad_mortar,
 int num_faces_mortar,
 int num_faces_side,
 int* nodes_mortar_quad,
 d4est_mesh_size_parameters_t size_params
)
{
  if (face_h_type == FACE_H_EQ_J_DIV_SJ_QUAD){
    /* D4EST_ABORT("H_EQ_J_DIV_SJ_QUAD is no longer supported"); */
    int stride = 0;
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = j_div_sj_quad_mortar[ks];
      }
      stride += nodes_mortar_quad[f];
    }
  }
  else if (face_h_type == FACE_H_EQ_TREE_H){
    int stride = 0;
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = elems_side[0]->dq/(double)P4EST_ROOT_LEN;
      }
      stride += nodes_mortar_quad[f];
    }
  }
  else if (face_h_type == FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO){

    /* when we try to compute h for the "plus" side and if this side exists on another processor, this will fail because it tries to receive
     * the j_div_sj_min from the local memory not the other processors memeory. This needs to be fixed before use, hence the ABORT */
    
    int stride = 0;
    double h [P4EST_HALF];

    for (int f = 0; f < num_faces_mortar; f++){
      int element_id = elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->id;
      int stride = (p4est->mpirank == elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->mpirank) ? 0 : p4est->local_num_quadrants;
      int element_index = element_id + stride;
      
      h[f] = size_params.j_div_sj_min[element_index*(P4EST_FACES) + face_side]; 
    } 
                 
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = h[f];
      }
      stride += nodes_mortar_quad[f];
    }
  }

  else if (face_h_type == FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO){
    int stride = 0;
    double h [P4EST_HALF];
    for (int f = 0; f < num_faces_mortar; f++){

      int element_id = elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->id;
      int stride = (p4est->mpirank == elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->mpirank) ? 0 : p4est->local_num_quadrants;
      int element_index = element_id + stride;
      h[f] = size_params.j_div_sj_mean[element_index*(P4EST_FACES) + face_side]; 
    }
                 
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = h[f];
      }
      stride += nodes_mortar_quad[f];
    }
  }
  
  else if (face_h_type == FACE_H_EQ_TOTAL_VOLUME_DIV_TOTAL_AREA){
    int stride = 0;
    double area_of_mortar = 0.;
    double volume_of_elements_touching_mortar = 0.;
    for (int f = 0; f < num_faces_side; f++){
      area_of_mortar += size_params.area[elems_side[f]->id*(P4EST_FACES) + face_side];
      volume_of_elements_touching_mortar += size_params.volume[elems_side[f]->id];
    }
    double h = volume_of_elements_touching_mortar/area_of_mortar;
                 
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = h;
      }
      stride += nodes_mortar_quad[f];
    }
  }

  else if (face_h_type == FACE_H_EQ_VOLUME_DIV_AREA){
    int stride = 0;
    double h [P4EST_HALF];    
    for (int f = 0; f < num_faces_mortar; f++){
      int element_id = elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->id;
      int stride = (p4est->mpirank == elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->mpirank) ? 0 : p4est->local_num_quadrants;
      int element_index = element_id + stride;

      double volume = size_params.volume[element_index];
      double area = size_params.area[element_index*(P4EST_FACES) + face_side];
      h[f] = volume/area;
    }
         
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = h[f];
      }
      stride += nodes_mortar_quad[f];
    }
  }  
  
  else if (face_h_type == FACE_H_EQ_FACE_DIAM){
    int stride = 0;
    double h [P4EST_HALF];
    for (int f = 0; f < num_faces_mortar; f++){
      h[f] = size_params.diam_face[elems_side[(num_faces_side == num_faces_mortar) ? f : 0]->id*(P4EST_FACES) + face_side]; 
    }
                 
    for (int f = 0; f < num_faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_quad[f]; k++){
        int ks = k + stride;
        h_quad_mortar[ks] = h[f];
      }
      stride += nodes_mortar_quad[f];
    }
  }
  
  else {
    D4EST_ABORT("this face_h_type is not supported");
  }
}

void
d4est_mesh_compute_mortar_quadrature_quantities_interface_callback
(
 p4est_t* p4est,
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

  double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, total_mortar_nodes_quad);
  double* j_div_sj_p_mortar_quad_porder = P4EST_ALLOC(double, total_mortar_nodes_quad); 
  double* j_div_sj_p_mortar_quad = P4EST_ALLOC(double, total_mortar_nodes_quad); 


  double* hm_mortar_quad = &d4est_factors->hm_mortar_quad[scalar_stride];
  double* hp_mortar_quad = &d4est_factors->hp_mortar_quad[scalar_stride];
  
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
     j_div_sj_m_mortar_quad,
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
     j_div_sj_p_mortar_quad_porder,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  
  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
      face_p = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
    }

    d4est_operators_reorient_face_data
        (
         d4est_ops,
         &j_div_sj_p_mortar_quad_porder[oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &j_div_sj_p_mortar_quad[face_mortar_stride]
        );
    
    face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
  }

  d4est_mesh_calculate_mortar_h
    (
     p4est,
     e_m,
     f_m,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     *(d4est_mesh_face_h_t*)params,
     j_div_sj_m_mortar_quad,
     hm_mortar_quad,
     faces_mortar,
     faces_m,
     mortar_nodes_quad,
     d4est_mesh_get_size_parameters(d4est_factors)
    );

  /* d4est_element_data_t* e_p_oriented [(P4EST_HALF)]; */
  /* d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented); */
  
  d4est_mesh_calculate_mortar_h
    (
     p4est,
     &e_p_oriented[0],
     f_p,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     *(d4est_mesh_face_h_t*)params,
     j_div_sj_p_mortar_quad,
     hp_mortar_quad,
     faces_mortar,
     faces_p,
     mortar_nodes_quad,
     d4est_mesh_get_size_parameters(d4est_factors)
    );

  P4EST_FREE(j_div_sj_p_mortar_quad);
  P4EST_FREE(j_div_sj_m_mortar_quad);
  P4EST_FREE(j_div_sj_p_mortar_quad_porder);  
}

void
d4est_mesh_compute_mortar_quadrature_quantities
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_face_h_t face_h_type
)
{
  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_mesh_compute_mortar_quadrature_quantities_interface_callback;
  flux_fcns.flux_boundary_fcn = d4est_mesh_compute_mortar_quadrature_quantities_boundary_callback;
  flux_fcns.user_ctx = (void*)&face_h_type;

  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     d4est_ghost,
     /* NULL, */
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
     //EXCHANGE_GHOST_DATA /* might not be needed because it's done for size computations */
    );

}


static void
d4est_mesh_compute_mortar_quadrature_sizes_boundary_callback
(
 p4est_t* p4est,
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
  d4est_mesh_local_sizes_t* sizes = params;  

  
  e_m->tree_that_touch_face [f_m][0] = e_m->tree;
  e_m->tree_quadid_that_touch_face [f_m][0] = e_m->tree_quadid;
  for (int i = 1; i < (P4EST_HALF); i++){
    e_m->tree_that_touch_face [f_m][i] = -1;
    e_m->tree_quadid_that_touch_face [f_m][i] = -1;
  }
   
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
 p4est_t* p4est,
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
  for (int m = 0; m < faces_m; m++){
    for (int i = 0; i < (P4EST_HALF); i++){
      e_m[m]->tree_that_touch_face [f_m][i] = -1;
      e_m[m]->tree_quadid_that_touch_face [f_m][i] = -1;
    }
    for (int i = 0; i < faces_p; i++){
      e_m[m]->tree_that_touch_face [f_m][i] = e_p[i]->tree;
      e_m[m]->tree_quadid_that_touch_face [f_m][i] = e_p[i]->tree_quadid;
    }
  }
  
  d4est_mesh_local_sizes_t* sizes = params;
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
 d4est_ghost_t* d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_local_sizes_t* local_sizes
){

  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_mesh_compute_mortar_quadrature_sizes_interface_callback;
  flux_fcns.flux_boundary_fcn = d4est_mesh_compute_mortar_quadrature_sizes_boundary_callback;
  flux_fcns.user_ctx = (void*)local_sizes;
  
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     d4est_ghost,
     /* NULL,/\* ghost_data, *\/ */
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
     /* EXCHANGE_GHOST_DATA */
    );
}

d4est_mesh_data_t*
d4est_mesh_data_init()
{
  d4est_mesh_data_t* d4est_factors = P4EST_ALLOC(d4est_mesh_data_t, 1);
  d4est_factors->J_quad = NULL;
  d4est_factors->xyz = NULL;
  d4est_factors->xyz_quad = NULL;
  d4est_factors->xyz_rst_quad = NULL;
  d4est_factors->rst_xyz_quad = NULL;

  d4est_factors->drst_dxyz_m_mortar_quad = NULL;
  d4est_factors->drst_dxyz_p_mortar_quad_porder = NULL;

  d4est_factors->element_touches_boundary = NULL;
  d4est_factors->element_data = NULL;
  d4est_factors->node_touches_boundary = NULL;
  d4est_factors->face_touches_element = NULL;
  d4est_factors->hm_mortar_quad = NULL;
  d4est_factors->hp_mortar_quad = NULL;

  d4est_factors->sj_m_mortar_quad = NULL;
  d4est_factors->n_m_mortar_quad = NULL;
  d4est_factors->xyz_m_mortar_quad = NULL;
  d4est_factors->xyz_m_mortar_lobatto = NULL;

  d4est_factors->j_div_sj_min = NULL;
  d4est_factors->j_div_sj_mean = NULL;
  d4est_factors->diam_face = NULL;
  d4est_factors->diam_volume = NULL;
  d4est_factors->area = NULL;
  d4est_factors->volume = NULL;
  
  return d4est_factors;
}

void
d4est_mesh_data_zero
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_local_sizes_t local_sizes
)
{
  int local_nodes = local_sizes.local_nodes;
  int local_nodes_quad = local_sizes.local_nodes_quad;
  int vector_nodes = local_nodes*(P4EST_DIM);
  int vector_nodes_quad = local_nodes_quad*(P4EST_DIM);
  int matrix_nodes_quad = local_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_matrix_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_vector_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM);
  int local_vector_boundary_nodes_quad = local_sizes.local_boundary_nodes_quad*(P4EST_DIM);

  d4est_util_zero_array(d4est_factors->xyz, vector_nodes);
  d4est_util_zero_array(d4est_factors->xyz_quad, vector_nodes_quad);
  d4est_util_zero_array(d4est_factors->J_quad, local_nodes_quad);
  d4est_util_zero_array(d4est_factors->xyz_rst_quad, matrix_nodes_quad);
  d4est_util_zero_array(d4est_factors->rst_xyz_quad, matrix_nodes_quad);   
  d4est_util_zero_array(d4est_factors->drst_dxyz_m_mortar_quad, local_matrix_mortar_nodes_quad);
  
  d4est_util_zero_array(d4est_factors->drst_dxyz_p_mortar_quad_porder, local_matrix_mortar_nodes_quad);   
  
  d4est_util_zero_array(d4est_factors->hm_mortar_quad, local_sizes.local_mortar_nodes_quad);  
  /* d4est_util_zero_array(d4est_factors->element_touches_boundary, p4est->local_num_quadrants);   */
  d4est_util_zero_array(d4est_factors->hp_mortar_quad, local_sizes.local_mortar_nodes_quad);   
  d4est_util_zero_array(d4est_factors->sj_m_mortar_quad, local_sizes.local_mortar_nodes_quad);   
  d4est_util_zero_array(d4est_factors->n_m_mortar_quad, local_vector_mortar_nodes_quad);
  d4est_util_zero_array(d4est_factors->xyz_m_mortar_quad, local_vector_boundary_nodes_quad);
  d4est_util_zero_array(d4est_factors->xyz_m_mortar_lobatto, local_vector_boundary_nodes_quad);
  
  d4est_util_zero_array(d4est_factors->j_div_sj_min, p4est->local_num_quadrants*(P4EST_FACES));
  d4est_util_zero_array(d4est_factors->j_div_sj_mean, p4est->local_num_quadrants*(P4EST_FACES));
  d4est_util_zero_array(d4est_factors->diam_face, p4est->local_num_quadrants*(P4EST_FACES));
  d4est_util_zero_array(d4est_factors->area, p4est->local_num_quadrants*(P4EST_FACES));
  d4est_util_zero_array(d4est_factors->volume,p4est->local_num_quadrants);
  d4est_util_zero_array(d4est_factors->diam_volume,p4est->local_num_quadrants);  
}


void
d4est_mesh_data_debug_print
(
 p4est_t* p4est,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_local_sizes_t local_sizes
)
{
  int local_nodes = local_sizes.local_nodes;
  int local_nodes_quad = local_sizes.local_nodes_quad;
  int vector_nodes = local_nodes*(P4EST_DIM);
  int vector_nodes_quad = local_nodes_quad*(P4EST_DIM);
  int matrix_nodes_quad = local_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_matrix_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_vector_mortar_nodes_quad = local_sizes.local_mortar_nodes_quad*(P4EST_DIM);
  int local_vector_boundary_nodes_quad = local_sizes.local_boundary_nodes_quad*(P4EST_DIM);

  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->xyz, vector_nodes);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->xyz_quad, vector_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->J_quad, local_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->xyz_rst_quad, matrix_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->rst_xyz_quad, matrix_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->drst_dxyz_m_mortar_quad, local_matrix_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->drst_dxyz_p_mortar_quad_porder, local_matrix_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->hm_mortar_quad, local_sizes.local_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->element_touches_boundary, p4est->local_num_quadrants);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->node_touches_boundary, local_nodes);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->face_touches_element, p4est->local_num_quadrants*(P4EST_FACES));
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->hp_mortar_quad, local_sizes.local_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->sj_m_mortar_quad, local_sizes.local_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->n_m_mortar_quad, local_vector_mortar_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->xyz_m_mortar_quad, local_vector_boundary_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->xyz_m_mortar_lobatto, local_vector_boundary_nodes_quad);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->j_div_sj_min, p4est->local_num_quadrants*(P4EST_FACES));
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->j_div_sj_mean, p4est->local_num_quadrants*(P4EST_FACES));
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->diam_face, p4est->local_num_quadrants*(P4EST_FACES));
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->area, p4est->local_num_quadrants*(P4EST_FACES));
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->volume,p4est->local_num_quadrants);
  DEBUG_PRINT_MPI_ARR_DBL_SUM(p4est->mpirank, d4est_factors->diam_volume,p4est->local_num_quadrants);  
}

void
d4est_mesh_data_realloc
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_local_sizes_t local_sizes
)
{
  d4est_factors->local_sizes = local_sizes;

  if (p4est->mpirank == 0){
    printf("[D4EST_INFO]: Reallocing storage for geometry data\n");
  }

  
  int local_nodes = local_sizes.local_nodes;
  int local_nodes_quad = local_sizes.local_nodes_quad;
  int ghost_elements = d4est_ghost->ghost->ghosts.elem_count;
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



  d4est_factors->element_touches_boundary = P4EST_REALLOC(d4est_factors->element_touches_boundary,
                                                      int,
                                                      p4est->local_num_quadrants);

  d4est_factors->element_data = P4EST_REALLOC(d4est_factors->element_data,
                                              d4est_element_data_t*,
                                              (p4est->local_num_quadrants + ghost_elements));

  
  d4est_factors->node_touches_boundary = P4EST_REALLOC(d4est_factors->node_touches_boundary,
                                                      int,
                                                      local_nodes);


  
  d4est_factors->face_touches_element = P4EST_REALLOC(d4est_factors->face_touches_element,
                                                      int,
                                                      p4est->local_num_quadrants*(P4EST_FACES));

  for (int i = 0; i < p4est->local_num_quadrants; i++){
    d4est_factors->element_touches_boundary[i] = 0;
  }

  for (int i = 0; i < local_nodes; i++){
    d4est_factors->node_touches_boundary[i] = 0;
  }
  
  
  d4est_factors->hm_mortar_quad = P4EST_REALLOC(d4est_factors->hm_mortar_quad,
                                                      double,
                                                      local_sizes.local_mortar_nodes_quad);

  d4est_factors->hp_mortar_quad = P4EST_REALLOC(d4est_factors->hp_mortar_quad,
                                                      double,
                                                      local_sizes.local_mortar_nodes_quad);
  
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



  
  d4est_factors->j_div_sj_min = P4EST_REALLOC(d4est_factors->j_div_sj_min,
                                              double,
                                              (p4est->local_num_quadrants + ghost_elements)*(P4EST_FACES));

  
  d4est_factors->j_div_sj_mean = P4EST_REALLOC(d4est_factors->j_div_sj_mean,
                                              double,
                                              (p4est->local_num_quadrants + ghost_elements)*(P4EST_FACES));


  d4est_factors->diam_face = P4EST_REALLOC(d4est_factors->diam_face,
                                           double,
                                           (p4est->local_num_quadrants + ghost_elements)*(P4EST_FACES));
  

  d4est_factors->area = P4EST_REALLOC(d4est_factors->area,
                                      double,
                                      (p4est->local_num_quadrants + ghost_elements)*(P4EST_FACES));

  d4est_factors->volume = P4EST_REALLOC(d4est_factors->volume,
                                        double,
                                        (p4est->local_num_quadrants + ghost_elements));

  d4est_factors->diam_volume = P4EST_REALLOC(d4est_factors->diam_volume,
                                             double,
                                             (p4est->local_num_quadrants + ghost_elements));
}


void
d4est_mesh_data_printout
(
 d4est_mesh_data_t* d4est_factors
)
{
  int local_matrix_mortar_nodes_quad = d4est_factors->local_sizes.local_mortar_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  int local_vector_mortar_nodes_quad = d4est_factors->local_sizes.local_mortar_nodes_quad*(P4EST_DIM);
  int local_vector_boundary_nodes_quad = d4est_factors->local_sizes.local_boundary_nodes_quad*(P4EST_DIM);
}


void
d4est_mesh_data_destroy
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
    P4EST_FREE(d4est_factors->hm_mortar_quad);
    P4EST_FREE(d4est_factors->element_touches_boundary);
    P4EST_FREE(d4est_factors->element_data);
    P4EST_FREE(d4est_factors->node_touches_boundary);
    P4EST_FREE(d4est_factors->face_touches_element);
    P4EST_FREE(d4est_factors->hp_mortar_quad);
    
    P4EST_FREE(d4est_factors->n_m_mortar_quad);
    P4EST_FREE(d4est_factors->xyz_m_mortar_quad);
    P4EST_FREE(d4est_factors->xyz_m_mortar_lobatto);
    P4EST_FREE(d4est_factors->sj_m_mortar_quad);
    P4EST_FREE(d4est_factors->diam_face);
    P4EST_FREE(d4est_factors->j_div_sj_min);
    P4EST_FREE(d4est_factors->j_div_sj_mean);
    P4EST_FREE(d4est_factors->diam_volume);
    P4EST_FREE(d4est_factors->area);
    P4EST_FREE(d4est_factors->volume);
    P4EST_FREE(d4est_factors);
  }
}


static void
d4est_mesh_init_element_size_parameters
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_mesh_data_t* d4est_factors,
 d4est_element_data_t* ed,
 d4est_mesh_volume_h_t volume_h_type,
 int stride
)
{
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t,1);
  d4est_quadrature_lobatto_new(d4est_quad, NULL, NULL);

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
  /* ed->diam_volume = diam_volume; */
  d4est_factors->diam_volume[ed->id + stride] = d4est_mesh_data_compute_volume_diam(
                                                                                    xyz_lobatto,
                                                                                    ed->deg,
                                                                                    volume_h_type
                                                                                   );
  
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

  /* ed->volume = d4est_quadrature_innerproduct */
  d4est_factors->volume[ed->id + stride]
    = d4est_quadrature_innerproduct
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

    d4est_factors->area[(ed->id+stride)*(P4EST_FACES) + f ]
      = d4est_quadrature_innerproduct
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
          
    d4est_factors->diam_face[(ed->id+stride)*(P4EST_FACES) + f] = diam_face;
    d4est_factors->j_div_sj_min[(ed->id+stride)*(P4EST_FACES) + f] = d4est_util_min_dbl_array(j_div_sj_on_f_lobatto, face_nodes_lobatto);
    d4est_factors->j_div_sj_mean[(ed->id+stride)*(P4EST_FACES) + f] = d4est_util_mean_dbl_array(j_div_sj_on_f_lobatto, face_nodes_lobatto); 
  }
     
  P4EST_FREE(ones);

  D4EST_FREE_DIM_VEC(xyz_lobatto);
  D4EST_FREE_DBYD_MAT(dxyz_drst_lobatto);
  D4EST_FREE_DIM_VEC(xyz_on_f_lobatto);
  D4EST_FREE(sj_on_f_lobatto);
  D4EST_FREE(j_div_sj_on_f_lobatto);
  D4EST_FREE(j_lobatto);


  P4EST_FREE(d4est_quad);
}

static void
d4est_mesh_init_element_size_parameters_ghost
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_volume_h_t volume_h_type
)
{
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
    ged->id = gid;
    d4est_factors->element_data[p4est->local_num_quadrants + gid] = ged;
    d4est_mesh_init_element_size_parameters
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       ged,
       volume_h_type,
       p4est->local_num_quadrants
      );    
  }
}

static void
d4est_mesh_init_element_size_parameters_local
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_volume_h_t volume_h_type
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree; tt++)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        d4est_mesh_init_element_size_parameters
          (
           p4est,
           d4est_ops,
           d4est_geom,
           d4est_factors,
           ed,
           volume_h_type,
           0
          );
      }
    }
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


int
d4est_mesh_get_local_degree_sum
(
 p4est_t* p4est
)
{
  int deg_sum = 0;
  int stride = 0;
  for
    (
     p4est_topidx_t tt = p4est->first_local_tree;
     tt <= p4est->last_local_tree;
     ++tt
    ){
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        deg_sum += ed->deg;
      }
    }
  return deg_sum;
}

int
d4est_mesh_get_local_max_degree
(
 p4est_t* p4est
)
{
  int max_deg = -1;
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
        if (ed->deg > max_deg){
          max_deg = ed->deg;
        };
      }
    }
  return max_deg;
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


/* void */
/* d4est_mesh_get_array_of_estimators */
/* ( */
/*  p4est_t* p4est, */
/*  double* eta2_array */
/* ) */
/* { */
/*   int stride = 0; */
/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*         d4est_element_data_t* ed = quad->p.user_data; */
/*         eta2_array[stride] = ed->local_estimator; */
/*         stride++; */
/*       } */
/*     } */
/* } */

void
d4est_mesh_compute_jacobian_on_lgl_grid
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double* jacobian_lgl,
 int max_deg
)
{
  int nodal_stride = 0;
  double* xyz_rst [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(xyz_rst, d4est_lgl_get_nodes((P4EST_DIM),max_deg));

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

double
d4est_mesh_compute_l2_norm_sqr
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double* nodal_vec,
 int local_nodes,
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


        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);

        
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
           J_quad,
           ed->deg_quad,
           &Mvec[ed->nodal_stride]
          );
      
        double norm2
          = d4est_linalg_vec_dot(&nodal_vec[ed->nodal_stride],
                                 &Mvec[ed->nodal_stride],
                                 volume_nodes);
        
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

d4est_mesh_size_parameters_t
d4est_mesh_get_size_parameters
(
 d4est_mesh_data_t* factors
)
{
  d4est_mesh_size_parameters_t size_params;
  size_params.j_div_sj_min = factors->j_div_sj_min;
  size_params.j_div_sj_mean = factors->j_div_sj_mean;
  size_params.diam_face = factors->diam_face;
  size_params.diam_volume = factors->diam_volume;
  size_params.area = factors->area;
  size_params.volume = factors->volume;
  return size_params;
}


d4est_mesh_local_sizes_t
d4est_mesh_init_element_data
(
 p4est_t* p4est,
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

        elem_data->mpirank = p4est->mpirank;
        elem_data->id = id_stride;
        elem_data->sqr_nodal_stride = sqr_nodal_stride;
        /* elem_data->sqr_mortar_stride = sqr_mortar_stride; */
        elem_data->nodal_stride = nodal_stride;
        elem_data->quad_stride = quad_stride;
        elem_data->region = d4est_geom->get_region(d4est_geom, elem_data->q, elem_data->dq, elem_data->tree);
        
        /* user_fcn should set degree,
           or the degree will be assumed to be set */
        if (user_fcn != NULL){
          user_fcn(elem_data, user_ctx);
        }

        
        D4EST_ASSERT(elem_data->deg > 0
                     &&
                     elem_data->deg_quad > 0
                  );

        if (elem_data->deg > d4est_ops->max_degree ||
            elem_data->deg_quad > d4est_ops->max_degree){
          printf("Element %d on processor %d has a degree that is too big\n", elem_data->id, p4est->mpirank);
          printf("Element deg, deg_quad, max_degree = %d, %d, %d\n",elem_data->deg, elem_data->deg_quad, d4est_ops->max_degree);
          D4EST_ABORT("Aborting...");
        }

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

  d4est_mesh_local_sizes_t local_sizes;
  local_sizes.local_nodes = local_nodes;
  local_sizes.local_sqr_nodes = local_sqr_nodes;
  local_sizes.local_nodes_quad = local_nodes_quad;
  local_sizes.local_mortar_nodes_quad = 0.;
  local_sizes.local_boundary_nodes_quad = 0.;

  /* local_sizes.local_sqr_nodes_invM = local_sqr_nodes_invM; */
  return local_sizes;
}


void
d4est_mesh_data_get_topological_coords
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double* a,
 double* b,
 double* c //NULL if 2-D
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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 2);
#endif
        
        for (int i = 0; i < volume_nodes; i++){
          a[ed->nodal_stride + i] = d4est_reference_rtox(rst_points_lobatto.r[i], (double)ed->q[0], (double)ed->dq)/(double)(P4EST_ROOT_LEN);
          b[ed->nodal_stride + i] = d4est_reference_rtox(rst_points_lobatto.s[i], (double)ed->q[1], (double)ed->dq)/(double)(P4EST_ROOT_LEN);
#if (P4EST_DIM)==3
          c[ed->nodal_stride + i] = d4est_reference_rtox(rst_points_lobatto.t[i], (double)ed->q[2], (double)ed->dq)/(double)(P4EST_ROOT_LEN);
#endif
        }
        
      }
    }        
}

void
d4est_mesh_data_compute
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_face_h_t face_h_type,
 d4est_mesh_volume_h_t volume_h_type
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

        d4est_factors->element_data[ed->id] = ed;
        
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


        d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                               (
                                                d4est_factors,
                                                ed
                                               );
        
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           tt,
           ed->deg,
           ed->q,
           ed->dq,
           md_on_e.xyz
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
           md_on_e.xyz_quad
          );

        
        if (d4est_geom->DX_compute_method == GEOM_COMPUTE_NUMERICAL){
            double* tmp = P4EST_ALLOC(double, volume_nodes);
            for (int d = 0; d < (P4EST_DIM); d++){
              for (int d1 = 0; d1 < (P4EST_DIM); d1++){
                d4est_operators_apply_dij(d4est_ops, &(md_on_e.xyz[d][0]), (P4EST_DIM), ed->deg, d1, tmp);
                d4est_quadrature_interpolate(d4est_ops, d4est_quad, d4est_geom, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, tmp, ed->deg, &(md_on_e.xyz_rst_quad[d][d1][0]), ed->deg_quad);
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
             md_on_e.xyz_rst_quad
            );
        }
        else {
          D4EST_ABORT("Not a supported compute method for DX");
        }
        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);
        d4est_geometry_compute_jacobian
          (
           md_on_e.xyz_rst_quad,
           J_quad,
           volume_nodes_quad
          );

        d4est_geometry_compute_drst_dxyz
          (
           md_on_e.xyz_rst_quad,
           J_quad,
           md_on_e.rst_xyz_quad,
           volume_nodes_quad
          );
      }
    }
      


  d4est_mesh_init_element_size_parameters_local
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     volume_h_type
    );


  d4est_mesh_init_element_size_parameters_ghost
    (
     p4est,
     ghost,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     volume_h_type
    );

  d4est_mesh_compute_mortar_quadrature_quantities
    (
     p4est,
     ghost,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     face_h_type
    ); 
  
}



/* void */
/* d4est_mesh_data_initialize_aliases */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_mesh_data_t* d4est_factors, */
/*  d4est_mesh_local_sizes_t local_sizes */
/* ) */
/* { */
/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int QQ = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int qq = 0; qq < QQ; ++qq) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq); */
/*         d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data); */

/*         /\* elem_data->J_quad = &d4est_factors->J_quad[elem_data->quad_stride]; *\/ */
/*         for (int i = 0; i < (P4EST_DIM); i++){ */
/*           /\* elem_data->xyz[i] = &d4est_factors->xyz[i*local_sizes.local_nodes + elem_data->nodal_stride]; *\/ */
/*           /\* elem_data->xyz_quad[i] = &d4est_factors->xyz_quad[i*local_sizes.local_nodes_quad + elem_data->quad_stride]; *\/ */
/*           for (int j = 0; j < (P4EST_DIM); j++){ */
/*             /\* elem_data->xyz_rst_quad[i][j] = &d4est_factors->xyz_rst_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride]; *\/ */
/*             /\* elem_data->rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride]; *\/ */
/*           } */
/*         } */

/*       } */
/*     } */
/* } */


d4est_mesh_data_on_element_t
d4est_mesh_data_on_element
(
 d4est_mesh_data_t* d4est_factors,
 d4est_element_data_t* ed
)
{
  d4est_mesh_local_sizes_t local_sizes = d4est_factors->local_sizes;
  d4est_mesh_data_on_element_t mesh_data_on_e;
  mesh_data_on_e.J_quad = &d4est_factors->J_quad[ed->quad_stride];
  for (int i = 0; i < (P4EST_DIM); i++){
    mesh_data_on_e.xyz[i] = &d4est_factors->xyz[i*local_sizes.local_nodes + ed->nodal_stride];
    mesh_data_on_e.xyz_quad[i] = &d4est_factors->xyz_quad[i*local_sizes.local_nodes_quad + ed->quad_stride];
    for (int j = 0; j < (P4EST_DIM); j++){
      mesh_data_on_e.xyz_rst_quad[i][j] = &d4est_factors->xyz_rst_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + ed->quad_stride];
      mesh_data_on_e.rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + ed->quad_stride];
    }
  }
  return mesh_data_on_e;
}


double*
d4est_mesh_get_jacobian_on_quadrature_points
(
 d4est_mesh_data_t* d4est_factors,
 d4est_element_data_t* ed
)
{
  return &d4est_factors->J_quad[ed->quad_stride];
}


d4est_mesh_local_sizes_t
d4est_mesh_update
(
 p4est_t* p4est,
 d4est_ghost_t** d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 d4est_mesh_ghost_init_option_t ghost_init_option,
 d4est_mesh_quadrature_data_init_option_t quad_init_option,
 d4est_mesh_geometry_data_init_option_t geom_init_option,
 d4est_mesh_geometry_aliases_init_option_t alias_init_option,
 void(*user_fcn)(d4est_element_data_t*, void*),
 void* user_ctx
)
{
  d4est_mesh_local_sizes_t local_sizes =
    d4est_mesh_init_element_data(p4est,
                                 d4est_ops,
                                 d4est_geom,
                                 d4est_quad,
                                 d4est_factors,
                                 user_fcn,//problem_set_degrees_donald_trump,
                                 user_ctx);



 if (ghost_init_option == INITIALIZE_GHOST){
    if (d4est_ghost != NULL && *d4est_ghost != NULL) {
      d4est_ghost_destroy(*d4est_ghost);
    }
    *d4est_ghost = d4est_ghost_init(p4est);
  }
  

  if (d4est_ghost != NULL && *d4est_ghost != NULL){
    d4est_mesh_compute_mortar_quadrature_sizes
      (
       p4est,
       *d4est_ghost,
       /* NULL, */
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       &local_sizes
      );
  }
  
  /* WILL BE DEPRECATED SOON */
  if (quad_init_option == INITIALIZE_QUADRATURE_DATA)
    {
      d4est_quadrature_reinit(
                              p4est,
                              d4est_ops,
                              d4est_geom,
                              d4est_quad
      );
    }

  if (geom_init_option == INITIALIZE_GEOMETRY_DATA)
    {
      d4est_mesh_data_realloc
        (
         p4est,
         *d4est_ghost,
         d4est_factors,
         local_sizes
        );
      d4est_mesh_data_compute
        (
         p4est,
         *d4est_ghost,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         d4est_factors,
         initial_extents->face_h_type,
         initial_extents->volume_h_type
        );
    }


  return local_sizes;
}

void
d4est_mesh_init_field
(
 p4est_t* p4est,
 double* node_vec,
 d4est_xyz_fcn_t init_fcn,
 d4est_operators_t* d4est_ops, // TODO: unused, remove?
 d4est_geometry_t* d4est_geom, // TODO: unused, remove?
 d4est_mesh_data_t* d4est_factors,
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

        d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                               (
                                                d4est_factors,
                                                ed
                                               );
        if (option == INIT_FIELD_ON_LOBATTO){
          int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

          
          for (int i = 0; i < volume_nodes; i++){
            node_vec[ed->nodal_stride + i] = init_fcn(md_on_e.xyz[0][i],
                                                      md_on_e.xyz[1][i],
#if (P4EST_DIM)==3
                                                      md_on_e.xyz[2][i],
#endif
                                                      user
                                                     );
          }
        }
        else if (option == INIT_FIELD_ON_QUAD){
          int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);


          for (int i = 0; i < volume_nodes_quad; i++){
            node_vec[ed->quad_stride + i] = init_fcn(md_on_e.xyz_quad[0][i],
                                                      md_on_e.xyz_quad[1][i],
#if (P4EST_DIM)==3
                                                      md_on_e.xyz_quad[2][i],
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
 

/* void */
/* d4est_mesh_init_field_ext */
/* ( */
/*  p4est_t* p4est, */
/*  double* node_vec, */
/*  d4est_xyz_fcn_ext_t xyz_fcn, */
/*  void* user, */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_geometry_t* d4est_geom, */
/*  int max_degree */
/* ) */
/* { */

/*   double* xyz_temp [(P4EST_DIM)]; */
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     xyz_temp[d] = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), max_degree)); */
/*   } */
  

/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*         d4est_element_data_t* ed = quad->p.user_data; */
/*         int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg); */


/*           d4est_quadrature_volume_t mesh_object; */
/*         mesh_object.dq = ed->dq; */
/*         mesh_object.tree = ed->tree; */
/*         mesh_object.element_id = ed->id; */
/*         mesh_object.q[0] = ed->q[0]; */
/*         mesh_object.q[1] = ed->q[1]; */
/* #if (P4EST_DIM)==3 */
/*         mesh_object.q[2] = ed->q[2]; */
/* #endif */
      
/*         d4est_rst_t rst_points_lobatto; */
/*         rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, */
/*                                                                 NULL, */
/*                                                                 NULL, */
/*                                                                 &mesh_object, */
/*                                                                 QUAD_OBJECT_VOLUME, */
/*                                                                 QUAD_INTEGRAND_UNKNOWN, */
/*                                                                 ed->deg, */
/*                                                                 0); */
/*         rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, */
/*                                                                 NULL, NULL, */
/*                                                                 &mesh_object, */
/*                                                                 QUAD_OBJECT_VOLUME, */
/*                                                                 QUAD_INTEGRAND_UNKNOWN, */
/*                                                                 ed->deg, 1); */
/*         rst_points_lobatto.t = NULL; */
/* #if (P4EST_DIM)==3 */
/*         rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, */
/*                                                                 NULL, */
/*                                                                 NULL, */
/*                                                                 &mesh_object, */
/*                                                                 QUAD_OBJECT_VOLUME, */
/*                                                                 QUAD_INTEGRAND_UNKNOWN, */
/*                                                                 ed->deg, 2); */
/* #endif */

        
/*         for (int i = 0; i < volume_nodes; i++){ */
/*           node_vec[ed->nodal_stride + i] = xyz_fcn(xyz_temp[0][i], */
/*                                                     xyz_temp[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                                     xyz_temp[2][i], */
/* #endif */
/*                                                     user, */
/*                                                d4est_geom, */
/*                                                ed */
/*                                               ); */

/*         } */
/*       } */
/*     } */


/*   for (int d = 0; d < (P4EST_DIM); d++) { */
/*     P4EST_FREE(xyz_temp[d]); */
/*   } */
  
/* } */

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
 d4est_ghost_t* d4est_ghost
)
{
  int ghost_nodes = 0;
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    ghost_nodes += d4est_lgl_get_nodes((P4EST_DIM),d4est_ghost->ghost_elements[gid].deg);
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
 d4est_mesh_data_t* d4est_factors,
 double* field1,
 double* field2,
 const char* msg,
 d4est_mesh_boundary_option_t boundary_option,
 d4est_mesh_print_option_t print_option,
 double eps
)
{
  printf("%s\n", msg);
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

        d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                               (
                                                d4est_factors,
                                                ed
                                               );

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
                         md_on_e.xyz[0][i],
                         md_on_e.xyz[1][i],
                         ((P4EST_DIM)==3) ? md_on_e.xyz[2][i] : 0.,
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
            /* printf("abc, xyz, f = %f, %f, %f, %f,%f,%f, %.15f\n", abc[0], abc[1], (P4EST_DIM)==3 ? abc[2] : 0, xyz[0], xyz[1], (P4EST_DIM)==3 ? xyz[2] : 0, data.f_at_xyz); */
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



int 
d4est_mesh_is_it_a_ghost_element
(
 p4est_t* p4est,
 d4est_element_data_t* ed
){
  return (p4est->mpirank != ed->mpirank);
}

double*
d4est_mesh_get_field_on_element
(
 p4est_t* p4est,
 d4est_element_data_t* ed,
 d4est_ghost_data_t* d4est_ghost_data,
 double* field,
 int local_nodes,
 int which_field
)
{
  int is_it_ghost = d4est_mesh_is_it_a_ghost_element(p4est,ed);
  D4EST_ASSERT(which_field >= 0);
  if (is_it_ghost && d4est_ghost_data == NULL){
    D4EST_ABORT("need d4est_ghost_data != NULL if element is ghost");
  }
  if(is_it_ghost){
    return d4est_ghost_data_get_field_on_element(ed,which_field,d4est_ghost_data);
  }
  else{
    return &field[which_field*local_nodes + ed->nodal_stride];
  }
}



double
d4est_mesh_data_compute_volume_diam
(
 double* xyz [(P4EST_DIM)],
 int deg,
 d4est_mesh_volume_h_t option
)
{
  double diam = 0.;
    
  /* Use an approximate method to calculate diam: iterate through corners of element*/
  /* if (option == DIAM_APPROX || option == DIAM_APPROX_CUBE){ */
  /*   for (int i = 0; i < (P4EST_CHILDREN); i++){ */
  /*     for (int j = 0; j < (P4EST_CHILDREN); j++){ */
  /*       int corner_node_i = d4est_reference_corner_to_node((P4EST_DIM), deg, i); */
  /*       int corner_node_j = d4est_reference_corner_to_node((P4EST_DIM), deg, j); */
  /*       double diam_temp = 0; */
  /*       for (int d = 0; d < (P4EST_DIM); d++){ */
  /*         double diam_dx = xyz[d][corner_node_i] - xyz[d][corner_node_j]; */
  /*         diam_temp += diam_dx*diam_dx; */
  /*       } */
  /*       diam_temp = sqrt(diam_temp);       */
  /*       diam = (diam_temp > diam) ? diam_temp : diam; */
  /*     } */
  /*   } */

  /*   if (option == DIAM_APPROX_CUBE){ */
  /*     diam *= 1./sqrt(3.); */
  /*   } */
    
  /* } */
  /* else { */
    int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), deg);
    for (int i = 0; i < volume_nodes; i++){
      for (int j = 0; j < volume_nodes; j++){
        double diam_temp = 0;
        /* printf("\n"); */
        for (int d = 0; d < (P4EST_DIM); d++){
          double diam_dx = xyz[d][i] - xyz[d][j];
          /* printf("diam_dx = %f\n", diam_dx); */
          diam_temp += diam_dx*diam_dx;
        }
        /* printf("diam_temp = %f\n",diam_temp);        */
        diam_temp = sqrt(diam_temp);
        diam = (diam_temp > diam) ? diam_temp : diam;
      }
    }


    if (option == VOL_H_EQ_CUBE_APPROX){
      diam *= 1./sqrt((P4EST_DIM));
    }

  return diam;
}





void
d4est_mesh_debug_boundary_elements
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 const char** field_names,
 double** fields,
 int local_nodes
)
{

  int num_fields = -1;
  while (field_names[++num_fields] != NULL)
    continue;


  
  double** fields_modal = P4EST_ALLOC(double*, num_fields);
  int max_degree = d4est_mesh_get_local_max_degree(p4est);
  for (int i = 0; i < num_fields; i++){
    fields_modal[i] = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), max_degree));
  }
          
  
  
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
        if (d4est_factors->element_touches_boundary[ed->id] == 1){
          int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

          printf("\n**************************************\n");
          printf("**************************************\n");
          printf("Element id, tree-id, tree = %d, %d, %d\n", ed->id, qq, tt);          
          printf("**************************************\n");
          printf("**************************************\n");
          for (int fi = 0; fi < num_fields; fi++){
            /* printf("field %d = %s\n", fi, field_names[fi]); */

            double* field_nodal = &fields[fi][ed->nodal_stride];
            double* field_modal = fields_modal[fi];
            
            d4est_operators_convert_nodal_to_modal(d4est_ops,
                                                   field_nodal,
                                                   (P4EST_DIM),
                                                   ed->deg,
                                                   field_modal
                                                  );

            /* DEBUG_PRINT_ARR_DBL(field_nodal, volume_nodes); */
            /* DEBUG_PRINT_ARR_DBL(field_modal, volume_nodes); */
            
          }
          


          for (int i_fields = 0; i_fields < num_fields; i_fields++){
            printf("%s ",field_names[i_fields]);
          }

         for (int i_fields = 0; i_fields < num_fields; i_fields++){            
            printf("%s_modal ",field_names[i_fields]);
          }
          printf("\n");
          
          double* xyz [3];


/*           double* x = &d4est_factors->xyz[ed->nodal_stride]; */
/*           xyz[0] = x; */
/*           double* y = &d4est_factors->xyz[ed->nodal_stride + local_nodes]; */
/*           xyz[1] = y; */
/* #if (P4EST_DIM)==3 */
/*           double* z = &d4est_factors->xyz[ed->nodal_stride + 2*local_nodes]; */
/*           xyz[2] = z; */
/* #endif */

          xyz[0] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), ed->deg, 0);
          xyz[1] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), ed->deg, 1);
#if (P4EST_DIM)==3
          xyz[2] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), ed->deg, 2);
#endif

          

          

          for (int i = 0; i < volume_nodes; i++){
            for (int i_fields = 0; i_fields < (P4EST_DIM); i_fields++){
              printf("%.15f  ", xyz[i_fields][i]);
            }
            for (int i_fields=0; i_fields < num_fields; i_fields++) {
              printf("%.15f  ",fields[i_fields][ed->nodal_stride + i]);
            }
            for (int i_fields=0; i_fields < num_fields; i_fields++) {
              printf("%.15f  ",fields_modal[i_fields][i]);
            }
            printf("\n");
          }
        }
      }
    }

 P4EST_FREE(fields_modal);
 for (int i = 0; i < num_fields; i++){
   P4EST_FREE(fields_modal[i]);
 }
}



void
d4est_mesh_apply_invM_on_field
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_mesh_data_t* d4est_factors,
 double* in,
 double* out
){

  
  
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
        d4est_quadrature_apply_inverse_mass_matrix
          (
           d4est_ops,
           &in[ed->nodal_stride],
           ed->deg,
           &d4est_factors->J_quad[ed->quad_stride],
           ed->deg_quad,
           (P4EST_DIM),
           &out[ed->nodal_stride]
        );
      }

    }
  
}

