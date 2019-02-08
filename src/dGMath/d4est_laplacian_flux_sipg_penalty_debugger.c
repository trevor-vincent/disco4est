#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_mortars.h>
#include <d4est_laplacian_flux.h>
#include <d4est_laplacian_flux_sipg_penalty_debugger.h>
#include <ini.h>

static void
d4est_laplacian_flux_sipg_penalty_debugger_boundary
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* flux_parameter_data
)
{ 
  d4est_laplacian_flux_sipg_penalty_debugger_t* penalty_data = flux_parameter_data;
  int deg_quad = e_m->deg_quad;
  int face_nodes_m_quad = d4est_lgl_get_nodes((P4EST_DIM) - 1, deg_quad);
  int mortar_quad_scalar_stride = d4est_factors->local_strides[e_m->id].mortar_quad_stride[f_m];
  double* sigma = P4EST_ALLOC(double, face_nodes_m_quad);
  double * restrict  h_quad = &d4est_factors->hm_mortar_quad[mortar_quad_scalar_stride];
  double multiplier = d4est_laplacian_flux_sipg_penalty_debugger_get_multiplier(penalty_data,e_m->region,e_m->region);
  
  for (int i = 0; i < face_nodes_m_quad; i++){
    sigma[i] = multiplier*penalty_data->sipg_penalty_fcn
               (
                e_m->deg,
                h_quad[i],
                e_m->deg,
                h_quad[i],
                penalty_data->sipg_penalty_prefactor
               );
  }

  double min = sigma[0];
  double max = sigma[0];
  double sum = 0.;
  
  for (int i = 0; i < face_nodes_m_quad; i++){  
    sum += sigma[i];
    max = (max > sigma[i]) ? max : sigma[i];
    min = (min < sigma[i]) ? min : sigma[i];
  }
  double mean = sum/(double)face_nodes_m_quad;
  penalty_data->average_mean_penalty_vtk[e_m->id] += mean/(double)(P4EST_FACES);
  penalty_data->average_max_penalty_vtk[e_m->id] += min/(double)(P4EST_FACES);
  penalty_data->average_min_penalty_vtk[e_m->id] += max/(double)(P4EST_FACES);
  penalty_data->mean_penalty_vtk_per_face[f_m][e_m->id] = mean;
  penalty_data->max_penalty_vtk_per_face[f_m][e_m->id] = min;
  penalty_data->min_penalty_vtk_per_face[f_m][e_m->id] = max;
  P4EST_FREE(sigma);
}

/* d4est_laplacian_with_opt_flux_sipg_params
/* d4est_laplacian_flux_sipg_penalty_debugger */

static void
d4est_laplacian_flux_sipg_penalty_debugger_interface
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
  d4est_laplacian_flux_sipg_penalty_debugger_t* penalty_data = params;  

  int deg_m_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int deg_p_quad [(P4EST_HALF)];
  int deg_p_lobatto [(P4EST_HALF)];
  int nodes_mortar_quad [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  for (int i = 0; i < faces_m; i++){
    deg_m_quad[i] = e_m[i]->deg_quad;
    deg_m_lobatto[i] = e_m[i]->deg;
  }
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_quad[i] = e_p_oriented[i]->deg_quad;
  }
  int total_nodes_mortar_quad = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = d4est_util_max_int
                             (
                              deg_m_quad[i],
                              deg_p_quad[j]
                             );
      nodes_mortar_quad[i+j]
        = d4est_lgl_get_nodes(
                              (P4EST_DIM) - 1,
                              deg_mortar_quad[i+j]
                             );
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
    }

  int mortar_quad_scalar_stride = d4est_factors->local_strides[e_m[0]->id].mortar_quad_stride[f_m];
  double* hm_mortar_quad = &d4est_factors->hm_mortar_quad[mortar_quad_scalar_stride];
  double* hp_mortar_quad = &d4est_factors->hp_mortar_quad[mortar_quad_scalar_stride];  

  double* sigma = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double multiplier = d4est_laplacian_flux_sipg_penalty_debugger_get_multiplier(penalty_data,e_m[0]->region,e_p[0]->region);
  
  int stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;      
      sigma[ks] = multiplier*penalty_data->sipg_penalty_fcn
                          (
                           (faces_m == faces_mortar) ? deg_m_lobatto[f] : deg_m_lobatto[0],
                           hm_mortar_quad[ks],
                           (faces_p == faces_mortar) ? deg_p_lobatto[f] : deg_p_lobatto[0],
                           hp_mortar_quad[ks],
                           penalty_data->sipg_penalty_prefactor
                          );
    }
    stride += nodes_mortar_quad[f];
  }

  double max [(P4EST_HALF)];
  double min [(P4EST_HALF)];
  double mean [(P4EST_HALF)];
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    max[f] = sigma[0];
    min[f] = sigma[0];
    mean[f] = 0;
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;      
      min[f] = (sigma[ks] < min[f]) ? sigma[ks] : min[f];
      max[f] = (sigma[ks] > max[f]) ? sigma[ks] : max[f];
      mean[f] += sigma[ks];
    }
    mean[f] /= (double)nodes_mortar_quad[f];
    stride += nodes_mortar_quad[f];
  }

  if (faces_m == faces_mortar){
    for (int f = 0; f < faces_m; f++){
      penalty_data->average_min_penalty_vtk[e_m[f]->id] += min[f]/(double)(P4EST_FACES);
      penalty_data->average_mean_penalty_vtk[e_m[f]->id] += mean[f]/(double)(P4EST_FACES);
      penalty_data->average_max_penalty_vtk[e_m[f]->id] += max[f]/(double)(P4EST_FACES);
      penalty_data->min_penalty_vtk_per_face[f_m][e_m[f]->id] = min[f];
      penalty_data->mean_penalty_vtk_per_face[f_m][e_m[f]->id] = mean[f];
      penalty_data->max_penalty_vtk_per_face[f_m][e_m[f]->id] = max[f];
    }
  }
  else {
    double min_over_mortar = min[0];
    double max_over_mortar = max[0];
    double mean_over_mortar = 0.;
    for (int f = 0; f < faces_mortar; f++){        
      mean_over_mortar += mean[f]/(double)(faces_mortar);
      min_over_mortar = (min[f] < min_over_mortar) ? min[f] : min_over_mortar;
      max_over_mortar = (max[f] > max_over_mortar) ? max[f] : max_over_mortar;
    }
    penalty_data->average_min_penalty_vtk[e_m[0]->id] += min_over_mortar/(double)(P4EST_FACES);
    penalty_data->average_mean_penalty_vtk[e_m[0]->id] += mean_over_mortar/(double)(P4EST_FACES);
    penalty_data->average_max_penalty_vtk[e_m[0]->id] += max_over_mortar/(double)(P4EST_FACES);
    penalty_data->min_penalty_vtk_per_face[f_m][e_m[0]->id] = min_over_mortar;
    penalty_data->mean_penalty_vtk_per_face[f_m][e_m[0]->id] = mean_over_mortar;
    penalty_data->max_penalty_vtk_per_face[f_m][e_m[0]->id] = max_over_mortar;
  }

  P4EST_FREE(sigma);
}
/* d4est_laplacian_flux_sipg_penalty_debugger_t* */

static
int d4est_laplacian_with_opt_flux_sipg_region_multiplier_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_flux_sipg_penalty_debugger_t* pconfig = (d4est_laplacian_flux_sipg_penalty_debugger_t*)user;

  int hit = 0;
  for (int i = 0; i < pconfig->sipg_number_of_regions; i++){
    char* mult_name;
    asprintf(&mult_name,"region%d_multiplier", i);

    if (d4est_util_match_couple(section,"flux",name,mult_name)) {
      D4EST_ASSERT(pconfig->sipg_region_multipliers[i] == -1);
      pconfig->sipg_region_multipliers[i] = atof(value);
      hit++;
    }
    free(mult_name);
  }
  if (hit)
    return 1;
  else
    return 0;
}

static
int d4est_laplacian_with_opt_flux_sipg_region_boundary_multiplier_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_flux_sipg_penalty_debugger_t* pconfig = (d4est_laplacian_flux_sipg_penalty_debugger_t*)user;

  int hit = 0;

  for (int i = 0; i < pconfig->sipg_number_of_region_boundaries; i++){
    char* mult_name;
    asprintf(&mult_name,"region_boundary%d_multiplier", i);

    if (d4est_util_match_couple(section,"flux",name,mult_name)) {
      D4EST_ASSERT(pconfig->sipg_region_boundary_multipliers[i] == -1);
      pconfig->sipg_region_boundary_multipliers[i] = atof(value);
      hit++;
    }
    free(mult_name);
  }
  
  if (hit)
    return 1;
  else
    return 0;
}


static
int d4est_laplacian_flux_sipg_penalty_debugger_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_laplacian_flux_sipg_penalty_debugger_t* pconfig = (d4est_laplacian_flux_sipg_penalty_debugger_t*)user;

  if (d4est_util_match_couple(section,"flux",name,"sipg_use_region_multiplier")) {
    D4EST_ASSERT(pconfig->sipg_use_region_multipliers == 0);
    pconfig->sipg_use_region_multipliers = atoi(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_use_region_boundary_multiplier")) {
    D4EST_ASSERT(pconfig->sipg_use_region_boundary_multipliers == 0);
    pconfig->sipg_use_region_boundary_multipliers = atoi(value);
  }
  else if (d4est_util_match_couple(section,"flux",name,"sipg_number_of_regions")) {
    D4EST_ASSERT(pconfig->sipg_number_of_regions == 0);
    pconfig->sipg_number_of_regions = atoi(value);
  } 
  else if (d4est_util_match_couple(section,"flux",name,"sipg_number_of_region_boundaries")) {
    D4EST_ASSERT(pconfig->sipg_number_of_region_boundaries == 0);
    pconfig->sipg_number_of_region_boundaries = atoi(value);
  }  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

double
d4est_laplacian_flux_sipg_penalty_debugger_get_multiplier
(
 d4est_laplacian_flux_sipg_penalty_debugger_t* input,
 int region_m,
 int region_p
){
  if (!input->sipg_use_region_boundary_multipliers &&
      !input->sipg_use_region_multipliers){
    return 1.0;
  }

  double multiplier = 1.0;
  if (region_m == region_p){
    if (input->sipg_use_region_multipliers){
      multiplier = input->sipg_region_multipliers[region_m];
    }
  }
  else {
    int smaller_region = (region_m < region_p) ? region_m : region_p;
    if (input->sipg_use_region_boundary_multipliers){
      multiplier = input->sipg_region_boundary_multipliers[smaller_region];
    }
    else if (input->sipg_use_region_multipliers){
      double mult_m = input->sipg_region_multipliers[region_m];
      double mult_p = input->sipg_region_multipliers[region_p];
      multiplier = .5*(mult_m + mult_p);
    }
    else {
      multiplier = 1.0;
    }
  }

  return multiplier;
}


void
d4est_laplacian_flux_sipg_penalty_debugger_input
(
 p4est_t* p4est,
 const char* input_file,
 d4est_laplacian_flux_sipg_penalty_debugger_t* input
)
{
  input->sipg_use_region_multipliers = 0;
  input->sipg_use_region_boundary_multipliers = 0;
  input->sipg_number_of_regions = -1;
  input->sipg_number_of_region_boundaries = -1;

  if (ini_parse(input_file, d4est_laplacian_flux_sipg_penalty_debugger_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  if(input->sipg_use_region_multipliers){
    
    D4EST_CHECK_INPUT("flux", input->sipg_number_of_regions, -1);

    if (input->sipg_number_of_regions > 10){
      D4EST_ABORT("sipg_number_of_regions is capped at 10 in the code for now");
    }

    for (int i = 0; i < 10; i++){
      input->sipg_region_multipliers[i] = -1;
    }
    
    if (ini_parse(
                  input_file,
                  d4est_laplacian_with_opt_flux_sipg_region_multiplier_input_handler,
                  input) < 0) {
      D4EST_ABORT("Can't load input file");
    }
    
    for (int i = 0; i < input->sipg_number_of_regions; i++){
      D4EST_CHECK_INPUT("flux", input->sipg_region_multipliers[i], -1);
    }
  }

  if(input->sipg_use_region_boundary_multipliers){
    
    D4EST_CHECK_INPUT("flux", input->sipg_number_of_region_boundaries, -1);

    if (input->sipg_number_of_region_boundaries < 10){
      D4EST_ABORT("sipg number_of_region_boundaries capped at 10 now");
    }
    
    for (int i = 0; i < 10; i++){
      input->sipg_region_boundary_multipliers[i] = -1;
    }
    
    if (ini_parse(
                  input_file,
                  d4est_laplacian_with_opt_flux_sipg_region_boundary_multiplier_input_handler,
                  input) < 0) {
      D4EST_ABORT("Can't load input file");
    }
    
    for (int i = 0; i < input->sipg_number_of_region_boundaries; i++){
      D4EST_CHECK_INPUT("flux", input->sipg_region_boundary_multipliers[i], -1);
    }
  }
  
}


d4est_laplacian_flux_sipg_penalty_debugger_t*
d4est_laplacian_flux_sipg_penalty_debugger_init
(
 p4est_t* p4est,
 const char* input_file,
 penalty_calc_t penalty_fcn,
 double penalty_prefactor
)
{
  d4est_laplacian_flux_sipg_penalty_debugger_t* debugger =
    P4EST_ALLOC_ZERO(d4est_laplacian_flux_sipg_penalty_debugger_t,1);
  debugger->sipg_penalty_fcn = penalty_fcn;
  debugger->sipg_penalty_prefactor = penalty_prefactor;

  d4est_laplacian_flux_sipg_penalty_debugger_input(p4est, input_file, debugger);

  
  debugger->average_min_penalty_vtk =
    P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);

  debugger->average_max_penalty_vtk =
    P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);

  debugger->average_mean_penalty_vtk =
    P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);

  for (int f = 0; f < (P4EST_FACES); f++){
    debugger->min_penalty_vtk_per_face[f] = P4EST_ALLOC_ZERO(double,p4est->local_num_quadrants);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    debugger->mean_penalty_vtk_per_face[f] = P4EST_ALLOC_ZERO(double,p4est->local_num_quadrants);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    debugger->max_penalty_vtk_per_face[f] = P4EST_ALLOC_ZERO(double,p4est->local_num_quadrants);
  }
  
  return debugger;
}


d4est_laplacian_flux_sipg_penalty_debugger_t*
d4est_laplacian_flux_sipg_penalty_debugger_get_vtk_data
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_laplacian_flux_sipg_penalty_debugger_t* debugger
)
{ 
  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.user_ctx = debugger;
  flux_fcns.flux_interface_fcn = d4est_laplacian_flux_sipg_penalty_debugger_interface;
  flux_fcns.flux_boundary_fcn = d4est_laplacian_flux_sipg_penalty_debugger_boundary;
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     d4est_ghost,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
    );
}

void
d4est_laplacian_flux_sipg_penalty_debugger_destroy
(
 d4est_laplacian_flux_sipg_penalty_debugger_t* debugger
)
{
  P4EST_FREE(debugger->average_mean_penalty_vtk);
  P4EST_FREE(debugger->average_min_penalty_vtk);
  P4EST_FREE(debugger->average_max_penalty_vtk);

  for (int f = 0; f < (P4EST_FACES); f++){
    P4EST_FREE(debugger->mean_penalty_vtk_per_face[f]);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    P4EST_FREE(debugger->min_penalty_vtk_per_face[f]);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    P4EST_FREE(debugger->max_penalty_vtk_per_face[f]);
  }
  
  P4EST_FREE(debugger);
}
