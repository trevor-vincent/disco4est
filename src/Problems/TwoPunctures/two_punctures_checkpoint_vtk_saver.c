#define _GNU_SOURCE
#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi_new.h>
#include <d4est_solver_cg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_brick.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_estimator_stats.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>
#include <d4est_solver_newton.h>
#include <d4est_solver_multigrid.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_solver_multigrid_matrix_operator.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_util.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <time.h>
#include "two_punctures_cactus_fcns_with_opt.h"
#include <p4est_vtk_ext.h>

int
keep_region_2_skipper
(
 d4est_element_data_t* ed
){
  if (ed->region != 2)
    return 1; /* skip */
  else
    return 0; /* do_not_skip */
}


int
keep_region_1_skipper
(
 d4est_element_data_t* ed
){
  if (ed->region != 1)
    return 1; /* skip */
  else
    return 0; /* do_not_skip */
}


int
keep_region_0_skipper
(
 d4est_element_data_t* ed
){
  if (ed->region != 0)
    return 1; /* skip */
  else
    return 0; /* do_not_skip */
}

static
double solve_for_c
(
 double c,
 void* user
)
{
  double* Rs = user;
  double R1 = Rs[0];
  double R2 = Rs[1];
  double Rc = Rs[2];
  double m = (2 - 1)/((1/R2) - (1/R1));
  double n = (1*R1 - 2*R2)/(R1 - R2);
  double R = m/(c - n);
  double pp = 2 - c;
  double q = R/sqrt(1 + 2*pp);
  double x = q;
  return x - Rc;  
}

static
double solve_for_c_outer
(
 double c,
 void* user
)
{
  double* Rs = user;
  double R1 = Rs[0];
  double R2 = Rs[1];
  double Rc = Rs[2];
  double m = (2 - 1)/((1/R2) - (1/R1));
  double n = (1*R1 - 2*R2)/(R1 - R2);
  double R = m/(c - n);
  /* double pp = 2 - c; */
  /* double q = R; */
  /* double x = R; */
  return R - Rc;  
}


double
get_inverted_outer_wedge_point(double R1, double R2, double Rc, int compactified){
  D4EST_ASSERT(Rc >= R1 && Rc <= R2);
  if (compactified){
    double c;
    if (Rc == R2){
      c = 2;
    }
    else {
      double Rs [] = {R1,R2,Rc};
      int success = d4est_util_bisection(solve_for_c_outer, 1, 2, DBL_EPSILON, 100000, &c, &Rs[0]);
      D4EST_ASSERT(!success);
    }
    return c - 1;
  }
  else{
    D4EST_ABORT("get_inverted_outer_wedge_point not accepted yet");
  }
}

double
get_inverted_inner_wedge_point(double R1, double R2, double Rc, int compactified){
  D4EST_ASSERT(Rc >= R1 && Rc <= R2);
  if (compactified){
    double c;
    if (Rc == R2){
      c = 2;
    }
    else {
      double Rs [] = {R1,R2,Rc};
      int success = d4est_util_bisection(solve_for_c, 1, 2, DBL_EPSILON, 100000, &c, &Rs[0]);
      D4EST_ASSERT(!success);
    }
    return c - 1;
  }
  else{
    return ((2*pow(R1,2) - 3*R1*R2 + pow(R2,2) - pow(Rc,2) + sqrt(pow(Rc,2)*(pow(R1,2) - 4*R1*R2 + 3*pow(R2,2) + pow(Rc,2))))/pow(R1 - R2,2)) - 1;
  }
}

double
get_inverted_box_point(double R0, double x){
  double a = R0/sqrt(3);
  D4EST_ASSERT(x <= a && x >= -a);
  return (x + a)/(2*a);
}


int
skip_curved_elements
(
 d4est_element_data_t* elem
)
{
  if (elem->tree == 6)
    return 0;
  else
    return 1;
}

typedef struct {
  
  int use_pointwise_estimator;
  int use_dirichlet;
  int interpolate_f;
  int amr_level_for_uniform_p;
  
} two_punctures_init_params_t;


static double
get_tree_coordinate(double R0, double R1, double R){
  double m = (2. - 1.)/((1./R1) - (1./R0));
  double t = (1.*R0 - 2.*R1)/(R0 - R1);
  return t + (m/R) - 1;
}

static
int two_punctures_init_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  two_punctures_init_params_t* pconfig = (two_punctures_init_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"use_pointwise_estimator")) {
    D4EST_ASSERT(pconfig->use_pointwise_estimator == -1);
    pconfig->use_pointwise_estimator = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"use_dirichlet")) {
    D4EST_ASSERT(pconfig->use_dirichlet == -1);
    pconfig->use_dirichlet = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"interpolate_f")) {
    D4EST_ASSERT(pconfig->interpolate_f == -1);
    pconfig->interpolate_f = atoi(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"amr_level_for_uniform_p")) {
    D4EST_ASSERT(pconfig->amr_level_for_uniform_p == -1);
    pconfig->amr_level_for_uniform_p = atoi(value);
  }
  else {
    return 0;
  }
  return 1;
}


static
two_punctures_init_params_t
two_punctures_init_params_input
(
 const char* input_file
)
{
  two_punctures_init_params_t input;
  input.use_dirichlet = -1;
  input.use_pointwise_estimator = -1;
  input.interpolate_f = -1;
  input.amr_level_for_uniform_p = 1000;
  
  if
    (
     ini_parse(input_file,
               two_punctures_init_params_handler,
               &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }
  
  D4EST_CHECK_INPUT("problem", input.use_dirichlet, -1);
  D4EST_CHECK_INPUT("problem", input.use_pointwise_estimator, -1);
  D4EST_CHECK_INPUT("problem", input.interpolate_f, -1);
  /* D4EST_CHECK_INPUT("amr", input.amr_level_for_uniform_p, -1); */
  
  return input;
}

static
int
amr_mark_element
(
 p4est_t* p4est,
 double eta2,
 d4est_estimator_stats_t* stats,
 d4est_element_data_t* elem_data,
 d4est_amr_smooth_pred_params_t* params,
 void* user
)
{
  problem_ctx_t* ctx = user;
  double eta2_percentile = stats->estimator_at_percentile;  
  return ((eta2 >= eta2_percentile)
          || fabs(eta2 - eta2_percentile)
          < eta2*1e-4);
}

static
gamma_params_t
amr_set_element_gamma
(
 p4est_t* p4est,
 d4est_estimator_stats_t* stats,
 d4est_element_data_t* elem_data,
 d4est_amr_smooth_pred_params_t* params,
 void* user
)
{
  problem_ctx_t* ctx = user;
  gamma_params_t gamma_hpn;
  gamma_hpn.gamma_h = params->gamma_h;
  gamma_hpn.gamma_p = params->gamma_p;
  gamma_hpn.gamma_n = params->gamma_n;
  return gamma_hpn;
}


void
problem_init
(
 p4est_t* p4est,
 d4est_ghost_t** d4est_ghost,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{ 
  int initial_nodes = initial_extents->initial_nodes;
  two_punctures_init_params_t init_params = two_punctures_init_params_input(input_file); 
  two_punctures_params_t two_punctures_params;
  init_two_punctures_data(&two_punctures_params);
  two_punctures_params.interpolate_f = init_params.interpolate_f;
  
  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_for_bi;
  bc_data_for_bi.dirichlet_fcn = zero_fcn;
  bc_data_for_bi.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_laplacian_with_opt_flux_data_t* flux_data_for_bi
    = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_bi);

  
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_jac = NULL;
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_res = NULL;

  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_dirichlet_for_jac;
  bc_data_dirichlet_for_jac.dirichlet_fcn = zero_fcn;
  bc_data_dirichlet_for_jac.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_dirichlet_for_res;
  bc_data_dirichlet_for_res.dirichlet_fcn = zero_fcn;
  bc_data_dirichlet_for_res.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_laplacian_with_opt_robin_bc_t bc_data_robin_for_jac;
  bc_data_robin_for_jac.robin_coeff = two_punctures_robin_coeff_sphere_fcn;
  bc_data_robin_for_jac.robin_rhs = two_punctures_robin_bc_rhs_fcn;

  d4est_laplacian_with_opt_robin_bc_t bc_data_robin_for_res;
  bc_data_robin_for_res.robin_coeff = two_punctures_robin_coeff_sphere_fcn;
  bc_data_robin_for_res.robin_rhs = two_punctures_robin_bc_rhs_fcn;  
  
  if(init_params.use_dirichlet){
    flux_data_for_jac
      = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_jac);
    flux_data_for_res
      = d4est_laplacian_with_opt_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_dirichlet_for_res);
  }
  else {  
    flux_data_for_jac = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_ROBIN, &bc_data_robin_for_jac);
    flux_data_for_res = d4est_laplacian_with_opt_flux_new(p4est, input_file,  BC_ROBIN, &bc_data_robin_for_res);
  }
  
  problem_ctx_t ctx;
  ctx.two_punctures_params = &two_punctures_params;
  ctx.flux_data_for_jac = flux_data_for_jac;
  ctx.flux_data_for_res = flux_data_for_res;
  
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = two_punctures_build_residual;
  prob_fcns.apply_lhs = two_punctures_apply_jac;
  prob_fcns.user = &ctx;
  
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);
  double* u_prev = P4EST_ALLOC(double, prob_vecs.local_nodes);
  
  d4est_laplacian_with_opt_flux_sipg_params_t* sipg_params = flux_data_for_jac->flux_data;
  
  d4est_estimator_bi_new_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_data.size_params = NULL;
  penalty_data.u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_data.gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  
  d4est_amr_smooth_pred_marker_t amr_marker;
  amr_marker.user = (void*)&ctx;
  amr_marker.mark_element_fcn = amr_mark_element;
  amr_marker.set_element_gamma_fcn = amr_set_element_gamma;

  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     &amr_marker
    );

  int initial_level = 0;
  if (initial_extents->load_from_checkpoint == 0 || initial_extents->checkpoint_prefix == NULL){
    d4est_mesh_init_field
      (
       p4est,
       prob_vecs.u,
       two_punctures_initial_guess,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );
  }
  else {


    d4est_checkpoint_read_dataset
      (
       p4est,
       initial_extents->checkpoint_prefix,
       "u",
       H5T_NATIVE_DOUBLE,
       prob_vecs.u,
       initial_extents->checkpoint_number
      );


    double sum = d4est_util_sum_array_dbl(prob_vecs.u, prob_vecs.local_nodes);

    d4est_checkpoint_check_dataset(p4est,
                                   initial_extents->checkpoint_prefix,
                                   "u",
                                   H5T_NATIVE_DOUBLE,
                                   (void*)&sum,
                                   initial_extents->checkpoint_number
                                  );

    initial_level = initial_extents->checkpoint_number + 1;

  }
  /* keep track of mg_levels on the previous amr step 
   * to calculate the mg_levels on the next amr step */
  int num_of_mg_levels_last_step = -1;
  
  // Norm function contexts  
  d4est_norms_fcn_L2_ctx_t L2_norm_ctx;
  L2_norm_ctx.p4est = p4est;
  L2_norm_ctx.d4est_ops = d4est_ops;
  L2_norm_ctx.d4est_geom = d4est_geom;
  L2_norm_ctx.d4est_quad = d4est_quad;
  L2_norm_ctx.d4est_factors = d4est_factors;
  
  d4est_norms_fcn_energy_ctx_t energy_norm_ctx;
  energy_norm_ctx.p4est = p4est;
  energy_norm_ctx.d4est_ops = d4est_ops;
  energy_norm_ctx.d4est_geom = d4est_geom;
  energy_norm_ctx.d4est_quad = d4est_quad;
  energy_norm_ctx.d4est_factors = d4est_factors;
  /* energy_norm_ctx.fit = NULL; */
  // These are updated later
  energy_norm_ctx.which_field = 0;
  energy_norm_ctx.energy_norm_data = NULL;
  energy_norm_ctx.energy_estimator_sq_local = -1.;

  /* if (p4est->mpirank == 0) */
  /*   d4est_norms_write_headers( */
  /*                             (const char * []){"u", NULL}, */
  /*                             (const char * []){"L_2", "L_infty", "energy_norm", "energy_estimator", NULL} */
  /*   ); */

  zlog_category_t* zlog_points = zlog_get_category("d4est_points");
  d4est_norms_linear_fit_t* point_3m_fit = d4est_norms_linear_fit_init(); 
  
  d4est_util_copy_1st_to_2nd(prob_vecs.u, u_prev, prob_vecs.local_nodes);

  double point [4][30];
  double point_diff [4][30];
  double point_spec_diff [4][30];
  double point_err [4];
  double point_dof [30];
  
  point[0][0] = 0;
  point_diff[0][0] = 0;
  point[1][0] = 0;
  point_diff[1][0] = 0;
  point[2][0] = 0;
  point_diff[2][0] = 0;
  point[3][0] = 0;
  point_diff[3][0] = 0;
  point_dof[0] = 0;
  point_spec_diff[0][0] = 0;
  point_spec_diff[1][0] = 0;
  point_spec_diff[2][0] = 0;
  point_spec_diff[3][0] = 0;

  
  int iterations = 1;

  zlog_category_t *c_geom = zlog_get_category("d4est_geometry_compactified");
  d4est_geometry_t* d4est_geom_compactified = d4est_geometry_new(p4est->mpirank, input_file,"compactified_geometry",c_geom);
  d4est_mesh_data_t* d4est_factors_compactified = d4est_mesh_data_init(p4est);

  d4est_norms_linear_fit_t* l2_linear_fit = d4est_norms_linear_fit_init();

  d4est_amr_t* d4est_amr_uniform_p = d4est_amr_init_uniform_p(p4est,d4est_amr->num_of_amr_steps);

  d4est_field_type_t field_type = NODAL;
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               *d4est_ghost,
                                                               &field_type,
                                                               1);
  // Extract mesh data  
  d4est_mesh_data_realloc
    (
     p4est,
     *d4est_ghost,
     d4est_factors_compactified,
     d4est_factors->local_sizes
    );
    
  d4est_mesh_data_compute
    (
     p4est,
     *d4est_ghost,
     d4est_ops,
     d4est_geom_compactified,
     d4est_quad,
     d4est_factors_compactified,
     initial_extents->face_h_type,
     initial_extents->volume_h_type
    );


  d4est_mesh_size_parameters_t size_params = d4est_mesh_get_size_parameters(d4est_factors_compactified);
  d4est_ip_energy_norm_data_t ip_norm_data;
  penalty_data.size_params = NULL;
  ip_norm_data.size_params = NULL;
  sipg_params->size_params = NULL;

  /* if (init_params.use_compactified_size_params){ */
  /* penalty_data.size_params = &size_params; */
  /* ip_norm_data.size_params = &size_params; */
  /* sipg_params->size_params = &size_params; */
  /* } */

  double* residual_pointwise_quad
    = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
    
  two_punctures_pointwise_residual
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     prob_vecs.u,
     residual_pointwise_quad,
     0,
     &ctx
    );
    
  double* estimator =
    d4est_estimator_bi_new_compute
    (
     p4est,
     &prob_vecs,
     residual_pointwise_quad,
     penalty_data,
     zero_fcn,
     NULL,
     *d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     d4est_geom_compactified,
     d4est_factors_compactified,
     d4est_quad,
     0,
     NULL,
     NULL,
     (init_params.use_pointwise_estimator) ? 1 : 0,
     &prob_fcns
    );

  P4EST_FREE(residual_pointwise_quad);

  d4est_amr_smooth_pred_params_t* sp_params = d4est_amr_smooth_pred_params_input
                                              (
                                               input_file
                                              );

  d4est_estimator_stats_t* stats = P4EST_ALLOC(d4est_estimator_stats_t,1);
  d4est_estimator_stats_compute(p4est, estimator, stats, sp_params->percentile, 1, 0);
  d4est_linalg_vec_fabsdiff(prob_vecs.u, u_prev, error, prob_vecs.local_nodes);
  double* error_l2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
  P4EST_FREE(sp_params);
    
  d4est_mesh_compute_l2_norm_sqr
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     error,
     prob_vecs.local_nodes,
     NULL,
     error_l2
    );

  /* if(init_params.use_error_l2_as_estimator){ */
  /* d4est_util_copy_1st_to_2nd(error_l2, estimator, p4est->local_num_quadrants); */
  /* }     */
    
  d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk",
     (const char * []){"u","u_prev","error", NULL},
     (double* []){prob_vecs.u, u_prev, error},
     (const char * []){"estimator","error_l2",NULL},
     (double* []){estimator,error_l2},
     NULL,
     NULL,
     0
    );

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk_compactified",
     (const char * []){"u","u_prev","error", NULL},
     (double* []){prob_vecs.u, u_prev, error},
     (const char * []){"estimator","error_l2",NULL},
     (double* []){estimator,error_l2},
     NULL,
     NULL,
     0
    );

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk_compactified_corner",
     (const char * []){"u","u_prev","error", NULL},
     (double* []){prob_vecs.u, u_prev, error},
     (const char * []){"estimator","error_l2",NULL},
     (double* []){estimator,error_l2},
     NULL,
     NULL,
     0
    );
      
  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk_corner",
     (const char * []){"u","u_prev","error", NULL},
     (double* []){prob_vecs.u, u_prev, error},
     (const char * []){"estimator","error_l2",NULL},
     (double* []){estimator,error_l2},
     NULL,
     NULL,
     0
    );    
      

    double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));

    d4est_element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       prob_vecs.u,
       u_vertex
      );



    char* folder = d4est_util_add_cwd("VTK_corner");
    d4est_util_make_directory(folder,0);
    /* if (sub_folder_number >= 0){ */
    asprintf(&folder,"%s%d/", folder, 0);
    /* } */
    d4est_util_make_directory(folder,0);

    
    /* double* u_vertex = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN)); */

    /* double* u_min_one_vertex */
      /* = P4EST_ALLOC */
      /* ( */
       /* double, */
       /* p4est->local_num_quadrants*(P4EST_CHILDREN) */
      /* ); */

 
    d4est_element_data_store_nodal_vec_in_vertex_array
      (
       p4est,
       prob_vecs.u,
       u_vertex
      );


   /* for (int i = 0; i < p4est->local_num_quadrants*(P4EST_CHILDREN); */
         /* i++){ */
      /* u_min_one_vertex[i] = u_vertex[i] - 1.0; */
    /* } */
    
    
    char* u_corner_file = P4EST_ALLOC(char, 100);
    sprintf(u_corner_file, "%s_%d", "u_corner", 0);
    /* sprintf(u_corner_file, "%s_%d",  "u_corner"); */
    p4est_vtk_ext_write_all
      (p4est,
       NULL,
       0.99,
       1,
       1,
       1,
       1,
       1,
       0,
       u_corner_file,
       "u",
       u_vertex
      );

    P4EST_FREE(u_corner_file);    
    /* P4EST_FREE(u_analytic); */
    P4EST_FREE(u_vertex);
    /* P4EST_FREE(error); */


  
  if (d4est_ghost_data != NULL){
    d4est_ghost_data_destroy(d4est_ghost_data);
    d4est_ghost_data = NULL;
  } 

  
  d4est_mesh_data_destroy(d4est_factors_compactified);
  d4est_geometry_destroy(d4est_geom_compactified);
  d4est_amr_destroy(d4est_amr);
  d4est_amr_destroy(d4est_amr_uniform_p);
  d4est_norms_linear_fit_destroy(l2_linear_fit);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_jac);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_res);
  P4EST_FREE(error);

  d4est_norms_linear_fit_destroy(point_3m_fit);
  
  P4EST_FREE(u_prev);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
}
