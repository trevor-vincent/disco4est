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
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_checkpoint.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>
#include <d4est_laplacian_flux_sipg_penalty_debugger.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <d4est_util.h>
#include <time.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_multigrid_logger_residual.h>
#include <d4est_krylov_pc_multigrid.h>
#include <d4est_solver_multigrid_element_data_updater.h>
#include <d4est_hessian.h>
#include "poisson_lorentzian_fcns_with_opt.h"
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
  double pp = 2 - c;
  double q = R;
  double x = q;
  return x - Rc;  
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
  /* = d4est_estimator_stats_get_percentile(stats,params->percentile); */
  
  return ((eta2 >= eta2_percentile) || fabs(eta2 - eta2_percentile) < eta2*1e-4);
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




typedef struct {
  
  int use_dirichlet;
  int use_pointwise_estimator;
  int analyze_matrix;
  int do_not_save_checkpoint;
  
} poisson_lorentzian_init_params_t;


static
int poisson_lorentzian_init_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  poisson_lorentzian_init_params_t* pconfig = (poisson_lorentzian_init_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"use_dirichlet")) {
    D4EST_ASSERT(pconfig->use_dirichlet == -1);
    pconfig->use_dirichlet = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"use_pointwise_estimator")) {
    D4EST_ASSERT(pconfig->use_pointwise_estimator == -1);
    pconfig->use_pointwise_estimator = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"do_not_save_checkpoint")) {
    D4EST_ASSERT(pconfig->do_not_save_checkpoint == 0);
    pconfig->do_not_save_checkpoint = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"analyze_matrix")) {
    D4EST_ASSERT(pconfig->analyze_matrix == 0);
    pconfig->analyze_matrix = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
poisson_lorentzian_init_params_t
poisson_lorentzian_init_params_input
(
 const char* input_file
)
{
  poisson_lorentzian_init_params_t input;
  input.use_dirichlet = -1;
  input.use_pointwise_estimator = -1;
  input.analyze_matrix = 0;
  input.do_not_save_checkpoint = 0;

  if (ini_parse(input_file, poisson_lorentzian_init_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.use_dirichlet, -1);
  D4EST_CHECK_INPUT("problem", input.use_pointwise_estimator, -1);
  /* D4EST_CHECK_INPUT("problem", input.analyze_matrix, -1); */
  
  return input;
}

int
problem_set_mortar_degree
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg;
}

void
problem_set_degrees_after_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;
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
  
  poisson_lorentzian_init_params_t init_params = poisson_lorentzian_init_params_input
                                                 (
                                                  input_file
                                                 ); 

  
  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;


  lorentzian_params_t lorentzian_params;
  lorentzian_params.R_surface = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user)->R2;
  // Setup boundary conditions

  d4est_laplacian_with_opt_robin_bc_t bc_data_robin_for_lhs;
  bc_data_robin_for_lhs.robin_coeff = poisson_lorentzian_robin_coeff_fcn;
  bc_data_robin_for_lhs.robin_rhs = poisson_lorentzian_robin_bc_rhs_fcn;
  
  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_dirichlet_for_lhs;
  bc_data_dirichlet_for_lhs.dirichlet_fcn = zero_fcn;
  bc_data_dirichlet_for_lhs.eval_method = eval_method;
  bc_data_dirichlet_for_lhs.user = &lorentzian_params;
  
  d4est_laplacian_with_opt_dirichlet_bc_t bc_data_dirichlet_for_rhs;
  bc_data_dirichlet_for_rhs.dirichlet_fcn = poisson_lorentzian_boundary_fcn;
  bc_data_dirichlet_for_rhs.eval_method = eval_method;
  bc_data_dirichlet_for_rhs.user = &lorentzian_params;
  
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_lhs = NULL; //d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, 
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_rhs = NULL; //d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_rhs);

  if(init_params.use_dirichlet){
    flux_data_for_lhs
      = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_dirichlet_for_lhs);
  
    flux_data_for_rhs
      = d4est_laplacian_with_opt_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_dirichlet_for_rhs);
  }
  else {  
    flux_data_for_lhs = d4est_laplacian_with_opt_flux_new(p4est, input_file, BC_ROBIN, &bc_data_robin_for_lhs);
    flux_data_for_rhs = d4est_laplacian_with_opt_flux_new(p4est, input_file,  BC_ROBIN, &bc_data_robin_for_lhs);
  }
  


  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_rhs;

                           
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_lorentzian_build_residual;
  prob_fcns.apply_lhs = poisson_lorentzian_apply_lhs;
  prob_fcns.user = &ctx;
  
  /* double* error = NULL; */
  /* double* u_analytic = NULL; */
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_laplacian_with_opt_flux_sipg_params_t* sipg_params = flux_data_for_lhs->flux_data;
  
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

  d4est_amr_t* d4est_amr_uniform = d4est_amr_init_uniform_h(p4est, d4est_amr->num_of_amr_steps);

  int initial_level = 0;
  if (initial_extents->load_from_checkpoint == 0 || initial_extents->checkpoint_prefix == NULL){
    d4est_mesh_init_field
      (
       p4est,
       prob_vecs.u,
       poisson_lorentzian_initial_guess,
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
  
    

  d4est_field_type_t field_type = NODAL;
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               *d4est_ghost,
                                                               &field_type,
                                                               1);
  
  d4est_laplacian_with_opt_build_rhs_with_strong_bc
    (
     p4est,
     *d4est_ghost,
     d4est_ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &prob_vecs,
     flux_data_for_rhs,
     prob_vecs.rhs,
     poisson_lorentzian_rhs_fcn,
     INIT_FIELD_ON_LOBATTO,
     &ctx,
     0
    );


  // Norm function contexts
    
  d4est_norms_fcn_L2_ctx_t L2_norm_ctx;
  L2_norm_ctx.p4est = p4est;
  L2_norm_ctx.d4est_ops = d4est_ops;
  L2_norm_ctx.d4est_geom = d4est_geom;
  L2_norm_ctx.d4est_quad = d4est_quad;
  L2_norm_ctx.d4est_factors = d4est_factors;
    
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry_compactified");
  d4est_geometry_t* d4est_geom_compactified = d4est_geometry_new(p4est->mpirank, input_file,"compactified_geometry", c_geom);
  d4est_mesh_data_t* d4est_factors_compactified = d4est_mesh_data_init(p4est);

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
     /* *ghost_data, */
     d4est_ops,
     d4est_geom_compactified,
     d4est_quad,
     d4est_factors_compactified,
     initial_extents->face_h_type,
     initial_extents->volume_h_type
    );

    double* estimator_vtk = P4EST_ALLOC(double, 4*p4est->local_num_quadrants);
    double* estimator_vtk_per_face = P4EST_ALLOC(double, 18*p4est->local_num_quadrants);

    double* residual_pointwise_quad = P4EST_ALLOC_ZERO(double, d4est_factors->local_sizes.local_nodes_quad);
    double* f_quad = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
    double* f = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes);

    d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       HESSIAN_ANALYTICAL,
       prob_vecs.u,
       residual_pointwise_quad
      );
    
    d4est_mesh_init_field
      (
       p4est,
       f,
       poisson_lorentzian_rhs_fcn,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );

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
      d4est_quadrature_interpolate
        (
         d4est_ops,
         d4est_quad,
         d4est_geom,
         &mesh_object,
         QUAD_OBJECT_VOLUME,
         QUAD_INTEGRAND_UNKNOWN,
         &f[ed->nodal_stride],
         ed->deg,
         &f_quad[ed->quad_stride],
         ed->deg_quad
        );

    }
  }
    
    
    for (int qnode = 0; qnode < d4est_factors->local_sizes.local_nodes_quad; qnode++){

      residual_pointwise_quad[qnode] =
        f_quad[qnode] - residual_pointwise_quad[qnode];
      
    }
    
    double* estimator =
      d4est_estimator_bi_new_compute
                        (
                         p4est,
                         &prob_vecs,
                         residual_pointwise_quad,
                         penalty_data,
                         poisson_lorentzian_boundary_fcn,
                         &lorentzian_params,
                         *d4est_ghost,
                         d4est_ghost_data,
                         d4est_ops,
                         d4est_geom,
                         d4est_factors,
                         d4est_geom_compactified,
                         d4est_factors_compactified,
                         d4est_quad,
                         0,
                         estimator_vtk,
                         estimator_vtk_per_face,
                         init_params.use_pointwise_estimator ? 1 : 0, 
                         &prob_fcns
                        );

    P4EST_FREE(residual_pointwise_quad);
    P4EST_FREE(f_quad);
    P4EST_FREE(f);

    d4est_amr_smooth_pred_params_t* sp_params = d4est_amr_smooth_pred_params_input
                                                (
                                                 input_file
                                                );
    
    d4est_estimator_stats_t* stats = P4EST_ALLOC(d4est_estimator_stats_t,1);
    d4est_estimator_stats_compute(p4est, estimator, stats, sp_params->percentile, 1, 0);
    d4est_estimator_stats_print(stats);

    double* error_l2 = P4EST_ALLOC(double, p4est->local_num_quadrants);
    double* u_analytic = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* error = P4EST_ALLOC(double, prob_vecs.local_nodes);

    d4est_mesh_init_field
      (
       p4est,
       u_analytic,
       poisson_lorentzian_analytic_solution,
       d4est_ops,
       d4est_geom,
       d4est_factors,
       INIT_FIELD_ON_LOBATTO,
       &ctx
      );

    for (int i = 0; i < prob_vecs.local_nodes; i++){
      error[i] = fabs(prob_vecs.u[i] - u_analytic[i]);
    }

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

    double* a = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* b = P4EST_ALLOC(double, prob_vecs.local_nodes);
    double* c = P4EST_ALLOC(double, prob_vecs.local_nodes);

    d4est_mesh_data_get_topological_coords
      (
       p4est,
       d4est_ops,
       a,
       b,
       c
      );

    d4est_laplacian_flux_sipg_penalty_debugger_t* debugger =
      d4est_laplacian_flux_sipg_penalty_debugger_init
      (
       p4est,
       input_file,
       sipg_params->sipg_penalty_fcn,
       sipg_params->sipg_penalty_prefactor
      );

    d4est_laplacian_flux_sipg_penalty_debugger_get_vtk_data
      (
       p4est,
       *d4est_ghost,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       debugger
      );
    
    d4est_vtk_save
      (
       p4est,
       d4est_ops,
       input_file,
       "d4est_vtk",
       (const char * []){"rhs", "u","u_analytic","error","a","b","c", NULL},
       (double* []){prob_vecs.rhs, prob_vecs.u, u_analytic, error,a,b,c},
       (const char * []){"estimator","error_l2","estimator_R", "estimator_Je1", "estimator_Je2", "estimator_Je3",
                           "estimator_Je1_face0",
                           "estimator_Je1_face1",
                           "estimator_Je1_face2",
                           "estimator_Je1_face3",
                           "estimator_Je1_face4",
                           "estimator_Je1_face5",
                           "estimator_Je2_face0",
                           "estimator_Je2_face1",
                           "estimator_Je2_face2",
                           "estimator_Je2_face3",
                           "estimator_Je2_face4",
                           "estimator_Je2_face5",
                           "estimator_Je3_face0",
                           "estimator_Je3_face1",
                           "estimator_Je3_face2",
                           "estimator_Je3_face3",
                           "estimator_Je3_face4",
                           "estimator_Je3_face5",
                           "average_min_penalty",
                           "average_max_penalty",
                           "average_mean_penalty",
                           NULL},
       (double* []){estimator,error_l2,&estimator_vtk[0], &estimator_vtk[1*p4est->local_num_quadrants],&estimator_vtk[2*p4est->local_num_quadrants],&estimator_vtk[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[0*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[1*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[2*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[4*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[5*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[6*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[7*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[8*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[9*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[10*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[11*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[12*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[13*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[14*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[15*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[16*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[17*p4est->local_num_quadrants],
                      debugger->average_min_penalty_vtk,
                      debugger->average_max_penalty_vtk,
                      debugger->average_mean_penalty_vtk
                      },
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
       (const char * []){"rhs", "u","u_analytic","error","a","b","c", NULL},
       (double* []){prob_vecs.rhs, prob_vecs.u, u_analytic, error,a,b,c},
       (const char * []){"estimator","error_l2","estimator_R", "estimator_Je1", "estimator_Je2", "estimator_Je3",
                           "estimator_Je1_face0",
                           "estimator_Je1_face1",
                           "estimator_Je1_face2",
                           "estimator_Je1_face3",
                           "estimator_Je1_face4",
                           "estimator_Je1_face5",
                           "estimator_Je2_face0",
                           "estimator_Je2_face1",
                           "estimator_Je2_face2",
                           "estimator_Je2_face3",
                           "estimator_Je2_face4",
                           "estimator_Je2_face5",
                           "estimator_Je3_face0",
                           "estimator_Je3_face1",
                           "estimator_Je3_face2",
                           "estimator_Je3_face3",
                           "estimator_Je3_face4",
                           "estimator_Je3_face5",
                           "average_min_penalty",
                           "average_max_penalty",
                           "average_mean_penalty",
                           NULL},
       (double* []){estimator,error_l2,&estimator_vtk[0], &estimator_vtk[1*p4est->local_num_quadrants],&estimator_vtk[2*p4est->local_num_quadrants],&estimator_vtk[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[0*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[1*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[2*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[4*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[5*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[6*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[7*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[8*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[9*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[10*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[11*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[12*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[13*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[14*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[15*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[16*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[17*p4est->local_num_quadrants],
                      debugger->average_min_penalty_vtk,
                      debugger->average_max_penalty_vtk,
                      debugger->average_mean_penalty_vtk
                      },
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
       (const char * []){"rhs", "u","u_analytic","error","a","b","c", NULL},
       (double* []){prob_vecs.rhs, prob_vecs.u, u_analytic, error,a,b,c},
       (const char * []){"estimator","error_l2","estimator_R", "estimator_Je1", "estimator_Je2", "estimator_Je3",
                           "estimator_Je1_face0",
                           "estimator_Je1_face1",
                           "estimator_Je1_face2",
                           "estimator_Je1_face3",
                           "estimator_Je1_face4",
                           "estimator_Je1_face5",
                           "estimator_Je2_face0",
                           "estimator_Je2_face1",
                           "estimator_Je2_face2",
                           "estimator_Je2_face3",
                           "estimator_Je2_face4",
                           "estimator_Je2_face5",
                           "estimator_Je3_face0",
                           "estimator_Je3_face1",
                           "estimator_Je3_face2",
                           "estimator_Je3_face3",
                           "estimator_Je3_face4",
                           "estimator_Je3_face5",
                           "average_min_penalty",
                           "average_max_penalty",
                           "average_mean_penalty",
                           NULL},
       (double* []){
         estimator,error_l2,&estimator_vtk[0],
           &estimator_vtk[1*p4est->local_num_quadrants],
           &estimator_vtk[2*p4est->local_num_quadrants],&estimator_vtk[3*p4est->local_num_quadrants],
           &estimator_vtk_per_face[0*p4est->local_num_quadrants],
           &estimator_vtk_per_face[1*p4est->local_num_quadrants],
           &estimator_vtk_per_face[2*p4est->local_num_quadrants],
           &estimator_vtk_per_face[3*p4est->local_num_quadrants],
           &estimator_vtk_per_face[4*p4est->local_num_quadrants],
           &estimator_vtk_per_face[5*p4est->local_num_quadrants],
           &estimator_vtk_per_face[6*p4est->local_num_quadrants],
           &estimator_vtk_per_face[7*p4est->local_num_quadrants],
           &estimator_vtk_per_face[8*p4est->local_num_quadrants],
           &estimator_vtk_per_face[9*p4est->local_num_quadrants],
           &estimator_vtk_per_face[10*p4est->local_num_quadrants],
           &estimator_vtk_per_face[11*p4est->local_num_quadrants],
           &estimator_vtk_per_face[12*p4est->local_num_quadrants],
           &estimator_vtk_per_face[13*p4est->local_num_quadrants],
           &estimator_vtk_per_face[14*p4est->local_num_quadrants],
           &estimator_vtk_per_face[15*p4est->local_num_quadrants],
           &estimator_vtk_per_face[16*p4est->local_num_quadrants],
           &estimator_vtk_per_face[17*p4est->local_num_quadrants],
           debugger->average_min_penalty_vtk,
           debugger->average_max_penalty_vtk,
           debugger->average_mean_penalty_vtk
           },
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
       (const char * []){"rhs", "u","u_analytic","error","a","b","c", NULL},
       (double* []){prob_vecs.rhs, prob_vecs.u, u_analytic, error,a,b,c},
       (const char * []){"estimator","error_l2","estimator_R", "estimator_Je1", "estimator_Je2", "estimator_Je3",
                           "estimator_Je1_face0",
                           "estimator_Je1_face1",
                           "estimator_Je1_face2",
                           "estimator_Je1_face3",
                           "estimator_Je1_face4",
                           "estimator_Je1_face5",
                           "estimator_Je2_face0",
                           "estimator_Je2_face1",
                           "estimator_Je2_face2",
                           "estimator_Je2_face3",
                           "estimator_Je2_face4",
                           "estimator_Je2_face5",
                           "estimator_Je3_face0",
                           "estimator_Je3_face1",
                           "estimator_Je3_face2",
                           "estimator_Je3_face3",
                           "estimator_Je3_face4",
                           "estimator_Je3_face5",
                           "average_min_penalty",
                           "average_max_penalty",
                           "average_mean_penalty",
                           NULL},
       (double* []){estimator,error_l2,&estimator_vtk[0], &estimator_vtk[1*p4est->local_num_quadrants],&estimator_vtk[2*p4est->local_num_quadrants],&estimator_vtk[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[0*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[1*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[2*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[3*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[4*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[5*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[6*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[7*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[8*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[9*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[10*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[11*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[12*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[13*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[14*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[15*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[16*p4est->local_num_quadrants],
                      &estimator_vtk_per_face[17*p4est->local_num_quadrants],
                      debugger->average_min_penalty_vtk,
                      debugger->average_max_penalty_vtk,
                      debugger->average_mean_penalty_vtk
                      },
       NULL,
       NULL,
       0
      );    

    d4est_laplacian_flux_sipg_penalty_debugger_destroy
      (
       debugger
      );


    /* d4est_vtk_save */
    /*   ( */
    /*    p4est, */
    /*    d4est_ops, */
    /*    input_file, */
    /*    "d4est_vtk", */
    /*    (const char * []){"u","u_analytic","error", NULL}, */
    /*    (double* []){prob_vecs.u, u_analytic, error}, */
    /*    (const char * []){NULL}, */
    /*    (double* []){NULL}, */
    /*    (const char * []){"degrees",NULL}, */
    /*    (int* []){deg_array,NULL}, */
    /*    0 */
    /*   ); */

    /* P4EST_FREE(deg_array); */
    
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
    P4EST_FREE(u_analytic);
    P4EST_FREE(u_vertex);
    P4EST_FREE(error);



     
    P4EST_FREE(u_analytic);
    P4EST_FREE(error_l2);
    P4EST_FREE(a);
    P4EST_FREE(b);
    P4EST_FREE(c);

      
    if (d4est_ghost_data != NULL){
      d4est_ghost_data_destroy(d4est_ghost_data);
      d4est_ghost_data = NULL;
    } 
    
        
    P4EST_FREE(estimator_vtk);
    P4EST_FREE(estimator_vtk_per_face);
    P4EST_FREE(estimator);

  
  d4est_mesh_data_destroy(d4est_factors_compactified);
  d4est_geometry_destroy(d4est_geom_compactified);  
  printf("[D4EST_INFO]: Starting garbage collection...\n");
  d4est_amr_destroy(d4est_amr);
  d4est_amr_destroy(d4est_amr_uniform);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_lhs);
  d4est_laplacian_with_opt_flux_destroy(flux_data_for_rhs);
  /* d4est_norms_linear_fit_destroy(point_3m_fit); */

  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
