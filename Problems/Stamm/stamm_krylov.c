#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_solver_cg.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_output.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <newton_petsc.h>
#include <krylov_petsc.h>
#include <d4est_util.h>
#include <time.h>
#include "stamm_fcns.h"

static
int
amr_mark_element
(
 p4est_t* p4est,
 double eta2,
 d4est_estimator_stats_t** stats,
 d4est_element_data_t* elem_data,
 void* user
)
{
  problem_ctx_t* ctx = user;
  d4est_amr_smooth_pred_params_t* params = ctx->smooth_pred_params;
  
  double eta2_avg = stats[0]->mean;
  printf("eta2_avg, eta2, params->sigma = %.25f,%.25f,%.25f\n",eta2_avg,eta2,params->sigma);
  return (eta2 >= params->sigma*eta2_avg);
}

static
gamma_params_t
amr_set_element_gamma
(
 p4est_t* p4est,
 double eta2,
 d4est_estimator_stats_t** stats,
 d4est_element_data_t* elem_data,
 void* user
)
{
  problem_ctx_t* ctx = user;
  d4est_amr_smooth_pred_params_t* params = ctx->smooth_pred_params;
  
  gamma_params_t gamma_hpn;
  gamma_hpn.gamma_h = params->gamma_h;
  gamma_hpn.gamma_p = params->gamma_p;
  gamma_hpn.gamma_n = params->gamma_n;

  return gamma_hpn;
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
  elem_data->deg_vol_quad = elem_data->deg;
}



void
problem_init
(
 p4est_t* p4est,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* geometric_factors,
 int initial_nodes,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{
  stamm_params_t stamm_params = stamm_params_input(input_file);
  d4est_amr_smooth_pred_params_t smooth_pred_params = d4est_amr_smooth_pred_params_input(input_file);


  d4est_poisson_dirichlet_bc_t bc_data_for_lhs;
  bc_data_for_lhs.dirichlet_fcn = zero_fcn;

  d4est_poisson_dirichlet_bc_t bc_data_for_rhs;
  bc_data_for_rhs.dirichlet_fcn = stamm_boundary_fcn;
  
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs = d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_lhs, problem_set_mortar_degree, NULL);
  
  d4est_poisson_flux_data_t* flux_data_for_build_rhs = d4est_poisson_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_for_rhs, problem_set_mortar_degree, NULL);

  problem_ctx_t ctx;
  ctx.stamm_params = &stamm_params;
  ctx.smooth_pred_params = &smooth_pred_params;
  ctx.flux_data_for_apply_lhs = flux_data_for_apply_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_build_rhs;
                           
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = stamm_build_residual;
  prob_fcns.apply_lhs = stamm_apply_lhs;
  prob_fcns.user = &ctx;
  
  double* error = NULL;
  double* u_analytic = NULL;
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

  d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_apply_lhs->flux_data;
  
  d4est_estimator_bi_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_data.u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_data.gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  penalty_data.sipg_flux_h = sipg_params->sipg_flux_h;
  
  d4est_amr_smooth_pred_marker_t amr_marker;
  amr_marker.user = (void*)&ctx;
  amr_marker.mark_element_fcn = amr_mark_element;
  amr_marker.set_element_gamma_fcn = amr_set_element_gamma;

  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     "[D4EST_AMR]:",
     &amr_marker
    );

  d4est_amr_t* d4est_amr_uniform = d4est_amr_init_uniform_h(p4est, d4est_amr->max_degree, d4est_amr->num_of_amr_steps);

  d4est_mesh_init_field
    (
     p4est,
     prob_vecs.u,
     stamm_initial_guess,
     d4est_ops,
     d4est_geom,
     NULL
    );
    

  d4est_poisson_build_rhs_with_strong_bc
    (
     p4est,
     *ghost,
     *ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &prob_vecs,
     flux_data_for_build_rhs,
     prob_vecs.rhs,
     stamm_rhs_fcn,
     &ctx
    );


  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){
    
    d4est_estimator_bi_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       penalty_data,
       zero_fcn,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       DIAM_APPROX_CUBE,
       problem_set_mortar_degree,
       NULL
      );

    d4est_estimator_stats_t* stats = P4EST_ALLOC(d4est_estimator_stats_t,1);
    d4est_estimator_stats_compute(p4est, stats);
    d4est_estimator_stats_print(stats);

    d4est_output_vtk_with_analytic_error
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       input_file,
       "uniform_stamm",
       stamm_analytic_solution,
       &ctx,
       1,
       level
      );

    d4est_output_vtk_degree_mesh
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       input_file,
       "uniform_stamm_degree_mesh",
       1,
       level
      );


    d4est_ip_energy_norm_data_t ip_norm_data;
    ip_norm_data.u_penalty_fcn = sipg_params->sipg_penalty_fcn;
    ip_norm_data.sipg_flux_h = sipg_params->sipg_flux_h;
    ip_norm_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
    /*  */
    printf("ip_norm_data.penalty_prefactor = %f\n", ip_norm_data.penalty_prefactor);
    
    d4est_output_norms_using_analytic_solution
      (
      p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       *ghost,
       *ghost_data,
       stats->total,
       &prob_vecs,
       &ip_norm_data,
       stamm_analytic_solution,
      &ctx,NULL,NULL);

    P4EST_FREE(stats);
    
    if (level != d4est_amr->num_of_amr_steps){

      if (p4est->mpirank == 0)
        printf("[D4EST_INFO]: AMR REFINEMENT LEVEL %d\n", level+1);

      d4est_amr_step
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         (level > 1) ? d4est_amr : d4est_amr_uniform,
         &prob_vecs.u,
         &stats
        );
      
    }


    prob_vecs.local_nodes = d4est_mesh_update
                  (
                   p4est,
                   *ghost,
                   *ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   geometric_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   problem_set_degrees_after_amr,
                   NULL
                  );

    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, prob_vecs.local_nodes);
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, prob_vecs.local_nodes);
    
    
    d4est_poisson_build_rhs_with_strong_bc
      (
       p4est,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       flux_data_for_build_rhs,
       prob_vecs.rhs,
       stamm_rhs_fcn,
       &ctx
      );

    krylov_petsc_params_t krylov_petsc_params;
    krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &krylov_petsc_params);

    krylov_petsc_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &krylov_petsc_params,
       NULL
      );

  }

  printf("[D4EST_INFO]: Starting garbage collection...\n");
  d4est_amr_destroy(d4est_amr);
  d4est_amr_destroy(d4est_amr_uniform);
  d4est_poisson_flux_destroy(flux_data_for_apply_lhs);  
  d4est_poisson_flux_destroy(flux_data_for_build_rhs);  
  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
  P4EST_FREE(prob_vecs.u);
  P4EST_FREE(prob_vecs.Au);
  P4EST_FREE(prob_vecs.rhs);
}
