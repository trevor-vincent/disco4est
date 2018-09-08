#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <problem.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_solver_newton.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_vtk.h>
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <d4est_solver_newton_petsc.h>
#include <d4est_solver_krylov_petsc.h>
#include <multigrid.h>
#include <d4est_krylov_pc_multigrid.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater.h>
#include <multigrid_matrix_operator.h>
#include <d4est_util.h>
#include <time.h>
#include "./okendon_fcns.h"


int
problem_set_mortar_degree
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg_quad;
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
 d4est_mesh_data_t* geometric_factors,
 d4est_mesh_initial_extents_t* initial_extents,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{
  okendon_params_t okendon_params = okendon_params_init(input_file);
  int initial_nodes = initial_extents->initial_nodes;

  
  d4est_poisson_dirichlet_bc_t bc_data_for_jac;
  bc_data_for_jac.dirichlet_fcn = zero_fcn;
  bc_data_for_jac.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;
  
  d4est_poisson_dirichlet_bc_t bc_data_for_res;
  bc_data_for_res.dirichlet_fcn = okendon_boundary_fcn;
  bc_data_for_res.eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;

  d4est_poisson_flux_data_t* flux_data_for_jac = d4est_poisson_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_jac, problem_set_mortar_degree, NULL);
  
  d4est_poisson_flux_data_t* flux_data_for_res = d4est_poisson_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_for_res, problem_set_mortar_degree, NULL);

  problem_ctx_t ctx;
  ctx.flux_data_for_jac = flux_data_for_jac;
  ctx.flux_data_for_res = flux_data_for_res;
  ctx.okendon_params = &okendon_params;

  bc_data_for_jac.user = &ctx;
  bc_data_for_res.user = &ctx;

  
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = okendon_build_residual;
  prob_fcns.apply_lhs = okendon_apply_jac;
  prob_fcns.user = &ctx;
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;

   d4est_poisson_flux_sipg_params_t* sipg_params = flux_data_for_jac->flux_data;
  
  d4est_estimator_bi_penalty_data_t penalty_data;
  penalty_data.u_penalty_fcn = houston_u_prefactor_maxp_minh;
  penalty_data.u_dirichlet_penalty_fcn = houston_u_dirichlet_prefactor_maxp_minh;
  penalty_data.gradu_penalty_fcn = houston_gradu_prefactor_maxp_minh;
  penalty_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
  penalty_data.sipg_flux_h = sipg_params->sipg_flux_h;
  penalty_data.user = &okendon_params;
  
  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     NULL
    );


  d4est_mesh_init_field
    (
     p4est,
     prob_vecs.u,
     okendon_initial_guess,
     d4est_ops,
     d4est_geom,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
  
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){

    d4est_estimator_bi_compute
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       penalty_data,
       okendon_boundary_fcn,
       *ghost,
       *ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       NO_DIAM_APPROX,
       problem_set_mortar_degree,
       NULL
      );
    
    d4est_estimator_stats_t stats;
    d4est_estimator_stats_compute(p4est, &stats);
    d4est_estimator_stats_print(&stats);

    
    d4est_output_vtk_with_analytic_error
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       input_file,
       "uniform_okendon",
       okendon_analytic_solution,
       &ctx,
       1,
       level
      );

    d4est_output_vtk_degree_mesh_with_analytic_error
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       okendon_analytic_solution,
       &ctx,
       prob_vecs.local_nodes,
       input_file,
       "okendon_degree_mesh",
       1,
       level
      );

    d4est_ip_energy_norm_data_t ip_norm_data;
    ip_norm_data.u_penalty_fcn = sipg_params->sipg_penalty_fcn;
    ip_norm_data.sipg_flux_h = sipg_params->sipg_flux_h;
    ip_norm_data.penalty_prefactor = sipg_params->sipg_penalty_prefactor;
    
    d4est_norms_norms_using_analytic_solution
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       *ghost,
       *ghost_data,
       stats.total,
       &prob_vecs,
       &ip_norm_data,
       okendon_analytic_solution,
       &ctx,NULL,NULL);


    if (level != d4est_amr->num_of_amr_steps){

      if (p4est->mpirank == 0)
        printf("[D4EST_INFO]: AMR REFINEMENT LEVEL %d\n", level+1);

      d4est_amr_step
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_amr,
         /* d4est_amr, */
         /* d4est_amr_uniform_p, */
         &prob_vecs.u,
         NULL
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
                   d4est_mesh_set_quadratures_after_amr,
                   initial_extents
                  );

    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, prob_vecs.local_nodes);

    
   int min_level, max_level;

    multigrid_get_level_range(p4est, &min_level, &max_level);
    if (p4est->mpirank == 0){
      printf("[min_level, max_level] = [%d,%d]\n", min_level, max_level);
    }
    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* mpi_assert(proc_size == 1); */
    int num_of_levels = max_level + 1;

 
    multigrid_logger_t* logger = multigrid_logger_residual_init
                                 (
                                 );
    
    multigrid_element_data_updater_t* updater = multigrid_element_data_updater_init
                                                (
                                                 num_of_levels,
                                                 ghost,
                                                 ghost_data,
                                                 geometric_factors,
                                                 d4est_mesh_set_quadratures_after_amr,
                                                 initial_extents
                                                );
    





    multigrid_user_callbacks_t* user_callbacks = multigrid_matrix_operator_init(p4est, num_of_levels);

    /* prob_vecs.u0 = prob_vecs.u; */
    
    /* multigrid_matrix_setup_fofufofvlilj_operator */
    /*   ( */
    /*    p4est, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    d4est_quad, */
    /*    prob_vecs.u0, */
    /*    NULL, */
    /*    neg_10pi_rho_up1_neg4, */
    /*    &ctx, */
    /*    NULL, */
    /*    NULL, */
    /*    user_callbacks->user */
    /*   ); */

    
    d4est_solver_multigrid_data_t* mg_data = multigrid_data_init(p4est,
                                                    d4est_ops,
                                                    d4est_geom,
                                                    d4est_quad,
                                                    num_of_levels,
                                                    logger,
                                                    user_callbacks,
                                                    updater,
                                                    input_file
                                                   );

    d4est_krylov_pc_t* pc = d4est_krylov_pc_multigrid_create(mg_data, okendond4est_krylov_pc_setup_fcn);
    ctx.use_matrix_operator = 1;
    ctx.mg_data = mg_data;

    d4est_solver_newton_petsc_params_t newton_params;
    d4est_solver_newton_petsc_input(p4est, input_file, "[NEWTON_PETSC]", &newton_params);

    d4est_solver_krylov_petsc_params_t krylov_params;
    d4est_solver_krylov_petsc_input(p4est, input_file, "d4est_solver_krylov_petsc", "[KRYLOV_PETSC]", &krylov_params);
      
      d4est_solver_newton_petsc_solve
        (
         p4est,
         &prob_vecs,
         &prob_fcns,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         &krylov_params,
         &newton_params,
         pc
        );

    d4est_krylov_pc_multigrid_destroy(pc);
    multigrid_logger_residual_destroy(logger);
    multigrid_element_data_updater_destroy(updater, num_of_levels);
    multigrid_data_destroy(mg_data);
    multigrid_matrix_operator_destroy(user_callbacks);
  }

  printf("[D4EST_INFO]: Starting garbage collection...\n");

  d4est_amr_destroy(d4est_amr);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_poisson_flux_destroy(flux_data_for_jac);
  d4est_poisson_flux_destroy(flux_data_for_res);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
}
