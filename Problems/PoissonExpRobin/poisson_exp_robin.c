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
#include <d4est_norms.h>
#include <d4est_mesh.h>
#include <ini.h>
#include <d4est_element_data.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux_sipg.h>
#include <newton_petsc.h>
#include <d4est_solver_test_symmetry.h>
#include <krylov_petsc.h>
#include <d4est_util.h>
#include <time.h>
#include "poisson_exp_robin_fcns.h"

void
problem_set_degrees_after_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_vol_quad = elem_data->deg;
}

int
problem_set_mortar_degree
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  return elem_data->deg_vol_quad;
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
 int initial_nodes,
 const char* input_file,
 sc_MPI_Comm mpicomm
)
{
  D4EST_ASSERT(d4est_geom->geom_type == GEOM_BRICK);

  d4est_poisson_robin_bc_t bc_data_for_lhs;
  bc_data_for_lhs.robin_coeff = poisson_exp_robin_coeff_fcn;
  bc_data_for_lhs.robin_rhs = poisson_exp_robin_bc_rhs_fcn;

  d4est_poisson_dirichlet_bc_t bc_data_for_rhs;
  bc_data_for_rhs.dirichlet_fcn = zero_fcn;
  
  d4est_poisson_flux_data_t* flux_data_for_apply_lhs = d4est_poisson_flux_new(p4est, input_file, BC_ROBIN, &bc_data_for_lhs, problem_set_mortar_degree, NULL);
  
  d4est_poisson_flux_data_t* flux_data_for_build_rhs = d4est_poisson_flux_new(p4est, input_file,  BC_DIRICHLET, &bc_data_for_rhs, problem_set_mortar_degree, NULL);

  problem_ctx_t ctx;
  ctx.flux_data_for_apply_lhs = flux_data_for_apply_lhs;
  ctx.flux_data_for_build_rhs = flux_data_for_build_rhs;
  
  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = poisson_exp_robin_build_residual;
  prob_fcns.apply_lhs = poisson_exp_robin_apply_lhs;
  prob_fcns.user = &ctx;
  
  double* error = NULL;
  double* u_analytic = NULL;
  
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.Au = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.u = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.rhs = P4EST_ALLOC(double, initial_nodes);
  prob_vecs.local_nodes = initial_nodes;



  d4est_amr_t* d4est_amr =
    d4est_amr_init
    (
     p4est,
     input_file,
     NULL
    );
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps + 1; ++level){

    int local_nodes = d4est_mesh_update
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

    if (level == 0){
      prob_vecs.u = P4EST_REALLOC(prob_vecs.u, double, local_nodes);
      d4est_mesh_init_field
        (
         p4est,
         prob_vecs.u,
         poisson_exp_robin_initial_guess,
         d4est_ops,
         d4est_geom,
         INIT_FIELD_ON_LOBATTO,
         NULL
        );
    }
    
    prob_vecs.Au = P4EST_REALLOC(prob_vecs.Au, double, local_nodes);
    prob_vecs.rhs = P4EST_REALLOC(prob_vecs.rhs, double, local_nodes);
    prob_vecs.local_nodes = local_nodes;
    
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
       poisson_exp_robin_build_rhs_fcn,
       INIT_FIELD_ON_LOBATTO,
       NULL
      );
    
    /* d4est_solver_test_symmetry */
    /*   ( */
    /*    p4est, */
    /*    *ghost, */
    /*    *ghost_data, */
    /*    local_nodes, */
    /*    &prob_fcns, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    d4est_quad, */
    /*    SYM_PRINT_UNEQUAL_PAIRS_AND_XYZ, */
    /*    1e-10 */
    /*   ); */

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

    d4est_output_vtk_with_analytic_error
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &prob_vecs,
       input_file,
       "uniform_poisson_exp_robin",
       poisson_exp_robin_analytic_solution,
       NULL,
       0,
       level
      );
    
    d4est_norms_norms_using_analytic_solution
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       *ghost,
       *ghost_data,
       -1,
       &prob_vecs,
       NULL,
       poisson_exp_robin_analytic_solution,
       NULL       ,
       NULL,
       NULL
      );


    
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
         &prob_vecs.u,
         NULL
        );
      
    }
    
  }

  printf("[D4EST_INFO]: Starting garbage collection...\n");


  d4est_amr_destroy(d4est_amr);
  d4est_poisson_flux_destroy(flux_data_for_apply_lhs);
  d4est_poisson_flux_destroy(flux_data_for_build_rhs);

  P4EST_FREE(error);
  P4EST_FREE(u_analytic);
}
/*  */
