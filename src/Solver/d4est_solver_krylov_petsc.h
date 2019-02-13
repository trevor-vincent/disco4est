#ifndef KSP_PETSC_SOLVE_H
#define KSP_PETSC_SOLVE_H

#include <petscsnes.h>
#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>

#include <d4est_solver_petsc.h>
#include <d4est_krylov_pc.h>

typedef struct {

  char input_section [50]; /* useful when using KSP in multiple contexts */
  
  char ksp_type [25];   /* the biggest type is "richardson" = 10 */
  char ksp_atol [25];
  char ksp_rtol [25];
  char ksp_max_it [25];

  int ksp_monitor;
  int ksp_view;
  int ksp_converged_reason;
  int ksp_initial_guess_nonzero;
  int ksp_monitor_singular_value;
  int ksp_do_not_use_preconditioner;
  
  char ksp_chebyshev_esteig_steps [25];
  char ksp_chebyshev_esteig [25];
  int ksp_chebyshev_esteig_random;

  int checkpoint_every_n_krylov_its;
    
} d4est_solver_krylov_petsc_params_t;


typedef struct {

  int total_krylov_iterations;
  double residual_norm;
  
} d4est_solver_krylov_petsc_info_t;


/* This file was automatically generated.  Do not edit! */
d4est_solver_krylov_petsc_info_t d4est_solver_krylov_petsc_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t **ghost,d4est_ghost_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_solver_krylov_petsc_params_t *d4est_solver_krylov_petsc_params,d4est_krylov_pc_t *d4est_krylov_pc,int amr_level);
PetscErrorCode d4est_solver_krylov_petsc_monitor(KSP ksp,PetscInt it,PetscReal norm,void *ctx);
void d4est_solver_krylov_petsc_input(p4est_t *p4est,const char *input_file,const char *input_section,d4est_solver_krylov_petsc_params_t *input);
void d4est_solver_krylov_petsc_set_options_database_from_params(d4est_solver_krylov_petsc_params_t *input);

#endif
