#ifndef KSP_PETSC_SOLVE_H
#define KSP_PETSC_SOLVE_H 

#include <petscsnes.h>
#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>

#include <d4est_petsc.h>
#include <krylov_pc.h>
#include <krylov_info.h>

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
  
    
} krylov_petsc_params_t;

/* This file was automatically generated.  Do not edit! */
krylov_info_t krylov_petsc_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,p4est_ghost_t **ghost,d4est_element_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,krylov_petsc_params_t *krylov_petsc_params,krylov_pc_t *krylov_pc);
void krylov_petsc_input(p4est_t *p4est,const char *input_file,const char *input_section,const char *printf_prefix,krylov_petsc_params_t *input);
void krylov_petsc_set_options_database_from_params(krylov_petsc_params_t *input);

#endif
