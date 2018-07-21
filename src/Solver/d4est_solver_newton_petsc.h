#ifndef NEWTON_PETSC_H
#define NEWTON_PETSC_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_solver_petsc.h>
#include <d4est_solver_krylov_petsc.h>

typedef struct {

  int snes_monitor;
  int snes_linesearch_monitor;
  int snes_view;
  int snes_converged_reason;
  int snes_ksp_ew;
  
  char snes_atol [25];
  char snes_type [25];
  char snes_rtol [25];
  char snes_stol [25];
  char snes_max_funcs [25];
  char snes_max_it [25];
  char snes_trtol [25];
  char snes_linesearch_order [25];  
  
} newton_petsc_params_t;


typedef struct {

  int total_newton_iterations;
  int total_krylov_iterations;
  double residual_norm;
  
} newton_petsc_info_t;

/* This file was automatically generated.  Do not edit! */
newton_petsc_info_t newton_petsc_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t **ghost,d4est_ghost_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,krylov_petsc_params_t *krylov_options,newton_petsc_params_t *newton_options,d4est_krylov_pc_t *d4est_krylov_pc);
void newton_petsc_set_options_database_from_params(newton_petsc_params_t *input);
void newton_petsc_input(p4est_t *p4est,const char *input_file,const char *printf_prefix,newton_petsc_params_t *input);

#endif
