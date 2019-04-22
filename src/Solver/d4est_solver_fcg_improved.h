#ifndef D4EST_SOLVER_FCG_IMPROVED_H
#define D4EST_SOLVER_FCG_IMPROVED_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_krylov_pc.h>

typedef struct {

  char input_section [50];  
  int imax;
  int do_not_use_preconditioner;
  int checkpoint_every_n_krylov_its;
  double rtol;
  double atol;
  
} d4est_solver_fcg_params_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_fcg_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t **ghost,d4est_ghost_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *fcg_params,d4est_krylov_pc_t *d4est_krylov_pc,int amr_level,int newton_iteration);
void d4est_solver_fcg_input(p4est_t *p4est,const char *input_file,const char *input_section,d4est_solver_fcg_params_t *input);

#endif
