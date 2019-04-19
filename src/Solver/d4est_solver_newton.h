#ifndef D4EST_SOLVER_NEWTON_H
#define D4EST_SOLVER_NEWTON_H 

#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_krylov_pc.h>
#include <d4est_solver_krylov.h>

typedef void
(*d4est_solver_newton_krylov_fcn_t)
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_ghost_t** ghost,
 d4est_ghost_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params,
 d4est_krylov_pc_t* d4est_krylov_pc,
 int amr_level
);


typedef struct {

  double rtol;
  double atol;
  int imax;
  int imin;
  int monitor;
  
} d4est_solver_newton_params_t;

/* This file was automatically generated.  Do not edit! */
int d4est_solver_newton_solve(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t **ghost,d4est_ghost_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,char *input_file,d4est_solver_newton_krylov_fcn_t krylov_fcn,void *krylov_fcn_params,d4est_krylov_pc_t *d4est_krylov_pc,int amr_level);
void d4est_solver_newton_solve_ksp(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t **ghost,d4est_ghost_data_t **ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params,d4est_krylov_pc_t *d4est_krylov_pc,int amr_level);
d4est_solver_newton_params_t d4est_solver_newton_input(p4est_t *p4est,const char *input_file);

#endif
