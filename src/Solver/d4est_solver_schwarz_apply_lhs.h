#ifndef D4EST_SOLVER_SCHWARZ_APPLY_LHS_H
#define D4EST_SOLVER_SCHWARZ_APPLY_LHS_H 

#include <pXest.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_elliptic_data.h>
/*  */
typedef
void
(*d4est_solver_schwarz_apply_lhs_fcn_t)
(
 p4est_t*,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_mesh_data_t*,
 d4est_ghost_t*,
 d4est_solver_schwarz_operators_t*,
 d4est_solver_schwarz_metadata_t*,
 d4est_solver_schwarz_geometric_data_t*,
 int,
 double*,
 double*,
 void*
);

typedef
void
(*d4est_solver_schwarz_apply_lhs_getctx_fcn_t)
(
 p4est_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_mesh_data_t*,
 d4est_ghost_t*,
 d4est_solver_schwarz_operators_t*,
 d4est_solver_schwarz_metadata_t*,
 d4est_solver_schwarz_geometric_data_t*,
 d4est_elliptic_data_t*,
 void*
);

typedef struct {
  
  d4est_solver_schwarz_apply_lhs_fcn_t apply_lhs_fcn;
  d4est_solver_schwarz_apply_lhs_getctx_fcn_t get_ctx_fcn;
  void* apply_lhs_ctx;

} d4est_solver_schwarz_apply_lhs_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_apply_lhs_destroy(d4est_solver_schwarz_apply_lhs_t *apply_lhs);
d4est_solver_schwarz_apply_lhs_t *d4est_solver_schwarz_apply_lhs_init(d4est_solver_schwarz_apply_lhs_fcn_t fcn,d4est_solver_schwarz_apply_lhs_getctx_fcn_t get_ctx_fcn,void *fcn_ctx);

#endif
