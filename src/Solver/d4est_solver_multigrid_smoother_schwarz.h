#ifndef D4EST_SOLVER_MULTIGRID_SMOOTHER_SCHWARZ_H
#define D4EST_SOLVER_MULTIGRID_SMOOTHER_SCHWARZ_H 

#include <d4est_solver_multigrid.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz.h>
#include <d4est_solver_schwarz_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <pXest.h>

typedef struct {

  d4est_solver_schwarz_t** schwarz_on_level;
  d4est_solver_schwarz_operators_t* schwarz_ops;
  d4est_solver_schwarz_apply_lhs_t* apply_lhs;

  int apply_lhs_is_set;
  char* input_file;
  int num_of_levels;
  int verbose;
  int iterations;
  
} d4est_solver_multigrid_smoother_schwarz_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_smoother_schwarz_destroy(d4est_solver_multigrid_smoother_t *smoother);
void d4est_solver_multigrid_smoother_schwarz_update(p4est_t *p4est,int level,d4est_elliptic_data_t *vecs);
d4est_solver_multigrid_smoother_t *d4est_solver_multigrid_smoother_schwarz_init(p4est_t *p4est,int num_of_levels,d4est_operators_t *d4est_ops,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,const char *input_file);
void d4est_solver_multigrid_smoother_schwarz_set_apply_lhs(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,d4est_solver_schwarz_apply_lhs_t *apply_lhs,d4est_solver_multigrid_t *mg_data,const char *input_file);

#endif
