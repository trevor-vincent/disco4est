#ifndef D4EST_SOLVER_MULTIGRID_SMOOTHER_SCHWARZ_H
#define D4EST_SOLVER_MULTIGRID_SMOOTHER_SCHWARZ_H 

#include <d4est_solver_multigrid.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <pXest.h>

typedef struct {

  d4est_solver_schwarz_metadata_t** schwarz_metadata_on_level;
  d4est_solver_schwarz_operators_t* schwarz_ops;
  char* input_file;
  int num_of_levels;

  d4est_laplacian_dirichlet_bc_t bc_data_dirichlet_for_lhs;
  d4est_laplacian_flux_data_t* flux_data_for_lhs;
  d4est_solver_schwarz_laplacian_mortar_data_t** laplacian_mortar_data_on_level;
  
} d4est_solver_multigrid_smoother_schwarz_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_smoother_schwarz_destroy(d4est_solver_multigrid_smoother_t *smoother);
void d4est_solver_multigrid_smoother_schwarz_update(p4est_t *p4est,int level,d4est_elliptic_data_t *vecs);
d4est_solver_multigrid_smoother_t *d4est_solver_multigrid_smoother_schwarz_init(p4est_t *p4est,int num_of_levels,d4est_operators_t *d4est_ops,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,const char *input_file);

#endif
