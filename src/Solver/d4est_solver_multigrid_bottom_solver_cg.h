#ifndef MULTIGRID_BOTTOM_SOLVER_CG_D4EST_H
#define MULTIGRID_BOTTOM_SOLVER_CG_D4EST_H 

#include <d4est_solver_multigrid.h>

typedef struct {

  double bottom_atol;
  double bottom_rtol;  
  int bottom_imax;
  int bottom_print_residual_norm;
  
} d4est_solver_multigrid_bottom_solver_cg_t;


/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_bottom_solver_cg_destroy(d4est_solver_multigrid_bottom_solver_t *bottom);
d4est_solver_multigrid_bottom_solver_t *d4est_solver_multigrid_bottom_solver_cg_init(p4est_t *p4est,const char *input_file);

#endif
