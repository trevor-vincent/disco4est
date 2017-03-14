#ifndef MULTIGRID_BOTTOM_SOLVER_CG_H
#define MULTIGRID_BOTTOM_SOLVER_CG_H 

#include <multigrid.h>

typedef struct {

  double bottom_atol;
  double bottom_rtol;  
  int bottom_imax;
  int bottom_print_residual_norm;
  
} multigrid_bottom_solver_cg_t;


/* This file was automatically generated.  Do not edit! */
void multigrid_bottom_solver_cg_destroy(multigrid_bottom_solver_t *bottom);
multigrid_bottom_solver_t *multigrid_bottom_solver_cg_init(p4est_t *p4est,const char *input_file);

#endif
