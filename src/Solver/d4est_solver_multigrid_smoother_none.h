#ifndef D4EST_SOLVER_MULTIGRID_SMOOTHER_NONE_H
#define D4EST_SOLVER_MULTIGRID_SMOOTHER_NONE_H 

#include <pXest.h>
#include <d4est_solver_multigrid.h>


/* This file was automatically generated.  Do not edit! */
d4est_solver_multigrid_smoother_t *d4est_solver_multigrid_smoother_none_init(p4est_t *p4est,int num_of_levels,const char *input_file);
void d4est_solver_multigrid_smoother_none_destroy(d4est_solver_multigrid_smoother_t *smoother);

#endif
