#ifndef MULTIGRID_BOTTOM_SOLVER_CHEBY_D4EST_H
#define MULTIGRID_BOTTOM_SOLVER_CHEBY_D4EST_H 

#include <multigrid.h>

typedef struct {

  double eig;
  int cheby_imax;
  int cheby_eigs_cg_imax;
  double cheby_eigs_lmax_lmin_ratio;
  double cheby_eigs_max_multiplier;
  int cheby_print_residual_norm;
  int cheby_print_eig;

} multigrid_bottom_solver_cheby_d4est_t;

/* This file was automatically generated.  Do not edit! */
multigrid_bottom_solver_t* multigrid_bottom_solver_cheby_d4est_init(p4est_t *p4est,int num_of_levels,const char *input_file);

void multigrid_bottom_solver_cheby_d4est_destroy(multigrid_bottom_solver_t* bottom_solver);

#endif
