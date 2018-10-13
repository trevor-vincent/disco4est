#ifndef MULTIGRID_BOTTOM_SOLVER_CHEBY_D4EST_H
#define MULTIGRID_BOTTOM_SOLVER_CHEBY_D4EST_H 

#include <d4est_solver_multigrid.h>

typedef struct {

  double eig;
  int cheby_imax;
  int cheby_eigs_cg_imax;
  double cheby_eigs_lmax_lmin_ratio;
  double cheby_eigs_max_multiplier;
  int cheby_print_residual_norm;
  int cheby_print_spectral_bound;
  int cheby_print_spectral_bound_iterations;
  int cheby_use_new_cg_eigs;
  
} d4est_solver_multigrid_bottom_solver_cheby_t;

/* This file was automatically generated.  Do not edit! */
d4est_solver_multigrid_bottom_solver_t* d4est_solver_multigrid_bottom_solver_cheby_init(p4est_t *p4est,int num_of_levels,const char *input_file);

void d4est_solver_multigrid_bottom_solver_cheby_destroy(d4est_solver_multigrid_bottom_solver_t* bottom_solver);

#endif
