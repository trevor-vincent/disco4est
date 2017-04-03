#ifndef MULTIGRID_BOTTOM_SOLVER_CHEBY_H
#define MULTIGRID_BOTTOM_SOLVER_CHEBY_H 

typedef struct {

  double eig;
  int cheby_imax;
  int cheby_eigs_cg_imax;
  double cheby_eigs_lmax_lmin_ratio;
  double cheby_eigs_max_multiplier;
  int cheby_print_residual_norm;
  int cheby_print_eig;

} multigrid_bottom_solver_cheby_t;

#endif
