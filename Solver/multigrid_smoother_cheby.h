#ifndef MULTIGRID_SMOOTHER_CHEBY_H
#define MULTIGRID_SMOOTHER_CHEBY_H 

#include <multigrid.h>

typedef struct {

  /* INTERNAL */
  int mpirank;
  int smoother_eigs_compute;

  /* SET IN INPUT FILE */
  int smoother_imax;
  int smoother_eigs_cg_imax;
  double smoother_eigs_lmax_lmin_ratio;
  double smoother_eigs_max_multiplier;
  int smoother_eigs_reuse_fromdownvcycle;
  int smoother_eigs_reuse_fromlastvcycle;
  int smoother_print_residual_norm;
  int smoother_print_eigs;

} multigrid_smoother_cheby_t;

/* This file was automatically generated.  Do not edit! */
void multigrid_smoother_cheby_destroy(multigrid_smoother_t *bottom);
multigrid_smoother_t *multigrid_smoother_cheby_init(p4est_t *p4est,int num_of_levels,const char *input_file);


#endif
