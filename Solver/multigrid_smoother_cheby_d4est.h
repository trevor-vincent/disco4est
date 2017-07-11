#ifndef MULTIGRID_SMOOTHER_CHEBY_D4EST_H
#define MULTIGRID_SMOOTHER_CHEBY_D4EST_H 

#include <multigrid.h>

typedef struct {

  /* INTERNAL */
  int mpirank;
  int cheby_eigs_compute;

  /* SET IN INPUT FILE */
  double* eigs;
  int cheby_imax;
  int cheby_eigs_cg_imax;
  double cheby_eigs_lmax_lmin_ratio;
  double cheby_eigs_max_multiplier;
  int cheby_eigs_reuse_fromdownvcycle;
  int cheby_eigs_reuse_fromlastvcycle;
  int cheby_print_residual_norm;
  int cheby_print_eigs;

} multigrid_smoother_cheby_d4est_t;

/* This file was automatically generated.  Do not edit! */
void multigrid_smoother_cheby_d4est_destroy(multigrid_smoother_t *smoother);
multigrid_smoother_t *multigrid_smoother_cheby_d4est_init(p4est_t *p4est,int num_of_levels,const char *input_file);

void 
multigrid_smoother_cheby_d4est_iterate
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int iter,
 double lmin,
 double lmax,
 int print_residual_norm
 /* multigrid_cheby_params_t* cheby_params */
);


#endif
