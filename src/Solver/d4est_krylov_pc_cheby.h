#ifndef D4EST_KRYLOV_PC_CHEBY_H
#define D4EST_KRYLOV_PC_CHEBY_H 

#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_solver_petsc.h>
#include <d4est_krylov_pc.h>

typedef struct {

  double eig;
  int cheby_imax;
  int cheby_eigs_cg_imax;
  double cheby_eigs_lmax_lmin_ratio;
  double cheby_eigs_max_multiplier;
  int cheby_print_residual_norm;
  int cheby_print_spectral_bound;
  int cheby_print_spectral_bound_iterations;
  int cheby_reuse_eig;
  int cheby_use_new_cg_eigs;
  int cheby_use_zero_guess_for_eigs;
  void(*user_setup_fcn)(d4est_krylov_pc_t* kpc);
  
} d4est_krylov_pc_cheby_data_t;


/* This file was automatically generated.  Do not edit! */
d4est_krylov_pc_t *d4est_krylov_pc_cheby_create(const char *input_file,void(*user_setup_fcn)(d4est_krylov_pc_t *kpc));
void d4est_krylov_pc_cheby_apply(d4est_krylov_pc_t *kpc,double *xp,double *yp);
void d4est_krylov_pc_cheby_destroy(d4est_krylov_pc_t *kpc);
void d4est_krylov_pc_cheby_setup(d4est_krylov_pc_t *kpc);

#endif
