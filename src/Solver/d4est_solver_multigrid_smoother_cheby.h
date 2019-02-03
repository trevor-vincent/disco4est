#ifndef MULTIGRID_SMOOTHER_CHEBY_D4EST_H
#define MULTIGRID_SMOOTHER_CHEBY_D4EST_H 

#include <pXest.h>
#include <d4est_solver_multigrid.h>

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
  int cheby_print_spectral_bound;
  int cheby_print_spectral_bound_iterations;
  int cheby_use_new_cg_eigs;
  int cheby_use_zero_guess_for_eigs;
  
} d4est_solver_multigrid_smoother_cheby_t;


/* This file was automatically generated.  Do not edit! */
d4est_solver_multigrid_smoother_t *d4est_solver_multigrid_smoother_cheby_init(p4est_t *p4est,int num_of_levels,const char *input_file);
void d4est_solver_multigrid_smoother_cheby_iterate(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,double *r,int iter,double lmin,double lmax,int print_residual_norm,int mg_level);
void d4est_solver_multigrid_smoother_cheby_iterate_aux(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *d4est_ghost,d4est_ghost_data_t *d4est_ghost_data,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,double *r,int iter,double lmin,double lmax,int print_residual_norm,int mg_level,int compute_residual_at_end);
void d4est_solver_multigrid_smoother_cheby_destroy(d4est_solver_multigrid_smoother_t *smoother);

#endif
