#ifndef D4EST_ESTIMATOR_BI_NEW_H
#define D4EST_ESTIMATOR_BI_NEW_H 


#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_estimator_bi.h>
#include <d4est_laplacian_flux_sipg.h>

typedef struct {

  double* estimator;
  int output_vtk;
  double* estimator_vtk; /* residual term, je1 term, je2 term, je2 boundary term */
  double* estimator_vtk_per_face; /* residual term, je1 term, je2 term, je2 boundary term */
  d4est_mesh_size_parameters_t* size_params;
  d4est_mesh_data_t* d4est_factors_compactified;
  penalty_calc_t u_penalty_fcn;
  penalty_calc_t u_dirichlet_penalty_fcn;
  penalty_calc_t gradu_penalty_fcn;
  double penalty_prefactor;
  void* user;
  
} d4est_estimator_bi_new_penalty_data_t;


/* This file was automatically generated.  Do not edit! */
double *d4est_estimator_bi_new_compute(p4est_t *p4est,d4est_elliptic_data_t *d4est_elliptic_data,double *pointwise_residual_on_physical_quadrature_points,d4est_estimator_bi_new_penalty_data_t penalty_data,d4est_xyz_fcn_t u_bndry_fcn,void *bndry_ctx,d4est_ghost_t *d4est_ghost,d4est_ghost_data_t *d4est_ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom_physical,d4est_mesh_data_t *d4est_factors_physical,d4est_geometry_t *d4est_geom_compactified,d4est_mesh_data_t *d4est_factors_compactified,d4est_quadrature_t *d4est_quad,int which_field,double *estimator_vtk,double *estimator_vtk_per_face,int use_pointwise_residual,d4est_elliptic_eqns_t *fcns);

#endif
