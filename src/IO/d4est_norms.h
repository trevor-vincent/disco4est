#ifndef D4EST_NORMS
#define D4EST_NORMS_H

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_ip_energy_norm.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_mesh.h>
#include <d4est_estimator_stats.h>


typedef struct {

  double* log_norm_data;
  double* dof_data;
  int num_of_data_entries;
  int stride;
  
} d4est_norms_linear_fit_t;


typedef struct {

  p4est_t *p4est;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_data_t* d4est_factors;

  
} d4est_norms_fcn_L2_ctx_t;

typedef struct {

  p4est_t *p4est;
  d4est_ghost_t* ghost;
  d4est_ghost_data_t* ghost_data;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_data_t* d4est_factors;
  int which_field;
  
  d4est_ip_energy_norm_data_t* energy_norm_data;
  double energy_estimator_sq_local;
  double* energy_estimator;
  
} d4est_norms_fcn_energy_ctx_t;


typedef double
(*d4est_norm_fcn_t)
(
  p4est_t* p4est,
  double *field_value_errors,
  int num_nodes_local,
  void *ctx,
  int (*skip_element_fcn)(d4est_element_data_t*)
);


/* This file was automatically generated.  Do not edit! */
void d4est_norms_save(p4est_t *p4est,d4est_mesh_data_t *d4est_factors,const char **field_names,double **field_values,double **field_values_compare,d4est_xyz_fcn_t *analytical_solutions,void **analytical_solution_ctxs,const char **norm_names,d4est_norm_fcn_t *norm_fcns,void **norm_fcn_ctxs,d4est_norms_linear_fit_t **linear_fits,const char *zlog_output_category_prefix,int(*skip_element_fcn)(d4est_element_data_t *));
void d4est_norms_write_headers(const char **field_names,const char **norm_names,const char *zlog_output_category_prefix);
double d4est_norms_fcn_energy_estimator(p4est_t *p4est,double *field_value_errors,int num_nodes_local,void *norm_fcn_ctx,int(*skip_element_fcn)(d4est_element_data_t *));
double d4est_norms_fcn_energy(p4est_t *p4est,double *field_value_errors,int num_nodes_local,void *norm_fcn_ctx,int(*skip_element_fcn)(d4est_element_data_t *));
void d4est_norms_linear_fit_destroy(d4est_norms_linear_fit_t *fit);
void d4est_norms_linear_fit_add_entry_and_fit(p4est_t *p4est,d4est_norms_linear_fit_t *fit,double global_norm_sqr,double global_dof);
void d4est_norms_linear_fit(p4est_t *p4est,d4est_norms_linear_fit_t *fit);
d4est_norms_linear_fit_t *d4est_norms_linear_fit_init();
double d4est_norms_fcn_Linfty(p4est_t *p4est,double *field_value_errors,int num_nodes_local,void *norm_fcn_ctx,int(*skip_element_fcn)(d4est_element_data_t *));
double d4est_norms_fcn_L2(p4est_t *p4est,double *field_value_errors,int num_nodes_local,void *norm_fcn_ctx,int(*skip_element_fcn)(d4est_element_data_t *));
 
#endif
 
