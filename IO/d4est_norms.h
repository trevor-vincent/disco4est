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

  double* log_energy_norm_data;
  double* dof_data;
  int num_of_data_entries;
  int stride;
  
} d4est_norms_energy_norm_fit_t;


typedef struct {

  p4est_t *p4est;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  
} d4est_norms_L2_ctx_t;

typedef struct {

  p4est_t *p4est;
  p4est_ghost_t* ghost;
  d4est_element_data_t* ghost_data;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_data_t* d4est_factors;

  d4est_ip_energy_norm_data_t* energy_norm_data;
  double energy_estimator_sq_local;
  d4est_norms_energy_norm_fit_t* fit;

} d4est_norms_energy_ctx_t;


typedef double
(*d4est_norm_fcn_t)
(
  double *field_value_errors,
  int num_nodes_local,
  void *ctx
);


typedef struct {

  int global_nodes;
  int global_num_quadrants;
  int avg_deg;
  int global_quad_nodes;
  int avg_deg_quad;
  double global_estimator;
  double global_l2_norm_sqr;
  double global_linf;
  
} d4est_norms_norms_t;


/* This file was automatically generated.  Do not edit! */
double d4est_norms_L2(double *field_value_errors, int num_nodes_local, void *norm_fcn_ctx);
double d4est_norms_Linfty(double *field_value_errors, int num_nodes_local, void *norm_fcn_ctx);
double d4est_norms_energy(double *field_value_errors, int num_nodes_local, void *norm_fcn_ctx);
double d4est_norms_energy_estimator(double *field_value_errors, int num_nodes_local, void *norm_fcn_ctx);

void d4est_norms_norms_using_analytic_solution(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,double local_estimator,d4est_elliptic_data_t *prob_vecs,d4est_ip_energy_norm_data_t *energy_norm_data,d4est_xyz_fcn_t analytic_solution,void *ctx,d4est_norms_energy_norm_fit_t *fit,int(*skip_element_fcn)(d4est_element_data_t *));
d4est_norms_norms_t d4est_norms_norms(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_ip_energy_norm_data_t *energy_norm_data,double estimator,double *error,d4est_norms_energy_norm_fit_t *fit,int(*skip_element_fcn)(d4est_element_data_t *));
void d4est_norms_destroy_energy_norm_fit(d4est_norms_energy_norm_fit_t *fit);
void d4est_norms_energy_norm_add_entry_and_fit(p4est_t *p4est,d4est_norms_energy_norm_fit_t *fit,double global_energy_norm_sqr,double global_dof);
void d4est_norms_energy_norm_fit(p4est_t *p4est,d4est_norms_energy_norm_fit_t *fit);
void d4est_norms_write_header();
void d4est_norms_write_headers(const char** field_names,const char** norm_names);
void d4est_norms_save(p4est_t *p4est,const char** field_names,double** field_values,double** field_values_compare,d4est_xyz_fcn_t *analytical_solution,void** analytical_solution_ctxs,const char** norm_names,d4est_norm_fcn_t *norm_fcns,void **norm_fcn_ctxs);

#endif
