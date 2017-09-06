#ifndef D4EST_OUTPUT_H
#define D4EST_OUTPUT_H 

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
  
} d4est_output_energy_norm_fit_t;

typedef struct {

  double* u;
  double* u_compare;
  double* error;
  double* eta2;
  double* jacobian;

} d4est_output_vtk_dg_fields_t;


typedef struct {

  double* deg;
  double* eta2;

} d4est_output_vtk_cell_fields_t;


/* This file was automatically generated.  Do not edit! */
void d4est_output_vtk_with_analytic_error(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_data_t *prob_vecs,const char *input_file,const char *save_as_prefix,d4est_xyz_fcn_t analytic_solution,void *ctx,int save_estimator,int level);
void d4est_output_vtk_degree_mesh(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,const char *input_file,const char *save_as_prefix,int save_estimator,int level);
void d4est_output_vtk(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *u,double *u_compare,double *error,const char *input_file,const char *save_as_prefix,int local_nodes,int level,int save_estimator);
void d4est_output_norms_using_analytic_solution(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,double local_estimator,d4est_elliptic_data_t *prob_vecs,d4est_ip_energy_norm_data_t *energy_norm_data,d4est_xyz_fcn_t analytic_solution,void *ctx,d4est_output_energy_norm_fit_t *fit);
void d4est_output_norms(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_ip_energy_norm_data_t *energy_norm_data,double estimator,double *error,d4est_output_energy_norm_fit_t *fit);
void d4est_output_destroy_energy_norm_fit(d4est_output_energy_norm_fit_t *fit);
void d4est_output_energy_norm_add_entry_and_fit(d4est_output_energy_norm_fit_t *fit,double global_energy_norm_sqr,double global_dof);
void d4est_output_energy_norm_fit(d4est_output_energy_norm_fit_t *fit);
d4est_output_energy_norm_fit_t *d4est_output_new_energy_norm_fit(int num_of_entries);

#endif
