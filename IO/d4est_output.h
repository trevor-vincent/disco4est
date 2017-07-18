#ifndef D4EST_OUTPUT_H
#define D4EST_OUTPUT_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_mesh.h>
#include <d4est_estimator_stats.h>


typedef struct {

  double* u;
  double* error;
  double* eta2;
  double* jacobian;

} d4est_output_vtk_fields_t;

typedef enum {OUTPUT_ESTIMATOR, DO_NOT_OUTPUT_ESTIMATOR} d4est_output_estimator_option_t;

/* This file was automatically generated.  Do not edit! */
void d4est_output_vtk_with_analytic_error(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_data_t *prob_vecs,const char *input_file,const char *save_as_prefix,d4est_xyz_fcn_t analytic_solution,void *ctx,d4est_output_estimator_option_t eta2_option,int level);
void d4est_output_vtk(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *u,double *error,int local_nodes,const char *input_file,const char *save_as_prefix,d4est_output_estimator_option_t eta2_option,int level);
void d4est_output_norms_using_analytic_solution(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_estimator_stats_t *stats,d4est_elliptic_data_t *prob_vecs,d4est_xyz_fcn_t analytic_solution,void *ctx,int level);
void d4est_output_norms(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_estimator_stats_t *stats,double *error,int local_nodes,int level);
void d4est_output_calculate_analytic_error(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_data_t *prob_vecs,d4est_xyz_fcn_t analytic_solution,void *analytic_solution_ctx,double *error);

#endif
