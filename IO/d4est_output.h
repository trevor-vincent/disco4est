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
void d4est_output_vtk_with_analytic_error(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_data_t *prob_vecs,const char *input_file,const char *save_as_prefix,d4est_xyz_fcn_t analytic_solution,void *ctx,int level);
void d4est_output_vtk_degree_mesh(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,const char *input_file,const char *save_as_prefix,int save_estimator,int level);
void d4est_output_vtk(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *u,double *u_compare,double *error,const char *input_file,const char *save_as_prefix,int local_nodes,int level,int save_estimator);
void d4est_output_norms_using_analytic_solution(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_estimator_stats_t *stats,d4est_elliptic_data_t *prob_vecs,d4est_xyz_fcn_t analytic_solution,void *ctx);

#endif
