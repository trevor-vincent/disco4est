#ifndef D4EST_OUTPUT_H
#define D4EST_OUTPUT_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_mesh.h>


typedef struct {

  double* u;
  double* error;
  double* eta2;
  double* jacobian;

} d4est_output_vtk_fields_t;

typedef enum {SAVE_ESTIMATOR, DO_NOT_SAVE_ESTIMATOR} d4est_output_estimator_option_t;

/* This file was automatically generated.  Do not edit! */
void d4est_output_vtk(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,int level,double *u,double *error,int local_nodes,const char *input_file,const char *save_as_prefix,d4est_output_estimator_option_t eta2_option);
void d4est_output_error(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,int level,double *error,double local_estimator,int local_nodes);

#endif
