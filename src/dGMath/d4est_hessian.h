#ifndef D4EST_HESSIAN_H
#define D4EST_HESSIAN_H 

#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_geometry.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_util.h>
#include <sc_reduce.h>
#include <d4est_linalg.h>

typedef enum {HESSIAN_ANALYTICAL, HESSIAN_NUMERICAL} d4est_hessian_compute_method_t;

/* This file was automatically generated.  Do not edit! */
void d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_hessian_compute_method_t compute_method,double *field_lobatto,double *del2field);
void d4est_hessian_compute_d2xdrdr_on_quadrature_points_of_element(d4est_element_data_t *ed,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_hessian_compute_method_t compute_method,double *d2xdrdr_quad[(P4EST_DIM)][(P4EST_DIM)][(P4EST_DIM)]);
void d4est_hessian(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,void *object,d4est_quadrature_object_type_t object_type,d4est_quadrature_integrand_type_t integrand_type,double *u,int deg_lobatto,double *hessian_u_quad[(P4EST_DIM)][(P4EST_DIM)],double *hessian_trace_u_quad,int deg_quad,double *drst_dxyz_quad[(P4EST_DIM)][(P4EST_DIM)],double *d2xyz_drstdrst_quad[(P4EST_DIM)][(P4EST_DIM)][P4EST_DIM]);

#endif
