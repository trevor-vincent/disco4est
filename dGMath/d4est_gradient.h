#ifndef D4EST_GRADIENT_H
#define D4EST_GRADIENT_H 

#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>


/* This file was automatically generated.  Do not edit! */
double d4est_gradient_l2_norm(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,void *object,d4est_quadrature_object_type_t object_type,d4est_quadrature_integrand_type_t integrand_type,double *u,int deg_lobatto,int deg_quad,double *rst_xyz_quad[(P4EST_DIM)][(P4EST_DIM)],double *jac_quad);
void d4est_gradient(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,void *object,d4est_quadrature_object_type_t object_type,d4est_quadrature_integrand_type_t integrand_type,double *u,int deg_lobatto,double *grad_u_quad[(P4EST_DIM)],int deg_quad,double *rst_xyz_quad[(P4EST_DIM)][(P4EST_DIM)]);



#endif
