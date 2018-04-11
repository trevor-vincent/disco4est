#ifndef D4EST_MORTARS_AUX_H
#define D4EST_MORTARS_AUX_H 

#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_geometry.h>
#include <d4est_operators.h>

/* This file was automatically generated.  Do not edit! */
void d4est_mortars_compute_dxyz_drst_face_isoparametric(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int which_tree,int face,int deg,int deg_quad,double *xyz_rst_face_quad[(P4EST_DIM)][(P4EST_DIM)]);

#endif
