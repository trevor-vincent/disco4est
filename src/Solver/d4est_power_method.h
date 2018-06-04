#ifndef D4EST_POWER_METHOD_H
#define D4EST_POWER_METHOD_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <krylov_pc.h>


/* This file was automatically generated.  Do not edit! */
double d4est_power_method(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,int imax);

#endif
