#ifndef D4EST_POISSON_H
#define D4EST_POISSON_H 

#include <pXest.h>
#include <d4est_util.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_xyz_functions.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>

/* This file was automatically generated.  Do not edit! */
void d4est_poisson_apply_aij(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_elliptic_data_t *prob_vecs,d4est_poisson_flux_data_t *flux_fcn_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);
void d4est_poisson_build_rhs_with_strong_bc(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_data_t *prob_vecs,d4est_poisson_flux_data_t *flux_fcn_data_for_build_rhs,double *rhs,d4est_xyz_fcn_t problem_rhs_fcn,void *ctx);

#endif
