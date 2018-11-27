#ifndef D4EST_SOLVER_SCHWARZ_LAPLACIAN_H
#define D4EST_SOLVER_SCHWARZ_LAPLACIAN_H 

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_laplacian_flux.h>

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_laplacian_apply_over_subdomain(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *d4est_ghost,d4est_solver_schwarz_metadata_t *schwarz_data,d4est_solver_schwarz_operators_t *schwarz_ops,d4est_laplacian_flux_data_t *flux_fcn_data,double *u_restricted_field_over_subdomain,double *Au_restricted_field_over_subdomain,int subdomain);

#endif
