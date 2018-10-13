#ifndef CG_EIGS_H
#define CG_EIGS_H 

#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>

/* This file was automatically generated.  Do not edit! */
void cg_eigs(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,d4est_ghost_t *ghost,d4est_ghost_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,int imax,int print_spectral_bound_iterations,int use_new,double *spectral_bound);

#endif
