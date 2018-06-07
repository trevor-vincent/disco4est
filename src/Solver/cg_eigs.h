#ifndef CG_EIGS_H
#define CG_EIGS_H 

#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>

void cg_eigs(p4est_t *p4est,d4est_elliptic_data_t *vecs,d4est_elliptic_eqns_t *fcns,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,int imax,int print_spectral_bound_iterations,double *spectral_bound);

#endif
