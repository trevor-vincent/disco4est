#ifndef CG_EIGS_H
#define CG_EIGS_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/d4est_elliptic_eqns.h"
#include "../dGMath/d4est_operators.h"

void cg_eigs(p4est_t *p4est,d4est_elliptic_problem_data_t *vecs,d4est_elliptic_eqns_t *fcns,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,int imax,double *eig_max);

#endif
