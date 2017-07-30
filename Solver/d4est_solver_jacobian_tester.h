#ifndef D4EST_ANALYTIC_JACOBIAN_TESTER_H
#define D4EST_ANAKYTIC_JACOBIAN_TESTER_H 

#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>


typedef enum {JAC_TEST_CENTRAL_DIFFERENCE, JAC_TEST_FORWARD_DIFFERENCE} d4est_analytic_jacobian_tester_type_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_jacobian_tester(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_elliptic_eqns_t *elliptic_eqns,int local_nodes,d4est_xyz_fcn_t initial_guess,void *initial_guess_ctx,double eps,d4est_analytic_jacobian_tester_type_t type,int num_vecs_to_try);

#endif
