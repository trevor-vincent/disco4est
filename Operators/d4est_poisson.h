#ifndef D4EST_POISSON_H
#define D4EST_POISSON_H 

#include <pXest.h>
#include <d4est_element_data.h>
#include <problem_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_poisson.h>
#include <grid_functions.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <util.h>


/* This file was automatically generated.  Do not edit! */
void d4est_poisson_apply_aij(p4est_t *p4est,p4est_ghost_t *ghost,void *ghost_data,d4est_elliptic_problem_data_t *prob_vecs,d4est_mortar_fcn_ptrs_t *flux_fcn_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);

#endif
