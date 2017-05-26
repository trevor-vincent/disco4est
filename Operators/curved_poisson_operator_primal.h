#ifndef CURVED_POISSON_OPERATOR_PRIMAL_H
#define CURVED_POISSON_OPERATOR_PRIMAL_H 

#include <pXest.h>
#include <curved_element_data.h>
#include <problem_data.h>
#include <d4est_operators.h>
#include <linalg.h>
#include <grid_functions.h>
#include <d4est_quadrature.h>
#include <util.h>

/* This file was automatically generated.  Do not edit! */
void curved_poisson_operator_primal_apply_aij(p4est_t *p4est,p4est_ghost_t *ghost,void *ghost_data,problem_data_t *prob_vecs,d4est_operators_t *d4est_ops,d4est_geometry_t *geom,d4est_quadrature_t *d4est_quad);
void curved_poisson_operator_primal_compute_stiffmatrixterm(p4est_iter_volume_info_t *info,void *user_data);
void curved_poisson_operator_primal_init_vecs(p4est_iter_volume_info_t *info,void *user_data);

#endif
