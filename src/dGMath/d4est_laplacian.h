#ifndef D4EST_POISSON_H
#define D4EST_POISSON_H 

#include <pXest.h>
#include <d4est_util.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_mesh.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_xyz_functions.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <d4est_ghost.h>
#include <d4est_ghost_data.h>

/* This file was automatically generated.  Do not edit! */
void d4est_laplacian_apply_mortar_matrices(p4est_t *p4est,d4est_ghost_t *ghost,d4est_laplacian_flux_data_t *flux_fcn_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors);
void d4est_laplacian_compute_dudr(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_ghost_data_t *d4est_ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,double *dudr_local[(P4EST_DIM)],double *dudr_ghost[(P4EST_DIM)],double *u,int local_nodes,int which_field);
void d4est_laplacian_apply_stiffness_matrix(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,double *u,double *Au,int local_nodes,int which_field);
void d4est_laplacian_apply_aij(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_ghost_data_t *d4est_ghost_data,d4est_elliptic_data_t *d4est_elliptic_data,d4est_laplacian_flux_data_t *flux_fcn_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,int which_field);
void d4est_laplacian_build_rhs_with_strong_bc(p4est_t *p4est,d4est_ghost_t *ghost,d4est_ghost_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_elliptic_data_t *prob_vecs,d4est_laplacian_flux_data_t *flux_fcn_data_for_build_rhs,double *rhs,d4est_xyz_fcn_t problem_rhs_fcn,d4est_mesh_init_field_option_t init_option,void *ctx,int which_field);


#endif
