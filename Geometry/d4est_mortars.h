#ifndef D4EST_MORTARS_H
#define D4EST_MORTARS_H 

#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

/* This file was automatically generated.  Do not edit! */
void d4est_mortars_project_side_onto_mortar_space(d4est_operators_t *d4est_ops,double *in_side,int faces_side,int *deg_side,double *out_mortar,int faces_mortar,int *deg_mortar);
void d4est_mortars_project_mass_mortar_onto_side(d4est_operators_t *dgmath,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_mortars_project_mortar_onto_side(d4est_operators_t *d4est_ops,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_geometry_compute_rst_on_mortar(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_quadrature_integrand_type_t integrand_type,d4est_rst_t *d4est_rst,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int face,int mortar_side_id,int *deg_mortar_quad);
void d4est_mortars_compute_geometric_data_on_mortar(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_quadrature_integrand_type_t integrand_type,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int mortar_side_id,int num_faces_side,int num_faces_mortar,int *deg_mortar_quad,int face,double *drst_dxyz_on_mortar_quad[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_quad,double *n_on_mortar_quad[(P4EST_DIM)],double *n_sj_on_mortar_quad[(P4EST_DIM)],double *j_div_sj_mortar_quad,normal_compute_method_t n_compute_method);
void d4est_geometry_compute_qcoords_on_mortar(p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int face,p4est_qcoord_t mortar_q0[(P4EST_HALF)][(P4EST_DIM)],p4est_qcoord_t *mortar_dq);
void d4est_mortars_compute_geometric_data_on_mortar_aux(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,d4est_rst_t *rst_quad_mortar,int num_faces_side,int num_faces_mortar,int *deg_mortar_quad,int face,double *drst_dxyz_on_mortar_quad[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_quad,double *n_on_mortar_quad[(P4EST_DIM)],double *n_sj_on_mortar_quad[(P4EST_DIM)],double *j_div_sj_mortar_quad,normal_compute_method_t n_compute_method);



#endif
