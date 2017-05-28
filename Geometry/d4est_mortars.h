#ifndef D4EST_MORTARS_H
#define D4EST_MORTARS_H 

#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

/* This file was automatically generated.  Do not edit! */
void d4est_geometry_compute_qcoords_on_mortar(p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int face,p4est_qcoord_t mortar_q0[(P4EST_HALF)][(P4EST_DIM)],p4est_qcoord_t *mortar_dq);
void d4est_mortars_compute_geometric_data_on_mortar(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int *deg_mortar_quad,int face,double *drst_dxyz_on_mortar_quad[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_quad,double *n_on_mortar_quad[(P4EST_DIM)],double *n_sj_on_mortar_quad[(P4EST_DIM)],double *j_div_sj_mortar_quad,normal_compute_method_t n_compute_method);


#endif
