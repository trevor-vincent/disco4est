#ifndef D4EST_MORTARS_H
#define D4EST_MORTARS_H 

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_mesh.h>

/* typedef enum */
/*   { */
/*     EXCHANGE_GHOST_DATA, */
/*     DO_NOT_EXCHANGE_GHOST_DATA */
/*   } d4est_mortars_exchange_data_option_t; */




/* No problem specific data needed */
typedef void (*d4est_mortars_interface_fcn_t)
(
 p4est_t* p4est,
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user_ctx
);

/* We need the problem specific data */
typedef void (*d4est_mortars_boundary_fcn_t)
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user_ctx
);

typedef struct {

  d4est_mortars_interface_fcn_t flux_interface_fcn;
  d4est_mortars_boundary_fcn_t flux_boundary_fcn;
  void* user_ctx;
  
} d4est_mortars_fcn_ptrs_t;

typedef struct {

  d4est_operators_t* d4est_ops;
  d4est_quadrature_t* d4est_quad;
  d4est_mortars_fcn_ptrs_t* flux_fcn_ptrs;
  d4est_mesh_data_t* d4est_factors;
  d4est_geometry_t* geom;
  int mortar_stride;

} d4est_mortars_compute_flux_user_data_t;


/* This file was automatically generated.  Do not edit! */
double *d4est_mortars_reorient_if_needed(d4est_operators_t *d4est_ops,d4est_element_data_t **e_m,d4est_element_data_t **e_p,int orientation,int f_m,int f_p,double *field_porder,int *deg_morder,int *deg_porder,int total_nodes,int faces);
int d4est_mortars_compute_flux_on_local_elements(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_mortars_fcn_ptrs_t *fcn_ptrs);
void d4est_mortars_compute_flux_on_local_elements_aux(p4est_iter_face_info_t *info,void *user_data);
void d4est_mortars_project_side_onto_mortar_space(d4est_operators_t *d4est_ops,double *in_side,int faces_side,int *deg_side,double *out_mortar,int faces_mortar,int *deg_mortar);
void d4est_mortars_project_mass_mortar_onto_side(d4est_operators_t *dgmath,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_mortars_project_mortar_onto_side(d4est_operators_t *d4est_ops,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_mortars_compute_rst_on_mortar(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_quadrature_integrand_type_t integrand_type,d4est_rst_t *d4est_rst,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int face,int mortar_side_id,int *deg_mortar_quad);
void d4est_mortars_compute_geometric_data_on_mortar(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_quadrature_integrand_type_t integrand_type,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int mortar_side_id,int num_faces_side,int num_faces_mortar,int *deg_mortar_lobatto,int *deg_mortar_quad,int face,double *xyz_on_mortar_quad[(P4EST_DIM)],double *drst_dxyz_on_mortar_quad[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_quad,double *n_on_mortar_quad[(P4EST_DIM)],double *n_sj_on_mortar_quad[(P4EST_DIM)],double *j_div_sj_mortar_quad,normal_compute_method_t n_compute_method);
void d4est_mortars_compute_qcoords_on_mortar(p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int face,p4est_qcoord_t mortar_q0[(P4EST_HALF)][(P4EST_DIM)],p4est_qcoord_t *mortar_dq);
void d4est_mortars_compute_geometric_data_on_mortar_aux(d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,d4est_rst_t *rst_quad_mortar,int num_faces_side,int num_faces_mortar,int *deg_mortar_lobatto,int *deg_mortar_quad,int face,double *xyz_on_mortar_quad[(P4EST_DIM)],double *drst_dxyz_on_mortar_quad[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_quad,double *n_on_mortar_quad[(P4EST_DIM)],double *n_sj_on_mortar_quad[(P4EST_DIM)],double *j_div_sj_mortar_quad,normal_compute_method_t n_compute_method);

#endif
