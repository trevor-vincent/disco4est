#ifndef D4EST_POISSON_SIPG_FLUX_H
#define D4EST_POISSON_SIPG_FLUX_H 

#include <pXest.h>
#include <d4est_mortars.h>
#include <d4est_poisson_flux.h>

typedef double
(*penalty_calc_t)
(
 int, //minus side degree
 double, //minus element size h
 int, //plus side degree
 double, //plus element size h
 double //penalty prefactor
);

typedef enum { H_EQ_J_DIV_SJ_QUAD, H_EQ_J_DIV_SJ_MIN_LOBATTO, H_EQ_J_DIV_SJ_MEAN_LOBATTO, H_EQ_TREE_H, H_EQ_VOLUME_DIV_AREA, H_EQ_FACE_DIAM, H_EQ_TOTAL_VOLUME_DIV_TOTAL_AREA, H_EQ_NOTSET } h_calc_method_t;

typedef struct {
  double sipg_penalty_prefactor;

  /* Given h and the order of the polynomial p, for the +,- sides sharing an interface
   * compute the penalty function */
  penalty_calc_t sipg_penalty_fcn;

  /* For curved elements, h is not well defined, therefore we try different methods  */
  h_calc_method_t sipg_flux_h;

  /* In case you want size_parameters from a different geometry (maybe good for stretched sphere grids?) */
  d4est_mesh_size_parameters_t* size_params;
  
} d4est_poisson_flux_sipg_params_t;


/* This file was automatically generated.  Do not edit! */
void d4est_poisson_flux_sipg_params_destroy(d4est_poisson_flux_data_t *data);
void d4est_poisson_flux_sipg_params_new(p4est_t *p4est,const char *print_prefix,const char *input_file,d4est_poisson_flux_data_t *d4est_poisson_flux_data);
void d4est_poisson_flux_sipg_params_input(p4est_t *p4est,const char *printf_prefix,const char *input_file,d4est_poisson_flux_sipg_params_t *input);
void d4est_poisson_flux_sipg_interface_aux(p4est_t *p4est,d4est_element_data_t **e_m,int faces_m,int f_m,int mortar_side_id_m,d4est_element_data_t **e_p,int faces_p,int f_p,int mortar_side_id_p,int *e_m_is_ghost,int orientation,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_poisson_flux_interface_data_t *mortar_data,void *params,double *lifted_proj_VT_w_term1_mortar_lobatto,double *proj_VT_w_term1_mortar_lobatto,double *VT_w_term1_mortar_lobatto,double *term1_mortar_quad,double *DT_lifted_proj_VT_w_term2_mortar_lobatto[(P4EST_DIM)],double *lifted_proj_VT_w_term2_mortar_lobatto[(P4EST_DIM)],double *proj_VT_w_term2_mortar_lobatto[(P4EST_DIM)],double *VT_w_term2_mortar_lobatto[(P4EST_DIM)],double *term2_mortar_quad[(P4EST_DIM)],double *lifted_proj_VT_w_term3_mortar_lobatto,double *proj_VT_w_term3_mortar_lobatto,double *VT_w_term3_mortar_lobatto,double *term3_mortar_quad);
void d4est_poisson_flux_sipg_robin_aux(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_poisson_flux_boundary_data_t *boundary_data,void *boundary_condition_fcn_data,void *flux_parameter_data,double *term1_quad,double *VT_w_term1_lobatto,double *lifted_VT_w_term1_lobatto);
void d4est_poisson_flux_sipg_dirichlet_aux(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_poisson_flux_boundary_data_t *boundary_data,void *boundary_condition_fcn_data,void *flux_parameter_data,double *term1_quad,double *VT_w_term1_lobatto,double *lifted_VT_w_term1_lobatto,double *term2_quad[P4EST_DIM],double *VT_w_term2_lobatto[P4EST_DIM],double *lifted_VT_w_term2_lobatto[P4EST_DIM],double *DT_lifted_VT_w_term2_lobatto[P4EST_DIM],double *term3_quad,double *VT_w_term3_lobatto,double *lifted_VT_w_term3_lobatto,double *sigma,double *u_at_bndry_lobatto,double *u_at_bndry_lobatto_to_quad);
void d4est_poisson_flux_sipg_calculate_h(p4est_t *p4est,d4est_element_data_t **elems_side,int face_side,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,h_calc_method_t sipg_flux_h,double *j_div_sj_quad_mortar,double *h_quad_mortar,int num_faces_mortar,int num_faces_side,int *nodes_mortar_quad,d4est_mesh_size_parameters_t size_params);

#endif
