#ifndef D4EST_LAPLACIAN_WITH_OPT_FLUX_SIPG_H
#define D4EST_LAPLACIAN_WITH_OPT_FLUX_SIPG_H 

#include <d4est_mesh.h>
#include <d4est_laplacian_with_opt_flux.h>
#include <d4est_laplacian_with_opt_flux_sipg.h>

typedef double
(*penalty_calc_t)
(
 int, //minus side degree
 double, //minus element size h
 int, //plus side degree
 double, //plus element size h
 double //penalty prefactor
);

typedef struct {
  double sipg_penalty_prefactor;

  int sipg_use_region_multipliers;
  int sipg_use_region_boundary_multipliers;
  
  int sipg_number_of_regions;
  int sipg_number_of_region_boundaries;
  /* For now cap at 10 */
  double sipg_region_multipliers [10];
  double sipg_region_boundary_multipliers [10];
  
  /* Given h and the order of the polynomial p, for the +,- sides sharing an interface
   * compute the penalty function */
  penalty_calc_t sipg_penalty_fcn;

  /* In case you want size_parameters from a different geometry (maybe good for stretched sphere grids?) */
  d4est_mesh_size_parameters_t* size_params;

  /* internal */
  int last_mortar_side_id_m;
  
} d4est_laplacian_with_opt_flux_sipg_params_t;


/* This file was automatically generated.  Do not edit! */
void d4est_laplacian_with_opt_flux_sipg_params_destroy(d4est_laplacian_with_opt_flux_data_t *data);
void d4est_laplacian_with_opt_flux_sipg_params_new(p4est_t *p4est,const char *print_prefix,const char *input_file,d4est_laplacian_with_opt_flux_data_t *d4est_laplacian_with_opt_flux_data);
void d4est_laplacian_with_opt_flux_sipg_params_input(p4est_t *p4est,const char *printf_prefix,const char *input_file,d4est_laplacian_with_opt_flux_sipg_params_t *input);
void d4est_laplacian_with_opt_flux_sipg_interface_aux(p4est_t *p4est,d4est_element_data_t **e_m,int faces_m,int f_m,int mortar_side_id_m,int *e_m_is_ghost,d4est_element_data_t **e_p,int faces_p,int f_p,int mortar_side_id_p,int *e_p_is_ghost,int orientation,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_laplacian_with_opt_flux_interface_data_t *mortar_data,void *params,double *restrict lifted_proj_VT_w_term1_mortar_lobatto_m,double *restrict proj_VT_w_term1_mortar_lobatto_m,double *restrict lifted_proj_VT_w_term1_mortar_lobatto_p,double *restrict proj_VT_w_term1_mortar_lobatto_p,double *restrict VT_w_term1_mortar_lobatto,double *restrict term1_mortar_quad,double *DT_lifted_proj_VT_w_term2_mortar_lobatto_m[(P4EST_DIM)],double *lifted_proj_VT_w_term2_mortar_lobatto_m[(P4EST_DIM)],double *proj_VT_w_term2_mortar_lobatto_m[(P4EST_DIM)],double *DT_lifted_proj_VT_w_term2_mortar_lobatto_p[(P4EST_DIM)],double *lifted_proj_VT_w_term2_mortar_lobatto_p[(P4EST_DIM)],double *proj_VT_w_term2_mortar_lobatto_p[(P4EST_DIM)],double *VT_w_term2_mortar_lobatto[(P4EST_DIM)],double *term2_mortar_quad[(P4EST_DIM)],double *restrict lifted_proj_VT_w_term3_mortar_lobatto_m,double *restrict proj_VT_w_term3_mortar_lobatto_m,double *restrict lifted_proj_VT_w_term3_mortar_lobatto_p,double *restrict proj_VT_w_term3_mortar_lobatto_p,double *restrict VT_w_term3_mortar_lobatto,double *restrict term3_mortar_quad);
void d4est_laplacian_with_opt_flux_sipg_robin_aux(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_laplacian_with_opt_flux_boundary_data_t *boundary_data,void *boundary_condition_fcn_data,void *flux_parameter_data,double *restrict term1_quad,double *restrict VT_w_term1_lobatto,double *restrict lifted_VT_w_term1_lobatto);
double d4est_laplacian_with_opt_flux_sipg_get_multiplier(d4est_laplacian_with_opt_flux_sipg_params_t *input,int region_m,int region_p);
void d4est_laplacian_with_opt_flux_sipg_dirichlet_aux(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_laplacian_with_opt_flux_boundary_data_t *boundary_data,void *boundary_condition_fcn_data,void *flux_parameter_data,double *restrict term1_quad,double *restrict VT_w_term1_lobatto,double *restrict lifted_VT_w_term1_lobatto,double *term2_quad[P4EST_DIM],double *VT_w_term2_lobatto[P4EST_DIM],double *lifted_VT_w_term2_lobatto[P4EST_DIM],double *DT_lifted_VT_w_term2_lobatto[P4EST_DIM],double *restrict term3_quad,double *restrict VT_w_term3_lobatto,double *restrict lifted_VT_w_term3_lobatto,double *restrict sigma,double *restrict u_at_bndry_lobatto,double *restrict u_m_on_f_m_min_u_at_bndry_quad);

#endif
