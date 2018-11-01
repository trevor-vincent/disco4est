#ifndef D4EST_LAPLACIAN_WITH_OPT_H
#define D4EST_LAPLACIAN_WITH_OPT_H 

#include <pXest.h>
#include <d4est_mortars_with_opt.h>
#include <d4est_laplacian_aux.h>

typedef struct {

  d4est_quadrature_mortar_t* mortar_face_object;
  d4est_quadrature_mortar_t* mortar_face_object_porder;
  
  int total_side_nodes_m_lobatto;
  int total_side_nodes_p_lobatto;
  int total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad;
  int faces_mortar;
  
  double* u_m_on_f_m_mortar_quad;
  double* u_m_on_f_m;
  double* u_p_on_f_p;
  double* sj_on_f_m_mortar_quad;
  /* double* j_div_sj_on_f_m_mortar_quad; */
  double* u_p_on_f_p_mortar_quad;
  /* double* j_div_sj_on_f_p_mortar_quad; */
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];

  double* Au_m [(P4EST_FACES)];
  double* Au_p [(P4EST_FACES)];

  double* hm_mortar_quad;
  double* hp_mortar_quad;
  
  int* face_nodes_m_lobatto;
  int* face_nodes_p_lobatto;
  int* deg_mortar_quad;
  int* deg_mortar_quad_porder;
  int* nodes_mortar_quad;
  int* nodes_mortar_quad_porder;
  int* deg_mortar_lobatto;
  int* deg_mortar_lobatto_porder;
  int* nodes_mortar_lobatto;
  int* nodes_mortar_lobatto_porder;
  int* deg_m_lobatto;
  int* deg_p_lobatto;
  int* deg_p_lobatto_porder;

} d4est_laplacian_with_opt_flux_interface_data_t;



typedef struct {

  d4est_quadrature_mortar_t* face_object;

  double* h_quad;
  int face_nodes_m_lobatto;
  int face_nodes_m_quad;

  int deg_mortar_quad;
  double* xyz_on_f_m_lobatto [(P4EST_DIM)];
  double* u_m_on_f_m_quad;  
  double* u_m_on_f_m;  
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];
  double* n_on_f_m_quad [(P4EST_DIM)];
  double* xyz_on_f_m_quad [(P4EST_DIM)];
  double* sj_on_f_m_quad;
  /* double* j_div_sj_quad; */

  double* Au_m;
  
} d4est_laplacian_with_opt_flux_boundary_data_t;


typedef
void (*d4est_laplacian_with_opt_flux_interface_fcn_t)
(
 p4est_t*,
 d4est_element_data_t**,
 int,
 int,
 int,
 int*,
 d4est_element_data_t**,
 int,
 int,
 int,
 int*,
 int,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_mesh_data_t*,
 d4est_laplacian_with_opt_flux_interface_data_t*,
 void*
);

typedef
void (*d4est_laplacian_with_opt_flux_boundary_fcn_t)
(
 p4est_t*,
 d4est_element_data_t*,
 int,
 int,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_mesh_data_t*,
 d4est_laplacian_with_opt_flux_boundary_data_t*,
 void*,
 void*
);

typedef struct d4est_laplacian_with_opt_flux_data d4est_laplacian_with_opt_flux_data_t;

struct d4est_laplacian_with_opt_flux_data{

 
  d4est_laplacian_flux_type_t flux_type; 
  d4est_laplacian_with_opt_flux_interface_fcn_t interface_fcn;
  d4est_laplacian_with_opt_flux_boundary_fcn_t boundary_fcn;
  void* flux_data;

  int (*get_deg_mortar_quad)(d4est_element_data_t*, void*);
  void* get_deg_mortar_quad_ctx;

  int skip_p_side;
  int last_mortar_side_id_m;
  d4est_laplacian_bc_t bc_type;
  void* bc_data;

  /* internally set */
  d4est_ghost_t* d4est_ghost;
  d4est_ghost_data_t* d4est_ghost_data;
  p4est_t* p4est;
  double* u;
  double* dudr_local [(P4EST_DIM)];
  double* dudr_ghost [(P4EST_DIM)];
  double* Au;
  int local_nodes;
  int which_field;
  
  void (*destroy)(d4est_laplacian_with_opt_flux_data_t*);
  
};

typedef 
double (*d4est_laplacian_with_opt_mortar_xyz_fcn_t)(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
 int mortar_node
);

typedef struct {

  d4est_laplacian_with_opt_mortar_xyz_fcn_t robin_rhs;
  d4est_laplacian_with_opt_mortar_xyz_fcn_t robin_coeff;
  void* user;

} d4est_laplacian_with_opt_robin_bc_t;


/* typedef enum {EVAL_BNDRY_FCN_NOT_SET, EVAL_BNDRY_FCN_ON_QUAD, EVAL_BNDRY_FCN_ON_LOBATTO} dirichlet_bndry_eval_method_t; */


typedef struct {

  d4est_xyz_fcn_t dirichlet_fcn;
  dirichlet_bndry_eval_method_t eval_method;
  void* user;

} d4est_laplacian_with_opt_dirichlet_bc_t;

/* This file was automatically generated.  Do not edit! */
void d4est_laplacian_with_opt_flux_destroy(d4est_laplacian_with_opt_flux_data_t *data);
d4est_laplacian_with_opt_flux_data_t *d4est_laplacian_with_opt_flux_new(p4est_t *p4est,const char *input_file,d4est_laplacian_bc_t bc_type,void *bc_data);
d4est_mortars_with_opt_fcn_ptrs_t d4est_laplacian_with_opt_flux_fetch_fcns(d4est_laplacian_with_opt_flux_data_t *data);


#endif
