#ifndef D4EST_SOLVER_SCHWARZ_GEOMETRIC_DATA_H
#define D4EST_SOLVER_SCHWARZ_GEOMETRIC_DATA_H 

#include <pXest.h>
#include <d4est_ghost_data_ext.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_mesh.h>


typedef struct {
  
  int mortar_side_id;
  int is_ghost;
  
  int faces_m;
  int faces_p; /* equals 0 on physical boundary */
  int f_p;
  int f_m;
  int tree_p;
  int tree_m;
  int orientation;

  d4est_element_data_t e_m [(P4EST_HALF)];
  d4est_element_data_t e_p [(P4EST_HALF)];
  
  int boundary_quad_stride;
  int mortar_quad_stride;
  int total_mortar_nodes_quad;
  int total_mortar_nodes_lobatto;
  
} d4est_mortar_side_data_t;

typedef struct {

  const char* input_section;
  /* local mortar data */
  d4est_mortar_side_data_t* mortar_side_data;
  int* mortar_which_touches_face;

  /* ghost mortar data */
  d4est_ghost_data_ext_t* mortar_side_ghost_data;

  /* saved aliases to mortars for quick computation */
  int* num_of_mortars_per_subdomain;
  int* mortar_strides_per_subdomain;
  d4est_mortar_side_data_t** subdomain_mortars;
  int* zero_and_skip_m;
  int* zero_and_skip_p;
  int* nodal_stride_m;
  int* nodal_stride_p;
  /* int* skip_p_sum; */

  /* geometry data for ghost mortars */
  d4est_mesh_face_h_t face_h_type;
  double* drst_dxyz_m_mortar_quad;
  double* drst_dxyz_p_mortar_quad_porder;
  double* sj_m_mortar_quad;
  double* n_m_mortar_quad;
  double* hm_mortar_quad;
  double* hp_mortar_quad;
  double* xyz_on_f_m_quad;
  double* xyz_on_f_m_lobatto;

  /* volume geometry data for ghosts */
  double* J_quad_ghost;
  double* rst_xyz_quad_ghost;
  int* volume_quad_strides_per_ghost;
  int total_ghost_volume_quad_size;
  
  /* internal temporary variables, do not use */
  int boundary_quad_stride;
  int mortar_quad_stride;
  int mortar_side_stride;  
  
}d4est_solver_schwarz_geometric_data_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_geometric_data_sum_test(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data);
void d4est_solver_schwarz_geometric_data_check_hp(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data);
d4est_solver_schwarz_geometric_data_t *d4est_solver_schwarz_geometric_data_init(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,d4est_solver_schwarz_metadata_t *schwarz_metadata,const char *input_file,const char *input_section);
void d4est_solver_schwarz_geometric_data_reduce_to_minimal_set(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data);
void d4est_solver_schwarz_geometric_data_destroy(d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data);
void d4est_solver_schwarz_geometric_data_interface_callback(p4est_t *p4est,d4est_element_data_t **e_m,int faces_m,int f_m,int mortar_side_id_m,d4est_element_data_t **e_p,int faces_p,int f_p,int mortar_side_id_p,int *e_m_is_ghost,int orientation,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);
void d4est_solver_schwarz_geometric_data_boundary_callback(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);
void d4est_solver_schwarz_geometric_data_input(p4est_t *p4est,const char *input_file,const char *input_section,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data);

#endif
