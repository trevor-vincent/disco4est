#ifndef D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H
#define D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H 

#include <pXest.h>
#include <d4est_ghost_data_ext.h>


typedef struct {
  
  int mortar_side_id;
  
  int faces_m;
  int faces_p; /* equals 0 on physical boundary */
  int f_p;
  int f_m;
  int tree_p;
  int tree_m;
  int orientation;
  
  int q_p [(P4EST_HALF)][(P4EST_DIM)];
  int dq_p [(P4EST_HALF)];
  int tree_quadid_p [(P4EST_HALF)];
  int deg_p [(P4EST_HALF)];
  int deg_quad_p [(P4EST_HALF)];

  int q_m [(P4EST_HALF)][(P4EST_DIM)];
  int dq_m [(P4EST_HALF)];
  int tree_quadid_m [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int deg_quad_m [(P4EST_HALF)];

  int total_mortar_nodes_quad;
  /* deprecated soon for the ones below */
  int boundary_quad_vector_stride;
  int mortar_quad_scalar_stride;
  int mortar_quad_vector_stride;
  int mortar_quad_matrix_stride;

  /* the ones above are redundant and can be reduced to */
  /* int boundary_quad_stride; */
  /* int mortar_quad_stride; */
  
} d4est_mortar_side_data_t;

typedef struct {

  double* drst_dxyz_m_mortar_quad;
  double* drst_dxyz_p_mortar_quad_porder;
  double* sj_m_mortar_quad;
  double* n_m_mortar_quad;
  double* hm_mortar_quad;
  double* hp_mortar_quad;
  double* xyz_on_f_m_quad;
  double* xyz_on_f_m_lobatto;

  int stride;  
  d4est_mortar_side_data_t* mortar_side_data;
  d4est_ghost_data_ext_t* mortar_side_ghost_data;

  int* which_mortar_face_belongs_to;
  
}d4est_solver_schwarz_mortar_data_t;


d4est_solver_schwarz_mortar_data_t *d4est_solver_schwarz_mortar_data_init(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *d4est_ghost,d4est_mesh_face_h_t face_h_type);
void d4est_solver_schwarz_mortar_data_destroy(d4est_solver_schwarz_mortar_data_t *schwarz_mortar_data);
void d4est_solver_schwarz_mortar_data_interface_callback(p4est_t *p4est,d4est_element_data_t **e_m,int faces_m,int f_m,int mortar_side_id_m,d4est_element_data_t **e_p,int faces_p,int f_p,int mortar_side_id_p,int *e_m_is_ghost,int orientation,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);
void d4est_solver_schwarz_mortar_data_boundary_callback(p4est_t *p4est,d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);

#endif
