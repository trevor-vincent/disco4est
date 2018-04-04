#ifndef D4EST_MESH_H
#define D4EST_MESH_H

#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

typedef enum {INITIALIZE_QUADRATURE_DATA, DO_NOT_INITIALIZE_QUADRATURE_DATA} d4est_mesh_quadrature_data_init_option_t;
typedef enum {INITIALIZE_GEOMETRY_DATA, DO_NOT_INITIALIZE_GEOMETRY_DATA} d4est_mesh_geometry_data_init_option_t;
typedef enum {INITIALIZE_GEOMETRY_ALIASES,              DO_NOT_INITIALIZE_GEOMETRY_ALIASES}d4est_mesh_geometry_aliases_init_option_t;
typedef enum {DISCARD_BOUNDARY, DISCARD_INTERIOR, DISCARD_NOTHING} d4est_mesh_boundary_option_t;
typedef enum {INIT_FIELD_NOT_SET, INIT_FIELD_ON_LOBATTO, INIT_FIELD_ON_QUAD} d4est_mesh_init_field_option_t;
typedef enum {PRINT, DO_NOT_PRINT, PRINT_ON_ERROR} d4est_mesh_print_option_t;

typedef struct {
  
  int local_nodes;
  int local_sqr_nodes;
  int local_mortar_nodes_quad;
  int local_boundary_nodes_quad;
  int local_nodes_quad;
  
} d4est_mesh_local_sizes_t;


typedef struct {

  /* error flag */
  int err;
  
  /* point data */
  double f_at_xyz;
  int nodal_stride;
  double xyz [3];
  double rst [3];
  double abc [3];
  
  /* element data */
  int id;
  p4est_qcoord_t q [3];
  p4est_qcoord_t dq;

} d4est_mesh_interpolate_data_t;

typedef struct {
  
  double* J_quad;
  double* xyz_quad [(P4EST_DIM)];
  double* xyz_rst_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* rst_xyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* xyz [(P4EST_DIM)];
  
} d4est_mesh_data_on_element_t;

typedef struct {

  double* j_div_sj_min;
  double* j_div_sj_mean;
  double* diam_face;
  double* diam_volume;
  double* area;
  double* volume;
  
} d4est_mesh_size_parameters_t;

typedef struct {

  d4est_mesh_local_sizes_t local_sizes;
  
  double* xyz;
  double* xyz_quad;
  double* xyz_rst_quad;
  double* J_quad;
  double* rst_xyz_quad;
  
  double* drst_dxyz_m_mortar_quad;
  double* drst_dxyz_p_mortar_quad_porder;

  double* j_div_sj_m_mortar_quad;
  double* j_div_sj_p_mortar_quad_porder;
  
  double* sj_m_mortar_quad;
  double* n_m_mortar_quad;
  double* xyz_m_mortar_quad;
  double* xyz_m_mortar_lobatto;

  double* j_div_sj_min;
  double* j_div_sj_mean;
  double* diam_face;
  double* diam_volume;
  double* area;
  double* volume;

} d4est_mesh_data_t;

typedef struct {

  int min_quadrants;
  int min_level;
  int fill_uniform;
  int initial_nodes;

  int number_of_regions;
  int* deg;
  int* deg_quad_inc;

  int load_from_checkpoint;
  int* checkpoint_deg_array;
  char* checkpoint_prefix;
  
} d4est_mesh_initial_extents_t;

/* This file was automatically generated.  Do not edit! */
d4est_mesh_interpolate_data_t d4est_mesh_interpolate_at_tree_coord(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double abc[(P4EST_DIM)],int tree_id,double *f,int print);
double d4est_mesh_volume_integral(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,int(*is_it_in_volume)(d4est_element_data_t *,void *),double(*compute_volume_integral)(d4est_element_data_t *,void *),void *user);
double d4est_mesh_surface_integral(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,int(*is_it_on_surface)(d4est_element_data_t *,int,void *),double(*compute_face_integral)(d4est_element_data_t *,int,double *,double *[(P4EST_DIM)],double *[(P4EST_DIM)],double *[(P4EST_DIM)][(P4EST_DIM)],void *),void *user);
double d4est_mesh_compare_two_fields(p4est_t *p4est,d4est_mesh_data_t *d4est_factors,double *field1,double *field2,const char *msg,d4est_mesh_boundary_option_t boundary_option,d4est_mesh_print_option_t print_option,double eps);
int d4est_mesh_get_local_nodes(p4est_t *p4est);
int d4est_mesh_get_ghost_nodes(p4est_ghost_t *ghost,d4est_element_data_t *ghost_data);
void d4est_mesh_get_local_nodes_callback(p4est_iter_volume_info_t *info,void *user_data);
void d4est_mesh_init_field_ext(p4est_t *p4est,double *node_vec,d4est_xyz_fcn_ext_t xyz_fcn,void *user,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom);
void d4est_mesh_compute_point_error(double *v1,double *v2,double *error,int local_nodes);
void d4est_mesh_init_field(p4est_t *p4est,double *node_vec,d4est_xyz_fcn_t init_fcn,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_mesh_data_t *d4est_factors,d4est_mesh_init_field_option_t option,void *user);
int d4est_mesh_update(p4est_t *p4est,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_mesh_quadrature_data_init_option_t quad_init_option,d4est_mesh_geometry_data_init_option_t geom_init_option,d4est_mesh_geometry_aliases_init_option_t alias_init_option,void(*user_fcn)(d4est_element_data_t *,void *),void *user_ctx);
void d4est_mesh_data_compute(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors);
d4est_mesh_local_sizes_t d4est_mesh_init_element_data(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void(*user_fcn)(d4est_element_data_t *,void *),void *user_ctx);
d4est_mesh_size_parameters_t d4est_mesh_get_size_parameters(d4est_mesh_data_t *factors);
double *d4est_mesh_get_jacobian_on_quadrature_points(d4est_mesh_data_t *d4est_factors,d4est_element_data_t *ed);
double d4est_mesh_compute_l2_norm_sqr(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,double *nodal_vec,int local_nodes,int(*skip_element_fcn)(d4est_element_data_t *),double *l2_array);
int d4est_mesh_debug_find_node(p4est_t *p4est,int node);
int d4est_mesh_global_node_to_local_node(p4est_t *p4est,int global_node);
d4est_element_data_t *d4est_mesh_get_element_data(p4est_t *p4est,int local_element_id);
void d4est_mesh_print_number_of_elements_per_tree(p4est_t *p4est);
int d4est_mesh_get_local_quad_nodes(p4est_t *p4est);
int d4est_mesh_get_local_matrix_nodes(p4est_t *p4est);
void d4est_mesh_compute_jacobian_on_lgl_grid(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *jacobian_lgl);
void d4est_mesh_get_array_of_quadrature_degrees(p4est_t *p4est,void *deg_array,d4est_builtin_t type);
void d4est_mesh_get_array_of_degrees(p4est_t *p4est,void *deg_array,d4est_builtin_t type);
void d4est_mesh_data_destroy(d4est_mesh_data_t *d4est_factors);
void d4est_mesh_data_printout(d4est_mesh_data_t *d4est_factors);
void d4est_mesh_data_realloc(p4est_t *p4est,d4est_mesh_data_t *d4est_factors,d4est_mesh_local_sizes_t local_sizes);
d4est_mesh_data_t *d4est_mesh_data_init();
void d4est_mesh_compute_mortar_quadrature_sizes(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_mesh_local_sizes_t *local_sizes);
void d4est_mesh_compute_mortar_quadrature_quantities(p4est_t *p4est,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors);
void d4est_mesh_compute_mortar_quadrature_quantities_interface_callback(d4est_element_data_t **e_m,int faces_m,int f_m,int mortar_side_id_m,d4est_element_data_t **e_p,int faces_p,int f_p,int mortar_side_id_p,int *e_m_is_ghost,int orientation,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);
d4est_mesh_data_on_element_t d4est_mesh_data_on_element(d4est_mesh_data_t *d4est_factors,d4est_element_data_t *ed);
void d4est_mesh_compute_mortar_quadrature_quantities_boundary_callback(d4est_element_data_t *e_m,int f_m,int mortar_side_id_m,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,void *params);
d4est_mesh_initial_extents_t *d4est_mesh_initial_extents_parse(const char *input_file,d4est_geometry_t *d4est_geom);
void d4est_mesh_initial_extents_destroy(d4est_mesh_initial_extents_t *initial_extents);
void d4est_mesh_set_quadratures_after_amr(d4est_element_data_t *elem_data,void *user_ctx);
void d4est_mesh_set_initial_extents(d4est_element_data_t *elem_data,void *user_ctx);

#endif



