#ifndef D4EST_MESH_H
#define D4EST_MESH_H 

#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_xyz_functions_ext.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

typedef enum { DO_NOT_STORE_LOCALLY, STORE_LOCALLY } norm_storage_option_t;
typedef enum {INITIALIZE_QUADRATURE_DATA, DO_NOT_INITIALIZE_QUADRATURE_DATA} d4est_mesh_quadrature_data_init_option_t;
typedef enum {INITIALIZE_GEOMETRY_DATA, DO_NOT_INITIALIZE_GEOMETRY_DATA} d4est_mesh_geometry_data_init_option_t;
typedef enum {INITIALIZE_GEOMETRY_ALIASES, DO_NOT_INITIALIZE_GEOMETRY_ALIASES}d4est_mesh_geometry_aliases_init_option_t;

typedef enum {DISCARD_BOUNDARY, DISCARD_INTERIOR, DISCARD_NOTHING} d4est_mesh_boundary_option_t;


typedef enum {PRINT, DO_NOT_PRINT, PRINT_ON_ERROR} d4est_mesh_print_option_t;

typedef struct {
  
  int local_nodes;
  int local_sqr_nodes;
  int local_sqr_mortar_nodes;
  int local_nodes_quad;
  
} d4est_local_sizes_t;

typedef struct {

  double* xyz;
  double* xyz_quad;
  double* xyz_rst_quad;
  double* J_quad;
  double* rst_xyz_quad;

} d4est_mesh_geometry_storage_t;

/* This file was automatically generated.  Do not edit! */
double d4est_mesh_compare_two_fields(p4est_t *p4est,double *field1,double *field2,const char *msg,d4est_mesh_boundary_option_t boundary_option,d4est_mesh_print_option_t print_option,double eps);
int d4est_mesh_get_local_nodes(p4est_t *p4est);
void d4est_mesh_get_local_nodes_callback(p4est_iter_volume_info_t *info,void *user_data);
void d4est_mesh_init_field_ext(p4est_t *p4est,double *node_vec,d4est_xyz_fcn_ext_t xyz_fcn,void *user,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom);
void d4est_mesh_compute_point_error(double *v1,double *v2,double *error,int local_nodes);
void d4est_mesh_init_field(p4est_t *p4est,double *node_vec,d4est_xyz_fcn_t init_fcn,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,void *user);
int d4est_mesh_update(p4est_t *p4est,p4est_ghost_t *ghost,void *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_geometry_storage_t *geometric_factors,d4est_mesh_quadrature_data_init_option_t quad_init_option,d4est_mesh_geometry_data_init_option_t geom_init_option,d4est_mesh_geometry_aliases_init_option_t alias_init_option,void(*user_fcn)(d4est_element_data_t *,void *),void *user_ctx);
void d4est_mesh_geometry_storage_initialize_aliases(p4est_t *p4est,d4est_mesh_geometry_storage_t *geometric_factors,d4est_local_sizes_t local_sizes);
d4est_local_sizes_t d4est_mesh_init_element_data(p4est_t *p4est,d4est_operators_t *d4est_ops,void(*user_fcn)(d4est_element_data_t *,void *),void *user_ctx);
double d4est_mesh_compute_l2_norm_sqr(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,double *nodal_vec,int local_nodes,norm_storage_option_t store_local);
void d4est_mesh_print_element_data_debug(p4est_t *p4est);
int d4est_mesh_debug_find_node(p4est_t *p4est,int node);
int d4est_mesh_global_node_to_local_node(p4est_t *p4est,int global_node);
d4est_element_data_t *d4est_mesh_get_element_data(p4est_t *p4est,int local_element_id);
void d4est_mesh_print_number_of_elements_per_tree(p4est_t *p4est);
int d4est_mesh_get_local_matrix_nodes(p4est_t *p4est);
void d4est_mesh_compute_jacobian_on_lgl_grid(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *jacobian_lgl);
void d4est_mesh_get_array_of_estimators(p4est_t *p4est,double *eta2_array);
void d4est_mesh_get_array_of_degrees(p4est_t *p4est,int *deg_array);
void d4est_mesh_geometry_storage_destroy(d4est_mesh_geometry_storage_t *geometric_factors);
d4est_mesh_geometry_storage_t *d4est_mesh_geometry_storage_init();

#endif
