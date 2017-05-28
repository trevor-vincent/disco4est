#ifndef CURVED_ELEMENT_DATA_H
#define CURVED_ELEMENT_DATA_H 

#include <ip_flux_params.h>
#include <grid_functions.h>
#include <d4est_operators.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>



typedef enum { DO_NOT_STORE_LOCALLY, STORE_LOCALLY } norm_storage_option_t;
typedef enum { DIAM_APPROX, NO_DIAM_APPROX, DIAM_APPROX_CUBE} diam_compute_option_t;



typedef struct {

  /* identification */
  int id;
  int mpi_rank;
  
  int sqr_nodal_stride;
  int sqr_trace_stride;
  int nodal_stride;
  int quad_stride;
  
  int tree;
  int tree_quadid;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;
  
  /* geometric factors */
  double* J_quad; /* Jacobian */
  double* xyz [(P4EST_DIM)]; /* points on lobatto grid */
  double* xyz_quad [(P4EST_DIM)]; /* points on quadrature grid */
  double* xyz_rst_quad[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  double* rst_xyz_quad[(P4EST_DIM)][(P4EST_DIM)]; /* inverse mapping derivatives */

  double diam; /* approximate value of element diameter*/
  double volume;
  double surface_area [(P4EST_FACES)];
  
  /* aposteriori/apriori error indicator for hp_amr or h_amr */
  double local_estimator;
  double local_predictor;

  double* Au_elem;  /* alias for Au */
  double u_elem[MAX_NODES];   /* storage for MPI transfers */
  double dudr_elem[(P4EST_DIM)][MAX_NODES];   /* storage for MPI transfers */
  
  int deg; /* nodal degree */
  int deg_quad; /* deg for quadrature */
  
#ifndef NDEBUG
  /* useful flag for debugging */
  int debug_flag;
  int on_bdry;
#endif
  
} d4est_element_data_t;

typedef void
(*d4est_element_data_user_fcn_t)
(
 void*,
 void*
);

/* This file was automatically generated.  Do not edit! */
void d4est_element_data_get_array_of_degrees(p4est_t *p4est,int *deg_array);
void d4est_element_data_compute_jacobian_on_lgl_grid(p4est_t *p4est,d4est_geometry_t *d4est_geometry,d4est_operators_t *d4est_ops,double *jacobian_lgl);
int d4est_element_data_get_local_matrix_nodes(p4est_t *p4est);
void d4est_element_data_print_number_of_elements_per_tree(p4est_t *p4est);
void d4est_element_data_print_local_estimator(p4est_t *p4est);
void d4est_element_data_print_element_data_debug(p4est_t *p4est);
int d4est_element_data_count_boundary_elements(p4est_t *p4est);
d4est_element_data_t *d4est_element_data_get_element_data(p4est_t *p4est,int local_element_id);
int d4est_element_data_global_node_to_local_node(p4est_t *p4est,int global_node);
int d4est_element_data_debug_find_node(p4est_t *p4est,int node);
void d4est_element_data_reorient_f_p_elements_to_f_m_order(d4est_element_data_t **e_p,int face_dim,int f_m,int f_p,int o,int faces_p,d4est_element_data_t *e_p_oriented[(P4EST_HALF)]);
void d4est_element_data_store_nodal_vec_in_vertex_array(p4est_t *p4est,double *nodal_vec,double *corner_vec);
void d4est_element_data_debug_spheresym(p4est_t *p4est,d4est_operators_t *d4est_ops,double *vec);
void d4est_element_data_store_element_scalar_in_vertex_array(p4est_t *p4est,double *vertex_array,double(*get_local_scalar_fcn)(d4est_element_data_t *));
void d4est_element_data_copy_from_storage_to_vec(p4est_t *p4est,double *vec);
void d4est_element_data_copy_from_vec_to_storage(p4est_t *p4est,double *vec);
int d4est_element_data_get_local_nodes(p4est_t *p4est);
void d4est_element_data_get_local_nodes_callback(p4est_iter_volume_info_t *info,void *user_data);
double d4est_element_data_compute_l2_norm_sqr(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,double *nodal_vec,int local_nodes,norm_storage_option_t store_local);
void d4est_element_data_init_new(p4est_t *p4est,d4est_geometry_storage_t *geometric_factors,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_element_data_user_fcn_t user_fcn,void *user_ctx,int compute_geometric_data,int set_geometric_aliases);
d4est_local_sizes_t d4est_element_data_compute_strides_and_sizes(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geometry,d4est_element_data_user_fcn_t user_fcn,void *user_ctx);
double d4est_element_data_compute_element_face_area(d4est_element_data_t *elem_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,int face,int deg);
double d4est_element_data_compute_element_volume(d4est_operators_t *d4est_ops,int deg_GL,double *jac_GL);
void d4est_element_data_compute_grid_volume_and_surface_area(p4est_t *p4est,double *volume,double *surface_area);
double d4est_element_data_compute_diam(double *xyz[(P4EST_DIM)],int deg,diam_compute_option_t option);
void d4est_element_data_debug_print_node_vecs(p4est_t *p4est,double **vecs,int num_vecs,int *elems,int num_elems);
void d4est_element_data_print_node_vec(p4est_t *p4est,double *vec);
void d4est_element_data_init_node_vec_ext(p4est_t *p4est,double *node_vec,grid_fcn_ext_t fofxyzv,double *v,double *fofxyzv_user,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom);
void d4est_element_data_init_node_vec(p4est_t *p4est,double *node_vec,grid_fcn_t init_fcn,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom);
void d4est_element_data_set_degrees(p4est_t *p4est,int(*set_deg_fcn)(int,int,int));


#endif
