#ifndef CURVED_ELEMENT_DATA_H
#define CURVED_ELEMENT_DATA_H 

#include <ip_flux_params.h>
#include "../dGMath/dgmath.h"
#include "../GridFunctions/grid_functions.h"
#include <d4est_geometry.h>



typedef enum { DO_NOT_STORE_LOCALLY, STORE_LOCALLY } norm_storage_option_t;
typedef enum { DIAM_APPROX, NO_DIAM_APPROX, DIAM_APPROX_CUBE} diam_compute_option_t;

typedef struct {

  /* double* J; */

  double* xyz;
  double* xyz_integ;
  double* xyz_rst_integ;
  double* J_integ;
  double* rst_xyz_integ;

  /* double* custom_abscissas; */
  /* double* custom_weights; */
  /* double* custom_interp; */
  
} geometric_factors_t;


typedef struct {

  /* identification */
  int id;
  int mpi_rank;
  
  int sqr_nodal_stride;
  int sqr_trace_stride;
  int nodal_stride;
  int integ_stride;
  
  int tree;
  int tree_quadid;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;
  
  /* auxiliary node vector aliases */
  double *Au_elem;

  /* Storage aliases for trace space calculations */
  double *M_ustar_min_u_n [(P4EST_DIM)];
  double *M_qstar_min_q_dot_n;

  /* geometric factors */
  /* double* J; /\* Jacobian *\/ */
  double* J_integ; /* Jacobian */
  /* double* invM; */
  /* double* invMface [P4EST_FACES]; */
  double* xyz [(P4EST_DIM)]; /* collocation points on physical grid */
  double* xyz_integ [(P4EST_DIM)]; /* collocation points on physical grid */
  /* double* xyz_rst[(P4EST_DIM)][(P4EST_DIM)]; /\* mapping derivatives *\/ */
  double* xyz_rst_integ[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  double* xyz_rst_Lobatto_integ[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  /* double* rst_xyz[(P4EST_DIM)][(P4EST_DIM)]; /\* inverse mapping derivatives *\/ */
  double* rst_xyz_integ[(P4EST_DIM)][(P4EST_DIM)]; /* inverse mapping derivatives */
  double diam; /* approximate value of element diameter*/

  double volume;
  double surface_area [(P4EST_FACES)];
  
  /* aposteriori/apriori error indicator for hp_amr or h_amr */
  double local_estimator;
  double local_predictor;
  
  /* storage for MPI transfers */
  double u_storage[MAX_NODES];
  /* double du_elem[(P4EST_DIM)][MAX_NODES]; */
  double dudr_elem[(P4EST_DIM)][MAX_NODES];
  /* double q_elem[(P4EST_DIM)][MAX_NODES]; */
  
  /* nodal degree */
  int deg;
  int deg_integ; /* for integration other than stiffness matrix */
  int deg_stiffness; /* for integration of the stiffness matrix */
  
#ifndef NDEBUG
  /* useful flag for debugging */
  int debug_flag;
  int on_bdry;
#endif
  
} curved_element_data_t;

typedef void
(*curved_element_data_user_fcn_t)
(
 void*,
 void*
);

typedef struct {
  
  int local_nodes;
  int local_sqr_nodes;
  int local_sqr_trace_nodes;
  int local_nodes_integ;
  /* int local_sqr_nodes_invM; */
  
} curved_element_data_local_sizes_t;

/* This file was automatically generated.  Do not edit! */
void curved_element_data_get_array_of_degrees(p4est_t *p4est,int *deg_array);
void curved_element_data_compute_jacobian_on_lgl_grid(p4est_t *p4est,d4est_geometry_t *d4est_geometry,dgmath_jit_dbase_t *dgmath_jit_dbase,double *jacobian_lgl);
void curved_element_data_apply_curved_stiffness_matrix(dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,curved_element_data_t *elem_data,quadrature_type_t quad_type,double *vec,double *stiff_vec);
int curved_element_data_get_local_matrix_nodes(p4est_t *p4est);
void curved_element_data_form_fofufofvlilj_matrix(dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,double *u,double *v,curved_element_data_t *elem_data,int deg_quad,quadrature_type_t quad_type,int dim,double *mat,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void curved_element_data_apply_fofufofvlj(dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,double *u,double *v,curved_element_data_t *elem_data,int deg_quad,quadrature_type_t quad_type,int dim,double *out,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void curved_element_data_apply_fofufofvlilj(dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,double *vec,double *u,double *v,curved_element_data_t *elem_data,int deg_quad,quadrature_type_t quad_type,int dim,double *Mvec,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void curved_element_data_print_number_of_elements_per_tree(p4est_t *p4est);
void curved_element_data_print_local_estimator(p4est_t *p4est);
void curved_element_data_print_element_data_debug(p4est_t *p4est);
int curved_element_data_count_boundary_elements(p4est_t *p4est);
curved_element_data_t *curved_element_data_get_element_data(p4est_t *p4est,int local_element_id);
int curved_element_data_global_node_to_local_node(p4est_t *p4est,int global_node);
int curved_element_data_debug_find_node(p4est_t *p4est,int node);
void curved_element_data_reorient_f_p_elements_to_f_m_order(curved_element_data_t **e_p,int face_dim,int f_m,int f_p,int o,int faces_p,curved_element_data_t *e_p_oriented[(P4EST_HALF)]);
void curved_element_data_store_nodal_vec_in_vertex_array(p4est_t *p4est,double *nodal_vec,double *corner_vec);
void curved_element_data_debug_spheresym(p4est_t *p4est,dgmath_jit_dbase_t *dgbase,double *vec);
void curved_element_data_store_element_scalar_in_vertex_array(p4est_t *p4est,double *vertex_array,double(*get_local_scalar_fcn)(curved_element_data_t *));
void curved_element_data_copy_from_storage_to_vec(p4est_t *p4est,double *vec);
int curved_element_data_get_local_nodes(p4est_t *p4est);
void curved_element_data_get_local_nodes_callback(p4est_iter_volume_info_t *info,void *user_data);
void curved_element_data_copy_from_vec_to_storage(p4est_t *p4est,double *vec);
double curved_element_data_compute_dg_norm_sqr(p4est_t *p4est,double *nodal_vec,int local_nodes,ip_flux_params_t *ip_flux_params,d4est_geometry_t *d4est_geom,p4est_ghost_t *ghost,void *ghost_data,dgmath_jit_dbase_t *dgmath_jit_dbase);
double curved_element_data_compute_l2_norm_sqr(p4est_t *p4est,double *nodal_vec,int local_nodes,dgmath_jit_dbase_t *dgmath_jit_dbase,quadrature_type_t quad_type,norm_storage_option_t store_local);
void curved_element_data_init_new(p4est_t *p4est,geometric_factors_t *geometric_factors,dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,curved_element_data_user_fcn_t user_fcn,void *user_ctx,int compute_geometric_data,int set_geometric_aliases);
curved_element_data_local_sizes_t curved_element_data_compute_strides_and_sizes(p4est_t *p4est,dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geometry,curved_element_data_user_fcn_t user_fcn,void *user_ctx);
double curved_element_data_compute_element_face_area(curved_element_data_t *elem_data,dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *geom,int face,int deg);
double curved_element_data_compute_element_volume(dgmath_jit_dbase_t *dgmath_jit_dbase,int deg_GL,double *jac_GL);
void curved_element_data_compute_grid_volume_and_surface_area(p4est_t *p4est,double *volume,double *surface_area);
double curved_element_data_compute_diam(double *xyz[(P4EST_DIM)],int deg,diam_compute_option_t option);
void curved_element_data_debug_print_node_vecs(p4est_t *p4est,double **vecs,int num_vecs,int *elems,int num_elems);
void curved_element_data_print_node_vec(p4est_t *p4est,double *vec);
void geometric_factors_destroy(geometric_factors_t *geometric_factors);
void geometric_factors_reinit(p4est_t *p4est,geometric_factors_t *geometric_factors,curved_element_data_local_sizes_t local_sizes);
geometric_factors_t *geometric_factors_init(p4est_t *p4est);
void curved_element_data_init_node_vec_ext(p4est_t *p4est,double *node_vec,grid_fcn_ext_t fofxyzv,double *v,double *fofxyzv_user,dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geom);
void curved_element_data_init_node_vec(p4est_t *p4est,double *node_vec,grid_fcn_t init_fcn,dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geom);
void curved_element_data_set_degrees(p4est_t *p4est,int(*set_deg_fcn)(int,int,int));

#endif
