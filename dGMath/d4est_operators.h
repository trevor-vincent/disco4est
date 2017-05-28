#ifndef D4EST_OPERATORS_H
#define D4EST_OPERATORS_H

#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"

#define MAX_DEGREE 20
#if (P4EST_DIM) == 3
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#else
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#endif

typedef struct {
  int d4est_operators_max_degree_used;
  int d4est_operators_max_storage;

  double** d4est_operators_ref_Mij_1d_table;
  double** d4est_operators_GLL_nodes_1d_table;
  double** d4est_operators_GLL_weights_1d_table;
  double** d4est_operators_GL_nodes_1d_table;
  double** d4est_operators_GL_weights_1d_table;
  double** d4est_operators_ref_invMij_1d_table;
  double** d4est_operators_ref_invVij_1d_table;
  double** d4est_operators_ref_GaussVij_1d_table;
  double** d4est_operators_ref_invGaussVij_1d_table;
  double** d4est_operators_vtk_interp_1d_table;
  
  double** d4est_operators_ref_GLL_to_GL_interp_1d_table;
  double** d4est_operators_ref_GLL_to_GL_interp_1d_inverse_table;
  double** d4est_operators_ref_GLL_to_GL_interp_trans_1d_table;
  double** d4est_operators_ref_GLL_to_GL_interp_trans_1d_inverse_table;
  
  double** d4est_operators_ref_Dij_1d_table;
  double** d4est_operators_ref_GaussDij_1d_table;
  double** d4est_operators_LIFT_1d_table;
  double** d4est_operators_ref_xyz_nd_table;
  double** d4est_operators_vtk_rst_2d_table;
  double** d4est_operators_vtk_rst_3d_table;
  double** d4est_operators_ref_Gauss_xyz_nd_table;
  
  double** d4est_operators_FLIP_1d_table;
  double** d4est_operators_hp_prolong_1d_table;
  double** d4est_operators_p_prolong_1d_table;
  double** d4est_operators_p_prolong_1d_inverse_table;
  double** d4est_operators_hp_restrict_1d_table;
  double** d4est_operators_p_restrict_1d_table;
  double** d4est_operators_hp_prolong_transpose_1d_table;
  double** d4est_operators_p_prolong_transpose_1d_table;  
  double** d4est_operators_p_prolong_transpose_1d_inverse_table;  
  double** d4est_operators_hp_restrict_interp_1d_table;
  double** d4est_operators_p_restrict_interp_1d_table;
  
} d4est_operators_t;

typedef void (*d4est_operators_hp_restrict_fcn_t)(d4est_operators_t *, double *, int *,
                                         int, int, double *);

typedef void (*d4est_operators_p_restrict_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int, double *);

typedef void (*d4est_operators_hp_prolong_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int *, double *);

typedef void (*d4est_operators_p_prolong_fcn_t)(d4est_operators_t *, double *, int, int,
                                       int, double *);

//         +----------------------+
//         |			  |
//         |			  |	     b
//         |			  |	     |
//         |			  |	     |
//         |	   face f	  |	     |
//         |	       	       	  |  	     |
//         |			  |	     /---------	a
//         |			  |	    /
//         |			  |	   /
//         |			  |	  c
//         |----------------------+

typedef struct {
  int a;      /* "x" coord on this face (z-ordering) */
  int b;      /* "y" coord on this face (z-ordering) */
  int c;      /* "normal" coord on this face (z-ordering) */
  double sgn; /* is the normal in the - or +  c-direction */
} d4est_operators_face_info_t;

typedef struct {

  double* r;
  double* s;
  double* t;
  
} d4est_rst_t;


typedef struct {
  
  int local_nodes;
  int local_sqr_nodes;
  int local_sqr_trace_nodes;
  int local_nodes_quad;
  
} d4est_local_sizes_t;

/* This file was automatically generated.  Do not edit! */
double *d4est_operators_fetch_vtk_rst(d4est_operators_t *d4est_ops,int deg,int dim);
void d4est_operators_convert_nodal_to_vtk(d4est_operators_t *d4est_ops,double *vec,int dim,int deg,double *vtk_vec);
void d4est_operators_apply_vtk_interp(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int c,double *out);
int d4est_operators_corner_to_node(int dim,int deg,int corner);
void d4est_operators_reorient_face_data(d4est_operators_t *d4est_ops,double *in,int face_dim,int deg,int o,int f_m,int f_p,double *out);
int d4est_operators_reorient_face_order(int face_dim,int f_m,int f_p,int o,int i);
void d4est_operators_apply_FLIP(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
void d4est_operators_compute_invM_vec_wo_aliasing(d4est_operators_t *d4est_ops,double *xyz_rst[(P4EST_DIM)][P4EST_DIM],double *vec,int degH,double *invMvec,int degh);
void d4est_operators_compute_M_vec_wo_aliasing(d4est_operators_t *d4est_ops,double *xyz_rst[(P4EST_DIM)][P4EST_DIM],double *vec,int degH,double *MJvec,int degh);
void d4est_operators_get_face_info(int f,d4est_operators_face_info_t *face_info);
void d4est_operators_apply_p_restrict_interp(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
void d4est_operators_apply_hp_restrict_interp(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
void d4est_operators_project_mortar_onto_side(d4est_operators_t *d4est_ops,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_operators_project_mass_mortar_onto_side(d4est_operators_t *dgmath,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void d4est_operators_project_side_onto_mortar_space(d4est_operators_t *d4est_ops,double *in_side,int faces_side,int *deg_side,double *out_mortar,int faces_mortar,int *deg_mortar);
void d4est_operators_apply_p_prolong_transpose(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
void d4est_operators_apply_hp_prolong_transpose(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
double *d4est_operators_fetch_p_prolong_transpose_1d_inverse(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_build_p_prolong_transpose_1d_inverse(d4est_operators_t *d4est_ops,double *p_prolong_transpose_1d_inverse,int degH,int degh);
double *d4est_operators_fetch_p_prolong_transpose_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_rtox_array(double *r,double xl,double h,double *x,int nodes);
double d4est_operators_rtox(double r,double xl,double h);
void d4est_operators_convert_nodal_to_modal(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
void d4est_operators_grad_vandermonde_1d(double *dr_v1d,double *lobatto_nodes,int degree);
void d4est_operators_vandermonde_1d(double *v1d,double *lobatto_nodes,int degree);
void d4est_operators_get_normal(int face,int dim,double *n);
void d4est_operators_apply_slicer(d4est_operators_t *d4est_ops,double *in,int dim,int face,int deg,double *out);
void d4est_operators_apply_LIFT(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int face,double *out);
double *d4est_operators_fetch_Gauss_xyz_nd(d4est_operators_t *d4est_ops,int dim,int deg,int dir);
double *d4est_operators_fetch_xyz_nd(d4est_operators_t *d4est_ops,int dim,int deg,int dir);
int d4est_operators_fetch_max_degree_used(d4est_operators_t *d4est_ops);
void d4est_operators_apply_hp_restrict(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
void d4est_operators_apply_p_restrict(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
double *d4est_operators_fetch_p_restrict_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_build_hp_restrict_1d_aux(int degh,int degH,double *hp_prolong_matrix_1d,double *mass_matrix_rs_degh,double *inv_mass_matrix_rs_degH,double *hp_restrict_matrix_1d);
void d4est_operators_build_p_prolong_1d_inverse(d4est_operators_t *d4est_ops,double *p_prolong_1d_inverse,int degH,int degh);
double *d4est_operators_fetch_p_prolong_1d_inverse(d4est_operators_t *d4est_ops,int degH,int degh);
double *d4est_operators_fetch_p_prolong_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_prolong_transpose_matrix_1d);
void d4est_operators_build_p_prolong_1d(d4est_operators_t *d4est_ops,double *p_prolong_1d,int degH,int degh);
void d4est_operators_apply_Dij_transpose(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
void d4est_operators_apply_Dij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
void d4est_operators_apply_GaussDij(d4est_operators_t *d4est_ops,double *in,int dim,int deg_Lobatto,int deg_Gauss,int dir,double *out);
double *d4est_operators_fetch_GaussDij_1d(d4est_operators_t *d4est_ops,int deg_Lobatto,int deg_Gauss);
void d4est_operators_apply_invMij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
void d4est_operators_apply_Mij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
double *d4est_operators_fetch_Mij_1d(d4est_operators_t *d4est_ops,int deg);
double *d4est_operators_fetch_GL_weights_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_GL_nodes_and_weights_1d(d4est_operators_t *d4est_ops,double *GL_nodes,double *GL_weights,int deg);
double *d4est_operators_fetch_GLL_weights_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_GLL_nodes_and_weights_1d(d4est_operators_t *d4est_ops,double *GLL_nodes,double *GLL_weights,int deg);
void d4est_operators_build_Mij_1d(d4est_operators_t *d4est_ops,double *Mij_1d,int deg);
double d4est_operators_gradjacobi(double r,double alpha,double beta,int N);
void d4est_operators_hp_apply_nd_restrict_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_restrict_matrix_1d);
void d4est_operators_compute_PT_mat_P(d4est_operators_t *d4est_ops,double *mat,int degH,int dim,int *degh,int children,double *PT_mat_P);
void d4est_operators_apply_p_prolong(d4est_operators_t *d4est_ops,double *in,int degH,int dim,int degh,double *out);
void d4est_operators_apply_hp_prolong(d4est_operators_t *d4est_ops,double *in,int degH,int dim,int *degh,double *out);
void d4est_operators_compute_prolong_matrix(d4est_operators_t *d4est_ops,int degH,int dim,int *degh,int children,double *prolong_mat);
double *d4est_operators_fetch_ref_GLL_to_GL_interp_trans_1d(d4est_operators_t *d4est_ops,int deg_Lobatto,int deg_Gauss);
void d4est_operators_build_ref_GLL_to_GL_interp_trans_1d(d4est_operators_t *d4est_ops,double *ref_GLL_to_GL_interp_trans_1d,int Lobatto_degree,int Gauss_degree);
void d4est_operators_build_ref_GLL_to_GL_interp_1d_inverse(d4est_operators_t *d4est_ops,double *ref_GLL_to_GL_interp_1d_inverse,int Lobatto_degree,int Gauss_degree);
double *d4est_operators_fetch_ref_GLL_to_GL_interp_1d_inverse(d4est_operators_t *d4est_ops,int deg_Lobatto,int deg_Gauss);
double *d4est_operators_fetch_ref_GLL_to_GL_interp_1d(d4est_operators_t *d4est_ops,int deg_Lobatto,int deg_Gauss);
void d4est_operators_build_custom_GL_interp_1d(d4est_operators_t *d4est_ops,double *custom_GL_interp_1d,int Lobatto_degree,int custom_degree,double *custom_points);
void d4est_operators_build_ref_GLL_to_GL_interp_1d(d4est_operators_t *d4est_ops,double *ref_GLL_to_GL_interp_1d,int Lobatto_degree,int Gauss_degree);
double *d4est_operators_fetch_invGaussVij_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_ref_invGaussVij_1d(d4est_operators_t *d4est_ops,double *ref_invGaussVij_1d,int Gauss_degree);
double *d4est_operators_fetch_ref_GaussVij_1d(d4est_operators_t *d4est_ops,int deg_Lobatto,int deg_Gauss);
double *d4est_operators_fetch_GL_nodes_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_ref_GaussVij_1d(d4est_operators_t *d4est_ops,double *ref_Vij_Gauss_1d,int Lobatto_degree,int Gauss_degree);
void d4est_operators_hp_apply_nd_prolong_with_ptr(double *Uh,int degh,double *UH,int degH,int dim,int c,double *hp_prolong_matrix_1d);
double *d4est_operators_fetch_invVij_1d(d4est_operators_t *d4est_ops,int deg);
double d4est_operators_jacobi(double r,double alpha,double beta,int N);
double *d4est_operators_fetch_GLL_nodes_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_Vij_1d(d4est_operators_t *d4est_ops,double *Vij_1d,int deg);
void d4est_operators_child_to_parent_ref_coords(int c,double *r);
int d4est_operators_is_child_left_or_right(int c,int dir);
int d4est_operators_get_nodes(int dim,int deg);
void d4est_ops_destroy(d4est_operators_t *d4est_ops);
d4est_operators_t *d4est_ops_init();

#endif
