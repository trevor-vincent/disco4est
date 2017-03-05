#ifndef DGMATH_H
#define DGMATH_H

#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"

#define MAX_DEGREE 8
#if (P4EST_DIM) == 3
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#else
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#endif

typedef struct {
  int dgmath_max_degree_used;
  int dgmath_max_storage;

  double** dgmath_ref_Mij_1d_table;
  double** dgmath_GLL_nodes_1d_table;
  double** dgmath_GLL_weights_1d_table;
  double** dgmath_GL_nodes_1d_table;
  double** dgmath_GL_weights_1d_table;
  double** dgmath_ref_invMij_1d_table;
  double** dgmath_ref_invVij_1d_table;
  double** dgmath_ref_GaussVij_1d_table;
  double** dgmath_ref_invGaussVij_1d_table;
  double** dgmath_vtk_interp_1d_table;
  double** dgmath_ref_GLL_to_GL_interp_1d_table;
  double** dgmath_ref_GLL_to_GL_interp_1d_inverse_table;
  double** dgmath_ref_GLL_to_GL_interp_trans_1d_table;
  double** dgmath_ref_GLL_to_GL_interp_trans_1d_inverse_table;
  
  double** dgmath_ref_Dij_1d_table;
  double** dgmath_ref_GaussDij_1d_table;
  double** dgmath_LIFT_1d_table;
  double** dgmath_ref_xyz_nd_table;
  double** dgmath_vtk_rst_2d_table;
  double** dgmath_vtk_rst_3d_table;
  double** dgmath_ref_Gauss_xyz_nd_table;
  
  double** dgmath_FLIP_1d_table;
  double** dgmath_hp_prolong_1d_table;
  double** dgmath_p_prolong_1d_table;
  double** dgmath_p_prolong_1d_inverse_table;
  double** dgmath_hp_restrict_1d_table;
  double** dgmath_p_restrict_1d_table;
  double** dgmath_hp_prolong_transpose_1d_table;
  double** dgmath_p_prolong_transpose_1d_table;  
  double** dgmath_p_prolong_transpose_1d_inverse_table;  
  double** dgmath_hp_restrict_interp_1d_table;
  double** dgmath_p_restrict_interp_1d_table;
  
} dgmath_jit_dbase_t;

typedef void (*dgmath_hp_restrict_fcn_t)(dgmath_jit_dbase_t *, double *, int *,
                                         int, int, double *);

typedef void (*dgmath_p_restrict_fcn_t)(dgmath_jit_dbase_t *, double *, int,
                                        int, int, double *);

typedef void (*dgmath_hp_prolong_fcn_t)(dgmath_jit_dbase_t *, double *, int,
                                        int, int *, double *);

typedef void (*dgmath_p_prolong_fcn_t)(dgmath_jit_dbase_t *, double *, int, int,
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
} dgmath_face_info_t;

typedef struct {

  double* r;
  double* s;
  double* t;
  
} dgmath_rst_t;

typedef enum {GAUSS, LOBATTO} quadrature_type_t; 

/* This file was automatically generated.  Do not edit! */
double *dgmath_fetch_vtk_rst(dgmath_jit_dbase_t *dgbase,int deg,int dim);
double dgmath_Gauss_quadrature(dgmath_jit_dbase_t *dgmath_jit_dbase,double *u,double *v,double *jac_GL,int deg_GL,int dim);
void dgmath_convert_nodal_to_vtk(dgmath_jit_dbase_t *dgbase,double *vec,int dim,int deg,double *vtk_vec);
void dgmath_apply_vtk_interp(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int dim,int deg,int c,double *out);
int dgmath_corner_to_node(int dim,int deg,int corner);
dgmath_rst_t dgmath_get_rst_points(dgmath_jit_dbase_t *dgmath_jit_dbase,int deg,int dim,quadrature_type_t type);
void dgmath_reorient_face_data(dgmath_jit_dbase_t *dgbase,double *in,int face_dim,int deg,int o,int f_m,int f_p,double *out);
int dgmath_reorient_face_order(int face_dim,int f_m,int f_p,int o,int i);
void dgmath_apply_FLIP(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,int dir,double *out);
void dgmath_compute_invM_vec_wo_aliasing(dgmath_jit_dbase_t *dgmath_jit_dbase,double *xyz_rst[(P4EST_DIM)][P4EST_DIM],double *vec,int degH,double *invMvec,int degh);
void dgmath_compute_M_vec_wo_aliasing(dgmath_jit_dbase_t *dgmath_jit_dbase,double *xyz_rst[(P4EST_DIM)][P4EST_DIM],double *vec,int degH,double *MJvec,int degh);
void dgmath_get_face_info(int f,dgmath_face_info_t *face_info);
void dgmath_apply_p_restrict_interp(dgmath_jit_dbase_t *dgbase,double *in,int degh,int dim,int degH,double *out);
void dgmath_apply_hp_restrict_interp(dgmath_jit_dbase_t *dgbase,double *in,int *degh,int dim,int degH,double *out);
void dgmath_project_mortar_onto_side(dgmath_jit_dbase_t *dgbase,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void dgmath_project_mass_mortar_onto_side(dgmath_jit_dbase_t *dgmath,double *in_mortar,int faces_mortar,int *deg_mortar,double *out_side,int faces_side,int *deg_side);
void dgmath_project_side_onto_mortar_space(dgmath_jit_dbase_t *dgbase,double *in_side,int faces_side,int *deg_side,double *out_mortar,int faces_mortar,int *deg_mortar);
void dgmath_apply_p_prolong_transpose(dgmath_jit_dbase_t *dgbase,double *in,int degh,int dim,int degH,double *out);
void dgmath_apply_hp_prolong_transpose(dgmath_jit_dbase_t *dgbase,double *in,int *degh,int dim,int degH,double *out);
void dgmath_build_p_prolong_transpose_1d_inverse(dgmath_jit_dbase_t *dgbase,double *p_prolong_transpose_1d_inverse,int degH,int degh);
void dgmath_rtox_array(double *r,double xl,double h,double *x,int nodes);
double dgmath_rtox(double r,double xl,double h);
void dgmath_convert_nodal_to_modal(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,double *out);
void dgmath_grad_vandermonde_1d(double *dr_v1d,double *lobatto_nodes,int degree);
void dgmath_vandermonde_1d(double *v1d,double *lobatto_nodes,int degree);
void dgmath_get_normal(int face,int dim,double *n);
void dgmath_apply_slicer(dgmath_jit_dbase_t *dgbase,double *in,int dim,int face,int deg,double *out);
void dgmath_apply_LIFT(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,int face,double *out);
double *dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase_t *dgbase,int dim,int deg,int dir);
double *dgmath_fetch_xyz_nd(dgmath_jit_dbase_t *dgbase,int dim,int deg,int dir);
int dgmath_fetch_max_degree_used(dgmath_jit_dbase_t *dgbase);
void dgmath_apply_hp_restrict(dgmath_jit_dbase_t *dgbase,double *in,int *degh,int dim,int degH,double *out);
void dgmath_apply_p_restrict(dgmath_jit_dbase_t *dgbase,double *in,int degh,int dim,int degH,double *out);
double *dgmath_fetch_p_restrict_1d(dgmath_jit_dbase_t *dgbase,int degH,int degh);
void dgmath_build_hp_restrict_1d_aux(int degh,int degH,double *hp_prolong_matrix_1d,double *mass_matrix_rs_degh,double *inv_mass_matrix_rs_degH,double *hp_restrict_matrix_1d);
void dgmath_build_p_prolong_1d_inverse(dgmath_jit_dbase_t *dgbase,double *p_prolong_1d_inverse,int degH,int degh);
void dgmath_hp_apply_nd_prolong_transpose_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_prolong_transpose_matrix_1d);
void dgmath_build_p_prolong_1d(dgmath_jit_dbase_t *dgbase,double *p_prolong_1d,int degH,int degh);
void dgmath_apply_GaussDij(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg_Lobatto,int deg_Gauss,int dir,double *out);
double *dgmath_fetch_GaussDij_1d(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
void dgmath_apply_invMij(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,double *out);
void dgmath_apply_Mij(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,double *out);
double *dgmath_fetch_Mij_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_build_GL_nodes_and_weights_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *GL_nodes,double *GL_weights,int deg);
void dgmath_build_GLL_nodes_and_weights_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *GLL_nodes,double *GLL_weights,int deg);
void dgmath_build_Mij_1d(dgmath_jit_dbase_t *dgbase,double *Mij_1d,int deg);
double dgmath_gradjacobi(double r,double alpha,double beta,int N);
void dgmath_hp_apply_nd_restrict_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_restrict_matrix_1d);
double *dgmath_fetch_p_prolong_1d_inverse(dgmath_jit_dbase_t *dgbase,int degH,int degh);
double *dgmath_fetch_p_prolong_transpose_1d_inverse(dgmath_jit_dbase_t *dgbase,int degH,int degh);
void dgmath_apply_curvedInverseLobattoMass(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int deg_Lobatto,double *jac_Lobatto_integ,int deg_Lobatto_integ,int dim,double *out);
void dgmath_apply_curvedInverseGaussMass(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int deg_Lobatto,double *jac_Gauss,int deg_Gauss,int dim,double *out);
void dgmath_apply_fofufofvlj_Gaussnodes(dgmath_jit_dbase_t *dgmath_jit_dbase,double *u,double *v,int deg_Lobatto,double *jac_Gauss,double *xyz_Gauss[(P4EST_DIM)],int deg_Gauss,int dim,double *out,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void dgmath_apply_fofufofvlilj_Gaussnodes(dgmath_jit_dbase_t *dgmath_jit_dbase,double *vec,double *u,double *v,int deg_Lobatto,double *jac_Gauss,double *xyz_Gauss[(P4EST_DIM)],int deg_Gauss,int dim,double *Mvec,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void dgmath_form_fofufofvlilj_matrix_Gaussnodes(dgmath_jit_dbase_t *dgmath_jit_dbase,double *u,double *v,int deg_Lobatto,double *xyz_Gauss[(P4EST_DIM)],double *jac_Gauss,int deg_Gauss,int dim,double *out,grid_fcn_ext_t fofu_fcn,void *fofu_ctx,grid_fcn_ext_t fofv_fcn,void *fofv_ctx);
void dgmath_compute_PT_mat_P(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in_mat,int degH,int dim,int *degh,int children,double *out_mat);
void dgmath_apply_hp_prolong(dgmath_jit_dbase_t *dgbase,double *in,int degH,int dim,int *degh,double *out);
void dgmath_compute_prolong_matrix(dgmath_jit_dbase_t *dgmath_jit_dbase,int degH,int dim,int *degh,int children,double *prolong_mat);
void dgmath_compute_curvedGaussMass(dgmath_jit_dbase_t *dgmath_jit_dbase,int deg_Lobatto,int deg_Gauss,int dim,double *jac_Gauss,double *M);
void dgmath_compute_curvedInverseGaussMass(dgmath_jit_dbase_t *dgmath_jit_dbase,int deg_Lobatto,int deg_Gauss,int dim,double *jac_Gauss,double *invM);
void dgmath_apply_curvedLobattoMass_onLobattoIntegNodeVec(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in_Lobatto_integ,int deg_Lobatto,double *jac_Lobatto_integ,int deg_Lobatto_integ,int dim,double *out);
void dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in_Gauss,int deg_Lobatto,double *jac_Gauss,int deg_Gauss,int dim,double *out);
void dgmath_apply_p_prolong(dgmath_jit_dbase_t *dgbase,double *in,int degH,int dim,int degh,double *out);
double *dgmath_fetch_p_prolong_transpose_1d(dgmath_jit_dbase_t *dgbase,int degH,int degh);
double *dgmath_fetch_GLL_weights_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_apply_curvedLobattoMass(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int deg_Lobatto,double *jac_Lobatto_integ,int deg_Lobatto_integ,int dim,double *out);
void dgmath_apply_Dij_transpose(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,int dir,double *out);
void dgmath_apply_Dij(dgmath_jit_dbase_t *dgbase,double *in,int dim,int deg,int dir,double *out);
void dgmath_apply_curvedGaussStiff(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int deg_Lobatto,double *jac_Gauss,double *rst_xyz[(P4EST_DIM)][(P4EST_DIM)],int deg_Gauss,int dim,double *out);
double *dgmath_fetch_GL_weights_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_apply_curvedGaussMass(dgmath_jit_dbase_t *dgmath_jit_dbase,double *in,int deg_Lobatto,double *jac_Gauss,int deg_Gauss,int dim,double *out);
void dgmath_interp_GLL_to_GL(dgmath_jit_dbase_t *dgmath_jit_dbase,double *u_Lobatto_in,int deg_Lobatto,int deg_Gauss,double *u_Gauss_out,int dim);
void dgmath_build_ref_GLL_to_GL_interp_trans_1d_inverse(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_GLL_to_GL_interp_trans_1d_inverse,int Lobatto_degree,int Gauss_degree);
double *dgmath_fetch_ref_GLL_to_GL_interp_trans_1d_inverse(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
void dgmath_build_ref_GLL_to_GL_interp_trans_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_GLL_to_GL_interp_trans_1d,int Lobatto_degree,int Gauss_degree);
double *dgmath_fetch_ref_GLL_to_GL_interp_trans_1d(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
void dgmath_build_ref_GLL_to_GL_interp_1d_inverse(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_GLL_to_GL_interp_1d_inverse,int Lobatto_degree,int Gauss_degree);
double *dgmath_fetch_ref_GLL_to_GL_interp_1d_inverse(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
double *dgmath_fetch_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
void dgmath_build_ref_GLL_to_GL_interp_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_GLL_to_GL_interp_1d,int Lobatto_degree,int Gauss_degree);
double *dgmath_fetch_invGaussVij_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_build_ref_invGaussVij_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_invGaussVij_1d,int Gauss_degree);
double *dgmath_fetch_ref_GaussVij_1d(dgmath_jit_dbase_t *dgbase,int deg_Lobatto,int deg_Gauss);
double *dgmath_fetch_GL_nodes_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_build_ref_GaussVij_1d(dgmath_jit_dbase_t *dgmath_jit_dbase,double *ref_Vij_Gauss_1d,int Lobatto_degree,int Gauss_degree);
void dgmath_hp_apply_nd_prolong_with_ptr(double *Uh,int degh,double *UH,int degH,int dim,int c,double *hp_prolong_matrix_1d);
double *dgmath_fetch_invVij_1d(dgmath_jit_dbase_t *dgbase,int deg);
double dgmath_jacobi(double r,double alpha,double beta,int N);
double *dgmath_fetch_GLL_nodes_1d(dgmath_jit_dbase_t *dgbase,int deg);
void dgmath_build_Vij_1d(dgmath_jit_dbase_t *dgbase,double *Vij_1d,int deg);
void dgmath_child_to_parent_ref_coords(int c,double *r);
int dgmath_is_child_left_or_right(int c,int dir);
int dgmath_get_nodes(int dim,int deg);
void dgmath_jit_dbase_destroy(dgmath_jit_dbase_t *dgbase);
dgmath_jit_dbase_t *dgmath_jit_dbase_init();

#endif
