#ifndef D4EST_OPERATORS_H
#define D4EST_OPERATORS_H

#include <pXest.h>
#include <grid_functions.h>
#include <d4est_lgl.h>
#include <d4est_reference.h>

typedef struct {

  int max_degree;
  
  double** lobatto_nodes_1d_table;
  double** lobatto_weights_1d_table;
  double** gauss_nodes_1d_table;
  double** gauss_weights_1d_table;

  double** mij_1d_table;
  double** dij_1d_table;
  double** invmij_1d_table;
  double** invvij_1d_table;
  double** lift_1d_table;
  double** flip_1d_table;
  double** vtk_interp_1d_table;
  
  double** vtk_rst_2d_table;
  double** vtk_rst_3d_table;
  double** gauss_rst_3d_table;
  double** gauss_rst_2d_table;
  double** lobatto_rst_3d_table;
  double** lobatto_rst_2d_table;
  
  double*** lobatto_to_gauss_interp_1d_table;
  double*** lobatto_to_gauss_interp_trans_1d_table;
  double*** hp_prolong_1d_table;
  double*** p_prolong_1d_table;
  double*** hp_restrict_1d_table;
  double*** p_restrict_1d_table;
  double*** hp_prolong_transpose_1d_table;
  double*** p_prolong_transpose_1d_table;   
  double*** hp_restrict_interp_1d_table;
  double*** p_restrict_interp_1d_table;
  
} d4est_operators_t;

typedef void (*d4est_operators_hp_restrict_fcn_t)(d4est_operators_t *, double *, int *,
                                         int, int, double *);

typedef void (*d4est_operators_p_restrict_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int, double *);

typedef void (*d4est_operators_hp_prolong_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int *, double *);

typedef void (*d4est_operators_p_prolong_fcn_t)(d4est_operators_t *, double *, int, int,
                                       int, double *);


/* This file was automatically generated.  Do not edit! */
double *d4est_operators_fetch_vtk_rst(d4est_operators_t *d4est_ops,int deg,int dim);
void d4est_operators_convert_nodal_to_vtk(d4est_operators_t *d4est_ops,double *vec,int dim,int deg,double *vtk_vec);
void d4est_operators_apply_vtk_interp(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int c,double *out);
void d4est_operators_reorient_face_data(d4est_operators_t *d4est_ops,double *in,int face_dim,int deg,int o,int f_m,int f_p,double *out);
void d4est_operators_apply_flip(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
void d4est_operators_apply_p_restrict_interp(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
void d4est_operators_apply_hp_restrict_interp(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
void d4est_operators_apply_p_prolong_transpose(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
void d4est_operators_apply_hp_prolong_transpose(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
void d4est_operators_build_p_prolong_transpose_1d_inverse(d4est_operators_t *d4est_ops,double *p_prolong_transpose_1d_inverse,int degH,int degh);
double *d4est_operators_fetch_p_prolong_transpose_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_build_p_prolong_1d_inverse(d4est_operators_t *d4est_ops,double *p_prolong_1d_inverse,int degH,int degh);
void d4est_operators_convert_nodal_to_modal(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
void d4est_operators_apply_slicer(d4est_operators_t *d4est_ops,double *in,int dim,int face,int deg,double *out);
void d4est_operators_apply_lift(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int face,double *out);
double *d4est_operators_fetch_gauss_rst_nd(d4est_operators_t *d4est_ops,int dim,int deg,int dir);
void d4est_operators_apply_dij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
double *d4est_operators_fetch_lobatto_rst_nd(d4est_operators_t *d4est_ops,int dim,int deg,int dir);
void d4est_operators_apply_hp_restrict(d4est_operators_t *d4est_ops,double *in,int *degh,int dim,int degH,double *out);
void d4est_operators_apply_p_restrict(d4est_operators_t *d4est_ops,double *in,int degh,int dim,int degH,double *out);
double *d4est_operators_fetch_p_restrict_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_build_hp_restrict_1d_aux(int degh,int degH,double *hp_prolong_matrix_1d,double *mass_matrix_rs_degh,double *inv_mass_matrix_rs_degH,double *hp_restrict_matrix_1d);
double *d4est_operators_fetch_p_prolong_1d(d4est_operators_t *d4est_ops,int degH,int degh);
void d4est_operators_hp_apply_nd_prolong_transpose_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_prolong_transpose_matrix_1d);
void d4est_operators_build_p_prolong_1d(d4est_operators_t *d4est_ops,double *p_prolong_1d,int degH,int degh);
void d4est_operators_apply_dij_transpose(d4est_operators_t *d4est_ops,double *in,int dim,int deg,int dir,double *out);
void d4est_operators_apply_invmij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
void d4est_operators_apply_mij(d4est_operators_t *d4est_ops,double *in,int dim,int deg,double *out);
double *d4est_operators_fetch_mij_1d(d4est_operators_t *d4est_ops,int deg);
double *d4est_operators_fetch_gauss_weights_1d(d4est_operators_t *d4est_ops,int deg);
double *d4est_operators_fetch_lobatto_weights_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_mij_1d(d4est_operators_t *d4est_ops,double *mij_1d,int deg);
void d4est_operators_hp_apply_nd_restrict_with_ptr(double *uH,int degH,double *uh,int degh,int dim,int c,double *hp_restrict_matrix_1d);
void d4est_operators_compute_PT_mat_P(d4est_operators_t *d4est_ops,double *mat,int degH,int dim,int *degh,int children,double *PT_mat_P);
void d4est_operators_apply_p_prolong(d4est_operators_t *d4est_ops,double *in,int degH,int dim,int degh,double *out);
void d4est_operators_apply_hp_prolong(d4est_operators_t *d4est_ops,double *in,int degH,int dim,int *degh,double *out);
void d4est_operators_compute_prolong_matrix(d4est_operators_t *d4est_ops,int degH,int dim,int *degh,int children,double *prolong_mat);
double *d4est_operators_fetch_lobatto_to_gauss_interp_trans_1d(d4est_operators_t *d4est_ops,int deg_lobatto,int deg_gauss);
void d4est_operators_build_lobatto_to_gauss_interp_trans_1d(d4est_operators_t *d4est_ops,double *ref_lobatto_to_gauss_interp_trans_1d,int lobatto_degree,int gauss_degree);
double *d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_operators_t *d4est_ops,int deg_lobatto,int deg_gauss);
void d4est_operators_build_custom_lobatto_interp_1d(d4est_operators_t *d4est_ops,double *custom_lobatto_interp_1d,int lobatto_degree,int custom_degree,double *custom_points);
double *d4est_operators_fetch_gauss_nodes_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_lobatto_to_gauss_interp_1d(d4est_operators_t *d4est_ops,double *ref_lobatto_to_gauss_interp_1d,int lobatto_degree,int gauss_degree);
void d4est_operators_hp_apply_nd_prolong_with_ptr(double *Uh,int degh,double *UH,int degH,int dim,int c,double *hp_prolong_matrix_1d);
double *d4est_operators_fetch_invvij_1d(d4est_operators_t *d4est_ops,int deg);
double *d4est_operators_fetch_lobatto_nodes_1d(d4est_operators_t *d4est_ops,int deg);
void d4est_operators_build_Vij_1d(d4est_operators_t *d4est_ops,double *Vij_1d,int deg);
double *d4est_operators_1index_2d_3d_fetch(d4est_operators_t *d4est_ops,int deg,int dim,int size,double **table_2d,double **table_3d,void(*build_fcn)(d4est_operators_t *,double *,int,int));
double *d4est_operators_2index_fetch(d4est_operators_t *d4est_ops,double ***table,int deg1,int deg2,int size,void(*build_fcn)(d4est_operators_t *,double *,int,int));
double *d4est_operators_1index_fetch(d4est_operators_t *d4est_ops,double **table,int deg,int size,void(*build_fcn)(d4est_operators_t *,double *,int));
void d4est_ops_destroy(d4est_operators_t *d4est_ops);
d4est_operators_t *d4est_ops_init(int deg_max);


#endif
