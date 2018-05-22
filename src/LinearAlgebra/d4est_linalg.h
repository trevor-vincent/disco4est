#ifndef LINALG_H
#define LINALG_H

#include <d4est_blas.h>

double d4est_linalg_vec1_trans_mat_vec2(double *vec1,double *mat,double *vec2,int N);
double d4est_linalg_vec_sum(double *vec,int N);
void d4est_linalg_cross_prod(double ax,double ay,double az,double bx,double by,double bz,double *axb_x,double *axb_y,double *axb_z);
void d4est_linalg_component_div(double *x,double *y,double *xdivy,int N);
void d4est_linalg_component_mult(double *x,double *y,double *xy,int N);
void d4est_linalg_vec_fabs(double *x,int N);
void d4est_linalg_vec_fabsdiff(double* x, double* y, double* result, int N);
void d4est_linalg_sym_eigvecs(double *A,double *eig_vals,int N);
void d4est_linalg_sym_eigvals(double *A,double *eig_vals,int N);
void d4est_util_fill_array(double *v,double val,int N);
void d4est_linalg_vec_xpby(double *x,double beta,double *y,int N);
void d4est_linalg_vec_axpyeqz(double alpha,double *x,double *y,double *z,int N);
void d4est_linalg_vec_axpy(double alpha,double *x,double *y,int N);
double d4est_linalg_vec_dot(double *x,double *y,int N);
void d4est_linalg_vec_normalize(double *x,int N);
void d4est_linalg_vec_scale(double alpha,double *x,int N);
void d4est_linalg_set_column_opt(double *A,double *column,int col,int N,int M);
void d4est_linalg_set_column(double *A,double *column,int col,int N,int M);
void d4est_linalg_mat_transpose(double *A,double *Atrans,int N);
void d4est_linalg_matvec_plus_vec(double alpha,double *A,double *v,double beta,double *b,int m,int n);
void d4est_linalg_mat_multiply(double *A,double *B,double *C,int m,int l,int n);
void d4est_linalg_mat_transpose_nonsqr(double *A,double *A_transpose,int A_rows,int A_cols);
void d4est_linalg_leftinverse(double *A,double *inv_A,int A_rows,int A_cols);
void d4est_linalg_invert_and_copy(double *A,double *inv_A,int m);
void d4est_linalg_invert(double *A,int m);

#endif
