#ifndef LINALG_H
#define LINALG_H

void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
            const int *K, double *ALPHA, double *A, const int *LDA, double *B,
            const int *LDB, double *BETA, double *C, const int *LDC);
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);
void daxpy_(const int *N, double *alpha, double *x, const int *incx, double *y,
            const int *incy);
void dgemv_(const char *TRANSA, const int *m, const int *n, double *alpha,
            double *A, const int *k, double *v, const int *incv, double *beta,
            double *b, const int *incb);
void dscal_(const int *N, double *alpha, double *x, const int *inc);
double ddot_(const int *N, double *x, const int *incx, double *y,
             const int *incy);
void dsyev_(const char *, const char *, const int *, double *, const int *,
            double *, double *, int *, int *);

/* This file was automatically generated.  Do not edit! */
double linalg_vec1_trans_mat_vec2(double *vec1,double *mat,double *vec2,int N);
double linalg_vec_sum(double *vec,int N);
void linalg_cross_prod(double ax,double ay,double az,double bx,double by,double bz,double *axb_x,double *axb_y,double *axb_z);
void linalg_component_div(double *x,double *y,double *xdivy,int N);
void linalg_component_mult(double *x,double *y,double *xy,int N);
void linalg_vec_fabs(double *x,int N);
int linalg_compare_vecs(double *v1,double *v2,double tol,int N);
void linalg_sym_eigvecs(double *A,double *eig_vals,int N);
void linalg_sym_eigvals(double *A,double *eig_vals,int N);
void linalg_fill_vec(double *v,double val,int N);
void linalg_vec_xpby(double *x,double beta,double *y,int N);
void linalg_vec_axpyeqz(double alpha,double *x,double *y,double *z,int N);
void linalg_vec_axpy(double alpha,double *x,double *y,int N);
double linalg_vec_dot(double *x,double *y,int N);
void linalg_vec_normalize(double *x,int N);
void linalg_vec_scale(double alpha,double *x,int N);
void linalg_set_column_opt(double *A,double *column,int col,int N,int M);
void linalg_set_column(double *A,double *column,int col,int N,int M);
void linalg_mat_transpose(double *A,double *Atrans,int N);
void linalg_copy_1st_to_2nd(double *v1,double *v2,int N);
void linalg_matvec_plus_vec(double alpha,double *A,double *v,double beta,double *b,int m,int n);
void linalg_leftinverse(double *A,double *inv_A,int A_rows,int A_cols);
void linalg_invert_and_copy(double *A,double *inv_A,int m);
void linalg_invert(double *A,int m);
void linalg_kron_IoIoVEC_TRANSx_SQR(double *IoIoVEC_TRANSx,double *VEC,double *X,int N);
void linalg_kron_IoIoVECx_SQR(double *IoIoVECx,double *VEC,double *X,int N);
void linalg_kron_IoVEC_TRANSoIx_SQR(double *IoVEC_TRANSoIx,double *VEC,double *X,int N);
void linalg_kron_IoVECoIx_SQR(double *IoVECoIx,double *VEC,double *X,int N);
void linalg_kron_VEC_TRANSoIoIx_SQR(double *VEC_TRANSoIoIx,double *VEC,double *X,int N);
void linalg_kron_VECoIoIx_SQR(double *VECoIoIx,double *VEC,double *X,int N);
void linalg_kron_IoIoMATx_SQR(double *IoIoMATx,double *MAT,double *X,int N);
void linalg_kron_IoMAToIx_SQR(double *IoMAToIx,double *MAT,double *X,int N);
void linalg_kron_MAToIoIx_SQR(double *MAToIoIx,double *MAT,double *X,int N);
void linalg_kron_A1A2A3x_NONSQR(double *A1A2A3x,double *A1,double *A2,double *A3,double *X,int a1_rows,int a1_cols,int a2_rows,int a2_cols,int a3_rows,int a3_cols);
void linalg_kron_IoVEC_TRANSx_SQR(double *IoVEC_TRANSx,double *VEC,double *X,int N);
void linalg_kron_IoVECx_SQR(double *IoVECx,double *VEC,double *X,int N);
void linalg_kron_VEC_TRANSoIx_SQR(double *VEC_TRANSoIx,double *VEC,double *X,int N);
void linalg_kron_VECoIx_SQR(double *VECoIx,double *VEC,double *X,int N);
void linalg_kron_IoMATx_SQR(double *IoMATx,double *MAT,double *X,int N);
void linalg_kron_MAToIx_SQR(double *MAToIx,double *MAT,double *X,int N);
void linalg_mat_multiply(double *A,double *B,double *C,int m,int l,int n);
void linalg_mat_transpose_nonsqr(double *A,double *A_transpose,int A_rows,int A_cols);
void linalg_kron_A1A2x_NONSQR(double *A1A2x,double *A1,double *A2,double *X,int a1_rows,int a1_cols,int a2_rows,int a2_cols);
void linalg_kron_AoBoC(double *A,double *B,double *C,double *D,int a_rows,int a_cols,int b_rows,int b_cols,int c_rows,int c_cols);
void linalg_kron_vec_o_vec_o_vec_dot_x(double *vec,double *x,int vec_size,double *vecvecvec_dot_x);
void linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vecvecvec_dot_xy);
void linalg_kron_vec_o_vec_o_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vecvecvec_dot_wxyz);
void linalg_kron_vec_o_vec_o_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vecvecvec_dot_xy);
void linalg_kron_oneover_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vec_dot_xy);
void linalg_kron_vec_dot_x(double *vec,double *x,int vec_size,double *vec_dot_x);
void linalg_kron_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vec_dot_wxyz);
void linalg_kron_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vec_dot_xy);
void linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vecvec_dot_xy);
void linalg_kron_vec_o_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vecvec_dot_wxyz);
void linalg_kron_vec_o_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vecvec_dot_xy);
void linalg_kron_vec1_o_vec2_o_vec3_dot_wxyz(double *vec1,double *vec2,double *vec3,double *w,double *x,double *y,double *z,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void linalg_kron_vec1_o_vec2_o_vec3_dot_xy(double *vec1,double *vec2,double *vec3,double *x,double *y,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void linalg_kron_vec1_o_vec2_o_vec3_dot_x(double *vec1,double *vec2,double *vec3,double *x,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void linalg_kron_vec1_o_vec2_dot_wxyz(double *vec1,double *vec2,double *w,double *x,double *y,double *z,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void linalg_kron_vec1_o_vec2_dot_xy(double *vec1,double *vec2,double *x,double *y,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void linalg_kron_vec1_o_vec2_dot_x(double *vec1,double *vec2,double *x,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void linalg_kron_vec_o_vec_dot_x(double *vec,double *x,int vec_size,double *vecvec_dot_x);
void linalg_kron_AoB(double *A,double *B,double *C,int a_rows,int a_cols,int b_rows,int b_cols);

#endif
