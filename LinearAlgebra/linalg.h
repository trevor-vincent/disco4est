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

void linalg_invert_and_copy(double *, double *, int);
void linalg_invert(double *, int);
void linalg_mat_multiply(double *, double *, double *, int, int, int);
void linalg_matvec_plus_vec(double, double *, double *, double, double *, int,
                            int);
void linalg_copy_1st_to_2nd(double *, double *, int);
void linalg_mat_transpose(double *, double *, int);
void linalg_vec_scale(double, double *, int);
void linalg_vec_axpy(double, double *, double *, int);
void linalg_vec_xpby(double *, double, double *, int);
double linalg_vec_dot(double *, double *, int);
void linalg_fill_vec(double *, double, int);
void linalg_sym_eigvals(double *, double *, int);
void linalg_sym_eigvecs(double *, double *, int);
void linalg_set_column(double *, double *, int, int, int);
void linalg_set_column_opt(double *A, double *column, int col, int N, int M);
void linalg_mat_transpose_nonsqr(double *, double *, int, int);
int linalg_compare_vecs(double *, double *, double, int);
void linalg_kron_A1A2x_NONSQR(double *A1A2x, double *A1, double *A2, double *X,
                              int a1_rows, int a1_cols, int a2_rows,
                              int a2_cols);
void linalg_kron_AoB(double *A, double *B, double *C, int a_rows, int a_cols,
                     int b_rows, int b_cols);
void linalg_kron_VECoIx_SQR(double *VECoIx, double *VEC, double *X, int N);
void linalg_kron_IoVECx_SQR(double *IoVECx, double *VEC, double *X, int N);
void linalg_kron_VECoIoIx_SQR(double *VECoIoIx, double *VEC, double *X, int N);
void linalg_kron_IoVECoIx_SQR(double *IoVECoIx, double *VEC, double *X, int N);
void linalg_kron_IoIoVECx_SQR(double *IoIoVECx, double *A, double *X, int N);
void linalg_kron_AoBoC(double *A, double *B, double *C, double *D, int a_rows,
                       int a_cols, int b_rows, int b_cols, int c_rows,
                       int c_cols);
void linalg_kron_IoMATx_SQR(double *IoMatx, double *MAT, double *X, int N);
void linalg_kron_MAToIx_SQR(double *MAToIx, double *MAT, double *X, int N);
void linalg_kron_MAToIoIx_SQR(double *MAToIoIx, double *MAT, double *X, int N);
void linalg_kron_IoMAToIx_SQR(double *IoMAToIx, double *MAT, double *X, int N);
void linalg_kron_IoIoMATx_SQR(double *IoIoMATx, double *MAT, double *X, int N);
void linalg_kron_A1A2A3x_NONSQR(double *A1A2A3x, double *A1, double *A2,
                                double *A3, double *X, int a1_rows, int a1_cols,
                                int a2_rows, int a2_cols, int a3_rows,
                                int a3_cols);
void linalg_kron_VEC_TRANSoIx_SQR(double *VEC_TRANSoIx, double *VEC, double *X,
                                  int N);
void linalg_kron_IoVEC_TRANSx_SQR(double *IoVEC_TRANSx, double *VEC, double *X,
                                  int N);
void linalg_kron_VEC_TRANSoIoIx_SQR(double *VEC_TRANSoIoIx, double *VEC,
                                    double *X, int N);
void linalg_kron_IoVEC_TRANSoIx_SQR(double *IoVEC_TRANSoIx, double *VEC,
                                    double *X, int N);
void linalg_kron_IoIoVEC_TRANSx_SQR(double *IoIoVEC_TRANSx, double *VEC,
                                    double *X, int N);
double linalg_Euc_norm(double *x, int local_nodes);
void linalg_vec_axpyeqz(double alpha, double *x, double *y, double *z, int N);

void linalg_component_mult(double *x, double *y, double *xy, int N);

void linalg_component_div(double *x, double *y, double *xdivy, int N);

void linalg_cross_prod(double ax, double ay, double az, double bx, double by,
                       double bz, double *axb_x, double *axb_y, double *axb_z);

double linalg_vec_sum(double *vec, int N);

void linalg_vec_fabs(double *x, int N);

void linalg_vec_normalize(double *x, int N);

void linalg_kron_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvec_dot_x);
void linalg_kron_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy);
void linalg_kron_vec_o_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy);
void linalg_kron_vec_o_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvecvec_dot_x);

void
linalg_kron_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy);

void
linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy);

void
linalg_kron_oneover_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy);
void
linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy);

void
linalg_kron_vec_dot_x (double *vec, double*x,int vec_size, double* vec_dot_x);

void
linalg_leftinverse(double *A, double* inv_A, int A_rows, int A_cols);

#endif
