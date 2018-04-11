#ifndef D4EST_BLAS_H
#define D4EST_BLAS_H 

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

#endif
