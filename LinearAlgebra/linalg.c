#include <linalg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// void Kron::A1oA2x(double *A1oA2x, double *A1, double *A2, double *X,
//                   int a1_rows, int a1_cols, int a2_rows, int a2_cols) {
  
//   double *aux = (double *)malloc(sizeof(double) * a1_cols * a2_rows);
//   // Compute A2*(X**T)
//   dgemm_('T', 'N', a1_cols, a2_rows, a2_cols, 1., X, a2_cols, A2, a2_cols, 0.,
//          aux, a1_cols);

//   // Compute A1*((A2*(X**T))**T)
//   dgemm_('T', 'N', a2_rows, a1_rows, a1_cols, 1., aux, a1_cols, A1, a1_cols, 0.,
//          A1oA2x, a2_rows);
//   free(aux);
// }

// void Kron::A1oA2oA3x(double *A1oA2oA3x, double *A1, double *A2, double *A3,
//                      double *X, int a1_rows, int a1_cols, int a2_rows,
//                      int a2_cols, int a3_rows, int a3_cols) {
  
//   double *tmp = (double *)malloc(sizeof(double) * a2_rows * a3_rows * a1_cols);
//   for (int i = 0; i < a1_cols; i++) {
//     Kron::A1oA2x(&tmp[i * a2_rows * a3_rows], A2, A3, &X[i * a2_cols * a3_cols],
//                  a2_rows, a2_cols, a3_rows, a3_cols);
//   }
//   int a2a3rows = a2_rows * a3_rows;
//   dgemm_('N', 'N', a2a3rows, a1_rows, a1_cols, 1., tmp, a2a3rows, A1, a1_cols,
//          0., A1oA2oA3x, a2a3rows);
//   free(tmp);
// }

/** 
 * Calculates C = A \otimes B. 
 * Only used for testing purposes
 *
 * @param A 
 * @param B 
 * @param C a_rows*b_rows by a_cols*b_cols
 * @param a_rows 
 * @param a_cols 
 * @param b_rows 
 * @param b_cols 
 */
void
linalg_kron_AoB (double *A, double *B, double *C, int a_rows, int a_cols,
                 int b_rows, int b_cols)
{
  int                 size2 = a_cols * b_cols;
  /* printf("size1, size2, a_rows, a_cols, b_rows, b_cols = %d,%d,%d,%d,%d,%d\n", size1, size2, a_rows, a_cols, b_rows, b_cols); */
  int                 i, j, k, p;
  for (i = 0; i < a_rows; i++)
    for (j = 0; j < a_cols; j++)
      for (k = 0; k < b_rows; k++)
        for (p = 0; p < b_cols; p++) {
          /* printf("(k + i*b_rows)*size2 + p + j*b_cols = %d\n", (k + i*b_rows)*size2 + p + j*b_cols); */
          C[(k + i * b_rows) * size2 + p + j * b_cols] =
            A[i * a_cols + j] * B[k * b_cols + p];
        }
}


void
linalg_kron_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvec_dot_x)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_x[stride] = vec[i] * vec[k] * x[stride];
    }
}
void
linalg_kron_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_xy[stride] = vec[i] * vec[k] * x[stride] * y[stride];
    }
}

void
linalg_kron_vec_o_vec_dot_wxyz (double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vecvec_dot_wxyz)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_wxyz[stride] = vec[i] * vec[k] * w[stride] * x[stride] * y[stride] * z[stride];
    }
}


void
linalg_kron_oneover_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_xy[stride] = (1./(vec[i] * vec[k] * x[stride])) * y[stride];
    }
}

void
linalg_kron_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_xy[i] = vec[i] * x[i] * y[i];
  }
}


void
linalg_kron_vec_dot_wxyz (double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vec_dot_wxyz)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_wxyz[i] = vec[i] * w[i] * x[i] * y[i] * z[i];
  }
}


void
linalg_kron_vec_dot_x (double *vec, double*x,int vec_size, double* vec_dot_x)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_x[i] = vec[i] * x[i];
  }
}

void
linalg_kron_oneover_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy)
{
  for (int i = 0; i < vec_size; i++){
    vec_dot_xy[i] = (1./(vec[i] * x[i])) * y[i];
  }
}

void
linalg_kron_vec_o_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy)
{
  int b_rows = vec_size*vec_size;
  
  for (int m = 0; m < vec_size; m++)
   for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int n = k + i * vec_size;
      vecvecvec_dot_xy[(n + m * b_rows)] = vec[m] * vec[i] * vec[k] * x[(n + m * b_rows)] * y[(n + m * b_rows)];
    }
}

void
linalg_kron_vec_o_vec_o_vec_dot_wxyz(double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vecvecvec_dot_wxyz)
{
  int b_rows = vec_size*vec_size;
  
  for (int m = 0; m < vec_size; m++)
   for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int n = k + i * vec_size;
      vecvecvec_dot_wxyz[(n + m * b_rows)] =
        vec[m]
        * vec[i]
        * vec[k]
        * w[(n + m * b_rows)]
        * x[(n + m * b_rows)]
        * y[(n + m * b_rows)]
        * z[(n + m * b_rows)]
        ;
    }
}


void
linalg_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy)
{
  int b_rows = vec_size*vec_size;
  
  for (int m = 0; m < vec_size; m++)
   for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int n = k + i * vec_size;
      vecvecvec_dot_xy[(n + m * b_rows)] = (1./(vec[m] * vec[i] * vec[k] * x[(n + m * b_rows)])) * y[(n + m * b_rows)];
    }
}

void
linalg_kron_vec_o_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvecvec_dot_x)
{
  int b_rows = vec_size*vec_size;
  
  for (int m = 0; m < vec_size; m++)
   for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int n = k + i * vec_size;
      vecvecvec_dot_x[(n + m * b_rows)] = vec[m] * vec[i] * vec[k] * x[(n + m * b_rows)];
    }
}


/** 
 * Calculates D = A \otimes B \otimes C.
 * Only used for testing purposes. 
 *
 * @param A 
 * @param B 
 * @param C 
 * @param D a_rows*b_rows*c_rows by a_cols*b_cols*c_cols
 * @param a_rows 
 * @param a_cols 
 * @param b_rows 
 * @param b_cols 
 * @param c_rows 
 * @param c_cols 
 */
void
linalg_kron_AoBoC (double *A, double *B, double *C, double *D, int a_rows,
                   int a_cols, int b_rows, int b_cols, int c_rows, int c_cols)
{
  double             *AoB =
    (double*) malloc (sizeof (double) * a_rows * a_cols * b_rows * b_cols);

  linalg_kron_AoB (A, B, AoB, a_rows, a_cols, b_rows, b_cols);
  linalg_kron_AoB (AoB, C, D, a_rows * b_rows, a_cols * b_cols, c_rows,
                   c_cols);

  free (AoB);
}


void
linalg_kron_A1A2x_NONSQR (double *A1A2x, double *A1, double *A2, double *X,
                          int a1_rows, int a1_cols, int a2_rows, int a2_cols)
{
  /* X should be a1_cols * a2_cols */
  /* A1A2x should be a1_rows * a2_rows */

  /* This can be optimized a lot, see dg-charm */
  
  int                 N =
    (a1_rows * a2_rows >
     a2_rows * a1_cols) ? a1_rows * a2_rows : a2_rows * a1_cols;
  N = (a1_cols * a2_cols > N) ? a1_cols * a2_cols : N;

  double             *tmp = (double *) malloc (sizeof (double) * N);
  double             *tmp1 = (double *) malloc (sizeof (double) * N);

  linalg_mat_transpose_nonsqr (X, tmp1, a1_cols, a2_cols);
  linalg_mat_multiply (A2, tmp1, tmp, a2_rows, a2_cols, a1_cols);
  linalg_mat_transpose_nonsqr (tmp, tmp1, a2_rows, a1_cols);
  linalg_mat_multiply (A1, tmp1, A1A2x, a1_rows, a1_cols, a2_rows);

  free (tmp);
  free (tmp1);
}

void
linalg_kron_MAToIx_SQR (double *MAToIx, double *MAT, double *X, int N)
{
  linalg_mat_multiply (MAT, X, MAToIx, N, N, N);
}

void
linalg_kron_IoMATx_SQR (double *IoMATx, double *MAT, double *X, int N)
{
  double             *tmp = (double *) malloc (sizeof (double) * N * N);
  double             *tmp1 = (double *) malloc (sizeof (double) * N * N);

  linalg_mat_transpose_nonsqr (X, tmp1, N, N);
  linalg_mat_multiply (MAT, tmp1, tmp, N, N, N);
  linalg_mat_transpose_nonsqr (tmp, IoMATx, N, N);

  free (tmp);
  free (tmp1);
}

void
linalg_kron_VECoIx_SQR (double *VECoIx, double *VEC, double *X, int N)
{

  linalg_mat_multiply (VEC, X, VECoIx, N, 1, N);

}

void
linalg_kron_VEC_TRANSoIx_SQR (double *VEC_TRANSoIx, double *VEC, double *X,
                              int N)
{
  linalg_mat_multiply (VEC, X, VEC_TRANSoIx, 1, N, N);
}

void
linalg_kron_IoVECx_SQR (double *IoVECx, double *VEC, double *X, int N)
{
  /* int N = (a1_rows*a2_rows > a2_rows*a1_cols) ? a1_rows*a2_rows : a2_rows*a1_cols; */

  double             *tmp = (double *) malloc (sizeof (double) * N * N);
  /* double* tmp1 = (double*)malloc(sizeof(double)*N); */

  linalg_mat_multiply (VEC, X, tmp, N, 1, N);
  linalg_mat_transpose_nonsqr (tmp, IoVECx, N, N);

  free (tmp);
  /* free(tmp1);   */
}

void
linalg_kron_IoVEC_TRANSx_SQR (double *IoVEC_TRANSx, double *VEC, double *X,
                              int N)
{
  double             *tmp1 = (double *) malloc (sizeof (double) * N * N);

  linalg_mat_transpose_nonsqr (X, tmp1, N, N);
  linalg_mat_multiply (VEC, tmp1, IoVEC_TRANSx, 1, N, N);

  free (tmp1);
}

void
linalg_kron_A1A2A3x_NONSQR (double *A1A2A3x, double *A1, double *A2,
                            double *A3, double *X, int a1_rows, int a1_cols,
                            int a2_rows, int a2_cols, int a3_rows,
                            int a3_cols)
{
  double             *tmp =
    (double *) malloc (sizeof (double) * a2_rows * a3_rows * a1_cols);
  int                 i;
  for (i = 0; i < a1_cols; i++) {
    linalg_kron_A1A2x_NONSQR (&tmp[i * a2_rows * a3_rows], A2, A3,
                              &X[i * a2_cols * a3_cols], a2_rows, a2_cols,
                              a3_rows, a3_cols);
  }

  linalg_mat_multiply (A1, tmp, A1A2A3x, a1_rows, a1_cols, a2_rows * a3_rows);
  free (tmp);
}

void
linalg_kron_MAToIoIx_SQR (double *MAToIoIx, double *MAT, double *X, int N)
{
  linalg_mat_multiply (MAT, X, MAToIoIx, N, N, N * N);
}

void
linalg_kron_IoMAToIx_SQR (double *IoMAToIx, double *MAT, double *X, int N)
{
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_MAToIx_SQR (&IoMAToIx[i * N * N], MAT, &X[i * N * N], N);
  }
}

void
linalg_kron_IoIoMATx_SQR (double *IoIoMATx, double *MAT, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_IoMATx_SQR (&IoIoMATx[i * N * N], MAT, &X[i * N * N], N);
  }
}

void
linalg_kron_VECoIoIx_SQR (double *VECoIoIx, double *VEC, double *X, int N)
{
  linalg_mat_multiply (VEC, X, VECoIoIx, N, 1, N * N);
}

void
linalg_kron_VEC_TRANSoIoIx_SQR (double *VEC_TRANSoIoIx, double *VEC,
                                double *X, int N)
{
  linalg_mat_multiply (VEC, X, VEC_TRANSoIoIx, 1, N, N * N);
}

void
linalg_kron_IoVECoIx_SQR (double *IoVECoIx, double *VEC, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_VECoIx_SQR (&IoVECoIx[i * N * N], VEC, &X[i * N], N);
  }
}

void
linalg_kron_IoVEC_TRANSoIx_SQR (double *IoVEC_TRANSoIx, double *VEC,
                                double *X, int N)
{
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_VEC_TRANSoIx_SQR (&IoVEC_TRANSoIx[i * N], VEC, &X[i * N * N],
                                  N);
  }
}

void
linalg_kron_IoIoVECx_SQR (double *IoIoVECx, double *VEC, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_IoVECx_SQR (&IoIoVECx[i * N * N], VEC, &X[i * N], N);
  }
}

void
linalg_kron_IoIoVEC_TRANSx_SQR (double *IoIoVEC_TRANSx, double *VEC,
                                double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    linalg_kron_IoVEC_TRANSx_SQR (&IoIoVEC_TRANSx[i * N], VEC, &X[i * N * N],
                                  N);
  }
}

void
linalg_invert (double *A, int m)
{
  int                 N = m;
  int                 lwork = N * N;

  int                *ipiv = (int *) malloc (sizeof (int) * (N + 1));
  double             *work = (double *) malloc (sizeof (double) * lwork);
  int                 info;

  dgetrf_ (&N, &N, A, &N, ipiv, &info);
  dgetri_ (&N, A, &N, ipiv, work, &lwork, &info);

  free (ipiv);
  free (work);
}

void
linalg_invert_and_copy (double *A, double *inv_A, int m)
{
  int                 i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      inv_A[i * m + j] = A[i * m + j];
  linalg_invert (inv_A, m);
}

void
linalg_leftinverse(double *A, double* inv_A, int A_rows, int A_cols)
{
  double* AT = (double *) malloc (sizeof (double) * A_rows * A_cols);
  double* ATA = (double *) malloc (sizeof (double) * A_cols * A_cols);
  double* inv_ATA = (double *) malloc (sizeof (double) * A_cols * A_cols);
  
  linalg_mat_transpose_nonsqr(A, AT, A_rows, A_cols);
  linalg_mat_multiply(AT, A, ATA, A_cols, A_rows, A_cols);
  linalg_invert_and_copy(ATA, inv_ATA, A_cols);
  linalg_mat_multiply(inv_ATA, AT, inv_A, A_cols, A_cols, A_rows);
  
  free(AT);
  free(ATA);
  free(inv_ATA);
}


/** 
 I'll bring  * C = A*B
 * 
 * @param A matrix of dimension m by l
 * @param B matrix of dimension l by n
 * @param C matrix of dimension m by n
 * @param m rows(A)
 * @param n cols(B)
 * @param l rows(B)
 */
void
linalg_mat_multiply (double *A, double *B, double *C, int m, int l, int n)
{
  char                ntran = 'N';
  double              one = 1.;
  double              zero = 0.;

  int                 dimM = m;
  int                 dimN = n;
  int                 dimL = l;

  dgemm_ (&ntran, &ntran, &dimN, &dimM, &dimL, &one, B, &dimN, A, &dimL,
          &zero, C, &dimN);
}

void
linalg_matvec_plus_vec (double alpha, double *A, double *v, double beta,
                        double *b, int m, int n)
{
  char                ytran = 'T';
  int                 dimM = m;
  int                 dimN = n;

  int                 incb = 1;
  int                 incv = 1;

  /* printf("\n\nmatvec plus vec\n\n"); */
  /* util_print_matrix(A, m, n, "A = ", 0); */
  /* util_print_matrix(v, n, 1, "v = ", 0); */

  dgemv_ (&ytran, &dimN, &dimM, &alpha, A, &dimN, v, &incv, &beta, b, &incb);
}

void
linalg_copy_1st_to_2nd (double *v1, double *v2, int N)
{
  memcpy (v2, v1, N * sizeof (double));
}

/*
void linalg_mat_transpose(double* A, double* Atrans, int N, int M)
{
  for(int n = 0; n<N*M; n++){
    int i = n/N;
    int j = n%N;
    Atrans[n] = A[M*j + i];
  }
}
*/

void
linalg_mat_transpose (double *A, double *Atrans, int N)
{
  int                 i, j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      Atrans[i * N + j] = A[j * N + i];
}

/** 
 * 
 * 
 * @param A 
 * @param A_transpose 
 * @param A_rows 
 * @param A_cols 
 */
void
linalg_mat_transpose_nonsqr (double *A, double *A_transpose, int A_rows,
                             int A_cols)
{
  int                 i, j;
  for (i = 0; i < A_rows; i++)
    for (j = 0; j < A_cols; j++)
      A_transpose[j * A_rows + i] = A[i * A_cols + j];
}

/** 
 * Only for debugging purposes, or small matrices.
 * This should not be called many times.
 * 
 * @param A 
 * @param column 
 * @param N rows
 * @param M columns
 */
void
linalg_set_column (double *A, double *column, int col, int N, int M)
{
  assert (col < M);
  int                 s = 0;
  int                 i, j, n;
  for (n = 0; n < N * M; n++) {
    i = n / M;
    j = n % M;
    if (j == col)
      A[M * i + j] = column[s++];
  }
}

void
linalg_set_column_opt (double *A, double *column, int col, int N, int M)
{
  assert (col < M);
  int                 s = 0;

  for (int i = 0; i < N; i++){
    A[M*i + col] = column[s++];
  }
}



void
linalg_vec_scale (double alpha, double *x, int N)
{
  int                 incx = 1;
  dscal_ (&N, &alpha, x, &incx);
}

void
linalg_vec_normalize(double*x, int N)
{
  double xdotx = linalg_vec_dot(x,x,N);
  double norm = sqrt(xdotx);
  linalg_vec_scale(norm, x, N); 
}

void
linalg_vec_axpy (double alpha, double *x, double *y, int N)
{
  int                 incx = 1;
  int                 incy = 1;
  daxpy_ (&N, &alpha, x, &incx, y, &incy);
}

void
linalg_vec_axpyeqz (double alpha, double *x, double *y, double *z, int N)
{
  linalg_copy_1st_to_2nd (y, z, N);
  linalg_vec_axpy (alpha, x, z, N);
}

void
linalg_vec_xpby (double *x, double beta, double *y, int N)
{
  linalg_vec_scale (beta, y, N);
  linalg_vec_axpy (1., x, y, N);
}

double
linalg_vec_dot (double *x, double *y, int N)
{
  int                 incx = 1;
  int                 incy = 1;
  return ddot_ (&N, x, &incx, y, &incy);
}

void
linalg_fill_vec (double *v, double val, int N)
{
  int                 i;
  for (i = 0; i < N; i++)
    v[i] = val;
}

/** 
 * This routine will work irrespective of the 
 * data order (COL or ROW major) b.c the matrix is
 * symmetric
 * 
 * @param A an NxN matrix
 * @param eig_vals vector of eigen values 
 * @param N 
 */
void
linalg_sym_eigvals (double *A, double *eig_vals, int N)
{
  const int           rows = N;
  char                no_vec = 'N';
  char                upper = 'U';
  int                 info = 0;
  int                 lwork = N * N;
  double             *work = (double *) malloc (sizeof (double) * lwork);
  dsyev_ (&no_vec, &upper, &rows, A, &rows, eig_vals, work, &lwork, &info);
  free (work);
}

void
linalg_sym_eigvecs (double *A, double *eig_vals, int N)
{
  const int           rows = N;
  char                jobz = 'V';
  char                upper = 'U';
  int                 info = 0;
  int                 lwork = N * N;
  double             *work = (double *) malloc (sizeof (double) * lwork);
  dsyev_ (&jobz, &upper, &rows, A, &rows, eig_vals, work, &lwork, &info);
  free (work);
}

int
linalg_compare_vecs (double *v1, double *v2, double tol, int N)
{
  int                 i;
  double              err = 0.;
  for (i = 0; i < N; i++) {
    err += fabs (v1[i] - v2[i]);
  }
  return (err < tol);
}

void
linalg_vec_fabs(double* x, int N)
{
  int i;
  for (i = 0; i < N; i++) {
    x[i] = fabs(x[i]);
  }
}

void
linalg_component_mult
(
 double* x,
 double* y,
 double* xy,
 int N
)
{
  int i;
  for (i = 0; i < N; i++) {
    xy[i] = x[i]*y[i];
  }
}


void
linalg_component_div
(
 double* x,
 double* y,
 double* xdivy,
 int N
)
{
  int i;
  for (i = 0; i < N; i++) {
    if (y[i] == 0)
      printf("[WARNING]: COMPONENT DIVISION BY ZERO\n");
    xdivy[i] = x[i]/y[i];
  }
}


/** 
 * return v1 \cross v2
 * 
 * @param v1
 * @param v2 
 */
void
linalg_cross_prod
(
 double ax,
 double ay,
 double az,
 double bx,
 double by,
 double bz,
 double* axb_x,
 double* axb_y,
 double* axb_z
)
{
  *axb_x = ay*bz - az*by;
  *axb_y = az*bx - ax*bz;
  *axb_z = ax*by - ay*bx;
}

double
linalg_vec_sum
(
 double* vec,
 int N
)
{
  int i;
  double sum = 0.;
  for (i = 0; i < N; i++) {
    sum += vec[i];
  }
  return sum;
}


double
linalg_vec1_trans_mat_vec2
(
 double* vec1,
 double* mat,
 double* vec2,
 int N
)
{

  double *mat_vec2 = (double *) malloc (sizeof (double) * N);
  linalg_matvec_plus_vec(1., mat, vec2, 0., mat_vec2, N, N);                      
  double dot = linalg_vec_dot(vec1, mat_vec2, N);
  free(mat_vec2);
  return dot;
}
