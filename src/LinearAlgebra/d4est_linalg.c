#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <d4est_util.h>
#include <d4est_linalg.h>


void
d4est_linalg_invert (double * restrict A, int m)
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
d4est_linalg_invert_and_copy (double * restrict A, double * restrict inv_A, int m)
{
  int                 i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      inv_A[i * m + j] = A[i * m + j];
  d4est_linalg_invert (inv_A, m);
}

void
d4est_linalg_leftinverse(double * restrict A, double * restrict  inv_A, int A_rows, int A_cols)
{
  double* AT = (double *) malloc (sizeof (double) * A_rows * A_cols);
  double* ATA = (double *) malloc (sizeof (double) * A_cols * A_cols);
  double* inv_ATA = (double *) malloc (sizeof (double) * A_cols * A_cols);
  
  d4est_linalg_mat_transpose_nonsqr(A, AT, A_rows, A_cols);
  d4est_linalg_mat_multiply(AT, A, ATA, A_cols, A_rows, A_cols);
  d4est_linalg_invert_and_copy(ATA, inv_ATA, A_cols);
  d4est_linalg_mat_multiply(inv_ATA, AT, inv_A, A_cols, A_cols, A_rows);
  
  free(AT);
  free(ATA);
  free(inv_ATA);
}


/**
 C = A*B
 *
 * @param A matrix of dimension m by l
 * @param B matrix of dimension l by n
 * @param C matrix of dimension m by n
 * @param m rows(A)
 * @param n cols(B)
 * @param l rows(B)
 */
void
d4est_linalg_mat_multiply (double * restrict A, double * restrict B, double * restrict C, int m, int l, int n)
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
d4est_linalg_matvec_plus_vec (double alpha, double * restrict A, double * restrict v, double beta,
                        double * restrict b, int m, int n)
{
  char                ytran = 'T';
  int                 dimM = m;
  int                 dimN = n;

  int                 incb = 1;
  int                 incv = 1;

  /* printf("\n\nmatvec plus vec\n\n"); */
  /* d4est_util_print_matrix(A, m, n, "A = ", 0); */
  /* d4est_util_print_matrix(v, n, 1, "v = ", 0); */

  dgemv_ (&ytran, &dimN, &dimM, &alpha, A, &dimN, v, &incv, &beta, b, &incb);
}

void
d4est_linalg_mat_transpose (double * restrict A, double * restrict Atrans, int N)
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
d4est_linalg_mat_transpose_nonsqr (double * restrict A, double * restrict A_transpose, int A_rows,
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
 d4est_linalg_set_column (double * restrict A, double * restrict column, int col, int N, int M)
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
d4est_linalg_set_column_opt (double * restrict A, double * restrict column, int col, int N, int M)
{
  assert (col < M);
  int                 s = 0;

  for (int i = 0; i < N; i++){
    A[M*i + col] = column[s++];
  }
}



void
d4est_linalg_vec_scale (double alpha, double * restrict x, int N)
{
  int                 incx = 1;
  dscal_ (&N, &alpha, x, &incx);
}

void
d4est_linalg_vec_normalize(double * restrict x, int N)
{
  double xdotx = d4est_linalg_vec_dot(x,x,N);
  double norm = sqrt(xdotx);
  d4est_linalg_vec_scale(norm, x, N);
}

void d4est_linalg_vec_gen_random(double * restrict  vec, int N, long int seed, double a, double b){
  int i;
  for (i = 0; i < N; i++) {
    vec[i] = d4est_util_uniform_rand(seed,a,b);
  }
}

void
d4est_linalg_vec_gen_normalized_random(double * restrict  vec, int N, long int seed, double a, double b)
{
  d4est_linalg_vec_gen_random(vec, N, seed, a, b);
  d4est_linalg_vec_normalize(vec, N);
}

void
d4est_linalg_vec_axpy (double alpha, double * restrict x, double * restrict y, int N)
{
  int                 incx = 1;
  int                 incy = 1;
  daxpy_ (&N, &alpha, x, &incx, y, &incy);
}

void
d4est_linalg_vec_axpyeqz (double alpha, double * restrict x, double * restrict y, double * restrict z, int N)
{
  d4est_util_copy_1st_to_2nd (y, z, N);
  d4est_linalg_vec_axpy (alpha, x, z, N);
}

void
d4est_linalg_vec_xpby (double * restrict x, double beta, double * restrict y, int N)
{
  d4est_linalg_vec_scale (beta, y, N);
  d4est_linalg_vec_axpy (1., x, y, N);
}

double
d4est_linalg_vec_dot (double * restrict x, double * restrict y, int N)
{
  int                 incx = 1;
  int                 incy = 1;
  return ddot_ (&N, x, &incx, y, &incy);
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
d4est_linalg_sym_eigvals (double * restrict A, double * restrict eig_vals, int N)
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
d4est_linalg_sym_eigvecs (double * restrict A, double * restrict eig_vals, int N)
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

void
d4est_linalg_vec_fabs(double * restrict  x, int N) 
{
  int i;
  for (i = 0; i < N; i++) {
    x[i] = fabs(x[i]);
  }
}

void
d4est_linalg_vec_fabsdiff(double * restrict  x, double * restrict  y, double * restrict  result, int N)
{
  d4est_linalg_vec_axpyeqz(-1., x, y, result, N);
  d4est_linalg_vec_fabs(result, N);
}

void
d4est_linalg_component_mult
(
 double * restrict  x,
 double * restrict  y,
 double * restrict  xy,
 int N
)
{
  int i;
  for (i = 0; i < N; i++) {
    xy[i] = x[i]*y[i];
  }
}


void
d4est_linalg_component_div
(
 double * restrict  x,
 double * restrict  y,
 double * restrict  xdivy,
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
d4est_linalg_cross_prod
(
 double ax,
 double ay,
 double az,
 double bx,
 double by,
 double bz,
 double * restrict  axb_x,
 double * restrict  axb_y,
 double * restrict  axb_z
)
{
  *axb_x = ay*bz - az*by;
  *axb_y = az*bx - ax*bz;
  *axb_z = ax*by - ay*bx;
}

double
d4est_linalg_vec_sum
(
 double * restrict  vec,
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
d4est_linalg_vec1_trans_mat_vec2
(
 double * restrict  vec1,
 double * restrict  mat,
 double * restrict  vec2,
 int N
)
{

  double *mat_vec2 = (double *) malloc (sizeof (double) * N);
  d4est_linalg_matvec_plus_vec(1., mat, vec2, 0., mat_vec2, N, N);
  double dot = d4est_linalg_vec_dot(vec1, mat_vec2, N);
  free(mat_vec2);
  return dot;
}
