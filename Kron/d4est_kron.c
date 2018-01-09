#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <d4est_linalg.h>
#include <d4est_kron.h>

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
d4est_kron_AoB (double *A, double *B, double *C, int a_rows, int a_cols,
                 int b_rows, int b_cols)
{
  int                 size2 = a_cols * b_cols;
  int                 i, j, k, p;
  for (i = 0; i < a_rows; i++)
    for (j = 0; j < a_cols; j++)
      for (k = 0; k < b_rows; k++)
        for (p = 0; p < b_cols; p++) {
          C[(k + i * b_rows) * size2 + p + j * b_cols] =
            A[i * a_cols + j] * B[k * b_cols + p];
        }
}


void
d4est_kron_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvec_dot_x)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_x[stride] = vec[i] * vec[k] * x[stride];
    }
}

void
d4est_kron_vec1_o_vec2_dot_x
(
 double *vec1,
 double* vec2,
 double*x,
 int vec1_size,
 int vec2_size,
 double* vec1vec2_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){
      int stride = k + i * vec2_size;
      vec1vec2_dot_x[stride] =
        vec1[i] * vec2[k] * x[stride];
    }
  }
}

double
d4est_kron_vec1_o_vec2_dot_x_sum
(
 double *vec1,
 double* vec2,
 double*x,
 int vec1_size,
 int vec2_size
)
{
  double sum = 0.;
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){
      int stride = k + i * vec2_size;
      sum +=
        vec1[i] * vec2[k] * x[stride];
    }
  }
  return sum;
}


double
d4est_kron_vec1_o_vec2_dot_ones_sum
(
 double *vec1,
 double* vec2,
 int vec1_size,
 int vec2_size
)
{
  double sum = 0.;
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){
      int stride = k + i * vec2_size;
      sum +=
        vec1[i] * vec2[k];
    }
  }
  return sum;
}


void
d4est_kron_vec1_o_vec2_dot_xy
(
 double *vec1,
 double* vec2,
 double*x,
 double*y,
 int vec1_size,
 int vec2_size,
 double* vec1vec2_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){
      int stride = k + i * vec2_size;
      vec1vec2_dot_x[stride] =
        vec1[i] * vec2[k] * x[stride] * y[stride];
    }
  }
}

void
d4est_kron_vec1_o_vec2_dot_wxyz
(
 double *vec1,
 double* vec2,
 double*w,
 double*x,
 double*y,
 double*z,
 int vec1_size,
 int vec2_size,
 double* vec1vec2_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){
      int stride = k + i * vec2_size;
      vec1vec2_dot_x[stride] =
        vec1[i] * vec2[k] * w[stride] * x[stride] * y[stride] * z[stride];
    }
  }
}

void
d4est_kron_vec1_o_vec2_o_vec3_dot_x
(
 double *vec1,
 double* vec2,
 double* vec3,
 double*x,
 int vec1_size,
 int vec2_size,
 int vec3_size,
 double* vec1vec2vec3_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){   
      for (int m = 0; m < vec3_size; m++){
        int stride = (m + (k + i * vec2_size) * vec3_size);
        vec1vec2vec3_dot_x[stride] =
          vec1[i] * vec2[k] * vec3[m] * x[stride];
      }
    }
  }  
}


double
d4est_kron_vec1_o_vec2_o_vec3_dot_x_sum
(
 double *vec1,
 double* vec2,
 double* vec3,
 double*x,
 int vec1_size,
 int vec2_size,
 int vec3_size
)
{
  double sum = 0;
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){   
      for (int m = 0; m < vec3_size; m++){
        int stride = (m + (k + i * vec2_size) * vec3_size);
        sum +=
          vec1[i] * vec2[k] * vec3[m] * x[stride];
      }
    }
  }
  return sum;
}


double
d4est_kron_vec1_o_vec2_o_vec3_dot_ones_sum
(
 double *vec1,
 double* vec2,
 double* vec3,
 int vec1_size,
 int vec2_size,
 int vec3_size
)
{
  double sum = 0;
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){   
      for (int m = 0; m < vec3_size; m++){
        sum +=
          vec1[i] * vec2[k] * vec3[m];
      }
    }
  }
  return sum;
}


void
d4est_kron_vec1_o_vec2_o_vec3_dot_xy
(
 double *vec1,
 double* vec2,
 double* vec3,
 double*x,
 double*y,
 int vec1_size,
 int vec2_size,
 int vec3_size,
 double* vec1vec2vec3_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){   
      for (int m = 0; m < vec3_size; m++){
        int stride = (m + (k + i * vec2_size) * vec3_size);
        vec1vec2vec3_dot_x[stride] =
          vec1[i] * vec2[k] * vec3[m] * x[stride] * y[stride];
      }
    }
  }  
}


void
d4est_kron_vec1_o_vec2_o_vec3_dot_wxyz
(
 double *vec1,
 double* vec2,
 double* vec3,
 double*w,
 double*x,
 double*y,
 double*z,
 int vec1_size,
 int vec2_size,
 int vec3_size,
 double* vec1vec2vec3_dot_x
)
{
  for (int i = 0; i < vec1_size; i++){
    for (int k = 0; k < vec2_size; k++){   
      for (int m = 0; m < vec3_size; m++){
        int stride = (m + (k + i * vec2_size) * vec3_size);
        vec1vec2vec3_dot_x[stride] =
          vec1[i] * vec2[k] * vec3[m] * w[stride] * x[stride] * y[stride] * z[stride];
      }
    }
  }  
}


void
d4est_kron_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_xy[stride] = vec[i] * vec[k] * x[stride] * y[stride];
    }
}

void
d4est_kron_vec_o_vec_dot_wxyz (double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vecvec_dot_wxyz)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_wxyz[stride] = vec[i] * vec[k] * w[stride] * x[stride] * y[stride] * z[stride];
    }
}


void
d4est_kron_oneover_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvec_dot_xy)
{
  for (int i = 0; i < vec_size; i++)
    for (int k = 0; k < vec_size; k++){
      int stride = k + i * vec_size;
      vecvec_dot_xy[stride] = (1./(vec[i] * vec[k] * x[stride])) * y[stride];
    }
}

void
d4est_kron_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_xy[i] = vec[i] * x[i] * y[i];
  }
}


void
d4est_kron_vec_dot_wxyz (double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vec_dot_wxyz)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_wxyz[i] = vec[i] * w[i] * x[i] * y[i] * z[i];
  }
}


void
d4est_kron_vec_dot_x (double *vec, double*x,int vec_size, double* vec_dot_x)
{
  for (int i = 0; i < vec_size; i++){
      vec_dot_x[i] = vec[i] * x[i];
  }
}

void
d4est_kron_oneover_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vec_dot_xy)
{
  for (int i = 0; i < vec_size; i++){
    vec_dot_xy[i] = (1./(vec[i] * x[i])) * y[i];
  }
}

void
d4est_kron_vec_o_vec_o_vec_dot_xy (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy)
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
d4est_kron_vec_o_vec_o_vec_dot_wxyz(double *vec, double* w, double*x, double* y, double* z, int vec_size, double* vecvecvec_dot_wxyz)
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
d4est_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y (double *vec, double*x, double* y, int vec_size, double* vecvecvec_dot_xy)
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
d4est_kron_vec_o_vec_o_vec_dot_x (double *vec, double*x, int vec_size, double* vecvecvec_dot_x)
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
d4est_kron_AoBoC (double *A, double *B, double *C, double *D, int a_rows,
                   int a_cols, int b_rows, int b_cols, int c_rows, int c_cols)
{
  double             *AoB =
    (double*) malloc (sizeof (double) * a_rows * a_cols * b_rows * b_cols);

  d4est_kron_AoB (A, B, AoB, a_rows, a_cols, b_rows, b_cols);
  d4est_kron_AoB (AoB, C, D, a_rows * b_rows, a_cols * b_cols, c_rows,
                   c_cols);

  free (AoB);
}


void
d4est_kron_A1A2x_nonsqr (double *A1A2x, double *A1, double *A2, double *X,
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

  d4est_linalg_mat_transpose_nonsqr (X, tmp1, a1_cols, a2_cols);
  d4est_linalg_mat_multiply (A2, tmp1, tmp, a2_rows, a2_cols, a1_cols);
  d4est_linalg_mat_transpose_nonsqr (tmp, tmp1, a2_rows, a1_cols);
  d4est_linalg_mat_multiply (A1, tmp1, A1A2x, a1_rows, a1_cols, a2_rows);

  free (tmp);
  free (tmp1);
}

void
d4est_kron_MAToIx_SQR (double *MAToIx, double *MAT, double *X, int N)
{
  d4est_linalg_mat_multiply (MAT, X, MAToIx, N, N, N);
}

void
d4est_kron_IoMATx_SQR (double *IoMATx, double *MAT, double *X, int N)
{
  double             *tmp = (double *) malloc (sizeof (double) * N * N);
  double             *tmp1 = (double *) malloc (sizeof (double) * N * N);

  d4est_linalg_mat_transpose_nonsqr (X, tmp1, N, N);
  d4est_linalg_mat_multiply (MAT, tmp1, tmp, N, N, N);
  d4est_linalg_mat_transpose_nonsqr (tmp, IoMATx, N, N);

  free (tmp);
  free (tmp1);
}

void
d4est_kron_VECoIx_SQR (double *VECoIx, double *VEC, double *X, int N)
{

  d4est_linalg_mat_multiply (VEC, X, VECoIx, N, 1, N);

}

void
d4est_kron_VEC_TRANSoIx_SQR (double *VEC_TRANSoIx, double *VEC, double *X,
                              int N)
{
  d4est_linalg_mat_multiply (VEC, X, VEC_TRANSoIx, 1, N, N);
}

void
d4est_kron_IoVECx_SQR (double *IoVECx, double *VEC, double *X, int N)
{
  /* int N = (a1_rows*a2_rows > a2_rows*a1_cols) ? a1_rows*a2_rows : a2_rows*a1_cols; */

  double             *tmp = (double *) malloc (sizeof (double) * N * N);
  /* double* tmp1 = (double*)malloc(sizeof(double)*N); */

  d4est_linalg_mat_multiply (VEC, X, tmp, N, 1, N);
  d4est_linalg_mat_transpose_nonsqr (tmp, IoVECx, N, N);

  free (tmp);
  /* free(tmp1);   */
}

void
d4est_kron_IoVEC_TRANSx_SQR (double *IoVEC_TRANSx, double *VEC, double *X,
                              int N)
{
  double             *tmp1 = (double *) malloc (sizeof (double) * N * N);

  d4est_linalg_mat_transpose_nonsqr (X, tmp1, N, N);
  d4est_linalg_mat_multiply (VEC, tmp1, IoVEC_TRANSx, 1, N, N);

  free (tmp1);
}

void
d4est_kron_A1A2A3x_nonsqr (double *A1A2A3x, double *A1, double *A2,
                            double *A3, double *X, int a1_rows, int a1_cols,
                            int a2_rows, int a2_cols, int a3_rows,
                            int a3_cols)
{
  double             *tmp =
    (double *) malloc (sizeof (double) * a2_rows * a3_rows * a1_cols);
  int                 i;
  for (i = 0; i < a1_cols; i++) {
    d4est_kron_A1A2x_nonsqr (&tmp[i * a2_rows * a3_rows], A2, A3,
                              &X[i * a2_cols * a3_cols], a2_rows, a2_cols,
                              a3_rows, a3_cols);
  }

  d4est_linalg_mat_multiply (A1, tmp, A1A2A3x, a1_rows, a1_cols, a2_rows * a3_rows);
  free (tmp);
}

void
d4est_kron_MAToIoIx_SQR (double *MAToIoIx, double *MAT, double *X, int N)
{
  d4est_linalg_mat_multiply (MAT, X, MAToIoIx, N, N, N * N);
}

void
d4est_kron_IoMAToIx_SQR (double *IoMAToIx, double *MAT, double *X, int N)
{
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_MAToIx_SQR (&IoMAToIx[i * N * N], MAT, &X[i * N * N], N);
  }
}

void
d4est_kron_IoIoMATx_SQR (double *IoIoMATx, double *MAT, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_IoMATx_SQR (&IoIoMATx[i * N * N], MAT, &X[i * N * N], N);
  }
}

void
d4est_kron_VECoIoIx_SQR (double *VECoIoIx, double *VEC, double *X, int N)
{
  d4est_linalg_mat_multiply (VEC, X, VECoIoIx, N, 1, N * N);
}

void
d4est_kron_VEC_TRANSoIoIx_SQR (double *VEC_TRANSoIoIx, double *VEC,
                                double *X, int N)
{
  d4est_linalg_mat_multiply (VEC, X, VEC_TRANSoIoIx, 1, N, N * N);
}

void
d4est_kron_IoVECoIx_SQR (double *IoVECoIx, double *VEC, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_VECoIx_SQR (&IoVECoIx[i * N * N], VEC, &X[i * N], N);
  }
}

void
d4est_kron_IoVEC_TRANSoIx_SQR (double *IoVEC_TRANSoIx, double *VEC,
                                double *X, int N)
{
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_VEC_TRANSoIx_SQR (&IoVEC_TRANSoIx[i * N], VEC, &X[i * N * N],
                                  N);
  }
}

void
d4est_kron_IoIoVECx_SQR (double *IoIoVECx, double *VEC, double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_IoVECx_SQR (&IoIoVECx[i * N * N], VEC, &X[i * N], N);
  }
}

void
d4est_kron_IoIoVEC_TRANSx_SQR (double *IoIoVEC_TRANSx, double *VEC,
                                double *X, int N)
{
  /* x should be N^2 by 1 */
  /* I should be N by N */
  /* VEC should be N by 1 */
  int                 i;
  for (i = 0; i < N; i++) {
    d4est_kron_IoVEC_TRANSx_SQR (&IoIoVEC_TRANSx[i * N], VEC, &X[i * N * N],
                                  N);
  }
}
