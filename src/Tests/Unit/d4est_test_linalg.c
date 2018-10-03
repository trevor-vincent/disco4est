#include <d4est_linalg.h>
#include <d4est_util.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double testd4est_linalg_point_error_max(double* u, double* u_analytic, int N){ 
  int i;
  double err_max, err_temp; 
  err_max = 0.;

  for (i = 0; i < N; i++){
    err_temp = fabs(u[i] - u_analytic[i]);
    if (err_temp > err_max)
      err_max = err_temp;
  }

  return err_max;
}

void testd4est_linalg_init_mat(double* A, int rows, int cols, double r){
  int i,j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      A[i*cols + j] = i*.23 + .5*j*j + (i*j) + 5 + r;
    } 
  }
}

void testd4est_linalg_init_vec(double* v, int rows, double r){
  int i;
  for (i=0; i<rows; i++)
    v[i] = i - 4 + cos(r*i);
}

/** 
 * Test Matrix Multiplication
 *
 * TODO: make dimensions different and random
 * TODO: initialize matrices randomly
 * 
 * @param testd4est_linalg_mat_multiply 
 * 
 * @return 
 */

void testd4est_linalg_transpose()
{
  const int A_rows = 3;
  const int A_cols = 4;

  double A [A_rows][A_cols];
  double A_transpose [A_cols][A_rows]; 
  double A_transpose_transpose [A_rows][A_cols];
  double err_max, err_temp;

  int i,j;

  testd4est_linalg_init_mat(&A[0][0], A_rows, A_cols, 4.323);

  d4est_linalg_mat_transpose_nonsqr(&A[0][0], &A_transpose[0][0], A_rows, A_cols);
  d4est_linalg_mat_transpose_nonsqr(&A_transpose[0][0], &A_transpose_transpose[0][0], A_cols, A_rows);
  err_max = 0.;

  for (i = 0; i < A_rows; i++)
    for (j = 0; j < A_cols; j++){
      err_temp = fabs(A[i][j] - A_transpose_transpose[i][j]);
      if (err_temp > err_max)
	err_max = err_temp;
    }

  if (err_max > .0000001){
    printf("testd4est_linalg_transpose failed\n");
    exit(1);
  }
}

void testd4est_linalg_mat_multiply()
{
  const int m = 3;
  const int l = 3;
  const int n = 3;

  double A [m][l];
  double B [l][n];
  double C [m][n];
  double C_ans [m][n];
  double err_max, err_temp;

  int i,j,k;

  testd4est_linalg_init_mat(&A[0][0], m, l, 4.323);
  testd4est_linalg_init_mat(&B[0][0], l, n, 53.23);
  
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++){
      C_ans[i][j] = 0.;
      for (k = 0; k < l; k++)
	C_ans[i][j] += A[i][k]*B[k][j];
    }

  d4est_linalg_mat_multiply(&A[0][0],&B[0][0],&C[0][0],m,n,l);
  err_max = 0.;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++){
      err_temp = fabs(C[i][j] - C_ans[i][j]);
      if (err_temp > err_max)
	err_max = err_temp;
    }
  if (err_max > .0000001){
    printf("testd4est_linalg_mat_multiply failed");
    exit(1);
  }
}



void testd4est_linalg_mat_vec()
{
  const int m = 10;
  const int l = 5;

  double A [m][l];
  double v [l];
  double Av [m];
  double Av_ans [m];
  double err_max, err_temp;

  int i,k;

  testd4est_linalg_init_mat(&A[0][0], m, l, 4.323);
  testd4est_linalg_init_vec(&v[0], l, 23.23);
  
  for (i = 0; i < m; i++){
      Av_ans[i] = 0.;
      for (k = 0; k < l; k++)
	Av_ans[i] += A[i][k]*v[k];
  }

  d4est_linalg_matvec_plus_vec(1.,&A[0][0],&v[0],0.0,&Av[0],m,l);
  err_max = 0.;

  for (i = 0; i < m; i++){
    err_temp = fabs(Av_ans[i] - Av[i]);
    /* printf("Ans = %f, Blas = %f\n", Av_ans[i], Av[i]); */
    if (err_temp > err_max)
      err_max = err_temp;
  }
  if (err_max > .0000001){
    printf("testd4est_linalg_mat_vec failed\n");
    exit(1);
  }
}


void testd4est_linalg_xpby()
{
  const int l = 20;
  int i;

  double x [l];
  double y [l];
  double b = -1.;
  double xpby_ans [l];

  double err_max, err_temp;

  testd4est_linalg_init_vec(&x[0], l, -22.22);
  testd4est_linalg_init_vec(&y[0], l, 21.23);
  
  for (i = 0; i < l; i++){
    xpby_ans[i] = x[i] + b*y[i];
  }

  d4est_linalg_vec_xpby(&x[0], b, &y[0], l);
  err_max = 0.;

  for (i = 0; i < l; i++){
    err_temp = fabs(xpby_ans[i] - y[i]);
    /* printf("Ans = %f, Blas = %f\n", xpby_ans[i], y[i]); */
    if (err_temp > err_max)
      err_max = err_temp;
  }
  if (err_max > .0000001){
    printf("testd4est_linalg_xpby failed\n");
    exit(1);
  }
}


void testd4est_linalg_sym_eigvals()
{
  int i;
  double A [10][10] =
    {
      {1., 11., 7., 9., 7., 11., 7., 9., 2., 11.} ,
      {11., 4., 10., 10., 6., 2., 9., 6., 10., 0.} ,
      {7., 10., 3., 5., 4., 4., 4., 4., 6., 10.} ,
      {9., 10., 5., 3., 8., 8., 3., 5., 1., 8.} ,
      {7., 6., 4., 8., 8., 10., 5., 6., 10., 0.} ,
      {11., 2., 4., 8., 10., 9., 4., 3., 5., 11.} ,
      {7., 9., 4., 3., 5., 4., 3., 10., 7., 2.} ,
      {9., 6., 4., 5., 6., 3., 10., 11., 1., 7.} ,
      {2., 10., 6., 1., 10., 5., 7., 1., 10., 5.} ,
      {11., 0., 10., 8., 0., 11., 2., 7., 5., 1.}
    };

  double eig_vals [10];
  double eig_vals_ans [10];

  eig_vals_ans[0] = -18.57327950993016;
  eig_vals_ans[1] = -10.83432690756231;
  eig_vals_ans[2] = -8.70591312578974;
  eig_vals_ans[3] = -5.75782552751678;
  eig_vals_ans[4] = -2.504218894968408;
  eig_vals_ans[5] = 3.484425996018023;
  eig_vals_ans[6] = 6.750775469094293;
  eig_vals_ans[7] = 12.02037231054202;
  eig_vals_ans[8] = 14.60708222029769;
  eig_vals_ans[9] = 62.51290796981536;

  d4est_linalg_sym_eigvals(&A[0][0], &eig_vals[0], 10);
  double err_max, err_temp;
  
  err_max = 0.;

  for (i = 0; i < 10; i++){
    err_temp = fabs(eig_vals_ans[i] - eig_vals[i]);
    /* printf("Ans = %f, LAPACK = %f\n", eig_vals_ans[i], eig_vals[i]); */
    if (err_temp > err_max)
      err_max = err_temp;
  }
  if (err_max > .0000001){
    printf("testd4est_linalg_sym_eigvals failed\n");
    exit(1);
  }
}


void testd4est_linalg_sym_eigvecs()
{
  int i;
  double A [10][10] =
    {
      {1., 11., 7., 9., 7., 11., 7., 9., 2., 11.} ,
      {11., 4., 10., 10., 6., 2., 9., 6., 10., 0.} ,
      {7., 10., 3., 5., 4., 4., 4., 4., 6., 10.} ,
      {9., 10., 5., 3., 8., 8., 3., 5., 1., 8.} ,
      {7., 6., 4., 8., 8., 10., 5., 6., 10., 0.} ,
      {11., 2., 4., 8., 10., 9., 4., 3., 5., 11.} ,
      {7., 9., 4., 3., 5., 4., 3., 10., 7., 2.} ,
      {9., 6., 4., 5., 6., 3., 10., 11., 1., 7.} ,
      {2., 10., 6., 1., 10., 5., 7., 1., 10., 5.} ,
      {11., 0., 10., 8., 0., 11., 2., 7., 5., 1.}
    };

  double A_copy [10][10];

  d4est_util_copy_1st_to_2nd(&A[0][0], &A_copy[0][0], 10*10);

  double eig_vals [10];

  d4est_linalg_sym_eigvecs(&A_copy[0][0], &eig_vals[0], 10);

  double evec [10];
  double err_max = 0.;

  for (i = 0; i < 10; i++){
    d4est_linalg_matvec_plus_vec(1.0, &A[0][0], &A_copy[i][0], 0.0, &evec[0], 10, 10);
    /* d4est_util_print_matrix(&evec[0], 10, 1, "evec before scale= ", 0); */
    d4est_linalg_vec_scale(1./eig_vals[i], &evec[0], 10);
    /* d4est_util_print_matrix(&A_copy[i][0], 10, 1, "A_copy = ", 0); */
    /* d4est_util_print_matrix(&evec[0], 10, 1, "evec = ", 0); */

    err_max += testd4est_linalg_point_error_max(&A_copy[i][0], &evec[0], 10);
  }

  if (err_max > .0000001){
    printf("testd4est_linalg_sym_eigvecs failed\n");
    exit(1);
  }
}

void testd4est_linalg_set_column()
{
  int i;
  double A [10][10] =
    {
      {1., 11., 7., 9., 7., 11., 7., 9., 2., 11.} ,
      {11., 4., 10., 10., 6., 2., 9., 6., 10., 0.} ,
      {7., 10., 3., 5., 4., 4., 4., 4., 6., 10.} ,
      {9., 10., 5., 3., 8., 8., 3., 5., 1., 8.} ,
      {7., 6., 4., 8., 8., 10., 5., 6., 10., 0.} ,
      {11., 2., 4., 8., 10., 9., 4., 3., 5., 11.} ,
      {7., 9., 4., 3., 5., 4., 3., 10., 7., 2.} ,
      {9., 6., 4., 5., 6., 3., 10., 11., 1., 7.} ,
      {2., 10., 6., 1., 10., 5., 7., 1., 10., 5.} ,
      {11., 0., 10., 8., 0., 11., 2., 7., 5., 1.}
    };

  double col [10];
  d4est_util_fill_array(col, 1., 10);
  d4est_linalg_set_column(&A[0][0], col, 3, 10, 10);

  double err_max, err_temp;
  err_max = 0.;

  for (i = 0; i < 10; i++){
    err_temp = fabs(A[i][3] - col[i]);
    if (err_temp > err_max)
      err_max = err_temp;
  }

  if (err_max > .0000001){
    printf("testd4est_linalg_set_column failed\n");
    exit(1);
  }
}

int main(int argc, char *argv[])
{
  testd4est_linalg_set_column();
  testd4est_linalg_sym_eigvecs();
  testd4est_linalg_sym_eigvals();
  testd4est_linalg_xpby();
  testd4est_linalg_mat_vec();
  testd4est_linalg_mat_multiply();
  testd4est_linalg_transpose();
  return 0;
}
