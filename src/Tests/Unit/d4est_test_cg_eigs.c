
#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <d4est_kron.h>
#include <d4est_util.h>

#define D4EST_TEST_ALLOC(t,n) (t *) malloc((n)*sizeof(t))

static void
fill_Tk
(
 double* Tk,
 double a0,
 double b0,
 double a1,
 double b1,
 int nodes,
 int i
)
{
  for (int j = 0; j < nodes; j++){
    Tk[i*nodes + j] = 0.;
  }
  
  if (i != 0 && i < nodes - 1) {
    Tk[i*nodes + i] = (1./a1 + b0/a0);
    Tk[i*nodes + i - 1] = -sqrt(b0)/a0;
    Tk[i*nodes + i + 1] = -sqrt(b1)/a1;
  }
  else if (i == 0) {
    Tk[i*nodes + i] = (1./a1);
    Tk[i*nodes + i + 1] = -sqrt(b1)/a1;
  }
  else {
    Tk[i*nodes + i] = (1./a1 + b0/a0);
    Tk[i*nodes + i - 1] = -sqrt(b0)/a0;
  }
}

static void
tridiag_gershgorin
(
 int i,
 int nodes,
 double a0,
 double b0,
 double a1,
 double b1,
 double* max,
 double* min
)
{
  double diag, offdiag_sum;
  if (i != 0 && i < nodes - 1) {
    diag = (1./a1 + b0/a0);
    offdiag_sum = fabs(sqrt(b1)/a1) + fabs(sqrt(b0)/a0);
  }
  else if (i == 0){
    diag = 1./a1;
    offdiag_sum = sqrt(b1)/a1;
  }
  else {
    diag = 1./a1 + b0/a0;
    offdiag_sum = fabs(sqrt(b0)/a0);
  }
  *max = diag + offdiag_sum;
  *min = diag - offdiag_sum;
}


void
apply_mat(double* A, double* x, double* Ax, int nodes){
  d4est_linalg_matvec_plus_vec (1.0, A, x, 0., Ax, nodes, nodes); 
}


double
cg_solve
(
 double* A,
 double* x0,
 double* b,
 int nodes,
 int iter,
 double* Tk
)
{

  double delta_new, delta_old, temp_max, temp_min, d_dot_Ad;
  double spectral_bound = -1;
  double alpha = -1.;
  double beta = -1.;
  double alpha_old = -1;
  double beta_old = -1;

  double* Ax = D4EST_TEST_ALLOC(double, nodes); 
  double* d = D4EST_TEST_ALLOC(double, nodes); 
  double* r = D4EST_TEST_ALLOC(double, nodes);

  apply_mat(A, x0, Ax, nodes);
  
  d4est_util_copy_1st_to_2nd(Ax, r, nodes);
  d4est_linalg_vec_xpby(b, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);
  
  for (int i = 0; i < iter; i++){

    apply_mat(A, d, Ax, nodes);
    d_dot_Ad = d4est_linalg_vec_dot(d, Ax, nodes);


    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, x0, nodes);
    d4est_linalg_vec_axpy(-alpha, Ax, r, nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, nodes);
    
    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);

    /* printf("d_dot_ad = %.15f\n", d_dot_Ad); */
    /* printf("beta = %.15f\n", beta); */
    /* printf("alpha = %.15f\n", alpha); */
    
    tridiag_gershgorin(
                       i,
                       nodes,
                       alpha_old,
                       beta_old,
                       alpha,
                       beta,
                       &temp_max,
                       &temp_min
                      );

    if (i > 0){
      spectral_bound = d4est_util_max( spectral_bound, temp_max );  
    }
    else{
      spectral_bound = temp_max;
    }
    
    fill_Tk(Tk, alpha_old, beta_old, alpha, beta, nodes, i);

    if (i >= 1){

      printf("\n*******i = %d*******\n", i);
      /* d4est_util_print_matrix(Tk, nodes, nodes, "Tk_sub =", 0); */
      double* Tk_sub = D4EST_TEST_ALLOC(double,(i+1)*(i+1));
      for (int m = 0; m < (i+1); m++){
        for (int n = 0; n < (i+1); n++){
          Tk_sub[m*(i+1) + n] = Tk[m*nodes + n];
        }
      }

      d4est_util_print_matrix(Tk_sub, (i+1), (i+1), "Tk_sub =", 0);

      double* Tk_eig_vals = D4EST_TEST_ALLOC(double, (i+1));
      d4est_linalg_sym_eigvals(Tk_sub, Tk_eig_vals, (i+1));

      d4est_util_print_matrix(Tk_eig_vals, (i+1), 1, "Tk_eigs =", 0);

      
    /*   printf("\n k = %d \n", i); */
      /* d4est_util_print_matrix(Tk_sub, nodes, nodes, "Tk_sub =", 0); */
      /* DEBUG_PRINT_ARR_DBL(Tk_eig_vals, (i+1)); */

      free(Tk_sub);
      free(Tk_eig_vals);
    }    
  }
  
  free(Ax);
  free(d);
  free(r);
  return spectral_bound;
}

int main(int argc, char *argv[])
{
  int nodes = 5;
  double* A = D4EST_TEST_ALLOC(double,nodes*nodes);
  double* Tk = D4EST_TEST_ALLOC(double,nodes*nodes);
  double* x0 = D4EST_TEST_ALLOC(double, nodes);
  double* b = D4EST_TEST_ALLOC(double, nodes);
    
  for (int i = 0; i < nodes; i++){

    x0[i] = (i!=0) ? 0. : 1.;
    b[i] = 1.;
    for (int j = 0; j < nodes; j++){
      Tk[i*nodes + j] = 0.;
      A[i*nodes + j] = 0.;
    }
    A[i*nodes + i] = 2;
    if (i != nodes - 1)
      A[i*nodes + i + 1] = -1;
    if (i != 0)
      A[i*nodes + i - 1] = -1;
  }

  double* A_inv = D4EST_TEST_ALLOC(double, nodes*nodes);
  double* A_inv_times_b = D4EST_TEST_ALLOC(double, nodes);
  double* A_eig_vals = D4EST_TEST_ALLOC(double, nodes);
  double* Tk_eig_vals = D4EST_TEST_ALLOC(double, nodes);

  d4est_util_print_matrix(A, nodes, nodes, "A =", 0);
  DEBUG_PRINT_ARR_DBL(x0, nodes);
  DEBUG_PRINT_ARR_DBL(b,nodes);

  double spectral_bound = cg_solve(A,x0,b,nodes,nodes,Tk);

  printf("\n****** POST SOLVE ******\n");
  
  d4est_linalg_invert_and_copy(A,A_inv,nodes);
  apply_mat(A_inv, b, A_inv_times_b, nodes);
  printf("spectral_bound = %.15f\n", spectral_bound);
  d4est_linalg_sym_eigvals(A, A_eig_vals, nodes);
  d4est_linalg_sym_eigvals(Tk, Tk_eig_vals, nodes);
  DEBUG_PRINT_ARR_DBL(x0, nodes);

  DEBUG_PRINT_ARR_DBL(A_inv_times_b,nodes);
  DEBUG_PRINT_ARR_DBL(A_eig_vals,nodes);

  d4est_util_print_matrix(Tk, nodes, nodes, "Tk =", 0);  
  DEBUG_PRINT_ARR_DBL(Tk_eig_vals,nodes);

  

  
  free(A);
  free(A_inv);
  free(A_inv_times_b);
  free(Tk);
  free(x0);
  free(b);
  free(A_eig_vals);
  free(Tk_eig_vals);
  return 0;
}


