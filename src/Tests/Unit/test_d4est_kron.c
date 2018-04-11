#include <d4est_kron.h>
#include <d4est_linalg.h>
#include <d4est_util.h>
#include <time.h>
#include <stdlib.h>

static void
gen_rand_vec(double* pVec, int N, double a, double b){
  int i;
  for (i = 0; i < N; i++) {
    double dev = (double) rand() / (double) RAND_MAX;
    pVec[i] = (b-a)*dev + a;
  }
}

void
test_d4est_kron_AoBx() 
{
  long int seed = time(NULL);
  srand(seed);
  
  int nRowsA = rand()%10 + 1;
  int nColsA = rand()%10 + 1;

  int nRowsB = rand()%10 + 1;
  int nColsB = rand()%10 + 1;
  
  int nRowsC = nRowsA*nRowsB;
  int nColsC = nColsA*nColsB;

  double* pA = (double*)malloc(sizeof(double)*nRowsA*nColsA);
  double* pB = (double*)malloc(sizeof(double)*nRowsB*nColsB);
  double* pC = (double*)malloc(sizeof(double)*nRowsC*nColsC);
  double* pX = (double*)malloc(sizeof(double)*nColsC);
  double* pY1 = (double*)malloc(sizeof(double)*nRowsC);
  double* pY2 = (double*)malloc(sizeof(double)*nRowsC);
    
  gen_rand_vec(pX, nColsC, 0., 1.);
  gen_rand_vec(pA, nRowsA*nColsA, 0., 1.);
  gen_rand_vec(pB, nRowsB*nColsB, 0., 1.);
  
  d4est_kron_AoB(pA, pB, pC, nRowsA, nColsA, nRowsB, nColsB);
  d4est_linalg_matvec_plus_vec(1.0, pC, pX, 0., pY1, nRowsC, nColsC);
  d4est_kron_A1A2x_nonsqr(pY2, pA, pB, pX, nRowsA, nColsA, nRowsB, nColsB);

  /* int test_check = ; */
  if(!d4est_util_compare_vecs(pY1, pY2, nRowsC, .000001)){
    printf("nRowsA = %d\n", nRowsA); 
    printf("nColsA = %d\n", nColsA); 
    printf("nRowsB = %d\n", nRowsB); 
    printf("nColsB = %d\n", nColsB);
    d4est_util_print_matrix(pC, nRowsC, nColsC, "C = ", 0);
    d4est_util_print_matrix(pX, nColsC, 1, "X = ", 0);
    d4est_util_print_matrices(pY1, pY2, nRowsC, 1, "Y1, Y2 = ");
    free(pA); 
   free(pB);
    free(pC);
    free(pX);
    free(pY1);
    free(pY2);
    D4EST_ABORT("[test_d4est_kron_AoBx failed]");
  }
  
  free(pA);
  free(pB);
  free(pC);
  free(pX);
  free(pY1);
  free(pY2);
}

static void
test_d4est_kron_AoBoCx()
{
  long int seed = time(NULL);
  srand(seed);
  
  int nRowsA = rand()%10 + 1;
  int nColsA = rand()%10 + 1;
  int nRowsB = rand()%10 + 1;
  int nColsB = rand()%10 + 1;
  int nRowsC = rand()%10 + 1;
  int nColsC = rand()%10 + 1;
  
  int nRowsD = nRowsA*nRowsB*nRowsC;
  int nColsD = nColsA*nColsB*nColsC;

  double* pA = (double*)malloc(sizeof(double)*nRowsA*nColsA);
  double* pB = (double*)malloc(sizeof(double)*nRowsB*nColsB);
  double* pAoB = (double*)malloc(sizeof(double)*nRowsA*nColsA*nRowsB*nColsB);
  
  double* pC = (double*)malloc(sizeof(double)*nRowsC*nColsC);
  double* pD = (double*)malloc(sizeof(double)*nRowsD*nColsD);
  
  double* pX = (double*)malloc(sizeof(double)*nColsD);
  double* pY1 = (double*)malloc(sizeof(double)*nRowsD);
  double* pY2 = (double*)malloc(sizeof(double)*nRowsD);
  
  gen_rand_vec(pX, nColsD, 0., 1.);
  gen_rand_vec(pA, nRowsA*nColsA, 0., 1.);
  gen_rand_vec(pB, nRowsB*nColsB, 0., 1.);
  gen_rand_vec(pC, nRowsC*nColsC, 0., 1.);
  
  d4est_kron_AoB(pA, pB, pAoB, nRowsA, nColsA, nRowsB, nColsB);
  d4est_kron_AoB(pAoB, pC, pD, nRowsA*nRowsB, nColsA*nColsB, nRowsC, nColsC);
  d4est_linalg_matvec_plus_vec(1.0, pD, pX, 0., pY1, nRowsD, nColsD);
  d4est_kron_A1A2A3x_nonsqr(pY2, pA, pB, pC, pX, nRowsA, nColsA, nRowsB, nColsB, nRowsC, nColsC);

  /* int test_check = ; */
  if(!d4est_util_compare_vecs(pY1, pY2, nRowsD, .000001)){
    printf("nRowsA = %d\n", nRowsA); 
    printf("nColsA = %d\n", nColsA); 
    printf("nRowsB = %d\n", nRowsB); 
    printf("nColsB = %d\n", nColsB); 
    printf("nRowsC = %d\n", nRowsC);   
    printf("nColsC = %d\n", nColsC);  
    free(pA);
    free(pAoB);
    free(pB);
    free(pC);
    free(pD);
    free(pX);
    free(pY1);
    free(pY2);
    D4EST_ABORT("[test_d4est_kron_AoBoCx failed]");
  }
    
  free(pA);
  free(pAoB);
  free(pB);
  free(pC);
  free(pD);
  free(pX);
  free(pY1);
  free(pY2);
}


int
main (int argc, char *argv[])
{
  test_d4est_kron_AoBx();
  test_d4est_kron_AoBoCx();
  return 0;
}
