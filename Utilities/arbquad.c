#include <arbquad.h>

void
arbquad_print_float_info()
{
  printf(" *******   SINGLE PRECISION INFO    ******* \n");
  printf("FLT_MIN      = %e\n", FLT_MIN);
  printf("FLT_MAX      = %e\n", FLT_MAX);
  printf("FLT_EPSILON  = %e\n", FLT_EPSILON);

  printf(" *******   DOUBLE PRECISION INFO    ******* \n");
  printf("DBL_MIN      = %e\n", DBL_MIN);
  printf("DBL_MAX      = %e\n", DBL_MAX);
  printf("DBL_EPSILON  = %e\n", DBL_EPSILON);

  printf(" ******* LONG DOUBLE PRECISION INFO ******* \n");
  printf("LDBL_MIN      = %Le\n", LDBL_MIN);
  printf("LDBL_MAX      = %Le\n", LDBL_MAX);
  printf("LDBL_EPSILON  = %Le\n", LDBL_EPSILON);    
}

void
arbquad_test_moments
(
 int n,
 long double* weights,
 long double* abscissas,
 arbquad_moment_fcn_t moment_fcn,
 arbquad_weight_choice_t weight_choice,
 arbquad_weight_fcn_t weight_fcn,
 void* user
)
{
  int problems = 0;
  long double moment_maxerror = 0.l;
  
  for (int i = 0; i < 2*n; i++){
    long double moment_analytic = moment_fcn(i, user);
    long double moment_numerical = 0.;
    long double moment_relerror;

    if(weight_choice == DO_NOT_DIVIDE_WEIGHTS_BY_WEIGHT_FCN){
      for (int j = 0; j < n; j++){
        moment_numerical += weights[j]*powl(abscissas[j], i);
      }
    }
    else if (weight_choice == DIVIDE_WEIGHTS_BY_WEIGHT_FCN){
      for (int j = 0; j < n; j++){
        moment_numerical += weights[j]*powl(abscissas[j], i)*weight_fcn(abscissas[j],user);
      }
    }
    else {
      D4EST_ABORT("[D4EST_ERROR]: weight choice not supported\n");
    }
    if (moment_analytic == 0.0l)
      moment_relerror = fabsl((moment_analytic - moment_numerical));
    else
      moment_relerror = fabsl((moment_analytic - moment_numerical)/moment_analytic);

    moment_maxerror = (moment_relerror > moment_maxerror) ? moment_relerror : moment_maxerror;
    
    printf(" analytic, numerical, rel error, dbl eps, error < eps = %.25Lf %.25Lf %Le %Le %d\n",
           moment_numerical,
           moment_analytic,
           moment_relerror,
           (long double)(DBL_EPSILON),
           (moment_relerror < (long double)(DBL_EPSILON))
          );


    if (moment_relerror > (long double)(DBL_EPSILON))
      problems++;
  }

  if (problems == 0){
    printf("\n");
    printf("[RESULT]: The weights and abscissas ARE at least DBL_PRECISION accurate \n");
    printf("[RESULT]: MAX_ERROR = %Le\n", moment_maxerror);
    printf("[RESULT]: DBL_PRECISION EPSILON = %Le\n", (long double)(DBL_EPSILON));
           
  }
  else {
    printf("[RESULT]: The weights and abscissas ARE NOT at least DBL_PRECISION accurate \n");
    printf("[RESULT]: MAX_ERROR = %Le\n", moment_maxerror);
    printf("[RESULT]: DBL_PRECISION EPSILON = %Le\n", (long double)(DBL_EPSILON));
  }
  
}

long double
arbquad_signl(const long double a, const long double b)
{
  return b >= 0.0l ? (a >= 0.0l ? a : -a) : (a >= 0.0l ? -a : a);
}

long double
arbquad_pythagl(const long double a, const long double b) {

  long double absa=fabsl(a);
  long double absb=fabsl(b);

  return (absa > absb ? absa*sqrtl(1.0l+(absb/absa)*(absb/absa)) :
          (absb == 0.0l ? 0.0l : absb*sqrtl(1.0l+(absa/absb)*(absa/absb))));
}

int
arbquad_ql_implicit
(
 int n,
 long double* d,
 long double* e,
 long double* z
)
{
  int m,l,iter,i,k;
  long double s,r,p,g,f,dd,c,b;
  const long double EPS = LDBL_EPSILON;
  for (i=1;i<n;i++) {
    e[i-1]=e[i];
  }
  e[n-1]=0.0l;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabsl(d[m])+fabsl(d[m+1]);
        if (fabsl(e[m]) <= EPS*dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("[ARBQUAD_ERROR]: Too many iterations in tqli, stopping...\n");
          return 0;
        }
        g=(d[l+1]-d[l])/(2.0l*e[l]);
        r=arbquad_pythagl(g,1.0l);
        g=d[m]-d[l]+e[l]/(g+arbquad_signl(r,g));
        s=c=1.0l;
        p=0.0l;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=arbquad_pythagl(f,g));
          if (r == 0.0l) {
            d[i+1] -= p;
            e[m]=0.0l;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0l*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          for (k=0;k<n;k++) {
            f=z[k*n+i+1];
            z[k*n+i+1]=s*z[k*n+i]+c*f;
            z[k*n+i]=c*z[k*n+i]-s*f;
          }
        }
        if (r == 0.0l && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0l;
      }
    } while (m != l);
  }
}


/** 
 * Obtain the coefficients of the symmetric
 * tridiagonal Jacobi matrix
 * The abscissas are the eigenvalues and the
 * weights are related to the eigenvectors
 * 
 * @param n 
 * @param moments 
 * @param aa 
 * @param bb 
 */
void
arbquad_get_jacobi_coefficients
(
 int n,
 long double* moments,
 long double* aa,
 long double* bb
)
{
  int max_size = n + 2;
  long double* p0 = (long double*)calloc(max_size, sizeof(long double));
  long double* p1 = (long double*)calloc(max_size, sizeof(long double));
  long double* p2 = (long double*)calloc(max_size, sizeof(long double));
  
  int i = 0;
  p0[0] = 1.0l;
  i++;
  long double pp0 = moments[0];
  long double xpp = moments[1];
  aa[0] = xpp/pp0;
  bb[0] = 0.l;
  p1[0] = -xpp/pp0;
  p1[1] = 1.l;

  for (i = 2; i <= n; i++){

    /* The dot product <p|p> */
    long double pp1 = 0.l;
    for (int j = 0; j <= i - 1; j++){
      for (int k = 0; k <= i - 1; k++){
        pp1 += p1[j]*p1[k]*moments[j+k];
      }
    }
    
    /* The dot product <xp|p> */
    xpp = 0.l;
    for (int j = 0; j <= i - 1; j++){
      for (int k = 0; k <= i - 1; k++){
        xpp += p1[j]*p1[k]*moments[j + k + 1];
      }
    }    

    /* compute Jacobi Coefficients */
    aa[i - 1] = xpp/pp1;
    bb[i - 1] = pp1/pp0;

    /* compute next polynomial using recurrence */
    for(int i = 0; i < max_size; i++){
      p2[i] = (-xpp/pp1)*p1[i] - p0[i]*(pp1/pp0);
      if (i != 0){
        p2[i] += p1[i-1];
      }
    }
    
    /* flip the pointers for the next iteration */
    long double* tmp = p0;
    p0 = p1;
    p1 = p2;
    p2 = tmp;

    /* pp0 = old <p|p> dot product */
    pp0 = pp1;
  }
  
  free(p0);
  free(p1);
  free(p2);
}

/** 
 * Main function of the library
 * Get the abscissas and weights
 * for order n gaussian quadrature
 * with the weight fcn defined through
 * the moment_fcn
 * 
 * @param n 
 * @param weights 
 * @param abscissas 
 * @param moment_fcn 
 * @param user 
 */
void
arbquad_get_abscissas_and_weights
(
 int n,
 long double* weights,
 long double* abscissas,
 arbquad_moment_fcn_t moment_fcn,
 void* user,
 arbquad_weight_choice_t weight_choice,
 arbquad_weight_fcn_t weight_fcn
)
{

  int num_moments = 2*n;
  long double* moments = (long double*)calloc(num_moments, sizeof(long double));
  long double* aa = (long double*)calloc(n, sizeof(long double));
  long double* bb = (long double*)calloc(n, sizeof(long double));
  long double* eig_vals = (long double*)malloc(sizeof(long double)*n);
  long double* eig_vecs = (long double*)malloc(sizeof(long double)*n*n);
  long double* jacobi_matrix = (long double*)malloc(sizeof(long double)*n*n);
  
  for (int j = 0; j < 2*n; j++){
    int i = 2*n - 1 - j;
    moments[i] = moment_fcn(i, user);
    /* printf("moments[%d] = %.25Lf\n", i, moments[i]); */
  }

  arbquad_get_jacobi_coefficients(n, moments, aa, bb);
  
  long double* z = (long double*)malloc(n*n*sizeof(long double));
  long double* e = (long double*)malloc(n*sizeof(long double));
  long double* d = (long double*)malloc(sizeof(long double)*n);


  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      z[i*n + j] = 0.l;

  for (int i = 0; i < n; i++){
    /* printf("aa[%d], bb[%d] = %.25Lf %.25Lf\n", i,i,aa[i], bb[i]); */
    d[i] = aa[i];
    if (i != 0)
      e[i] = sqrtl(bb[i]);
    z[i*n + i] = 1.0l;
  }

  arbquad_ql_implicit(n, d,e,z);

  if (weight_choice == DO_NOT_DIVIDE_WEIGHTS_BY_WEIGHT_FCN){
    for (int k = 0; k < n; k++){
      long double vv = 0.l;
      for (int i = 0; i < n; i++){
        long double vi = z[i*n + k];
        vv += vi*vi;
      }
      abscissas[k] = d[k];
      weights[k] = moments[0]*z[k]*z[k]/vv;
    }
  }
  else if (weight_choice == DIVIDE_WEIGHTS_BY_WEIGHT_FCN) {
    D4EST_ASSERT(weight_fcn != NULL);
    for (int k = 0; k < n; k++){
      long double vv = 0.l;
      for (int i = 0; i < n; i++){
        long double vi = z[i*n + k];
        vv += vi*vi;
      }
      abscissas[k] = d[k];
      weights[k] = moments[0]*z[k]*z[k]/(vv*weight_fcn(abscissas[k],user));
    }    

  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: This weight choice is not supported\n");
  }
  
  free(z);
  free(e);
  free(d);
  free(jacobi_matrix);
  free(aa);
  free(bb);
  free(moments);
  free(eig_vals);
  free(eig_vecs);
}


void
arbquad_get_abscissas_and_weights_use_aa_and_bb
(
 int n,
 long double* weights,
 long double* abscissas,
 arbquad_moment_fcn_t moment_fcn,
 arbquad_aa_and_bb_fcn_t aa_and_bb_fcn,
 void* user,
 arbquad_weight_choice_t weight_choice,
 arbquad_weight_fcn_t weight_fcn
)
{

  int num_moments = 2*n;
  long double* moments = (long double*)calloc(num_moments, sizeof(long double));
  long double* aa = (long double*)calloc(n, sizeof(long double));
  long double* bb = (long double*)calloc(n, sizeof(long double));
  long double* eig_vals = (long double*)malloc(sizeof(long double)*n);
  long double* eig_vecs = (long double*)malloc(sizeof(long double)*n*n);
  long double* jacobi_matrix = (long double*)malloc(sizeof(long double)*n*n);
  
  for (int j = 0; j < 2*n; j++){
    int i = 2*n - 1 - j;
    moments[i] = moment_fcn(i, user);
    /* printf("moments[%d] = %.25Lf\n", i, moments[i]); */
  }

  aa_and_bb_fcn(n, aa, bb, user);
  
  long double* z = (long double*)malloc(n*n*sizeof(long double));
  long double* e = (long double*)malloc(n*sizeof(long double));
  long double* d = (long double*)malloc(sizeof(long double)*n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      z[i*n + j] = 0.l;

  for (int i = 0; i < n; i++){
    /* printf("aa[%d], bb[%d] = %.25Lf %.25Lf\n", i,i,aa[i], bb[i]); */
    d[i] = aa[i];
    if (i != 0)
      e[i] = sqrtl(bb[i]);
    z[i*n + i] = 1.0l;
  }

  arbquad_ql_implicit(n, d,e,z);

  if (weight_choice == DO_NOT_DIVIDE_WEIGHTS_BY_WEIGHT_FCN){
    for (int k = 0; k < n; k++){
      long double vv = 0.l;
      for (int i = 0; i < n; i++){
        long double vi = z[i*n + k];
        vv += vi*vi;
      }
      abscissas[k] = d[k];
      weights[k] = moments[0]*z[k]*z[k]/vv;
    }
  }
  else if (weight_choice == DIVIDE_WEIGHTS_BY_WEIGHT_FCN) {
    D4EST_ASSERT(weight_fcn != NULL);
    for (int k = 0; k < n; k++){
      long double vv = 0.l;
      for (int i = 0; i < n; i++){
        long double vi = z[i*n + k];
        vv += vi*vi;
      }
      abscissas[k] = d[k];
      weights[k] = moments[0]*z[k]*z[k]/(vv*weight_fcn(abscissas[k],user));
    }    

  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: This weight choice is not supported\n");
  }
  
  free(z);
  free(e);
  free(d);
  free(jacobi_matrix);
  free(aa);
  free(bb);
  free(moments);
  free(eig_vals);
  free(eig_vecs);

}
