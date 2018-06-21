#include "d4est_util.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <mpsort.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <sc_reduce.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


double
d4est_util_secant_fcn(double x){
  return 1./cos(x);
}


char* d4est_util_add_cwd(const char* dir)
{
  char* full_dir;
  char *cwd_path =(char *)malloc(10000*sizeof(char));
  getcwd(cwd_path,10000);
  asprintf(&full_dir, "%s/%s/", cwd_path, dir);
  free(cwd_path);
  return full_dir;
}

void d4est_util_make_directory(const char* dir, int add_cwd_to_dir)
{
  char* full_dir = NULL;

  if (add_cwd_to_dir){
    full_dir = d4est_util_add_cwd(dir);
  }
  else {
    asprintf(&full_dir, "%s", dir);
  }
  
  struct stat st = {0};
  if (stat(full_dir, &st) == -1) {
    mkdir(full_dir, 0700);
  }

  free(full_dir);
}

void
d4est_util_fill_array (double * restrict v, double val, int N)
{
  int                 i;
  for (i = 0; i < N; i++)
    v[i] = val;
}

void
d4est_util_zero_array (double * restrict v, int N)
{
  d4est_util_fill_array(v, 0., N);
}

int d4est_util_compact(double * restrict array, int size)
{
  int i;
  int last = 0;
  assert(size >= 0);
  if (size <= 0)
    return size;
  for (i = 1; i < size; i++)
    {
      if (array[i] != array[last])
        array[++last] = array[i];
    }
  return(last + 1);
}

int d4est_util_compact_1st_alongwith_2nd(double * restrict array, double * restrict  restrict  array2, int size)
{
  int i;
  int last = 0;
  assert(size >= 0);
  if (size <= 0)
    return size;
  for (i = 1; i < size; i++)
    {
      if (array[i] != array[last]){
        array[++last] = array[i];
        array2[last] = array2[i];
      }
    }
  return(last + 1);
}

double
d4est_util_min(double a, double b)
{
  return a < b ? a : b;
}

double
d4est_util_max(double a, double b)
{
  return a > b ? a : b; 
}

int
d4est_util_min_int(int a, int b)
{
  return a < b ? a : b;
}

int
d4est_util_max_int(int a, int b)
{
  return a > b ? a : b; 
}

void
d4est_util_eye(double * restrict  eye, int N){
  memset(eye, 0., sizeof(double)*N*N);
  int i;
  for (i = 0; i < N; i++) {
    eye[i*N + i] = 1.;
  }
}

double
d4est_util_dbl_pow_int(double a, int b){
   if (b == 0)
    {
        /* base case: anything to the 0 power is 1 */
        return 1.0;
    }
    else if (a == 0.0)
    {
        /* save us some time, 0 to any power other than 0 is 0 */
        return 0.0;
    }
    else if (b < 0)
    {
        /* b is negative, take the reciprocal of the positive version */
        return 1.0 / d4est_util_dbl_pow_int(a, -b);
    }
    else
    {
        /* b is positive, normal recursion */
        double result = d4est_util_dbl_pow_int(a, b / 2);
        result *= result;
        if (b % 2 != 0)
        {
            /* account for the truncation of b / 2 due to integer division */
            result *= a;
        }
        return result;
    }
}

int
d4est_util_int_pow_int(int base, int exp){
  int result = 1;
  while (exp)
    {
      if (exp & 1)
        result *= base;
      exp >>= 1;
      base *= base;
    }

  return result;
}

int
d4est_util_compare_double(double a, double b, double eps){
  if (fabs(a-b) < eps)
    return 1;
  return 0;
}

int
d4est_util_compare_vecs(double * restrict  a, double * restrict b, int N, double eps){
  int i;
  for (i = 0; i < N; i++) {
    if(!d4est_util_compare_double(a[i],b[i],eps)){
      /* printf("a[%d],b[%d],eps = %.25f,%.25f,%.25f\n", i,i, a[i], b[i], eps); */
      return 0;
    }
  }
  return 1;
}


void
d4est_util_find_biggest_error
(
 double * restrict  a,
 double * restrict  b,
 int N,
 double* biggest_err,
 int* biggest_id
){
  *biggest_err = -1.;
  *biggest_id = -1;
  for (int i = 0; i < N; i++) {
    double err = fabs(a[i] - b[i]);
    if (err > *biggest_err){
      *biggest_err = err;
      *biggest_id = i;
    }
  }
}


void 
d4est_util_print_3d_matrix(double * restrict  mat_3d, int n, int m, int l, char message [], int print_rank){
  int mpiret,mpirank,i,j,slice;
  if (print_rank){
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);
    printf("Printing from %d \n",mpirank);
  }
  printf("%s\n",message);

  for (slice = 0; slice < l; slice++)
    for (i = 0; i < n; i++){
      printf("\n");
      for (j = 0; j < m; j++)
	printf("%f ", mat_3d[slice*n*m + i*m + j]);
    }
  printf("\n");
}

void
d4est_util_print_matrix(double * restrict  mat, int n, int m, char message [], int print_rank){
  int mpiret,mpirank,i,j;
  if (print_rank){
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);
    printf("Printing from %d \n",mpirank);
  }
  printf("%s\n",message);
  for (i = 0; i < n; i++){
    printf("\n");
    for (j = 0; j < m; j++)
      printf("%.15f ",mat[i*m+j]);
  }
  printf("\n");
}

void
d4est_util_print_matrix_int(int* mat, int n, int m, char message [], int print_rank){
  int mpiret,mpirank,i,j;
  if (print_rank){
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);
    printf("Printing from %d \n",mpirank);
  }
  printf("%s\n",message);
  for (i = 0; i < n; i++){
    printf("\n");
    for (j = 0; j < m; j++)
      printf("%d ",mat[i*m+j]);
  }
  printf("\n");
}


void
d4est_util_print_matrices(double * restrict  mat, double * restrict  mat2, int n, int m, char message [])
{
  int i,j;
  printf("%s\n",message);
  for (i = 0; i < n; i++){
    printf("\n");
    for (j = 0; j < m; j++)
      printf("%.15f %.15f",mat[i*m+j], mat2[i*m+j]);
  }
  printf("\n");  
}

void
d4est_util_print_matrix_for_code(double * restrict  mat, int n, int m, char name []){
  int i;
  printf("\n\n/* BEGIN %s */\n", name); 
  printf("/* %s is %d by %d */\n\n", name, n, m);
  printf("double %s [%d] = {\n",name,n*m);
  for (i = 0; i < n*m; i++){
    printf("\n");
    if (i != n*m -1)
      printf("%.16f,",mat[i]);
    else{
      printf("%.16f\n",mat[i]);
      printf("};\n");
    }
  }
  printf("\n\n/* END %s */\n\n", name);
}


void
d4est_util_print_matrix_for_mathematica(double * restrict  mat, int n, int m, const char* name){
  printf("\n\n%s = {", name);
  for (int i = 0; i < n; i++){
    printf("{");
    for (int j = 0; j < m; j++){
      if (j != m-1)
        printf("%.16f, ",mat[i*m+j]);
      else
        printf("%.16f ",mat[i*m+j]);
    }
    if (i != n-1)
      printf("},");
    else
      printf("}");
  }
  printf("}\n\n");
}

double 
d4est_util_uniform_rand(long int seed, double a, double b){
  static long int seed_store = -1;
  /* static int i = 0; */

  /* if first time or seed is new */
  if (seed != seed_store){
    srand(seed);
    seed_store = seed;
  }
  
  double dev = (double) rand() / (double) RAND_MAX;
  return (b-a)*dev + a;
}

int
d4est_util_uniform_rand_int(long int seed, int a, int b){ 
  static long int seed_store = -1;
  /* static int i = 0; */

  /* if first time or seed is new */
  if (seed != seed_store){
    srand(seed);
    seed_store = seed;
  }

  return rand()%b + a;
}

double
d4est_util_max_error(double * restrict  u, double * restrict  u_sol, int N){
  int i;
  double temp;
  double e_max = -1.;
  for(i=0;i<N;i++){
    if (temp = fabs(u[i] - u_sol[i]),temp > e_max)
      e_max = temp;
  }
  return e_max;
}


int d4est_util_sum_array_int(int* array, int N)
{
  int i;
  int sum = 0;
  for (i = 0; i < N; i++) {
    sum += array[i];
  }
  return sum;
}


double d4est_util_sum_array_dbl(double * restrict array, int N)
{
  double sum = 0.;
  for (int i = 0; i < N; i++) {
    sum += array[i];
  }
  return sum;
}

#include <stdlib.h>

/* Comparison function. Receives two generic (void) pointers. */
static int
d4est_util_sort_double_callback(const void *p, const void *q)
{
    int ret;
    double x = *(const double *)p;
    double y = *(const double *)q;

    /* Avoid return x - y, which can cause undefined behaviour
       because of signed integer overflow. */
    if (x == y)
        ret = 0;
    else if (x < y)
        ret = -1;
    else
        ret = 1;

    return ret;
}

/* Sort an array of n integers, pointed to by a. */
void d4est_util_sort_double(double * restrict a, size_t n)
{
    qsort(a, n, sizeof(double), d4est_util_sort_double_callback);
}

typedef struct {

  int stride;
  double* a;
  double* b;

} d4est_util_double_pair_t;


double
d4est_util_parallel_checksum_dbl
(
 double * restrict  vec,
 int N
)
{
  int i = 0;
  double local_chksum = 0.;
  double global_chksum = 0.;
  for (i = 0; i < N; i++)
    local_chksum += vec[i];
  
  sc_reduce
    (
     &local_chksum,
     &global_chksum,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     0,
     sc_MPI_COMM_WORLD
    );
  
  return global_chksum;
}

double
d4est_util_parallel_reduce_dbl_scalar
(
 double local
)
{
  double global = 0.;
  sc_reduce
    (
     &local,
     &global,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     0,
     sc_MPI_COMM_WORLD
    );
  
  return global;
}

void
d4est_util_linear_regression
(
 double * restrict  y,
 double * restrict  x,
 double * restrict  m,
 double * restrict  b,
 int n
)
{
  double x1 = 0.;
  double y1 = 0.;
  double x1y1 = 0.;
  double x12 = 0.;
  int i;
  for(i=0;i<n;i++) { 
    x1+=x[i]; 
    y1+=y[i]; 
    x1y1+=(x[i]*y[i]); 
    x12+=(x[i]*x[i]); 
  } 

  (*m)=(n*x1y1-x1*y1)/(n*x12-x1*x1); 
  (*b)=y1/n-(*m)*x1/n; 
  /* e=y[0]-(*b)-(*m)*x[0];  */
  /* (*b)+=e;  */
}


/** 
 * Bisection
 * 
 * @param f 
 * @param a 
 * @param b 
 * @param m 
 */
int d4est_util_bisection
(
 double funk(double, void*),
 double x1,
 double x2,
 double xacc,
 int iter,
 double * restrict  root,
 void* user
)
{
  double dx,xmid;
  double f = funk(x1, user);
  double fmid = funk(x2, user);
  if(f*fmid > 0.0)
    return 1;
  
  (*root) = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);

  int j;
  for (j=1;j<=iter;j++) {

    /* printf(" root = %.25f\n ", *root); */
    fmid=funk(xmid=(*root)+(dx *= 0.5), user); 
    if (fmid <= 0.0) (*root)=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return 0;
  }
  
  return 2; 
}


double
d4est_util_compute_median(double * restrict x, int n){
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

double
d4est_util_normal_deviate (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


int d4est_util_does_file_exist(const char *filename) {
  struct stat st;
  int result = stat(filename, &st);
  return result == 0;
}

int d4est_util_match_couple
(
 const char* section1,
 const char* section2,
 const char* name1,
 const char* name2
)
{
  return (strcmp(section1, section2) == 0 && strcmp(name1, name2) == 0);
}

int d4est_util_match
(
 const char* str1,
 const char* str2
)
{
  return (strcmp(str1, str2) == 0);
}

void
d4est_util_compute_error_array(double * restrict  arr1, double * restrict  arr2, double * restrict  err, int N){
  for (int i = 0; i < N; i++) {
    err[i] = fabs(arr2[i] - arr1[i]);
  }
}

double
d4est_util_min_dbl_array(double * restrict  arr, int N){

  double min = arr[0];
  for (int i = 1; i < N; i++) {
    if (arr[i] < min){
      min = arr[i];
    }
  }
  return min;
}

double
d4est_util_mean_dbl_array(double * restrict  arr, int N){

  double mean = 0.;
  for (int i = 1; i < N; i++) {
      mean += arr[i];
  }
  return mean/(double)N;
}


double
d4est_util_max_dbl_array(double * restrict  arr, int N){

  double max = arr[0];
  for (int i = 1; i < N; i++) {
    if (arr[i] > max){
      max = arr[i];
    }
  }
  return max;
}

void
d4est_util_copy_1st_to_2nd (double * restrict v1, double * restrict v2, int N)
{
  memcpy (v2, v1, N * sizeof (double));
}

void d4est_util_gen_rand_vec(double * restrict  vec, int N, int seed, double a, double b){
  int i;
  for (i = 0; i < N; i++) {
    vec[i] = d4est_util_uniform_rand(seed,a,b);
  }
}

static void radix_double(const void * ptr, void * radix, void * arg) {
    *(double*)radix = *(const double*) ptr;
}

void d4est_util_parallel_sort
(
 sc_MPI_Comm mpicomm,
 double * restrict  array,
 int local_size
)
{
  mpsort_mpi(array, local_size, sizeof(double),
             radix_double, sizeof(double),
             NULL, mpicomm);
}
