#ifndef UTIL_H
#define UTIL_H

#define _GNU_SOURCE
#include <stdio.h>


#ifndef NDEBUG

#define DEBUG_PRINT_INT(x) printf("%s = %d\n", #x, x)  
#define DEBUG_PRINT_DBL(x) printf("%s = %.16f\n",#x, x)

#define DEBUG_PRINT_2DBL(a, b) do {            \
    printf("%s, %s = \n",#a, #b);               \
    printf("%.16f %.16f\n",a, b);         \
  } while(0)

#define DEBUG_PRINT_4DBL(a, b, c, d) do {                              \
    printf("%s, %s, %s, %s= \n",#a, #b, #c, #d);                        \
      printf("%.16f %.16f %.16f %.16f\n",a, b, c, d);       \
  } while(0)

#define DEBUG_PRINT_8DBL(a, b, c, d, e, f, g, h) do {             \
    printf("%s, %s, %s, %s %s %s %s %s= \n",#a, #b, #c, #d, #e, #f, #g, #h);    \
      printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n",a, b, c, d, e, f, g, h); \
  } while(0)


#define DEBUG_PRINT_ARR_DBL(a, n) do {                \
    printf("%s = \n",#a);                       \
    for (int i = 0; i < n; i++) {               \
      printf("%.16f\n",a[i]);                   \
    }                                           \
  } while(0)

#define DEBUG_PRINT_ARR_DBL_SUM(a, n) do {                \
    double sum = 0.;                                  \
    for (int i = 0; i < n; i++) {                     \
      sum += a[i];                                    \
    }                                                 \
    printf("%s = %.25f\n",#a, sum);                   \
  } while(0)


#define DEBUG_PRINT_ARR_INT(a, n) do {                \
    printf("%s = \n",#a);                       \
    for (int i = 0; i < n; i++) {               \
      printf("%d\n",a[i]);                   \
    }                                           \
  } while(0)


#define DEBUG_PRINT_MAT_DBL(a, m, n) do {             \
    printf("%s = ",#a);                         \
    for (int i = 0; i < m; i++) {               \
      printf("\n");                             \
      for (int j = 0; j < n; j++)               \
        printf("%.16f",a[i*m + j]);             \
    }                                           \
  } while(0)

#define DEBUG_PRINT_MAT_INT(a, m, n) do {             \
    printf("%s = ",#a);                         \
    for (int i = 0; i < m; i++) {               \
      printf("\n");                             \
      for (int j = 0; j < n; j++)               \
        printf("%d",a[i*m + j]);                \
    }                                           \
    printf("\n");                               \
  } while(0)

#define DEBUG_PRINT_2ARR_DBL(a, b, n) do {            \
    printf("%s, %s = \n",#a, #b);               \
    for (int i = 0; i < n; i++) {               \
      printf("%.25f %.25f\n",a[i], b[i]);       \
    }                                           \
  } while(0)

#define DEBUG_PRINT_3ARR_DBL(a, b, c, n) do {                 \
    printf("%s, %s, %s = \n",#a, #b, #c);               \
    for (int i = 0; i < n; i++) {                       \
      printf("%.16f %.16f %.16f\n",a[i], b[i], c[i]);   \
    }                                                   \
  } while(0)

#define DEBUG_PRINT_4ARR_DBL(a, b, c, d, n) do {                              \
    printf("%s, %s, %s, %s= \n",#a, #b, #c, #d);                        \
    for (int i = 0; i < n; i++) {                                       \
      printf("%.16f %.16f %.16f %.16f\n",a[i], b[i], c[i], d[i]);       \
    }                                                                   \
  } while(0)


#define DEBUG_PRINT_6ARR_DBL(a, b, c, d, e, f, n) do {                        \
    printf("%s, %s, %s, %s %s %s = \n",#a, #b, #c, #d, #e, #f);                \
    for (int i = 0; i < n; i++) {                                       \
      printf("%.16f %.16f %.16f %.16f %.16f %.16f\n",a[i], b[i], c[i], d[i], e[i], f[i]); \
    }                                                                   \
  } while(0)

#define DEBUG_PRINT_8ARR_DBL(a, b, c, d, e, f, g, h, n) do {             \
    printf("%s, %s, %s, %s %s %s %s %s= \n",#a, #b, #c, #d, #e, #f, #g, #h);    \
    for (int i = 0; i < n; i++) {                                       \
      printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n",a[i], b[i], c[i], d[i], e[i], f[i], g[i], h[i]); \
    }                                                                   \
  } while(0)



#define DEBUG_PRINT_9ARR_DBL(a, b, c, d, e, f, g, h, j, n) do {          \
    printf("%s, %s, %s, %s %s %s %s %s %s= \n",#a, #b, #c, #d, #e, #f, #g, #h, #j); \
    for (int i = 0; i < n; i++) {                                       \
      printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n",a[i], b[i], c[i], d[i], e[i], f[i], g[i], h[i], j[i]); \
    }                                                                   \
  } while(0)


#endif



#define D4EST_NOOP()                                                            \
  do                                                                           \
  {                                                                            \
  } while (0)



#ifndef NDEBUG
#define mpi_assert(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define mpi_assert(c) D4EST_NOOP()
#endif

/* #define mpi_assert(expression)  \ */
/* do { \ */
/*   if (!(expression)) { \ */
/*      fprintf(stderr, "Failed assertion at %d in %s",__LINE__, __FILE__); \ */
/*      MPI_Abort(MPI_COMM_WORLD, 1); \ */
/*  } \ */
/* } while (0) */

#define mpi_abort(c) SC_ABORT(c)

#define mpi_abort_if(q, ...) ((q) ? mpi_abort(__VA_ARGS__) : (void)0)
#define mpi_abort_if_not(q, ...) mpi_abort_if(!(q), __VA_ARGS__)


/* #define mpi_abort(expression) \ */
/* do { \ */
/*     fprintf(stderr, "%s Failed assertion at %d in %s",expression,__LINE__, __FILE__); \ */
/*      MPI_Abort(MPI_COMM_WORLD, 1); \ */
/* } while (0) */


#define D4EST_ASPRINTF(write_to, ...) {              \
    char *tmp_string_for_extend = (write_to);        \
    if(asprintf(&(write_to), __VA_ARGS__) <0){       \
      mpi_abort("ASPRINTF_ERROR");                   \
    }                                                \
    free(tmp_string_for_extend);                     \
  }

#define D4EST_FREE_MAT(a, n1, n2) do {                                  \
    for (int i = 0; i < n1; i++) {                                      \
      for (int j = 0; j < n2; j++) {                                    \
        P4EST_FREE(a[i][j]);                                            \
      }                                                                 \
    }                                                                   \
  } while(0)

#define D4EST_FREE_VEC(a, n1) do {                                      \
    for (int i = 0; i < n1; i++) {                                      \
      P4EST_FREE(a[i]);                                                 \
    }                                                                   \
  } while(0)

#define D4EST_ALLOC_MAT(a, n1, n2, size) do {                           \
    for (int i = 0; i < n1; i++) {                                      \
      for (int j = 0; j < n2; j++) {                                    \
        a[i][j] = P4EST_ALLOC(double, size);                            \
      }                                                                 \
    }                                                                   \
  } while(0)

#define D4EST_ALLOC_VEC(a, n1, size) do {                               \
    for (int i = 0; i < n1; i++) {                                      \
      a[i] = P4EST_ALLOC(double, size);                                 \
    }                                                                   \
  } while(0)

#define D4EST_ALLOC_DBYD_MAT(a, size) D4EST_ALLOC_MAT(a, (P4EST_DIM), (P4EST_DIM), size);
#define D4EST_FREE_DBYD_MAT(a) D4EST_FREE_MAT(a, (P4EST_DIM), (P4EST_DIM));
#define D4EST_ALLOC_DIM_VEC(a, size) D4EST_ALLOC_VEC(a, (P4EST_DIM), size);
#define D4EST_FREE_DIM_VEC(a) D4EST_FREE_VEC(a, (P4EST_DIM));

/* #define D4EST_INPUT_CHECK(section, param, param_default) do {           \ */
/*   if (param == param_default){                                          \ */
/*     printf("Please set %s in section %s\n", section, #param);           \ */
/*     mpi_abort("Parameter not set");                                     \ */
/*   }                                                              */
/* } while(0) */


#define D4EST_CHECK_INPUT(section, param, param_default) do {           \
    if (param == param_default) {                                       \
      printf("Please set %s in input file section %s\n",#param, section); \
      mpi_abort("");                                                    \
    }                                                                   \
  } while(0)
    


void util_print_matrix(double*,int,int,char[],int);
void util_print_matrix_int(int*,int,int,char[],int);
void util_print_matrices(double*,double*,int,int,char[]);
void util_print_3d_matrix(double*,int,int,int,char[],int);
void util_print_matrix_for_code(double* mat, int n, int m, char name []);
double util_uniform_rand(long int, double, double);
double util_max_error(double*, double*, int);
int util_compare_double(double,double,double);
int util_compare_vecs(double*,double*,int,double);
void util_gen_rand_vec(double*, int, long int, double, double);
int util_uniform_rand_int(long int seed, int a, int b);
void util_eye(double* eye, int N);
int util_int_pow_int(int x, int N);
double util_dbl_pow_int(double x, int N);
int util_sum_array_int(int* array, int N);
double util_min(double a, double b);
double util_max(double a, double b);
int util_levi_civita(int a, int b, int c);
double util_min(double a, double b);
double util_max(double a, double b);
int util_min_int(int a, int b);
int util_max_int(int a, int b);
int util_sort_double_callback(const void *p, const void *q);
void util_sort_double(double *a, size_t n);
void util_linear_regression(double*y,double* x,double* m,double* b,int n);
int util_compact_1st_alongwith_2nd(double *array, double* array2, int size);
int util_compact(double *array, int size);
double
util_parallel_checksum_dbl
(
 double* vec,
 int N
);
double
util_parallel_reduce_dbl_scalar
(
 double local
);

int util_bisection
(
 double funk(double, void*),
 double x1,
 double x2,
 double xacc,
 int iter,
 double* root,
 void* user
);

void
util_find_biggest_error
(
 double* a,
 double*b,
 int N,
 double* biggest_err,
 int* biggest_id
);

double
util_compute_median(double *x, int n);

double
util_normal_deviate (double mu, double sigma);
int util_does_file_exist(const char *filename);

void
util_print_matrix_for_mathematica(double* mat, int n, int m, const char* name);

int util_match_couple
(
 const char* section1,
 const char* section2,
 const char* name1,
 const char* name2
);

int util_match
(
 const char* str1,
 const char* str2
);

void
util_compute_error_array(double* arr1, double* arr2, double* err, int N);

#endif
