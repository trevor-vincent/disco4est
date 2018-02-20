#ifndef UTIL_H
#define UTIL_H

#define _GNU_SOURCE
#include <stdio.h>
#include <sc.h>

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
    printf("%s sum = %.25f\n",#a, sum);                   \
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




#define D4EST_NOOP()                                                            \
  do                                                                           \
  {                                                                            \
  } while (0)

#ifndef NDEBUG
#define D4EST_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define D4EST_ASSERT(c) D4EST_NOOP()
#endif

#define D4EST_ALLOC(a,b) P4EST_ALLOC(a,b)
#define D4EST_FREE(a) P4EST_FREE(a)
#define D4EST_ALLOC_ZERO(a,b) P4EST_ALLOC_ZERO(a,b)
#define D4EST_REALLOC(a,b,c) P4EST_REALLOC(a,b,c)

#define D4EST_ABORT(c) SC_ABORT(c)
#define D4EST_ABORT_IF(q, ...) ((q) ? D4EST_ABORT(__VA_ARGS__) : (void)0)
#define D4EST_ABORT_IF_NOT(q, ...) D4EST_ABORT_IF(!(q), __VA_ARGS__)


/* #define D4EST_ASPRINTF(write_to, ...) {         \ */
/*     char *tmp_string_for_extend = (write_to);   \ */
/*     if(asprintf(&(write_to), __VA_ARGS__) <0){  \ */
/*       D4EST_ABORT("ASPRINTF_ERROR");              \ */
/*     }                                           \ */
/*     free(tmp_string_for_extend);                \ */
/*   } */

#define D4EST_FREE_MAT(a, n1, n2) do {          \
    for (int i = 0; i < n1; i++) {              \
      for (int j = 0; j < n2; j++) {            \
        P4EST_FREE(a[i][j]);                    \
      }                                         \
    }                                           \
  } while(0)

#define D4EST_FREE_VEC(a, n1) do {              \
    for (int i = 0; i < n1; i++) {              \
      P4EST_FREE(a[i]);                         \
    }                                           \
  } while(0)

#define D4EST_ALLOC_MAT(a, n1, n2, size) do {   \
    for (int i = 0; i < n1; i++) {              \
      for (int j = 0; j < n2; j++) {            \
        a[i][j] = P4EST_ALLOC(double, size);    \
      }                                         \
    }                                           \
  } while(0)


#define D4EST_COPY_MAT(a, b, n1, n2) do {       \
    for (int i = 0; i < n1; i++) {              \
      for (int j = 0; j < n2; j++) {            \
        b[i][j] = a[i][j];                      \
      }                                         \
    }                                           \
  } while(0)

#define D4EST_COPY_VEC(a, b, n1) do {           \
    for (int i = 0; i < n1; i++) {              \
      b[i] = a[i];                              \
    }                                           \
  } while(0)

/*  */
#define D4EST_ALLOC_VEC(a, n1, size) do {       \
    for (int i = 0; i < n1; i++) {              \
      a[i] = P4EST_ALLOC(double, size);         \
    }                                           \
  } while(0)

#define D4EST_ALLOC_DBYD_MAT(a, size) D4EST_ALLOC_MAT(a, (P4EST_DIM), (P4EST_DIM), size);
#define D4EST_COPY_DBYD_MAT(a, b) D4EST_COPY_MAT(a, b, (P4EST_DIM), (P4EST_DIM));
#define D4EST_FREE_DBYD_MAT(a) D4EST_FREE_MAT(a, (P4EST_DIM), (P4EST_DIM));
#define D4EST_ALLOC_DIM_VEC(a, size) D4EST_ALLOC_VEC(a, (P4EST_DIM), size);
#define D4EST_COPY_DIM_VEC(a, b) D4EST_COPY_VEC(a, b, (P4EST_DIM));
#define D4EST_FREE_DIM_VEC(a) D4EST_FREE_VEC(a, (P4EST_DIM));

/* #define D4EST_INPUT_CHECK(section, param, param_default) do {           \ */
/*   if (param == param_default){                                          \ */
/*     printf("Please set %s in section %s\n", section, #param);           \ */
/*     D4EST_ABORT("Parameter not set");                                     \ */
/*   }                                                              */
/* } while(0) */


#define D4EST_CHECK_INPUT(section, param, param_default) do {           \
    if (param == param_default) {                                       \
      printf("Please set %s in input file section %s\n",#param, section); \
      D4EST_ABORT("");                                                    \
    }                                                                   \
  } while(0)


typedef enum {D4EST_INT, D4EST_DOUBLE} d4est_builtin_t;


/* This file was automatically generated.  Do not edit! */
void d4est_util_gen_rand_vec(double *vec,int N,int seed,double a,double b);
double d4est_util_max_dbl_array(double *arr,int N);
double d4est_util_min_dbl_array(double *arr,int N);
void d4est_util_compute_error_array(double *arr1,double *arr2,double *err,int N);
int d4est_util_match(const char *str1,const char *str2);
int d4est_util_match_couple(const char *section1,const char *section2,const char *name1,const char *name2);
int d4est_util_does_file_exist(const char *filename);
double d4est_util_normal_deviate(double mu,double sigma);
double d4est_util_compute_median(double *x,int n);
int d4est_util_bisection(double funk(double,void *),double x1,double x2,double xacc,int iter,double *root,void *user);
void d4est_util_linear_regression(double *y,double *x,double *m,double *b,int n);
double d4est_util_parallel_reduce_dbl_scalar(double local);
double d4est_util_parallel_checksum_dbl(double *vec,int N);
void d4est_util_sort_double(double *a,size_t n);
double d4est_util_sum_array_dbl(double *array,int N);
int d4est_util_sum_array_int(int *array,int N);
double d4est_util_max_error(double *u,double *u_sol,int N);
int d4est_util_uniform_rand_int(long int seed,int a,int b);
double d4est_util_uniform_rand(long int seed,double a,double b);
void d4est_util_print_matrix_for_mathematica(double *mat,int n,int m,const char *name);
void d4est_util_print_matrix_for_code(double *mat,int n,int m,char name[]);
void d4est_util_print_matrices(double *mat,double *mat2,int n,int m,char message[]);
void d4est_util_print_matrix_int(int *mat,int n,int m,char message[],int print_rank);
void d4est_util_print_matrix(double *mat,int n,int m,char message[],int print_rank);
void d4est_util_print_3d_matrix(double *mat_3d,int n,int m,int l,char message[],int print_rank);
void d4est_util_find_biggest_error(double *a,double *b,int N,double *biggest_err,int *biggest_id);
int d4est_util_compare_vecs(double *a,double *b,int N,double eps);
int d4est_util_compare_double(double a,double b,double eps);
int d4est_util_int_pow_int(int base,int exp);
double d4est_util_dbl_pow_int(double a,int b);
void d4est_util_eye(double *eye,int N);
int d4est_util_max_int(int a,int b);
int d4est_util_min_int(int a,int b);
double d4est_util_max(double a,double b);
double d4est_util_min(double a,double b);
void d4est_util_copy_1st_to_2nd(double *v1,double *v2,int N);
int d4est_util_compact_1st_alongwith_2nd(double *array,double *array2,int size);
int d4est_util_compact(double *array,int size);
void d4est_util_make_directory(const char *dir,int add_cwd_to_dir);
char *d4est_util_add_cwd(const char *dir);
double d4est_util_secant_fcn(double x);

#endif
