#ifndef D4EST_SOLVER_SCHWARZ_OPERATORS_H
#define D4EST_SOLVER_SCHWARZ_OPERATORS_H 

#include <pXest.h>
#include <d4est_operators.h>

typedef enum {SCHWARZ_WEIGHT_NONE,
              SCHWARZ_WEIGHT_HAT_POLY_QUINTIC} d4est_solver_schwarz_weight_function_t;

typedef struct {

  d4est_operators_t* d4est_ops;

  double*** schwarz_restrictor_1d_table;
  double*** schwarz_nodes_1d_table;
  double*** schwarz_weights_1d_table;
  double*** schwarz_restrictor_transpose_1d_table;

} d4est_solver_schwarz_operators_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_operators_destroy(d4est_solver_schwarz_operators_t *d4est_schwarz_ops);
d4est_solver_schwarz_operators_t *d4est_solver_schwarz_operators_init(d4est_operators_t *d4est_ops);
void d4est_solver_schwarz_operators_apply_schwarz_weights(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,double *in,int dim,int *faces,int deg,int restricted_size,double *out);
void d4est_solver_schwarz_operators_apply_schwarz_weights_test(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,double *in,int dim,int deg,int restricted_size,double *out);
void d4est_solver_schwarz_operators_apply_schwarz_restrictor(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,double *in,int dim,int *faces,int deg,int restricted_size,d4est_ops_transpose_t transpose,double *out);
double *d4est_solver_schwarz_operators_fetch_schwarz_restrictor_transpose_1d(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,int deg,int res_size);
void d4est_solver_schwarz_operators_build_schwarz_restrictor_transpose_1d(d4est_operators_t *d4est_ops,double *restrict restrictor_transpose_1d,int deg,int res_size,void *ctx);
double *d4est_solver_schwarz_operators_fetch_schwarz_weights_1d(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,int deg,int restricted_size);
void d4est_solver_schwarz_operators_build_schwarz_weights_1d(d4est_operators_t *d4est_ops,double *restrict schwarz_weights_1d,int deg,int restricted_size);
double *d4est_solver_schwarz_operators_fetch_schwarz_restrictor_1d(d4est_solver_schwarz_operators_t *d4est_schwarz_ops,int deg,int restricted_size);
void d4est_solver_schwarz_operators_build_schwarz_restrictor_1d(d4est_operators_t *d4est_ops,double *restrict restrictor_1d,int deg,int restricted_size);

#endif
