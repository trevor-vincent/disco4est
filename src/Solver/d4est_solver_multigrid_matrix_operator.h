#ifndef D4EST_SOLVER_MULTIGRID_MATRIX_OPERATOR_H
#define D4EST_SOLVER_MULTIGRID_MATRIX_OPERATOR_H 

#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_solver_multigrid.h>

typedef struct d4est_solver_multigrid_matrix_op_t d4est_solver_multigrid_matrix_op_t;

struct d4est_solver_multigrid_matrix_op_t{  
  /* if d4est_solver_multigrid is used as a preconditioner, this changes everytime a new newton iteration
   * occurs, which would mean the operator should change, which means we should update the 
   * the d4est_solver_multigrid matrix operator */
  int newton_iteration;
  
  int* matrix_nodes_on_level;
  int total_matrix_nodes_on_d4est_solver_multigrid;
  
  int fine_matrix_nodes;
  int coarse_matrix_nodes;
  int fine_matrix_stride;
  int coarse_matrix_stride;

  int stride_to_fine_matrix_data;
  int completed_alloc; /* completed_alloc flag (1/0) specifies whether we have already allocated memory for every mg level, this is to stop reallocation and total node calculation everytime we do a v-cycle, or use MG as a preconditioner. This flag gets zeroed everytime a matrix_op is inited and stays non-zero after the first v-cycle, so make sure the d4est_solver_multigrid hierachy does not change (in elements or element deg) without the matrix op being destroyed first */
  
  double* matrix_at0; //points to data at beginning (fine level)
  double* matrix; //points to data at current level
  
};


/* This file was automatically generated.  Do not edit! */
void d4est_solver_multigrid_matrix_operator_destroy(d4est_solver_multigrid_user_callbacks_t *user_callbacks);
void d4est_solver_multigrid_matrix_setup_fofufofvlilj_operator(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,double *u,double *v,d4est_xyzu_fcn_t fofu_fcn,void *fofu_ctx,d4est_xyzu_fcn_t fofv_fcn,void *fofv_ctx,d4est_solver_multigrid_matrix_op_t *matrix_op);
d4est_solver_multigrid_user_callbacks_t *d4est_solver_multigrid_matrix_operator_init(p4est_t *p4est,int num_of_levels);

#endif
