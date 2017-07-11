#ifndef MULTIGRID_MATRIX_OPERATOR_H
#define MULTIGRID_MATRIX_OPERATOR_H 

#include <problem_data.h>
#include <multigrid.h>

typedef struct multigrid_matrix_op_t multigrid_matrix_op_t;

struct multigrid_matrix_op_t{  
  int(*get_local_matrix_nodes)(p4est_t*);
       
  int* matrix_nodes_on_level;
  int total_matrix_nodes_on_multigrid;
  
  int fine_matrix_nodes;
  int coarse_matrix_nodes;
  int fine_matrix_stride;
  int coarse_matrix_stride;

  int stride_to_fine_matrix_data;
  double* matrix_at0; //points to data at beginning (fine level)
  double* matrix; //points to data at current level
  void* user;
  
};


/* This file was automatically generated.  Do not edit! */
void multigrid_matrix_operator_destroy(multigrid_user_callbacks_t *user_callbacks);
void multigrid_matrix_curved_fofu_fofv_mass_operator_setup_deg_quad_eq_deg(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,double *u,double *v,d4est_xyz_fcn_ext_t fofu_fcn,void *fofu_ctx,d4est_xyz_fcn_ext_t fofv_fcn,void *fofv_ctx,multigrid_matrix_op_t *matrix_op,int(*set_deg_gauss)(void *,void *),void *set_deg_gauss_ctx);
void multigrid_matrix_fofu_fofv_mass_operator_setup_deg_quad_eq_deg(p4est_t *p4est,d4est_operators_t *d4est_ops,double *u,double *v,d4est_xyz_fcn_ext_t fofu_fcn,void *fofu_ctx,d4est_xyz_fcn_ext_t fofv_fcn,void *fofv_ctx,multigrid_matrix_op_t *matrix_op);
multigrid_user_callbacks_t *multigrid_matrix_operator_init(p4est_t *p4est,int num_of_levels,d4est_operators_t *d4est_ops,int(*get_local_matrix_nodes)(p4est_t *),void *user);


#endif
