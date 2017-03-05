#ifndef MULTIGRID_MATRIX_OPERATOR_H
#define MULTIGRID_MATRIX_OPERATOR_H 

#include <problem_data.h>
#include <multigrid.h>

typedef struct {

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
  
} multigrid_matrix_op_t;

void multigrid_matrix_operator_update_callback(p4est_t *p4est,int level, problem_data_t *vecs);
void multigrid_matrix_operator_restriction_callback(p4est_iter_volume_info_t *info,void *user_data,int *degh,int degH,int num_children_h);
void multigrid_matrix_operator_destroy(multigrid_matrix_op_t *matrix_op);

void
multigrid_matrix_fofu_fofv_mass_operator_setup_deg_integ_eq_deg
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u,
 double* v,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx,
 multigrid_matrix_op_t* matrix_op
);

multigrid_matrix_op_t*
multigrid_matrix_operator_init
(
 multigrid_data_t* mg_data,
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 void* user
);

#endif
