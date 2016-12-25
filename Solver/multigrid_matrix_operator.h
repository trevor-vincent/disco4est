#ifndef MULTIGRID_MATRIX_OPERATOR_H
#define MULTIGRID_MATRIX_OPERATOR_H 

typedef struct {

  int* matrix_nodes_on_level;
  int total_matrix_nodes_on_multigrid;
  
  int fine_matrix_nodes;
  int coarse_matrix_nodes;
  
  int fine_matrix_stride;
  int coarse_matrix_stride;

  double* matrix_at0;
  double* matrix;
  
} multigrid_matrix_op_t;

multigrid_matrix_opt_t*
multigrid_matrix_operator_init
(
 multigrid_data_t* mg_data,
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
)
{
  multigrid_matrix_op_t* matrix_op = P4EST_ALLOC(multigrid_matrix_op_t, 1);
  int local_matrix_nodes = element_get_local_matrix_nodes(p4est);
  matrix_op->matrix_at0 = P4EST_ALLOC(double, local_matrix_nodes);
  matrix_op->matrix_nodes_on_level = P4EST_ALLOC(int, mg_data->num_of_levels);
  matrix_op->total_matrix_nodes_on_multigrid = 0;

  matrix_op->fine_matrix_nodes = -1;
  matrix_op->coarse_matrix_nodes = -1;
  matrix_op->fine_matrix_stride = -1;
  matrix_op->coarse_matrix_stride = -1;
  matrix_op->matrix = NULL;
  
  return matrix_op
}

void
multigrid_matrix_fofu_fofv_mass_operator_setup
(
 p4est_t*,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* u,
 double* v,
 int deg_Lobatto,
 double* jac_Gauss,
 int deg_Gauss,
 int dim,
 double (*fofu_fcn)(double, void*),
 void* fofu_ctx,
 double (*fofv_fcn)(double, void*),
 void* fofv_ctx,
 multigrid_matrix_op_t* matrix_op
)
{
  int nodal_stride = 0;
  int matrix_nodal_stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        int matrix_volume_nodes = volume_nodes*volume_nodes;

        dgmath_form_fofufofvlilj_matrix_Gaussnodes
          (
           dgmath_jit_dbase,
           (u == NULL) ? NULL : &u[nodal_stride],
           (v == NULL) ? NULL : &v[nodal_stride],
           deg_Lobatto,
           jac_Gauss,
           deg_Gauss,
           (P4EST_DIM),
           &matrix_op->matrix_at0[nodal_matrix_stride],
           fofu_fcn,
           fofu_ctx,
           fofv_fcn,
           fofv_ctx
          );

        nodal_stride += volume_nodes;
        matrix_nodal_stride += matrix_volume_nodes;
      }
    }
}

multigrid_matrix_operator_destroy(multigrid_matrix_opt_t* matrix_op)
{
  P4EST_FREE(matrix_op->matrix_at0);
  P4EST_FREE(matrix_op->matrix_nodes_on_level);
  P4EST_FREE(matrix_op);
}

void
multigrid_matrix_operator_update_on_restriction
(
 p4est_iter_volume_info_t* info,
 void* user_data,
 int* degh,
 int degH,
 int num_children_h,
)
{
  multigrid_data_t* mg_data = (multigrid_data_t*) info->p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;
  multigrid_matrix_op_t* matrix_op = mg_data->user_ctx;

  int* fine_matrix_stride = &matrix_op->fine_matrix_stride;
  int* coarse_matrix_stride = &matrix_op->coarse_matrix_stride;
  int fine_matrix_nodes = matrix_op->fine_matrix_nodes;
  
  double* matrix_children = &matrix[*fine_matrix_stride];
  double* matrix_parent = &matrix[fine_matrix_nodes + *coarse_matrix_stride];

  dgmath_compute_PT_mat_P
    (
     dgmath_jit_dbase,
     matrix_children,
     degH,
     (P4EST_DIM),
     degh,
     num_children_h,
     matrix_parent
    );

  for (int i = 0; i < num_children_h; i++){
    int fine_volume_nodes = dgmath_get_nodes( (P4EST_DIM) , degh[i] );
    (*fine_matrix_stride) += fine_volume_nodes*fine_volume_nodes;
  }
  
  int coarse_volume_nodes = dgmath_get_nodes( (P4EST_DIM) , degH)
  (*coarse_matrix_stride) += coarse_volume_nodes*coarse_volume_nodes;
}


void
multigrid_matrix_operator_update_size_or_stride
(
 p4est_t* p4est,
 int level,
 problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_matrix_op_t* matrix_op = mg_data->user_ctx;
  
  /**
   * Assumes V starts at a fine level of num_of_levels - 1 and drops to level 0 
   * at the bottom of the V and then rides back up to level = num_of_levels - 1
   * 
   */

  
  if(
      mg_data->mg_state == DOWNV_POST_BALANCE
    )
    {
      /* update after coarsening and balancing*/
      matrix_op->matrix_nodes_on_level[level] = element_data_get_local_matrix_nodes(p4est);
      matrix_op->total_matrix_nodes_on_multigrid += matrix_op->matrix_nodes_on_level[level];
      matrix_op->matrix_at0 = P4EST_ALLOC
                              (
                               matrix_op->matrix_at0,
                               double,
                               matrix_op->total_matrix_nodes_on_multigrid
                              );

      /* setup for restriction */
      matrix_op->fine_matrix_nodes = matrix_op->matrix_nodes_on_level[level+1];
      matrix_op->coarse_matrix_nodes = matrix_op->matrix_nodes_on_level[level];
      matrix_op->matrix = &matrix_at0[total_nodes_on_multigrids - mg_data->fine_matrix_nodes - mg_data->coarse_matrix_nodes];
      matrix_op->coarse_matrix_stride = 0;
      matrix_op->fine_matrix_stride = 0;
    }
  else if
    (
     mg_data->mg_state == PREV_PRE_SMOOTH ||
     mg_data->mg_state == DOWNV_PRE_SMOOTH ||
     mg_data->mg_state == COARSE_PRE_SOLVE ||
     mg_data->mg_state == UPV_PRE_SMOOTH ||
     mg_data->mg_state == POSTV_PRE_SMOOTH
    )
    {
      int stride = 0;
      for (int i = num_of_levels - 1; i > level; i--){
        stride += matrix_op->matrix_nodes_on_level[level];
      }
      matrix_op->matrix = &matrix_at0[stride];
      vecs->user_ctx = matrix_op;
    }
  else {
    return;
  }
}


#endif
