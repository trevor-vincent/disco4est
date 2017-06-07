#include <pXest.h>
#include <d4est_linalg.h>
#include <grid_functions.h>
#include <element_data.h>
#include <d4est_element_data.h>
#include <multigrid_matrix_operator.h>

static void
multigrid_matrix_operator_restriction_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data,
 int* degh,
 int degH,
 int num_children_h
)
{
  multigrid_data_t* mg_data = (multigrid_data_t*) info->p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;
  d4est_operators_t* d4est_ops = mg_data->d4est_ops;
  multigrid_matrix_op_t* matrix_op = mg_data->user_callbacks->user;

  int* fine_matrix_stride = &matrix_op->fine_matrix_stride;
  int* coarse_matrix_stride = &matrix_op->coarse_matrix_stride;
  int fine_matrix_nodes = matrix_op->fine_matrix_nodes;
  
  double* matrix_children = &matrix_op->matrix[*fine_matrix_stride];
  double* matrix_parent = &matrix_op->matrix[fine_matrix_nodes + *coarse_matrix_stride];

  d4est_operators_compute_PT_mat_P
    (
     d4est_ops,
     matrix_children,
     degH,
     (P4EST_DIM),
     degh,
     num_children_h,
     matrix_parent
    );

  for (int i = 0; i < num_children_h; i++){
    int fine_volume_nodes = d4est_operators_get_nodes( (P4EST_DIM) , degh[i] );
    (*fine_matrix_stride) += fine_volume_nodes*fine_volume_nodes;
  }
  
  int coarse_volume_nodes = d4est_operators_get_nodes( (P4EST_DIM) , degH);
  (*coarse_matrix_stride) += coarse_volume_nodes*coarse_volume_nodes;
}

static void
multigrid_matrix_operator_update_callback
(
 p4est_t* p4est,
 int level,
 problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_matrix_op_t* matrix_op = mg_data->user_callbacks->user;
  int vcycle = mg_data->vcycle_num_finished;

  if(mg_data->mg_state == PRE_V){
    if (vcycle == 0){
      matrix_op->matrix_nodes_on_level[level] = matrix_op->get_local_matrix_nodes(p4est);
      matrix_op->total_matrix_nodes_on_multigrid = matrix_op->matrix_nodes_on_level[level];
    }
    matrix_op->stride_to_fine_matrix_data = 0;
    matrix_op->coarse_matrix_stride = 0;
    matrix_op->fine_matrix_stride = 0;
    matrix_op->fine_matrix_nodes = 0;
    matrix_op->coarse_matrix_nodes = 0;
  }

  else if(
          mg_data->mg_state == UPV_PRE_SMOOTH ||
          mg_data->mg_state == DOWNV_PRE_SMOOTH ||
          mg_data->mg_state == COARSE_PRE_SOLVE
  ){
    matrix_op->matrix = &(matrix_op->matrix_at0)[matrix_op->stride_to_fine_matrix_data];
    vecs->user = matrix_op;
  }
  else if(mg_data->mg_state == DOWNV_POST_BALANCE){
    if (vcycle == 0){
      matrix_op->matrix_nodes_on_level[level-1] = matrix_op->get_local_matrix_nodes(p4est);
      matrix_op->total_matrix_nodes_on_multigrid += matrix_op->matrix_nodes_on_level[level-1];
      matrix_op->matrix_at0 = P4EST_REALLOC
                              (
                               matrix_op->matrix_at0,
                               double,
                               matrix_op->total_matrix_nodes_on_multigrid
                              );
    }
    /* setup for restriction */
    matrix_op->fine_matrix_nodes = matrix_op->matrix_nodes_on_level[level];
    matrix_op->coarse_matrix_nodes = matrix_op->matrix_nodes_on_level[level - 1];
  }
  else if(mg_data->mg_state == DOWNV_PRE_RESTRICTION){
      matrix_op->matrix = &((matrix_op->matrix_at0)[matrix_op->stride_to_fine_matrix_data]);      
      matrix_op->coarse_matrix_stride = 0;
      matrix_op->fine_matrix_stride = 0;
  }
  else if(mg_data->mg_state == DOWNV_POST_RESTRICTION){
    matrix_op->stride_to_fine_matrix_data += matrix_op->matrix_nodes_on_level[level];
  }
  else if(mg_data->mg_state == UPV_PRE_REFINE){
    matrix_op->stride_to_fine_matrix_data -= matrix_op->matrix_nodes_on_level[level+1];
  }  
  else {
    return;
  }
}

multigrid_user_callbacks_t*
multigrid_matrix_operator_init
(
 p4est_t* p4est,
 int num_of_levels,
 d4est_operators_t* d4est_ops,
 int(*get_local_matrix_nodes)(p4est_t*),
 void* user
)
{
  multigrid_matrix_op_t* matrix_op = P4EST_ALLOC(multigrid_matrix_op_t, 1);
  int local_matrix_nodes = get_local_matrix_nodes(p4est);
  matrix_op->matrix_at0 = P4EST_ALLOC(double, local_matrix_nodes);
  matrix_op->matrix_nodes_on_level = P4EST_ALLOC(int, num_of_levels);
  matrix_op->get_local_matrix_nodes = get_local_matrix_nodes;
  matrix_op->total_matrix_nodes_on_multigrid = -1;
  matrix_op->stride_to_fine_matrix_data = -1;
  
  matrix_op->fine_matrix_nodes = -1;
  matrix_op->coarse_matrix_nodes = -1;
  matrix_op->fine_matrix_stride = -1;
  matrix_op->coarse_matrix_stride = -1;
  matrix_op->matrix = matrix_op->matrix_at0;
  matrix_op->user = user;

  multigrid_user_callbacks_t* user_callbacks = P4EST_ALLOC(multigrid_user_callbacks_t, 1);

  user_callbacks->mg_prolong_user_callback = NULL;
  user_callbacks->mg_restrict_user_callback = multigrid_matrix_operator_restriction_callback;
  user_callbacks->update = multigrid_matrix_operator_update_callback;
  user_callbacks->user = matrix_op;
  
  return user_callbacks;
}

void
multigrid_matrix_fofu_fofv_mass_operator_setup_deg_quad_eq_deg
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double* u,
 double* v,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
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
        int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), ed->deg);
        int matrix_volume_nodes = volume_nodes*volume_nodes;
        int volume_nodes_Gauss = d4est_operators_get_nodes((P4EST_DIM), ed->deg);

        double* jac_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
        d4est_linalg_fill_vec(jac_Gauss, ed->jacobian, volume_nodes_Gauss);

        double* r_GL = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), ed->deg, 0);
        double* s_GL = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), ed->deg, 1);
        double* x_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        d4est_operators_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_Gauss);
        d4est_operators_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_Gauss);


#if (P4EST_DIM)==3        
        double* t_GL = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), ed->deg, 2);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        d4est_operators_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_Gauss);
#endif
        
        double* xyz_Gauss [(P4EST_DIM)] = {x_GL, y_GL
#if(P4EST_DIM)==3
                                           ,z_GL
#endif
                                          };
          
        d4est_operators_form_fofufofvlilj_matrix
          (
           d4est_ops,
           (u == NULL) ? NULL : &u[nodal_stride],
           (v == NULL) ? NULL : &v[nodal_stride],
           ed->deg,
           xyz_Gauss,
           jac_Gauss,
           ed->deg,
           QUAD_GAUSS,
           (P4EST_DIM),
           &matrix_op->matrix_at0[matrix_nodal_stride],
           fofu_fcn,
           fofu_ctx,
           fofv_fcn,
           fofv_ctx
          );

        P4EST_FREE(jac_Gauss);
        nodal_stride += volume_nodes;
        matrix_nodal_stride += matrix_volume_nodes;
      }
    }
}


void
multigrid_matrix_curved_fofu_fofv_mass_operator_setup_deg_quad_eq_deg
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double* u,
 double* v,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx,
 multigrid_matrix_op_t* matrix_op,
 int(*set_deg_quad)(void*, void*),
 void* set_deg_quad_ctx
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
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), ed->deg);
        int matrix_volume_nodes = volume_nodes*volume_nodes;

        d4est_element_data_form_fofufofvlilj_matrix
          (
           d4est_ops,
           d4est_geom,
           (u == NULL) ? NULL : &u[nodal_stride],
           (v == NULL) ? NULL : &v[nodal_stride],
           ed,
           set_deg_quad(ed, set_deg_quad_ctx),
           d4est_geom->geom_quad_type,
           (P4EST_DIM),
           &matrix_op->matrix_at0[matrix_nodal_stride],
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



void
multigrid_matrix_operator_destroy(multigrid_user_callbacks_t* user_callbacks)
{
  multigrid_matrix_op_t* matrix_op = user_callbacks->user;
  P4EST_FREE(matrix_op->matrix);
  P4EST_FREE(matrix_op->matrix_nodes_on_level);
  P4EST_FREE(matrix_op);
  P4EST_FREE(user_callbacks);
}
