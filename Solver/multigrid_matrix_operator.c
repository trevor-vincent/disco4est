#include <pXest.h>
#include <linalg.h>
#include <grid_functions.h>
#include <element_data.h>
#include <multigrid_matrix_operator.h>


multigrid_matrix_op_t*
multigrid_matrix_operator_init
(
 multigrid_data_t* mg_data,
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 void* user
)
{
  multigrid_matrix_op_t* matrix_op = P4EST_ALLOC(multigrid_matrix_op_t, 1);
  int local_matrix_nodes = element_data_get_local_matrix_nodes(p4est);
  matrix_op->matrix_at0 = P4EST_ALLOC(double, local_matrix_nodes);
  matrix_op->matrix_nodes_on_level = P4EST_ALLOC(int, mg_data->num_of_levels);
  
  matrix_op->total_matrix_nodes_on_multigrid = -1;
  matrix_op->stride_to_fine_matrix_data = -1;
  
  matrix_op->fine_matrix_nodes = -1;
  matrix_op->coarse_matrix_nodes = -1;
  matrix_op->fine_matrix_stride = -1;
  matrix_op->coarse_matrix_stride = -1;
  matrix_op->matrix = matrix_op->matrix_at0;
  matrix_op->user = user;

  mg_data->user_defined_fields = 1;
  mg_data->mg_prolong_user_callback = NULL;
  mg_data->mg_restrict_user_callback = multigrid_matrix_operator_restriction_callback;
  mg_data->mg_update_user_callback = multigrid_matrix_operator_update_callback;
  mg_data->user_ctx = matrix_op;
  
  return matrix_op;
}

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
        int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), ed->deg);

        double* jac_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
        linalg_fill_vec(jac_Gauss, ed->jacobian, volume_nodes_Gauss);

        double* r_GL = dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), ed->deg, 0);
        double* s_GL = dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), ed->deg, 1);
        double* x_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        double* y_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        dgmath_rtox_array(r_GL, ed->xyz_corner[0], ed->h, x_GL, volume_nodes_Gauss);
        dgmath_rtox_array(s_GL, ed->xyz_corner[1], ed->h, y_GL, volume_nodes_Gauss);


#if (P4EST_DIM)==3        
        double* t_GL = dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), ed->deg, 2);
        double* z_GL = P4EST_ALLOC(double, volume_nodes_Gauss);
        dgmath_rtox_array(t_GL, ed->xyz_corner[2], ed->h, z_GL, volume_nodes_Gauss);
#endif
        
        double* xyz_Gauss [(P4EST_DIM)] = {x_GL, y_GL
#if(P4EST_DIM)==3
                                           ,z_GL
#endif
                                          };
          
        dgmath_form_fofufofvlilj_matrix_Gaussnodes
          (
           dgmath_jit_dbase,
           (u == NULL) ? NULL : &u[nodal_stride],
           (v == NULL) ? NULL : &v[nodal_stride],
           ed->deg,
           xyz_Gauss,
           jac_Gauss,
           ed->deg,
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
multigrid_matrix_operator_destroy(multigrid_matrix_op_t* matrix_op)
{
  P4EST_FREE(matrix_op->matrix);
  P4EST_FREE(matrix_op->matrix_nodes_on_level);
  P4EST_FREE(matrix_op);
}

void
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
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;
  multigrid_matrix_op_t* matrix_op = mg_data->user_ctx;

  int* fine_matrix_stride = &matrix_op->fine_matrix_stride;
  int* coarse_matrix_stride = &matrix_op->coarse_matrix_stride;
  int fine_matrix_nodes = matrix_op->fine_matrix_nodes;
  
  double* matrix_children = &matrix_op->matrix[*fine_matrix_stride];
  double* matrix_parent = &matrix_op->matrix[fine_matrix_nodes + *coarse_matrix_stride];

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
  
  int coarse_volume_nodes = dgmath_get_nodes( (P4EST_DIM) , degH);
  (*coarse_matrix_stride) += coarse_volume_nodes*coarse_volume_nodes;
}


void
multigrid_matrix_operator_update_callback
(
 p4est_t* p4est,
 int level,
 problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_matrix_op_t* matrix_op = mg_data->user_ctx;
  
  if(mg_data->mg_state == PRE_V){
    matrix_op->matrix_nodes_on_level[level] = element_data_get_local_matrix_nodes(p4est);
    matrix_op->total_matrix_nodes_on_multigrid = matrix_op->matrix_nodes_on_level[level];
    matrix_op->stride_to_fine_matrix_data = 0;
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
    matrix_op->matrix_nodes_on_level[level-1] = element_data_get_local_matrix_nodes(p4est);
    matrix_op->total_matrix_nodes_on_multigrid += matrix_op->matrix_nodes_on_level[level-1];
    matrix_op->matrix_at0 = P4EST_REALLOC
                            (
                             matrix_op->matrix_at0,
                             double,
                             matrix_op->total_matrix_nodes_on_multigrid
                            );

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
