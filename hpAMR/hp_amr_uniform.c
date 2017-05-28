#include "../hpAMR/hp_amr.h"
#include "../hpAMR/hp_amr_uniform.h"
#include "../ElementData/element_data.h"
#include <linalg.h>

static void
hp_amr_uniform_refine_replace_callback ( 
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
#ifdef SAFETY  
  mpi_assert(num_outgoing == 1);
#endif
  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*) p4est->user_pointer;
  d4est_operators_t* d4est_ops = hp_amr_data->d4est_ops;
  /* hp_amr_smooth_pred_data_t* smooth_pred_data = (hp_amr_smooth_pred_data_t*) (hp_amr_data->hp_amr_scheme_data); */

  
  element_data_t* parent_data = (element_data_t*) outgoing[0]->p.user_data;
  element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degH);
    
  /* int h_pow = -1; */
  /* if (smooth_pred_data->norm_type == l2_norm_type) */
  /*   h_pow = parent_data->deg+1; */
  /* else if (smooth_pred_data->norm_type == dg_norm_type) */
  /*   h_pow = parent_data->deg; */
  /* else{ */
  /*   mpi_abort("[ERROR]: hp_amr_smooth_pred norm_type not supported"); */
  /*   h_pow *= -1.; /\* remove warnings *\/ */
  /* } */
  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  d4est_operators_apply_hp_prolong
    (
     d4est_ops,
     &(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    /* printf("child_data->deg = %d\n", child_data->deg); */
    child_data->local_predictor = parent_data->local_predictor;
    linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_elem[0], volume_nodes);
  }

  P4EST_FREE(temp_data);
}



static void
hp_amr_uniform_balance_replace_callback (
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
#ifdef SAFETY  
  mpi_assert(num_outgoing == 1);
#endif
  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*) p4est->user_pointer;
  d4est_operators_t* d4est_ops = hp_amr_data->d4est_ops;
  
  element_data_t* parent_data = (element_data_t*) outgoing[0]->p.user_data;
  element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degH);
    
  int h_pow = -1;
  /* if (smooth_pred_data->norm_type == l2_norm_type) */
    /* h_pow = parent_data->deg+1; */
  /* else if (smooth_pred_data->norm_type == dg_norm_type) */
    h_pow = parent_data->deg;
  /* else */
    /* mpi_abort("[ERROR]: hp_amr_smooth_pred norm_type not supported"); */

  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  d4est_operators_apply_hp_prolong
    (
     d4est_ops,
     &(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_elem[0], volume_nodes);
  }

  P4EST_FREE(temp_data);
}


void
hp_amr_uniform_set_refinement
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  hp_amr_data_t* amr_data = (hp_amr_data_t*) info->p4est->user_pointer;
  element_data_t* elem_data = (element_data_t*) info->quad->p.user_data;
  amr_data->refinement_log[info->quadid] = -elem_data->deg;
}

void
hp_amr_uniform_destroy
(
 hp_amr_scheme_t* scheme
)
{
  P4EST_FREE(scheme);
}

hp_amr_scheme_t*
hp_amr_uniform
()
{
  hp_amr_scheme_t* scheme = P4EST_ALLOC(hp_amr_scheme_t, 1);
  scheme->iter_volume = hp_amr_uniform_set_refinement;
  scheme->balance_replace_callback_fcn_ptr = hp_amr_uniform_balance_replace_callback;
  scheme->refine_replace_callback_fcn_ptr = hp_amr_uniform_refine_replace_callback;;
  scheme->post_balance_callback = NULL;
  scheme->pre_refine_callback = NULL;
  scheme->hp_amr_scheme_data = NULL;
  return scheme;
}
