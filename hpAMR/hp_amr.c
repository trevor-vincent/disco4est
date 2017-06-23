#define SAFETY
#define NDEBUG

#include "../hpAMR/hp_amr.h"
#include "../GridFunctions/grid_functions.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../ElementData/element_data.h"
#include "../ElementData/d4est_element_data.h"


static int
hp_amr_refine_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{
  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*) p4est->user_pointer;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) quadrant->p.user_data;
  d4est_operators_t* d4est_ops = hp_amr_data->d4est_ops;
  
  int* refinement_log = hp_amr_data->refinement_log;
  int* refinement_log_stride = &hp_amr_data->refinement_log_stride;
  
  /* h-refine */
  if (refinement_log[*refinement_log_stride] < 0){
    (*refinement_log_stride)++;
    hp_amr_data->elements_marked_for_hrefine += 1;
    return 1;
  }

  /* p-refine, p-coarsen or do nothing */
  else {
    int degH = ed->deg; 
    int degh = refinement_log[*refinement_log_stride];
    int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degh);

    double* stored_data = hp_amr_data->get_storage(elem_data);
    double* temp_data = P4EST_ALLOC(double, volume_nodes);
    
    if (elem_data->deg < degh){
      hp_amr_data->elements_marked_for_prefine += 1;
      d4est_operators_apply_p_prolong
        (
         d4est_ops,
         stored_data,//&elem_data->u_elem[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data        
        );
      d4est_linalg_copy_1st_to_2nd(temp_data, stored_data, volume_nodes);
    }
    else if (elem_data->deg > degh){
      d4est_operators_apply_p_restrict
        (
         d4est_ops,
         stored_data,//&elem_data->u_elem[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data
        );
      d4est_linalg_copy_1st_to_2nd(temp_data, stored_data, volume_nodes);
    }     

    elem_data->deg = degh;
    P4EST_FREE(temp_data);
    (*refinement_log_stride)++;
    return 0;
  }  
}



static void
hp_amr_refine_replace_callback (
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
  /* hp_amr_curved_smooth_pred_data_t* smooth_pred_data = (hp_amr_curved_smooth_pred_data_t*) (hp_amr_data->hp_amr_scheme_data); */
  
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degH);

  double* stored_parent_data = hp_amr_data->get_storage(parent_data);
  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  d4est_operators_apply_hp_prolong
    (
     d4est_ops,
     stored_parent_data,
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    double* stored_child_data = hp_amr_data->get_storage(child_data)
    child_data->deg = parent_data->deg;
    d4est_linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], stored_child_data, volume_nodes);
  }

  P4EST_FREE(temp_data);

  if(hp_amr_data->refine_replace_callback_fcn_ptr != NULL)
    hp_amr_data->refine_replace_callback_fcn_ptr(p4est, which_tree, num_outgoing, outgoing, num_incoming, incoming);
}



static void
hp_amr_balance_replace_callback (
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
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;

  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degH);    

double* stored_parent_data = hp_amr_data->get_storage(parent_data);
 double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  d4est_operators_apply_hp_prolong
    (
     d4est_ops,
     stored_parent_data,//&(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    double* stored_child_data = hp_amr_data->get_storage(child_data)
    child_data->deg = parent_data->deg;
    /* d4est_linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], */
    /*                              &child_data->u_elem[0], */
    /*                              volume_nodes); */
    d4est_linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i],
                                 stored_child_data,
                                 volume_nodes);
  }

  P4EST_FREE(temp_data);

  if(hp_amr_data->balance_replace_callback_fcn_ptr != NULL)
    hp_amr_data->balance_replace_callback_fcn_ptr(p4est, which_tree, num_outgoing, outgoing, num_incoming, incoming);
}

/** 
 * Assumes iter_volume will impose a refinement criterion
 * 
 * @param p4est 
 * @param data_to_hp_refine 
 * @param iter_volume 
 * @param hp_amr_refine_store_fcn_ptr if NULL, we use default
 * @param hp_amr_balance_store_fcn_ptr if NULL, we use default
 */
void
hp_amr
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double** data_to_hp_refine,
 estimator_stats_t** stats,
 hp_amr_scheme_t* scheme,
)
{
  hp_amr_data_t hp_amr_data;
  hp_amr_data.refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
  hp_amr_data.refinement_log_stride = 0;
  hp_amr_data.data = *data_to_hp_refine;
  hp_amr_data.hp_amr_scheme_data = scheme->hp_amr_scheme_data;
  hp_amr_data.d4est_ops = d4est_ops;
  hp_amr_data.estimator_stats = stats;
  hp_amr_data.elements_marked_for_hrefine = 0;
  hp_amr_data.elements_marked_for_prefine = 0;
  hp_amr_data.refine_replace_callback_fcn_ptr
    = scheme->refine_replace_callback_fcn_ptr;
  hp_amr_data.balance_replace_callback_fcn_ptr
    = scheme->balance_replace_callback_fcn_ptr;
  
  if(scheme->pre_refine_callback != NULL){
    scheme->pre_refine_callback(p4est, scheme->hp_amr_scheme_data);
  }
  

  /* copy field to interpolate to storage */
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
        double* storage = hp_amr_data->get_storage(ed);
        d4est_linalg_copy_1st_to_2nd
          (
           &((*data_to_hp_refine)[ed->nodal_stride]),
           storage,
           volume_nodes
          );
      }
    }

  /* store p4est user pointer */
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = &hp_amr_data;

  /* iterate and store refinement boolean */
  p4est_iterate(p4est,
                NULL,
                NULL,
                scheme->iter_volume,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  
  /* interpolate solution and refine mesh */
  p4est_refine_ext
    (
     p4est,
     0,
     -1,
     hp_amr_refine_callback,
     NULL,
     hp_amr_refine_replace_callback_fcn_ptr
    );
 
  /* 2-1 balance mesh */
  p4est_balance_ext
    (
     p4est,
     P4EST_CONNECT_FULL,
     NULL,
     hp_amr_balance_replace_callback_fcn_ptr
    );


  
  /* restore pointer back to original state */
  p4est->user_pointer = tmp;

  int post_balance_local_nodes = 0;
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
        local_nodes += volume_nodes;
     }
    }  
  
  /* realloc room for new refined mesh data */
  *data_to_hp_refine = P4EST_REALLOC(*data_to_hp_refine,
                                     double,
                                     post_balance_local_nodes
                                    );


  /* copy storage back to field */
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
        double* storage = hp_amr_data->get_storage(ed);
        d4est_linalg_copy_1st_to_2nd
          (
           storage,
           &((*data_to_hp_refine)[ed->nodal_stride]),
           volume_nodes
          );
      }
    }

  if(scheme->post_balance_callback != NULL){
    scheme->post_balance_callback(p4est, scheme->hp_amr_scheme_data);
  }
  
  P4EST_FREE(hp_amr_data.refinement_log);
}

