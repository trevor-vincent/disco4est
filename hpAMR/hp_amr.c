#define SAFETY
#define NDEBUG

#include "../hpAMR/hp_amr.h"
#include "../GridFunctions/grid_functions.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/linalg.h"
#include "../ElementData/element_data.h"
#include "../ElementData/curved_element_data.h"
#include <hacked_p4est_vtk.h>



static int
hp_amr_refine_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{
  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*) p4est->user_pointer;
  element_data_t* elem_data = (element_data_t*) quadrant->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = hp_amr_data->dgmath_jit_dbase;
  
  int* refinement_log = hp_amr_data->refinement_log;
  int* refinement_log_stride = &hp_amr_data->refinement_log_stride;
  
  /* h-refine */
  if (refinement_log[*refinement_log_stride] < 0){
    (*refinement_log_stride)++;
    hp_amr_data->elements_marked_for_hrefine += 1;
    /* printf("Element will be refined, rank=%d, id=%d, estimator=%.25f\n", p4est->mpirank, elem_data->id, elem_data->local_estimator); */
    return 1;
  }

  /* p-refine, p-coarsen or do nothing */
  else {
    int degH = elem_data->deg;
    int degh = refinement_log[*refinement_log_stride];
    int volume_nodes = dgmath_get_nodes((P4EST_DIM), degh);
    
    double* temp_data = P4EST_ALLOC(double, volume_nodes);
    
    if (elem_data->deg < degh){
      hp_amr_data->elements_marked_for_prefine += 1;
      dgmath_apply_p_prolong
        (
         dgmath_jit_dbase,
         &elem_data->u_storage[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data        
        );
      linalg_copy_1st_to_2nd(temp_data, &elem_data->u_storage[0], volume_nodes);
    }
    else if (elem_data->deg > degh){
      dgmath_apply_p_restrict
        (
         dgmath_jit_dbase,
         &elem_data->u_storage[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data
        );
      linalg_copy_1st_to_2nd(temp_data, &elem_data->u_storage[0], volume_nodes);
    }     

    elem_data->deg = degh;
    P4EST_FREE(temp_data);
    (*refinement_log_stride)++;
    return 0;
  }  
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
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double** data_to_hp_refine,
 estimator_stats_t* stats,
 hp_amr_scheme_t* scheme,
 int curved
)
{
  hp_amr_data_t hp_amr_data;
  hp_amr_data.refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
  hp_amr_data.refinement_log_stride = 0;
  hp_amr_data.data = *data_to_hp_refine;
  hp_amr_data.hp_amr_scheme_data = scheme->hp_amr_scheme_data;
  hp_amr_data.dgmath_jit_dbase = dgmath_jit_dbase;
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
  
  /* iterate and store refinement boolean */
  if (!curved){
  element_data_copy_from_vec_to_storage
    (
     p4est,
     *data_to_hp_refine
    );
  }
  else {
    curved_element_data_copy_from_vec_to_storage
      (
       p4est,
       *data_to_hp_refine
      );
  }

  /* store p4est user pointer */
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = &hp_amr_data;
  
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
     hp_amr_data.refine_replace_callback_fcn_ptr
    );

  /* 2-1 balance mesh */
  p4est_balance_ext
    (
     p4est,
     P4EST_CONNECT_FULL,
     NULL,
     hp_amr_data.balance_replace_callback_fcn_ptr
    );
  
  /* restore pointer back to original state */
  p4est->user_pointer = tmp;

  int new_nodes;
  if (!curved){
    new_nodes = element_data_get_local_nodes(p4est);
  }
  else {
    new_nodes = curved_element_data_get_local_nodes(p4est);
  }
  
  /* realloc room for new refined mesh data */
  *data_to_hp_refine = P4EST_REALLOC(*data_to_hp_refine,
                                     double,
                                     new_nodes
                                    );

  if (!curved){
    element_data_copy_from_storage_to_vec
      (
       p4est,
       *data_to_hp_refine
      );
  }
  else {
    curved_element_data_copy_from_storage_to_vec
      (
       p4est,
       *data_to_hp_refine
      );    
  }
  
  if(scheme->post_balance_callback != NULL){
    scheme->post_balance_callback(p4est, scheme->hp_amr_scheme_data);
  }
  
  P4EST_FREE(hp_amr_data.refinement_log);
}

