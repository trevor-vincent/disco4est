#ifndef HP_AMR_SMOOTH_PRED_H
#define HP_AMR_SMOOTH_PRED_H 

#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"
#include "../ElementData/element_data.h"

typedef struct {

  /* parameters per tree */
  double* sigma;
  double* gamma_h;
  double* gamma_p;
  double* gamma_n;
  
  int num_trees;
  int max_degree;
  norm_t norm_type;

  int pred_storage_stride;
  double pred_storage_size;
  double* pred_storage; 
  int* percentiles;
  
} curved_hp_amr_smooth_pred_data_t;

curved_hp_amr_smooth_pred_data_t*
curved_hp_amr_smooth_pred_init
(
 p4est_t* p4est,
 double* sigma,
 double* gamma_h,
 double* gamma_p,
 double* gamma_n,
 int num_trees,
 double max_degree,
 norm_t norm_type,
 int* percentiles
);

void
curved_hp_amr_smooth_pred_destroy
(
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data
);

void
curved_hp_amr_smooth_pred_set_refinement_normal
(
 p4est_iter_volume_info_t* info,
 void* user_data
);

void
curved_hp_amr_smooth_pred_set_refinement_percentiles_no_prefine
(
 p4est_iter_volume_info_t* info,
 void* user_data
);


void
curved_hp_amr_smooth_pred_set_refinement_percentiles
(
 p4est_iter_volume_info_t* info,
 void* user_data
);

void
curved_hp_amr_smooth_pred_refine_replace_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
);

void
curved_hp_amr_smooth_pred_balance_replace_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
);

void
curved_hp_amr_smooth_pred_load_predictor
(
 p4est_t* p4est,
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data
);

void
curved_hp_amr_smooth_pred_save_predictor
(
 p4est_t* p4est,
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data
);

void
curved_hp_amr_smooth_pred_zero_predictor
(
 p4est_t* p4est
);

void
curved_hp_amr_smooth_pred_set_local_eta2_avg
(
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data,
 double local_eta2_avg
);


#endif
