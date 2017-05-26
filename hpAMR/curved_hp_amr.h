#ifndef CURVED_HP_AMR_H
#define CURVED_HP_AMR_H 

#include "../pXest/pXest.h"
#include "../dGMath/d4est_operators.h"
#include "../Estimators/estimator_stats.h"

typedef struct
{  
  double* data;
  int* refinement_log;
  int refinement_log_stride;
  void* curved_hp_amr_scheme_data;
  d4est_operators_t* d4est_ops;  
  p4est_replace_t refine_replace_callback_fcn_ptr;
  p4est_replace_t balance_replace_callback_fcn_ptr;
  estimator_stats_t* estimator_stats;
  
} curved_hp_amr_data_t;

void
curved_hp_amr
(
 p4est_t* p4est,
 double** data_to_hp_refine,
 p4est_iter_volume_t iter_volume,
 p4est_replace_t refine_replace_callback_fcn_ptr, 
 p4est_replace_t balance_replace_callback_fcn_ptr, 
 void* curved_hp_amr_scheme_data,
 estimator_stats_t* stats,
 d4est_operators_t* d4est_ops
);

void
curved_hp_amr_store_balance_changes
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
);

void
curved_hp_amr_set_deg_on_refine
(
 p4est_t* p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t* outgoing[],
 int num_incoming,
 p4est_quadrant_t* incoming[]
);

/* void */
/* curved_hp_amr_save_to_vtk */
/* ( */
/*  p4est_t* p4est, */
/*  p4est_geometry_t* geom, */
/*  char* filename, */
/*  int save_local_est */
/* ); */

#endif
