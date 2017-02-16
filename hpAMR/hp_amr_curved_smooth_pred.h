#ifndef HP_AMR_CURVED_SMOOTH_PRED_H
#define HP_AMR_CURVED_SMOOTH_PRED_H 

#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"
#include "../hpAMR/hp_amr.h"
#include "../ElementData/element_data.h"


typedef struct {

  double gamma_h;
  double gamma_p;
  double gamma_n;

} gamma_params_t;

typedef struct {


  int (*mark_element_fcn)
  (
   double, /* eta2 */
   estimator_stats_t*,
   curved_element_data_t*,
   void* /* user ptr */
  );
  
  gamma_params_t (*set_element_gamma_fcn)
  (
   curved_element_data_t*,
   void* /* user ptr */
  );
  
  const char* name;
  void* user;
  
} curved_smooth_pred_marker_t;

typedef struct {
  
  int max_degree;
  curved_smooth_pred_marker_t marker;
  double* predictors; 
  
} hp_amr_curved_smooth_pred_data_t;



void hp_amr_curved_smooth_pred_set_marker(hp_amr_scheme_t *scheme, curved_smooth_pred_marker_t marker);
void hp_amr_curved_smooth_pred_destroy(hp_amr_scheme_t *scheme);
void hp_amr_curved_smooth_pred_post_balance_callback(p4est_t *p4est,void *user);
void hp_amr_curved_smooth_pred_set_refinement(p4est_iter_volume_info_t *info,void *user_data);
void hp_amr_curved_smooth_pred_refine_replace_callback(p4est_t *p4est,p4est_topidx_t which_tree,int num_outgoing,p4est_quadrant_t *outgoing[],int num_incoming,p4est_quadrant_t *incoming[]);
void hp_amr_curved_smooth_pred_balance_replace_callback(p4est_t *p4est,p4est_topidx_t which_tree,int num_outgoing,p4est_quadrant_t *outgoing[],int num_incoming,p4est_quadrant_t *incoming[]);
void hp_amr_curved_smooth_pred_pre_refine_callback(p4est_t *p4est,void *user);
void hp_amr_curved_smooth_pred_print(p4est_t *p4est);
hp_amr_scheme_t* hp_amr_curved_smooth_pred_init
(
 p4est_t* p4est,
 int max_degree,
 curved_smooth_pred_marker_t marker
);


#endif
