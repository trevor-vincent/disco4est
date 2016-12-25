#ifndef HP_AMR_SMOOTH_PRED_H
#define HP_AMR_SMOOTH_PRED_H 

#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"
#include "../hpAMR/hp_amr.h"
#include "../ElementData/element_data.h"

typedef struct {

  void* user;
  int (*mark_element_fcn)(double, estimator_stats_t*, void*);
  const char* type;
  
} smooth_pred_marker_t;

typedef struct {

  double gamma_h;
  double gamma_p;
  double gamma_n;

  int max_degree;

  smooth_pred_marker_t marker;

  double* predictors; 
  
} hp_amr_smooth_pred_data_t;
/* This file was automatically generated.  Do not edit! */
smooth_pred_marker_t hp_amr_smooth_pred_get_sigaverage_marker(double *sigma);
smooth_pred_marker_t hp_amr_smooth_pred_get_NULL_marker();
smooth_pred_marker_t hp_amr_smooth_pred_get_percentile_marker(int *percentile);
void hp_amr_smooth_pred_set_gammah(hp_amr_scheme_t *scheme,double gamma_h);
void hp_amr_smooth_pred_set_marker(hp_amr_scheme_t *scheme,smooth_pred_marker_t marker);
void hp_amr_smooth_pred_destroy(hp_amr_scheme_t *scheme);
void hp_amr_smooth_pred_post_balance_callback(p4est_t *p4est,void *user);
void hp_amr_smooth_pred_set_refinement(p4est_iter_volume_info_t *info,void *user_data);
void hp_amr_smooth_pred_refine_replace_callback(p4est_t *p4est,p4est_topidx_t which_tree,int num_outgoing,p4est_quadrant_t *outgoing[],int num_incoming,p4est_quadrant_t *incoming[]);
void hp_amr_smooth_pred_balance_replace_callback(p4est_t *p4est,p4est_topidx_t which_tree,int num_outgoing,p4est_quadrant_t *outgoing[],int num_incoming,p4est_quadrant_t *incoming[]);
void hp_amr_smooth_pred_pre_refine_callback(p4est_t *p4est,void *user);
hp_amr_scheme_t *hp_amr_smooth_pred_init(p4est_t *p4est,double gamma_h,double gamma_p,double gamma_n,int max_degree,smooth_pred_marker_t marker);
void hp_amr_smooth_pred_print(p4est_t *p4est);

#endif
