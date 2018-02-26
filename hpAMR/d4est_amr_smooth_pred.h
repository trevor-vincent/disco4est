#ifndef D4EST_AMR_SMOOTH_PRED_H
#define D4EST_AMR_SMOOTH_PRED_H 

#include <pXest.h>
#include <d4est_xyz_functions.h>
#include <d4est_amr.h>
#include <d4est_element_data.h>


typedef struct {

  double gamma_h;
  double gamma_p;
  double gamma_n;

} gamma_params_t;

typedef struct {

  int (*mark_element_fcn)
  (
   p4est_t*,
   double, /* eta2 */
   d4est_estimator_stats_t**,
   d4est_element_data_t*,
   void* /* user ptr */
  );
  
  gamma_params_t (*set_element_gamma_fcn)
  (
   p4est_t*,
   d4est_estimator_stats_t**,
   d4est_element_data_t*,
   void* /* user ptr */
  );
  
  const char* name;
  void* user;
  
} d4est_amr_smooth_pred_marker_t;

typedef struct {
  
  double gamma_h;
  double gamma_p;
  double gamma_n;
  double sigma;
  int percentile;
  int inflation_size;
  
} d4est_amr_smooth_pred_params_t;

typedef struct {
  
  d4est_amr_smooth_pred_marker_t marker;
  double* predictor;
  
} d4est_amr_smooth_pred_data_t;

/* This file was automatically generated.  Do not edit! */
d4est_amr_scheme_t *d4est_amr_smooth_pred_init(p4est_t *p4est,const char *input_file,d4est_amr_scheme_t *scheme,void *marker);
d4est_amr_smooth_pred_params_t d4est_amr_smooth_pred_params_input(const char *input_file);
void d4est_amr_smooth_pred_destroy(d4est_amr_scheme_t *scheme);
void d4est_amr_smooth_pred_print(p4est_t *p4est);



#endif
