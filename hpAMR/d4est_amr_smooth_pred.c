#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_linalg.h>
#include <d4est_util.h>
#include <d4est_amr.h>
#include <ini.h>
#include <d4est_amr_smooth_pred.h>
#include <math.h>

#if (P4EST_DIM)==3
#define ONE_OVER_CHILDREN 0.125
#else
#define ONE_OVER_CHILDREN 0.25
#endif

static void
d4est_amr_smooth_pred_pre_refine_callback
(
 p4est_t* p4est,
 void* user
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) user;
  d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);
  
  if (smooth_pred_data->predictor == NULL){
    smooth_pred_data->predictor = P4EST_REALLOC
                                   (
                                    smooth_pred_data->predictor,
                                    double,
                                    p4est->local_num_quadrants
                                   );
    d4est_linalg_fill_vec(smooth_pred_data->predictor, smooth_pred_data->smooth_pred_params->initial_pred, p4est->local_num_quadrants);
    
  }
}

static void
d4est_amr_smooth_pred_compute_post_balance_predictor
(
 p4est_t* p4est,
 d4est_amr_t* d4est_amr,
 d4est_amr_smooth_pred_data_t* smooth_pred_data,
 double** predictors
)
{
  int post_balance_elements = p4est->local_num_quadrants;
  double *aux = D4EST_ALLOC(double, post_balance_elements);

  /* interpolate from initial to auxiliary (unbalanced) grid */
  int pre_refine_stride = 0;
  int post_refine_stride = 0;
  for (int i = 0; i < d4est_amr->initial_log_size; i++){
    int pre_volume_nodes = 1;
    int post_volume_nodes = 1;
        
    if(d4est_amr->refinement_log[i] < 0)
      post_volume_nodes *= (P4EST_CHILDREN);

    for (int j = 0; j < post_volume_nodes; j++){
      aux[post_refine_stride + j] = (*predictors)[pre_refine_stride];
    }

    pre_refine_stride += pre_volume_nodes;
    post_refine_stride += post_volume_nodes;
  }

  *predictors = D4EST_REALLOC(*predictors, double, post_balance_elements);
  
  /* interpolate from auxiliary to balanced grid */
  int pre_balance_stride = 0;
  int post_balance_stride = 0;
  for (int i = 0; i < d4est_amr->balance_log_size; i++){
    int pre_balance_volume_nodes = 1;
    int post_balance_volume_nodes = 1;
    if(d4est_amr->balance_log[i] < 0){
      post_balance_volume_nodes *= (P4EST_CHILDREN);
      for (int j = 0; j < post_balance_volume_nodes; j++){
        gamma_params_t gamma_hpn
          = smooth_pred_data->marker.set_element_gamma_fcn(p4est, NULL, NULL, smooth_pred_data->smooth_pred_params,smooth_pred_data->marker.user);
        int h_pow = abs(d4est_amr->balance_log[i]);
        (*predictors)[post_balance_stride + j] = (ONE_OVER_CHILDREN)*gamma_hpn.gamma_h*d4est_util_dbl_pow_int(.5, 2*(h_pow))*aux[pre_balance_stride];        
      }
    }
    else {
      /* only one element */
      (*predictors)[post_balance_stride] = aux[pre_balance_stride];
    }
    pre_balance_stride += pre_balance_volume_nodes;
    post_balance_stride += post_balance_volume_nodes;
  }
    
  D4EST_FREE(aux);
}
 



static void
d4est_amr_smooth_pred_post_balance_callback
(
 p4est_t* p4est,
 void* user
)
{

  d4est_amr_t* d4est_amr = (d4est_amr_t*) user;
  d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);

  d4est_amr_smooth_pred_compute_post_balance_predictor
    (
     p4est,
     d4est_amr,
     smooth_pred_data,
     &smooth_pred_data->predictor
    );
}

static void
d4est_amr_smooth_pred_mark_elements
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  d4est_estimator_stats_t* stats = d4est_amr->d4est_estimator_stats;
  
  double eta2 = d4est_amr->d4est_estimator[elem_data->id];
  double eta2_pred = smooth_pred_data->predictor[elem_data->id];//elem_data->local_predictor;

  gamma_params_t gamma_hpn =
    smooth_pred_data->marker.set_element_gamma_fcn
    (
     info->p4est,
     stats,
     elem_data,
     smooth_pred_data->smooth_pred_params,
     smooth_pred_data->marker.user
    );
  
  int is_marked = 
    smooth_pred_data->marker.mark_element_fcn
    (
     info->p4est,
     eta2,
     stats,
     elem_data,
     smooth_pred_data->smooth_pred_params,
     smooth_pred_data->marker.user
    );

  int max_degree = d4est_amr->max_degree;

  if (is_marked){
    if (eta2 <= eta2_pred && elem_data->deg < max_degree){
      d4est_amr->refinement_log[elem_data->id] = d4est_util_min_int(elem_data->deg + 1, max_degree);
      eta2_pred = gamma_hpn.gamma_p*eta2;
    }
    else {
      d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
      eta2_pred = gamma_hpn.gamma_h*eta2*d4est_util_dbl_pow_int(.5, 2*(elem_data->deg))*(ONE_OVER_CHILDREN);
    }
  }
  else {
    eta2_pred = gamma_hpn.gamma_n*eta2_pred;
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
  }  
  smooth_pred_data->predictor[elem_data->id] = eta2_pred;
}

void
d4est_amr_smooth_pred_destroy(d4est_amr_scheme_t* scheme){

  d4est_amr_smooth_pred_data_t* smooth_pred_data =
    (d4est_amr_smooth_pred_data_t*)scheme->amr_scheme_data;  
  P4EST_FREE(smooth_pred_data->predictor);
  P4EST_FREE(smooth_pred_data->smooth_pred_params);
  P4EST_FREE(smooth_pred_data);
  P4EST_FREE(scheme);
}


static
int d4est_amr_smooth_pred_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_amr_smooth_pred_params_t* pconfig = (d4est_amr_smooth_pred_params_t*)user;
   if (d4est_util_match_couple(section,"amr",name,"sigma")) {
    D4EST_ASSERT(pconfig->sigma == -1);
    pconfig->sigma = atof(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"gamma_h")) {
    D4EST_ASSERT(pconfig->gamma_h == -1);
    pconfig->gamma_h = atof(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"initial_pred")) {
    D4EST_ASSERT(pconfig->initial_pred == 0);
    pconfig->initial_pred = atof(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"gamma_p")) {
    D4EST_ASSERT(pconfig->gamma_p == -1);
    pconfig->gamma_p = atof(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"gamma_n")) {
    D4EST_ASSERT(pconfig->gamma_n == -1);
    pconfig->gamma_n = atof(value);
  }   
  else if (d4est_util_match_couple(section,"amr",name,"percentile")) {
    D4EST_ASSERT(pconfig->percentile == -1);
    pconfig->percentile = atoi(value);
  }
  else if (d4est_util_match_couple(section,"amr",name,"inflation_size")) {
    D4EST_ASSERT(pconfig->inflation_size == -1);
    pconfig->inflation_size = atoi(value);
  }
  else {
    return 0; 
  }
  return 1;
}

d4est_amr_smooth_pred_params_t*
d4est_amr_smooth_pred_params_input
(
 const char* input_file
)
{
  d4est_amr_smooth_pred_params_t* input = D4EST_ALLOC(d4est_amr_smooth_pred_params_t, 1);
  input->initial_pred;
  input->gamma_h = -1;
  input->gamma_n = -1;
  input->gamma_p = -1;
  input->sigma = -1;
  input->percentile = -1;
  input->inflation_size = -1;
  input->initial_pred = 0;
  
  if (ini_parse(input_file, d4est_amr_smooth_pred_params_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("amr", input->gamma_h, -1);
  D4EST_CHECK_INPUT("amr", input->gamma_p, -1);
  D4EST_CHECK_INPUT("amr", input->gamma_n, -1);
  D4EST_ASSERT((input->sigma != -1) || (input->percentile != -1));
  
  return input;
}



d4est_amr_scheme_t*
d4est_amr_smooth_pred_init
(
 p4est_t* p4est,
 const char* input_file,
 d4est_amr_scheme_t* scheme,
 void* marker
)
{  
  d4est_amr_smooth_pred_data_t* smooth_pred_data;
  smooth_pred_data = P4EST_ALLOC(d4est_amr_smooth_pred_data_t, 1);
  smooth_pred_data->marker = *((d4est_amr_smooth_pred_marker_t*)marker);
  smooth_pred_data->predictor = NULL;

  d4est_amr_smooth_pred_params_t* smooth_pred_params = d4est_amr_smooth_pred_params_input(input_file);

  smooth_pred_data->smooth_pred_params = smooth_pred_params;
  
  scheme->pre_refine_callback
    = d4est_amr_smooth_pred_pre_refine_callback;
  
  scheme->balance_replace_callback_fcn_ptr
    = NULL;

  scheme->refine_replace_callback_fcn_ptr
    = NULL;

  scheme->mark_elements
    = d4est_amr_smooth_pred_mark_elements;

  scheme->amr_scheme_data
    = smooth_pred_data;

  scheme->post_balance_callback
    = d4est_amr_smooth_pred_post_balance_callback;

  scheme->destroy
    = d4est_amr_smooth_pred_destroy;
  
  return scheme;
}
