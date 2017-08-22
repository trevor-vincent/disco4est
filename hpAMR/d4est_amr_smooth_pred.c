#include <pXest.h>
#include <d4est_element_data.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_linalg.h>
#include <d4est_util.h>
#include <d4est_amr.h>
#include <d4est_amr_smooth_pred.h>

#if (P4EST_DIM)==3
#define ONE_OVER_CHILDREN 0.125
#else
#define ONE_OVER_CHILDREN 0.25
#endif

void
d4est_amr_smooth_pred_print
(
 p4est_t* p4est
)
{
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
        printf("Element %d predictor = %.25f\n", ed->id, ed->local_predictor);
      }
    }
}



static void
d4est_amr_smooth_pred_pre_refine_callback
(
 p4est_t* p4est,
 void* user
)
{
  d4est_amr_smooth_pred_data_t* smooth_pred_data =
    (d4est_amr_smooth_pred_data_t*)user;

  if (smooth_pred_data->predictors == NULL){
    smooth_pred_data->predictors = P4EST_REALLOC
                                   (
                                    smooth_pred_data->predictors,
                                    double,
                                    p4est->local_num_quadrants
                                   );
    d4est_linalg_fill_vec(smooth_pred_data->predictors, 0., p4est->local_num_quadrants);
  }
  
  int pred_stride = 0;
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
        ed->local_predictor = smooth_pred_data->predictors[pred_stride];
        pred_stride++;
      }
    }
}

static void
d4est_amr_smooth_pred_post_balance_callback
(
 p4est_t* p4est,
 void* user
)
{
  d4est_amr_smooth_pred_data_t* smooth_pred_data =
    (d4est_amr_smooth_pred_data_t*)user;
  smooth_pred_data->predictors = P4EST_REALLOC
                                 (
                                  smooth_pred_data->predictors,
                                  double,
                                  p4est->local_num_quadrants
                                 );  
  int pred_stride = 0;
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
        smooth_pred_data->predictors[pred_stride] = ed->local_predictor;
        pred_stride++;
      }
    }
 
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
  d4est_estimator_stats_t** stats = d4est_amr->d4est_estimator_stats;
  
  double eta2 = elem_data->local_estimator;
  double eta2_pred = elem_data->local_predictor;

  gamma_params_t gamma_hpn =
    smooth_pred_data->marker.set_element_gamma_fcn
    (
     info->p4est,
     eta2,
     stats,
     elem_data,
     smooth_pred_data->marker
    );
  
  int is_marked = 
    smooth_pred_data->marker.mark_element_fcn
    (
     info->p4est,
     eta2,
     stats,
     elem_data,
     smooth_pred_data->marker
    );

  int max_degree = d4est_amr->max_degree;

  /* printf("\n **ELEMENT %d\n", elem_data->id); */
  /* printf("local_predictor = %.15f\n", elem_data->local_predictor); */
  /* printf("eta2 = %.15f\n", elem_data->local_estimator); */
  /* printf("max_degree = %d\n", max_degree); */
  
  if (is_marked){
    if (eta2 <= elem_data->local_predictor && elem_data->deg < max_degree){
      /* printf("ELEMENT %d P-REFINED\n", elem_data->id); */
      d4est_amr->refinement_log[elem_data->id] = d4est_util_min_int(elem_data->deg + 1, max_degree);
      eta2_pred = gamma_hpn.gamma_p*eta2;
    }
    else {
      /* printf("ELEMENT %d H-REFINED\n", elem_data->id); */
      d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
      eta2_pred = gamma_hpn.gamma_h*eta2*d4est_util_dbl_pow_int(.5, 2*(elem_data->deg))*(ONE_OVER_CHILDREN);
    }
  }
  else {
    eta2_pred = gamma_hpn.gamma_n*eta2_pred;
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
  }
  
  elem_data->local_predictor = eta2_pred;
}

static void
d4est_amr_smooth_pred_balance_replace_callback (
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
  D4EST_ASSERT(num_outgoing == 1);
  d4est_amr_t* d4est_amr = (d4est_amr_t*) p4est->user_pointer;
  d4est_operators_t* d4est_ops = d4est_amr->d4est_ops;
  d4est_amr_smooth_pred_data_t* smooth_pred_data = (d4est_amr_smooth_pred_data_t*) (d4est_amr->scheme->amr_scheme_data);
  d4est_estimator_stats_t** stats = d4est_amr->d4est_estimator_stats;
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;

  gamma_params_t gamma_hpn
    = smooth_pred_data->marker.set_element_gamma_fcn(p4est,parent_data->local_estimator, stats, parent_data, smooth_pred_data->marker.user);
  
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), degH);  
  int h_pow = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;   
    child_data->local_predictor = (ONE_OVER_CHILDREN)*gamma_hpn.gamma_h*d4est_util_dbl_pow_int(.5, 2*(h_pow))*parent_data->local_predictor;
  }

}


static void
d4est_amr_smooth_pred_refine_replace_callback (
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
  D4EST_ASSERT(num_outgoing == 1);
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;

  for (int i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    child_data->local_predictor = parent_data->local_predictor;
  }
}


void
d4est_amr_smooth_pred_destroy(d4est_amr_scheme_t* scheme){

  d4est_amr_smooth_pred_data_t* smooth_pred_data =
    (d4est_amr_smooth_pred_data_t*)scheme->amr_scheme_data;  
  P4EST_FREE(smooth_pred_data->predictors);
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
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_amr_smooth_pred_params_t
d4est_amr_smooth_pred_params_input
(
 const char* input_file
)
{
  d4est_amr_smooth_pred_params_t input;
  input.gamma_h = -1;
  input.gamma_n = -1;
  input.gamma_p = -1;
  input.sigma = -1;
  input.percentile = -1;
  input.inflation_size = -1;
  
  if (ini_parse(input_file, d4est_amr_smooth_pred_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("amr", input.gamma_h, -1);
  D4EST_CHECK_INPUT("amr", input.gamma_p, -1);
  D4EST_CHECK_INPUT("amr", input.gamma_n, -1);
  D4EST_ASSERT((input.sigma != -1) || (input.percentile != -1));
  
  return input;
}



d4est_amr_scheme_t*
d4est_amr_smooth_pred_init
(
 p4est_t* p4est,
 const char* input_file,
 d4est_amr_scheme_t* scheme,
 void* marker //smooth_pred_marker_t marker
)
{  

  d4est_amr_smooth_pred_data_t* smooth_pred_data;
  smooth_pred_data = P4EST_ALLOC(d4est_amr_smooth_pred_data_t, 1);
  smooth_pred_data->marker = *((smooth_pred_marker_t*)marker);
  smooth_pred_data->predictors = NULL;

  scheme->pre_refine_callback
    = d4est_amr_smooth_pred_pre_refine_callback;
  
  scheme->balance_replace_callback_fcn_ptr
    = d4est_amr_smooth_pred_balance_replace_callback;

  scheme->refine_replace_callback_fcn_ptr
    = d4est_amr_smooth_pred_refine_replace_callback;

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
