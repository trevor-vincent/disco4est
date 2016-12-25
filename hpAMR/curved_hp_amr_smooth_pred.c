#define NDEBUG
#define SAFETY

#include "../Debug/debug.h"
#include "../pXest/pXest.h"
#include "../ElementData/element_data.h"
#include "../Problem/Blueprints/problem_data.h"
#include "../Problem/Blueprints/problem_weakeqn_ptrs.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../hpAMR/curved_hp_amr.h"
#include "../hpAMR/curved_hp_amr_smooth_pred.h"

#if (P4EST_DIM)==3
#define ONE_OVER_CHILDREN 0.125
#else
#define ONE_OVER_CHILDREN 0.25
#endif

static
void
curved_hp_amr_smooth_pred_init_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;

  /* initialize it such that initial step is h-refinement */
  elem_data->local_predictor = 0.;
}

static
void
curved_hp_amr_smooth_pred_save_predictor_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data =
    (curved_hp_amr_smooth_pred_data_t*) user_data;
  
  smooth_pred_data->pred_storage[smooth_pred_data->pred_storage_stride] = elem_data->local_predictor;
  smooth_pred_data->pred_storage_stride += 1;
}

static
void
curved_hp_amr_smooth_pred_load_predictor_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) user_data;
  elem_data->local_predictor = smooth_pred_data->pred_storage[smooth_pred_data->pred_storage_stride];
  smooth_pred_data->pred_storage_stride += 1;
}

void
curved_hp_amr_smooth_pred_zero_predictor
(
 p4est_t* p4est
)
{
  /* 
     zero out prediction if this is the first init call, 
     otherwise we want to keep it around for the next
     round of refinement. 
  */

  p4est_iterate(p4est,
                NULL,
                NULL,
                curved_hp_amr_smooth_pred_init_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);
      
}

void
curved_hp_amr_smooth_pred_save_predictor
(
 p4est_t* p4est,
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data
)
{
  smooth_pred_data->pred_storage_size = p4est->local_num_quadrants;
  smooth_pred_data->pred_storage_stride = 0;
  smooth_pred_data->pred_storage = P4EST_REALLOC
    (
     smooth_pred_data->pred_storage,
     double,
     smooth_pred_data->pred_storage_size
     );
  
  p4est_iterate
    (
     p4est,
     NULL,
     (void*)smooth_pred_data,
     curved_hp_amr_smooth_pred_save_predictor_callback,
     NULL,
#if (P4EST_DIM)==3
     NULL,       
#endif                
     NULL
    );
}

void
curved_hp_amr_smooth_pred_load_predictor
(
 p4est_t* p4est,
 curved_hp_amr_smooth_pred_data_t* smooth_pred_data
)
{
  /* only load the data if it matches the size of local mesh */
  mpi_assert(smooth_pred_data->pred_storage_size == p4est->local_num_quadrants);
  smooth_pred_data->pred_storage_stride = 0;
  p4est_iterate
    (
     p4est,
     NULL,
     (void*)smooth_pred_data,
     curved_hp_amr_smooth_pred_load_predictor_callback,
     NULL,
#if (P4EST_DIM)==3
     NULL,       
#endif                
     NULL
    );
}

/* void */
/* curved_hp_amr_smooth_pred_set_local_eta2_avg */
/* ( */
/*  curved_hp_amr_smooth_pred_data_t* smooth_pred_data, */
/*  double local_eta2_avg */
/* ) */
/* { */
/*   smooth_pred_data->local_eta2_avg = local_eta2_avg; */
/* } */

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
)
{
#ifdef SAFETY
  mpi_assert( (ONE_OVER_CHILDREN) == (1./(double)(P4EST_CHILDREN)) );
#endif
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data;
  smooth_pred_data = P4EST_ALLOC(curved_hp_amr_smooth_pred_data_t, 1);
  
  smooth_pred_data->sigma = sigma;
  smooth_pred_data->gamma_h = gamma_h;
  smooth_pred_data->gamma_p = gamma_p;
  smooth_pred_data->gamma_n = gamma_n;
  smooth_pred_data->percentiles = percentiles;
  
  smooth_pred_data->max_degree = max_degree;
  smooth_pred_data->norm_type = norm_type;
  smooth_pred_data->num_trees = num_trees;

  smooth_pred_data->pred_storage_size = 0;
  smooth_pred_data->pred_storage = NULL; 
  
  return smooth_pred_data;
}

void
curved_hp_amr_smooth_pred_destroy(curved_hp_amr_smooth_pred_data_t* smooth_pred_data){
  if(smooth_pred_data->pred_storage != NULL)
    P4EST_FREE(smooth_pred_data->pred_storage);
  P4EST_FREE(smooth_pred_data);
}

void
curved_hp_amr_smooth_pred_set_refinement_normal
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) info->p4est->user_pointer;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) (curved_hp_amr_data->curved_hp_amr_scheme_data);
  curved_element_data_t* elem_data = (curved_element_data_t*) info->quad->p.user_data;

  estimator_stats_t* stats = curved_hp_amr_data->estimator_stats;
  
  double eta2 = elem_data->local_estimator;
  double eta2_pred = elem_data->local_predictor;
  double eta2_avg = stats->mean;
  
  int which_tree = info->treeid;
  double gamma_h = smooth_pred_data->gamma_h[which_tree];
  double gamma_p = smooth_pred_data->gamma_p[which_tree];
  double gamma_n = smooth_pred_data->gamma_n[which_tree];
  
  if (eta2 >= smooth_pred_data->sigma[which_tree] * eta2_avg){
    /* it's smooth! -> p-refine*/
    if (eta2 <= elem_data->local_predictor && elem_data->deg < smooth_pred_data->max_degree){
      /* printf("SHOULD NOT HAPPEN\n"); */
      curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = util_min_int(elem_data->deg + 1, smooth_pred_data->max_degree);
      eta2_pred = gamma_p*eta2;
    }

    else {
      curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = -1*elem_data->deg;
      if (smooth_pred_data->norm_type == dg_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg))*(ONE_OVER_CHILDREN);
      else if (smooth_pred_data->norm_type == l2_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg + 1))*(ONE_OVER_CHILDREN);
      else
        mpi_abort("ERROR: curved_hp_amr_smooth_pred_data->set_refinement");
    }
  }

  /* do not refine */
  else {
    eta2_pred = gamma_n*eta2_pred;
    curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = elem_data->deg;
  }

  elem_data->local_predictor = eta2_pred;
  curved_hp_amr_data->refinement_log_stride += 1;
}


void
curved_hp_amr_smooth_pred_set_refinement_percentiles
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) info->p4est->user_pointer;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) (curved_hp_amr_data->curved_hp_amr_scheme_data);
  curved_element_data_t* elem_data = (curved_element_data_t*) info->quad->p.user_data;

  estimator_stats_t* stats = curved_hp_amr_data->estimator_stats;
  
  double eta2 = elem_data->local_estimator;
  double eta2_pred = elem_data->local_predictor;
  /* double eta2_avg = stats->mean; */
  /* double eta2_min = stats->min; */
  
  int which_tree = info->treeid;
  double gamma_h = smooth_pred_data->gamma_h[which_tree];
  double gamma_p = smooth_pred_data->gamma_p[which_tree];
  double gamma_n = smooth_pred_data->gamma_n[which_tree];


  double percentile = -1;
  if (smooth_pred_data->percentiles[which_tree] == 5)
    percentile = stats->p5;
  else if (smooth_pred_data->percentiles[which_tree] == 10)
    percentile = stats->p10;
  else if (smooth_pred_data->percentiles[which_tree] == 15)
    percentile = stats->p15;
  else if (smooth_pred_data->percentiles[which_tree] == 20)
    percentile = stats->p20;
  else if (smooth_pred_data->percentiles[which_tree] == 25)
    percentile = stats->Q3;
  else
    mpi_abort("This percentile is not supported");
  
  if (eta2 >= percentile){
    /* it's smooth! -> p-refine*/
    if (eta2 <= elem_data->local_predictor && elem_data->deg < smooth_pred_data->max_degree){
      /* printf("SHOULD NOT HAPPEN\n"); */
      curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = util_min_int(elem_data->deg + 1, smooth_pred_data->max_degree);
      eta2_pred = gamma_p*eta2;
    }

    else {
      curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = -1*elem_data->deg;
      if (smooth_pred_data->norm_type == dg_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg))*(ONE_OVER_CHILDREN);
      else if (smooth_pred_data->norm_type == l2_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg + 1))*(ONE_OVER_CHILDREN);
      else
        mpi_abort("ERROR: curved_hp_amr_smooth_pred_data->set_refinement");
    }
  }

  /* do not refine */
  else {
    eta2_pred = gamma_n*eta2_pred;
    curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = elem_data->deg;
  }

  elem_data->local_predictor = eta2_pred;
  curved_hp_amr_data->refinement_log_stride += 1;
}

void
curved_hp_amr_smooth_pred_set_refinement_percentiles_no_prefine
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) info->p4est->user_pointer;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) (curved_hp_amr_data->curved_hp_amr_scheme_data);
  curved_element_data_t* elem_data = (curved_element_data_t*) info->quad->p.user_data;

  estimator_stats_t* stats = curved_hp_amr_data->estimator_stats;
  
  double eta2 = elem_data->local_estimator;
  double eta2_pred = elem_data->local_predictor;
  /* double eta2_avg = stats->mean; */
  /* double eta2_min = stats->min; */
  
  int which_tree = info->treeid;
  double gamma_h = smooth_pred_data->gamma_h[which_tree];
  /* double gamma_p = smooth_pred_data->gamma_p[which_tree]; */
  double gamma_n = smooth_pred_data->gamma_n[which_tree];


  double percentile = -1;
  if (smooth_pred_data->percentiles[which_tree] == 5)
    percentile = stats->p5;
  else if (smooth_pred_data->percentiles[which_tree] == 10)
    percentile = stats->p10;
  else if (smooth_pred_data->percentiles[which_tree] == 15)
    percentile = stats->p15;
  else if (smooth_pred_data->percentiles[which_tree] == 20)
    percentile = stats->p20;
  else if (smooth_pred_data->percentiles[which_tree] == 25)
    percentile = stats->Q3;
  else
    mpi_abort("This percentile is not supported");
  
  if (eta2 >= percentile){
      curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = -1*elem_data->deg;
      if (smooth_pred_data->norm_type == dg_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg))*(ONE_OVER_CHILDREN);
      else if (smooth_pred_data->norm_type == l2_norm_type)
        eta2_pred = gamma_h*eta2*util_dbl_pow_int(.5, 2*(elem_data->deg + 1))*(ONE_OVER_CHILDREN);
      else
        mpi_abort("ERROR: curved_hp_amr_smooth_pred_data->set_refinement");
  }

  /* do not refine */
  else {
    eta2_pred = gamma_n*eta2_pred;
    curved_hp_amr_data->refinement_log[curved_hp_amr_data->refinement_log_stride] = elem_data->deg;
  }

  elem_data->local_predictor = eta2_pred;
  curved_hp_amr_data->refinement_log_stride += 1;
}


void
curved_hp_amr_smooth_pred_balance_replace_callback (
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
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = curved_hp_amr_data->dgmath_jit_dbase;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) (curved_hp_amr_data->curved_hp_amr_scheme_data);
  
  curved_element_data_t* parent_data = (curved_element_data_t*) outgoing[0]->p.user_data;
  curved_element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = dgmath_get_nodes((P4EST_DIM), degH);
  double gamma_h = smooth_pred_data->gamma_h[which_tree];
    
  int h_pow = -1;
  if (smooth_pred_data->norm_type == l2_norm_type)
    h_pow = parent_data->deg+1;
  else if (smooth_pred_data->norm_type == dg_norm_type)
    h_pow = parent_data->deg;
  else
    mpi_abort("[ERROR]: curved_hp_amr_smooth_pred norm_type not supported");

  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  dgmath_apply_hp_prolong
    (
     dgmath_jit_dbase,
     &(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (curved_element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    child_data->local_predictor = (ONE_OVER_CHILDREN)*gamma_h*util_dbl_pow_int(.5, 2*(h_pow))*parent_data->local_predictor;
    linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_elem[0], volume_nodes);
  }

  P4EST_FREE(temp_data);
}


void
curved_hp_amr_smooth_pred_refine_replace_callback (
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
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = curved_hp_amr_data->dgmath_jit_dbase;
  curved_hp_amr_smooth_pred_data_t* smooth_pred_data = (curved_hp_amr_smooth_pred_data_t*) (curved_hp_amr_data->curved_hp_amr_scheme_data);
  
  curved_element_data_t* parent_data = (curved_element_data_t*) outgoing[0]->p.user_data;
  curved_element_data_t* child_data;
  int i;

  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;
    
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;

  int volume_nodes = dgmath_get_nodes((P4EST_DIM), degH);
    
  int h_pow = -1;
  if (smooth_pred_data->norm_type == l2_norm_type)
    h_pow = parent_data->deg+1;
  else if (smooth_pred_data->norm_type == dg_norm_type)
    h_pow = parent_data->deg;
  else{
    mpi_abort("[ERROR]: curved_hp_amr_smooth_pred norm_type not supported");
    h_pow *= -1.; /* remove warnings */
  }
  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  dgmath_apply_hp_prolong
    (
     dgmath_jit_dbase,
     &(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (curved_element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    child_data->local_predictor = parent_data->local_predictor;
    linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_elem[0], volume_nodes);
  }

  P4EST_FREE(temp_data);
}
