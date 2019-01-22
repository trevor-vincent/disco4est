#include <pXest.h>
#include <d4est_amr.h>
#include <d4est_amr_random.h>
#include <d4est_element_data.h>
#include <stdlib.h>

void
d4est_amr_random_mark_elements
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  double x = (double) rand() / (double) RAND_MAX;

  if (d4est_amr->scheme->amr_scheme_type == AMR_RANDOM_HP){
    if (x < 0.3){
      d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
    }
    else if (x < 0.6){
      d4est_amr->refinement_log[elem_data->id] = elem_data->deg + 1;
      /* printf("elem_data->deg+1 = %d\n", elem_data->deg + 1); */
    }
    else {
      d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
    }
  }
  else if (d4est_amr->scheme->amr_scheme_type == AMR_RANDOM_H){
    if (x < 0.5){
      d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
    }
    else {
      d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
    }
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: must be AMR_RANDOM_HP or AMR_RANDOM_H");
  }
}


d4est_amr_scheme_t*
d4est_amr_random_init
(
 p4est_t* p4est,
 const char* input_file,
 d4est_amr_scheme_t* scheme,
 void* user
)
{
  /* TODO: Put in input file */
  srand(14234331);
  
  scheme->pre_refine_callback
    = NULL;
  
  scheme->balance_replace_callback_fcn_ptr
    = NULL;

  scheme->refine_replace_callback_fcn_ptr
    = NULL;

  scheme->amr_scheme_data
    = NULL;

  scheme->post_h_balance_callback
    = NULL;

  scheme->post_p_balance_callback
    = NULL;
  
  scheme->mark_elements
    = d4est_amr_random_mark_elements;
  
  scheme->destroy
    = d4est_amr_random_destroy;

  D4EST_ASSERT(scheme->amr_scheme_type == AMR_RANDOM_H ||
               scheme->amr_scheme_type == AMR_RANDOM_HP);
  
  return scheme;
}

void
d4est_amr_random_destroy(d4est_amr_scheme_t* scheme)
{
  P4EST_FREE(scheme);
}
