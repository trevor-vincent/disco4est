#include <pXest.h>
#include <d4est_amr.h>
#include <d4est_amr_uniform.h>

d4est_amr_scheme_t*
d4est_amr_uniform_init
(
 p4est_t* p4est,
 const char* input_file,
 d4est_amr_scheme_t* scheme,
 void* user 
)
{  
  scheme->pre_refine_callback
    = NULL;
  
  scheme->balance_replace_callback_fcn_ptr
    = NULL;

  scheme->refine_replace_callback_fcn_ptr
    = NULL;

  scheme->amr_scheme_data
    = NULL;

  scheme->post_balance_callback
    = NULL;

  scheme->mark_elements
    = d4est_amr_uniform_mark_elements;
  
  scheme->destroy
    = d4est_amr_smooth_pred_destroy;

  D4EST_ASSERT(scheme->amr_scheme_type == AMR_UNIFORM_H ||
               scheme->amr_scheme_type == AMR_UNIFORM_P);
  
  return scheme;
}

void
d4est_amr_uniform_destroy(d4est_amr_scheme_t* scheme)
{
  P4EST_FREE(scheme);
}

static void
d4est_amr_uniform_mark_elements
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  
  if (d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_H){
    d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
  }
  else if (d4est_amr->scheme->amr_scheme_type == AMR_UNIFORM_P){
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg + 1;
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: Using uniform amr, but amr_scheme_type != AMR_UNIFORM_H OR AMR_UNIFORM_P");
  }
  
}
