#include "../hpAMR/hp_amr.h"
#include "../hpAMR/hp_amr_uniform.h"
#include "../ElementData/element_data.h"

hp_amr_scheme_t*
hp_amr_uniform
()
{
  hp_amr_scheme_t* scheme = P4EST_ALLOC(hp_amr_scheme_t, 1);
  scheme->iter_volume = hp_amr_uniform_set_refinement;
  scheme->balance_replace_callback_fcn_ptr = NULL;
  scheme->refine_replace_callback_fcn_ptr = NULL;
  scheme->post_balance_callback = NULL;
  scheme->pre_refine_callback = NULL;
  scheme->hp_amr_scheme_data = NULL;
  return scheme;
}

void
hp_amr_uniform_set_refinement
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  hp_amr_data_t* amr_data = (hp_amr_data_t*) info->p4est->user_pointer;
  element_data_t* elem_data = (element_data_t*) info->quad->p.user_data;
  amr_data->refinement_log[info->quadid] = -elem_data->deg;
}

void
hp_amr_uniform_destroy
(
 hp_amr_scheme_t* scheme
)
{
  P4EST_FREE(scheme);
}
