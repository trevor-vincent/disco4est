#ifndef CURVED_DG_NORM_H
#define CURVED_DG_NORM_H 

#include <d4est_mortars_compute_flux.h>
#include <ip_flux.h>

typedef struct{
  d4est_poisson_flux_sipg_params_t* ip_flux_params;
  double dg_norm_face_term;
}
curved_dg_norm_params_t;

d4est_mortars_fcn_ptrs_t curved_dg_norm_fetch_fcns(curved_dg_norm_params_t *curved_dg_params);


#endif
