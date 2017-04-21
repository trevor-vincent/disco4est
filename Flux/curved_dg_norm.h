#ifndef CURVED_DG_NORM_H
#define CURVED_DG_NORM_H 

#include <curved_compute_flux.h>
#include <ip_flux.h>

typedef struct{
  ip_flux_params_t* ip_flux_params;
  double dg_norm_face_term;
}
curved_dg_norm_params_t;

curved_flux_fcn_ptrs_t curved_dg_norm_fetch_fcns(curved_dg_norm_params_t *curved_dg_params);


#endif
