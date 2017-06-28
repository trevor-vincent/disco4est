#ifndef CURVED_GAUSS_SIPG_FLUX_HESTHAVEN_H
#define CURVED_GAUSS_SIPG_FLUX_HESTHAVEN_H 

#include "../Flux/curved_compute_flux.h"
#include <ip_flux_params.h>

curved_flux_fcn_ptrs_t curved_gauss_primal_sipg_hesthaven_flux_dirichlet_fetch_fcns(grid_fcn_t bndry_fcn,ip_flux_params_t *curved_gauss_sipg_params);

#endif
