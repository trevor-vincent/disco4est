#ifndef CURVED_GAUSS_SIPG_FLUX_VECTOR_FCNS_H
#define CURVED_GAUSS_SIPG_FLUX_VECTOR_FCNS_H 

#include <curved_compute_flux.h>
#include <ip_flux_aux.h>

curved_flux_fcn_ptrs_t
curved_Gauss_sipg_flux_vector_dirichlet_fetch_fcns(
                                                   grid_fcn_t bndry_fcn,
                                                   ip_flux_params_t *curved_Gauss_sipg_params
                                                  );

#endif
