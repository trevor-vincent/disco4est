#ifndef CURVED_GAUSS_CENTRAL_FLUX_VECTOR_FCNS_H
#define CURVED_GAUSS_CENTRAL_FLUX_VECTOR_FCNS_H 

#include "../Flux/curved_compute_flux.h"
#include "../Flux/central_flux_params.h"

curved_flux_fcn_ptrs_t
curved_gauss_central_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 central_flux_params_t* curved_gauss_central_params
);

#endif
