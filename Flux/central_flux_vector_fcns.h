#ifndef CENTRAL_FLUX_VECTOR_FCNS_H
#define CENTRAL_FLUX_VECTOR_FCNS_H 

#include "../GridFunctions/grid_functions.h"
#include "../Flux/compute_flux.h"
#include "../Flux/central_flux_params.h"

flux_fcn_ptrs_t
central_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 central_flux_params_t* central_params
);

#endif
