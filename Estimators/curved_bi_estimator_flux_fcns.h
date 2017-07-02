#ifndef CURVED_BI_ESTIMATOR_FLUX_FCNS_H
#define CURVED_BI_ESTIMATOR_FLUX_FCNS_H

#include "../Flux/d4est_mortar_compute_flux.h"
#include "../Flux/ip_flux_params.h"
#include "../Utilities/util.h"

d4est_mortar_fcn_ptrs_t
curved_bi_est_dirichlet_fetch_fcns
(
 d4est_grid_fcn_t bndry_fcn,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 double penalty_prefactor
);



#endif
