#ifndef CURVED_BI_ESTIMATOR_FLUX_FCNS_H
#define CURVED_BI_ESTIMATOR_FLUX_FCNS_H

#include "../Flux/curved_compute_flux.h"
#include "../Flux/ip_flux_aux.h"
#include "../Utilities/util.h"

curved_flux_fcn_ptrs_t
curved_bi_est_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 double penalty_prefactor
);



#endif
