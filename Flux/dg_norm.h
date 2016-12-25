#ifndef DG_NORM_H
#define DG_NORM_H 

flux_fcn_ptrs_t
dg_norm_ip_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* ip_params
);

#endif
