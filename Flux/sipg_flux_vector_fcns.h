#ifndef SIPG_FLUX_VECTOR_FCNS_H
#define SIPG_FLUX_VECTOR_FCNS_H 

#include "../GridFunctions/grid_functions.h"
#include "../Flux/compute_flux.h"
#include "../Flux/ip_flux_params.h"

inline static double
sipg_flux_vector_calc_penalty_meanp2_over_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double mean_p2 = .5*( (deg_m)*(deg_m) + (deg_p)*(deg_p) );
  double min_h = (h_m < h_p) ? h_m : h_p;
  return penalty_prefactor*mean_p2/min_h;
}

inline static double
sipg_flux_vector_calc_penalty_maxp2_over_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_deg = (deg_m > deg_p) ? deg_m : deg_p;
  double min_h = (h_m < h_p) ? h_m : h_p;
  return (penalty_prefactor*(max_deg)*(max_deg))/min_h;
}

inline static double
sipg_flux_vector_calc_penalty_maxpp12_over_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_deg = (deg_m > deg_p) ? deg_m : deg_p;
  double min_h = (h_m < h_p) ? h_m : h_p;
  return (penalty_prefactor*(max_deg+1)*(max_deg+1))/min_h;
}

inline static double
sipg_flux_vector_calc_penalty_max_pp12_over_h
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double minus_pp12_over_h = (deg_m + 1)*(deg_m + 1)/h_m;
  double plus_pp12_over_h = (deg_p + 1)*(deg_p + 1)/h_p;
  double max_pp12_over_h = (minus_pp12_over_h > plus_pp12_over_h)
                           ? minus_pp12_over_h : plus_pp12_over_h;    
  return (penalty_prefactor*max_pp12_over_h);
}

flux_fcn_ptrs_t
sipg_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* sipg_params
);

#endif
