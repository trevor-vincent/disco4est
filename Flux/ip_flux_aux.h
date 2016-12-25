#ifndef IP_FLUX_AUX_H
#define IP_FLUX_AUX_H

typedef double
(*penalty_calc_t)
(
 int, //minus side degree
 double, //minus element size h
 int, //plus side degree
 double, //plus element size h
 double //penalty prefactor
);

typedef struct {
  double ip_flux_penalty_prefactor;
  penalty_calc_t ip_flux_penalty_calculate_fcn;
} ip_flux_params_t;

#endif
