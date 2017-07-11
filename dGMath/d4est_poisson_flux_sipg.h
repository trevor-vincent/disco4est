#ifndef D4EST_POISSON_SIPG_FLUX_H
#define D4EST_POISSON_SIPG_FLUX_H 

#include <pXest.h>
#include <d4est_mortars.h>
#include <d4est_poisson_flux.h>

typedef double
(*penalty_calc_t)
(
 int, //minus side degree
 double, //minus element size h
 int, //plus side degree
 double, //plus element size h
 double //penalty prefactor
);

typedef enum { H_EQ_J_DIV_SJ, H_EQ_J_DIV_SJ_MIN, H_EQ_VOLUME_DIV_AREA, H_EQ_NOTSET } h_calc_method_t;

typedef struct {
  double sipg_penalty_prefactor;

  /* Given h and the order of the polynomial p, for the +,- sides sharing an interface
   * compute the penalty function */
  penalty_calc_t sipg_penalty_fcn;

  /* For curved elements, h is not well defined, therefore we try two methods  */
  /* 1. H_EQ_J_DIV_SJ ----->  h(+,-) ~ jacobian(+,-)/surface_jacobian */
  /* 2. H_EQ_J_DIV_SJ_MIN */
  
  h_calc_method_t sipg_flux_h;

} d4est_poisson_flux_sipg_params_t;

/* This file was automatically generated.  Do not edit! */
void d4est_poisson_flux_sipg_params_destroy(d4est_poisson_flux_data_t *data);
void d4est_poisson_flux_sipg_params_new(p4est_t *p4est,d4est_xyz_fcn_t boundary_condition,const char *print_prefix,const char *input_file,d4est_poisson_flux_data_t *d4est_poisson_flux_data);

#endif
