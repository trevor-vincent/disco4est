#ifndef IP_FLUX_PARAMS_H
#define IP_FLUX_PARAMS_H

#include <pXest.h>

typedef double
(*penalty_calc_t)
(
 int, //minus side degree
 double, //minus element size h
 int, //plus side degree
 double, //plus element size h
 double //penalty prefactor
);

typedef enum { BC_EVAL_ON_GAUSS_POINTS, BC_EVAL_ON_LOBATTO_POINTS, BC_EVAL_NOTSET} bc_eval_t;
typedef enum { H_EQ_J_DIV_SJ, H_EQ_VOLUME_DIV_AREA, H_EQ_NOTSET } h_calc_method_t;

typedef struct {
  double ip_flux_penalty_prefactor;

  /* Given h and the order of the polynomial p, for the +,- sides sharing an interface
   * compute the penalty function */
  penalty_calc_t ip_flux_penalty_calculate_fcn;

  /* For curved elements, h is not well defined, therefore we try two methods  */
  /* 1. H_EQ_J_DIV_SJ ----->  h(+,-) ~ jacobian(+,-)/surface_jacobian */
  /* 2. H_EQ_VOLUME_DIV_AREA ----> h(+,-) ~ volume_of_element(+,-)/surface_area_of_face */
  /* here (+,-) refers to the plus or minus side of a shared interface */
  h_calc_method_t ip_flux_h_calc;

  /* Evaluate boundary function on */
  /* 1. BC_EVAL_ON_GAUSS_POINTS -----> Gauss integration points */
  /* 2. BC_EVAL_ON_LOBATTO_POINTS -----> Lobatto points and then interpolate to Gauss points */
  bc_eval_t ip_flux_bc_eval;

  char name [50];

} ip_flux_params_t;


/* This file was automatically generated.  Do not edit! */
void ip_flux_params_destroy(ip_flux_params_t *params);
ip_flux_params_t *ip_flux_params_new(p4est_t *p4est,const char *print_prefix,const char *input_file);

#endif
