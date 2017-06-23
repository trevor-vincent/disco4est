#ifndef ARBQUAD_H
#define ARBQUAD_H 

#include <pXest.h>
#include <math.h>
#include <util.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

/**
 * @file   arbquad.h
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Sun May 14 15:00:37 2017
 * 
 * @brief  Find the abscissas and weights to at least
 * double precision for arbitrary weight functions and intervals.
 * This is based off the paper
 * 
 * { fukuda2005gaussian, 
 *   title={Gaussian quadrature rule for arbitrary weight function and interval}, 
 *   author={Fukuda, H and Katuya, M and Alt, EO and Matveenko, AV}, 
 *   journal={Computer physics communications}, 
 *   volume={167}, 
 *   number={2}, 
 *   pages={143--150}, 
 *   year={2005}, 
 *   publisher={Elsevier} 
 * } 
 *
 * 1) The paper above provides a mathematica code to compute the 
 * abscissas and weights to high precision (much higher than double).
 * We created a C-version of this code so that it could be embedded 
 * in a C-library and the abscissas and weights could be computed 
 * on demand (for example, if the weighting function has a parameter
 * which isn't known until runtime).
 *
 *
 * 2) We use long double precision in the calculations 
 * to obtain at least double precision results. We do not obtain
 * precision as high as the mathematica code, to do this we would need
 * something like the GMP library, but the precision we get should be more than
 * enough for double precision codes. The arbquad_test_precision function tests
 * if the moments can be calculated to at least DBL_EPSILON with the computed
 * quadrature values.
 *
 */


typedef enum {DO_NOT_DIVIDE_WEIGHTS_BY_WEIGHT_FCN, DIVIDE_WEIGHTS_BY_WEIGHT_FCN} arbquad_weight_choice_t;

typedef long double
(*arbquad_moment_fcn_t)
(
 int,
 void*
);

typedef long double
(*arbquad_weight_fcn_t)
(
 long double,
 void*
);

typedef void
(*arbquad_aa_and_bb_fcn_t)
(
 int n,
 long double*,
 long double*,
 void*
);

/* This file was automatically generated.  Do not edit! */
void arbquad_get_abscissas_and_weights_use_aa_and_bb(int n,long double *weights,long double *abscissas,arbquad_moment_fcn_t moment_fcn,arbquad_aa_and_bb_fcn_t aa_and_bb_fcn,void *user,arbquad_weight_choice_t weight_choice,arbquad_weight_fcn_t weight_fcn);
void arbquad_get_abscissas_and_weights(int n,long double *weights,long double *abscissas,arbquad_moment_fcn_t moment_fcn,void *user,arbquad_weight_choice_t weight_choice,arbquad_weight_fcn_t weight_fcn);
void arbquad_get_jacobi_coefficients(int n,long double *moments,long double *aa,long double *bb);
int arbquad_ql_implicit(int n,long double *d,long double *e,long double *z);
long double arbquad_pythagl(const long double a,const long double b);
long double arbquad_signl(const long double a,const long double b);
void arbquad_test_moments(int n,long double *weights,long double *abscissas,arbquad_moment_fcn_t moment_fcn,arbquad_weight_choice_t weight_choice,arbquad_weight_fcn_t weight_fcn,void *user);
void arbquad_print_float_info();

#endif
