#ifndef CURVED_TEST_MORTARJACOBIANTERMS_FLUX_FCNS_H
#define CURVED_TEST_MORTARJACOBIANTERMS_FLUX_FCNS_H 

#include <curved_compute_flux.h>

typedef struct {
  
  double nx_err;
  double ny_err;
  double nz_err;

  double sj_err;
  
  double dudx_err;
  double dudy_err;
  double dudz_err;

} test_mortarjacobianterms_data_t;

curved_flux_fcn_ptrs_t
curved_test_mortarjacobianterms_fetch_fcns
(
 test_mortarjacobianterms_data_t* data
);

#endif
