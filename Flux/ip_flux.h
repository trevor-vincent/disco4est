#ifndef IP_FLUX_H
#define IP_FLUX_H 

#include <curved_compute_flux.h>
#include <ip_flux_params.h>

typedef struct {

  ip_flux_params_t* ip_flux_params;
  curved_flux_fcn_ptrs_t curved_flux_fcn_ptrs;

} ip_flux_t;

/* This file was automatically generated.  Do not edit! */
void ip_flux_dirichlet_destroy(ip_flux_t *ip_flux);
ip_flux_t *ip_flux_dirichlet_new(p4est_t *p4est,const char *print_prefix,const char *input_file,grid_fcn_t bndry_fcn);



#endif
