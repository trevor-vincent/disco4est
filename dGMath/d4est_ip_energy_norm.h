#ifndef D4EST_IP_ENERGY_NORM_H
#define D4EST_IP_ENERGY_NORM_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_poisson_flux_sipg.h>

typedef struct {
  
  penalty_calc_t u_penalty_fcn;
  h_calc_method_t sipg_flux_h;
  double penalty_prefactor;
  void* user;
  double ip_energy_norm_sqr_volume_term;
  double ip_energy_norm_sqr_boundary_term;
  double ip_energy_norm_sqr_interface_term;
  
} d4est_ip_energy_norm_data_t;

/* This file was automatically generated.  Do not edit! */
double d4est_ip_energy_norm_compute(p4est_t *p4est,double *u,d4est_ip_energy_norm_data_t *energy_norm_data,p4est_ghost_t *ghost,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_geometry_storage_t *d4est_factors);


#endif
