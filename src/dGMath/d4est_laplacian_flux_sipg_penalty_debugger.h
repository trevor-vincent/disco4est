#ifndef D4EST_LAPLACIAN_FLUX_SIPG_PENALTY_DEBUGGER_H
#define D4EST_LAPLACIAN_FLUX_SIPG_PENALTY_DEBUGGER_H 

#include <pXest.h>
#include <d4est_mortars.h>
#include <d4est_laplacian_flux.h>
#include <d4est_laplacian_flux_sipg.h>

typedef struct {

  penalty_calc_t penalty_fcn;
  double penalty_prefactor;

  double* average_min_penalty_vtk;
  double* average_max_penalty_vtk;
  double* average_mean_penalty_vtk;

  double* min_penalty_vtk_per_face [(P4EST_FACES)];
  double* mean_penalty_vtk_per_face [(P4EST_FACES)];
  double* max_penalty_vtk_per_face [(P4EST_FACES)];
  
} d4est_laplacian_flux_sipg_penalty_debugger_t;

/* This file was automatically generated.  Do not edit! */
void d4est_laplacian_flux_sipg_penalty_debugger_destroy(d4est_laplacian_flux_sipg_penalty_debugger_t *debugger);
d4est_laplacian_flux_sipg_penalty_debugger_t *d4est_laplacian_flux_sipg_penalty_debugger_get_vtk_data(p4est_t *p4est,d4est_ghost_t *d4est_ghost,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_laplacian_flux_sipg_penalty_debugger_t *debugger);
d4est_laplacian_flux_sipg_penalty_debugger_t *d4est_laplacian_flux_sipg_penalty_debugger_init(p4est_t *p4est,penalty_calc_t penalty_fcn,double penalty_prefactor);


#endif
