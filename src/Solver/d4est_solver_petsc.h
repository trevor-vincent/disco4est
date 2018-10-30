#ifndef D4EST_PETSC_H
#define D4EST_PETSC_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>
#include <d4est_quadrature.h>
#include <d4est_geometry.h>
#include <petscsnes.h>
#include <time.h>

typedef struct {

  p4est_t* p4est;
  d4est_elliptic_data_t* vecs;
  d4est_elliptic_eqns_t* fcns;
  d4est_ghost_t** ghost;
  d4est_ghost_data_t** ghost_data;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_data_t* d4est_factors;
  KSP* ksp;

  int last_krylov_checkpoint_it;
  int last_newton_checkpoint_it;
  int checkpoint_every_n_krylov_its;
  int checkpoint_every_n_newton_its;
  int amr_level;

  clock_t time_start; /* time at which krylov solver or newton solver (if nonlinear eqn)  started */
  
} krylov_ctx_t;

#endif
