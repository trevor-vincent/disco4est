#ifndef D4EST_PETSC_H
#define D4EST_PETSC_H 

#include <pXest.h>
#include <d4est_elliptic_data.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_operators.h>
#include <d4est_quadrature.h>
#include <d4est_geometry.h>
#include <petscsnes.h>

typedef struct {

  p4est_t* p4est;
  d4est_elliptic_data_t* vecs;
  d4est_elliptic_eqns_t* fcns;
  p4est_ghost_t** ghost;
  d4est_element_data_t** ghost_data;
  d4est_operators_t* d4est_ops;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_geometry_storage_t* d4est_factors;
  
} petsc_ctx_t;

#endif
