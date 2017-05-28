#ifndef CURVED_COMPUTE_FLUX_H
#define CURVED_COMPUTE_FLUX_H 



/**
 * @file   compute_flux.h
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Fri Sep 25 01:53:54 2015
 * 
 * @brief  
 * 
 * 
 */

#include "../ElementData/d4est_element_data.h"
#include <d4est_geometry.h>
#include <d4est_quadrature.h>

/* No problem specific data needed */
typedef void (*curved_flux_interface_fcn_t)
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 d4est_quadrature_t* d4est_quad,
 void* params
);

/* We need the problem specific data */
typedef void (*curved_flux_boundary_fcn_t)
(
 d4est_element_data_t*,
 int,
 grid_fcn_t,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 d4est_quadrature_t* d4est_quad,
 void* params
);

typedef struct {

  curved_flux_interface_fcn_t flux_interface_fcn;
  curved_flux_boundary_fcn_t flux_boundary_fcn;
  grid_fcn_t bndry_fcn;
  void* params;
  
} curved_flux_fcn_ptrs_t;

typedef struct {

  d4est_operators_t* d4est_ops;
  d4est_quadrature_t* d4est_quad;
  curved_flux_fcn_ptrs_t* flux_fcn_ptrs;
  d4est_geometry_t* geom;

} curved_compute_flux_user_data_t;

void curved_compute_flux_on_local_elements
(
 p4est_iter_face_info_t * info,
 void *user_data
);

#endif
