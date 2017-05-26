/**
 * @file   compute_flux.h
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Fri Sep 25 01:53:54 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef COMPUTE_FLUX_H
#define COMPUTE_FLUX_H 

#include "../ElementData/element_data.h"

/* No problem specific data needed */
typedef void (*flux_interface_fcn_t)
(
 element_data_t** e_m,
 int faces_m,
 int f_m,
 element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 d4est_operators_t* d4est_ops,
 void*
);


/* We need the problem specific data */
typedef void (*flux_boundary_fcn_t)
(
 element_data_t*,
 int,
 grid_fcn_t,
 d4est_operators_t* d4est_ops,
 void*
);

typedef struct {

  flux_interface_fcn_t flux_interface_fcn;
  flux_boundary_fcn_t flux_boundary_fcn;
  grid_fcn_t bndry_fcn;
  void* params;
  
} flux_fcn_ptrs_t;

typedef struct {

  d4est_operators_t* d4est_ops;
  flux_fcn_ptrs_t* flux_fcn_ptrs;

} compute_flux_user_data_t;

void compute_flux_on_local_elements(p4est_iter_face_info_t * info, void *user_data);

#endif
