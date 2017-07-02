#ifndef D4EST_ELLIPTIC_EQNS_H
#define D4EST_ELLIPTIC_EQNS_H 

#include <problem_data.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_mortars.h>
#include <d4est_operators.h>

typedef
void (*d4est_apply_operator_fcn_t)
(
 p4est_t*,
 p4est_ghost_t*,
 d4est_element_data_t*, /* element_data */
 d4est_elliptic_problem_data_t*,
 d4est_mortar_fcn_ptrs_t*,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 void*
);

typedef struct {

  /* function pointer to apply LHS of weak equations */
  d4est_apply_operator_fcn_t apply_lhs;

  /* function pointer to build residual */
  d4est_apply_operator_fcn_t build_residual;

  /* flux fcn data */
  d4est_mortar_fcn_ptrs_t* flux_fcn_data;

  /* user-specified data */
  void* user;
  
} d4est_elliptic_eqns_t;


void d4est_elliptic_eqns_build_residual(p4est_t *p4est,p4est_ghost_t *ghost,d4est_elliptic_eqns_t *eqns,d4est_elliptic_problem_data_t *vecs,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);
void d4est_elliptic_eqns_apply_lhs(p4est_t *p4est,p4est_ghost_t *ghost,d4est_elliptic_eqns_t *eqns,d4est_elliptic_problem_data_t *vecs,d4est_element_data_t *ghost_data,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad);


#endif
