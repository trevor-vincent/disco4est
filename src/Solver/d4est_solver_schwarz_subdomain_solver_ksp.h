#ifndef D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_KSP_H
#define D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_KSP_H 

#include <pXest.h>
#include <d4est_solver_schwarz.h>
#include <d4est_mesh.h>
#include <petscsnes.h>

typedef struct {

  p4est_t* p4est;
  d4est_geometry_t* d4est_geom;
  d4est_quadrature_t* d4est_quad;
  d4est_mesh_data_t* d4est_factors;
  d4est_ghost_t* ghost;
  d4est_solver_schwarz_operators_t* schwarz_ops;
  d4est_solver_schwarz_metadata_t* schwarz_metadata;
  d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data;
  d4est_solver_schwarz_apply_lhs_t* apply_lhs;
  int subdomain;
  KSP* ksp;
  
} d4est_solver_schwarz_subdomain_solver_ksp_ctx_t;  


typedef struct {
  
  const char* input_section;
  int subdomain_iter;
  int subdomain_monitor;
  double subdomain_atol;
  double subdomain_rtol;
  
} d4est_solver_schwarz_subdomain_solver_ksp_data_t;  


/* This file was automatically generated.  Do not edit! */
d4est_solver_schwarz_subdomain_solver_info_t d4est_solver_schwarz_subdomain_solver_ksp(p4est_t *p4est,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *ghost,d4est_solver_schwarz_operators_t *schwarz_ops,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data,d4est_solver_schwarz_apply_lhs_t *apply_lhs,double *du_restricted_field_over_subdomain,double *rhs_restricted_field_over_subdomain,int subdomain,void *params);
d4est_solver_schwarz_subdomain_solver_ksp_data_t *d4est_solver_schwarz_subdomain_solver_ksp_init(p4est_t *p4est,const char *input_file,const char *input_section);
void d4est_solver_schwarz_subdomain_solver_ksp_destroy(void *params);


#endif
