#ifndef D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_GMRES_H
#define D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_GMRES_H 


#include <pXest.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_apply_lhs.h>
#include <d4est_solver_schwarz_subdomain_solver_info.h>

typedef struct {

  int subdomain_iter;
  int subdomain_inner_iter;
  double subdomain_atol;
  double subdomain_rtol;
  int verbose;

  const char* input_section;
  
} d4est_solver_schwarz_subdomain_solver_gmres_t;

/* This file was automatically generated.  Do not edit! */
d4est_solver_schwarz_subdomain_solver_gmres_t *d4est_solver_schwarz_subdomain_solver_gmres_init(p4est_t *p4est,const char *input_file,const char *input_section);
void d4est_solver_schwarz_subdomain_solver_gmres_destroy(void *gmres_params);
d4est_solver_schwarz_subdomain_solver_info_t d4est_solver_schwarz_subdomain_solver_gmres(p4est_t *p4est,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *ghost,d4est_solver_schwarz_operators_t *schwarz_ops,d4est_solver_schwarz_metadata_t *schwarz_metadata,d4est_solver_schwarz_geometric_data_t *schwarz_geometric_data,d4est_solver_schwarz_apply_lhs_t *apply_lhs,double *du_restricted_field_over_subdomain,double *rhs_restricted_field_over_subdomain,int subdomain,void *params);

#endif
