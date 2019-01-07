#ifndef D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_H
#define D4EST_SOLVER_SCHWARZ_SUBDOMAIN_SOLVER_H 

typedef enum {SUBDOMAIN_SOLVER_NOT_SET, SUBDOMAIN_SOLVER_CG, SUBDOMAIN_SOLVER_GMRES} d4est_solver_schwarz_subdomain_solver_type_t;

typedef
void
(*d4est_solver_schwarz_subdomain_solver_destroy_fcn_t)
(
 void*
);


typedef
void
(*d4est_solver_schwarz_subdomain_solver_fcn_t)
(
 p4est_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_mesh_data_t*,
 d4est_ghost_t*,
 d4est_solver_schwarz_operators_t*,
 d4est_solver_schwarz_metadata_t*,
 d4est_solver_schwarz_geometric_data_t*,
 d4est_solver_schwarz_apply_lhs_t*,
 double*,
 double*,
 int,
 void*
);


typedef struct {

  d4est_solver_schwarz_subdomain_solver_type_t solver_type;
  d4est_solver_schwarz_subdomain_solver_fcn_t solver_fcn;
  d4est_solver_schwarz_subdomain_solver_destroy_fcn_t destroy_fcn;
  void* solver_ctx;

  const char* input_section;
  
} d4est_solver_schwarz_subdomain_solver_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_subdomain_solver_destroy(d4est_solver_schwarz_subdomain_solver_t *schwarz_subdomain_solver);
d4est_solver_schwarz_subdomain_solver_t *d4est_solver_schwarz_subdomain_solver_init(p4est_t *p4est,const char *input_file,const char *input_section);

#endif
