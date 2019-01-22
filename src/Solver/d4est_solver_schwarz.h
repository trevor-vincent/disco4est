#ifndef D4EST_SOLVER_SCHWARZ_H
#define D4EST_SOLVER_SCHWARZ_H 

#include <pXest.h>
#include <d4est_ghost_data_ext.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_apply_lhs.h>
#include <d4est_solver_schwarz_subdomain_solver.h>
#include <d4est_operators.h>

typedef struct {

  /* flag for if operators was externally/internally (1/0) allocated */
  int operators_external_alloc;
  d4est_solver_schwarz_operators_t* operators;
  
  d4est_solver_schwarz_metadata_t* metadata;
  d4est_solver_schwarz_geometric_data_t* geometric_data;
  d4est_solver_schwarz_subdomain_solver_t* subdomain_solver;
  d4est_solver_schwarz_apply_lhs_t* apply_lhs;
  
  d4est_ghost_data_ext_t* correction_ghost_data;
  d4est_ghost_data_t* residual_ghost_data;

  double* subdomain_solve_residuals;
  int* subdomain_solve_iterations;

  /* temporary for debugging */
  int used_as_smoother;
  int debug_output_ksp_level;
  int debug_output_mg_level;
  int debug_output_amr_level;
  
} d4est_solver_schwarz_t;


/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_debug_vtk(p4est_t *p4est,d4est_solver_schwarz_t *schwarz,char *input_file,char *input_section,char *save_as,char *folder,int sub_folder_number,const char **dg_field_names,double **dg_fields);
void d4est_solver_schwarz_iterate(p4est_t *p4est,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_mesh_data_t *d4est_factors,d4est_ghost_t *ghost,d4est_solver_schwarz_t *schwarz,double *u,double *r);
void d4est_solver_schwarz_destroy(d4est_solver_schwarz_t *schwarz);
d4est_solver_schwarz_t *d4est_solver_schwarz_init(p4est_t *p4est,d4est_operators_t *d4est_ops,d4est_geometry_t *d4est_geom,d4est_quadrature_t *d4est_quad,d4est_ghost_t *d4est_ghost,d4est_mesh_data_t *d4est_factors,d4est_solver_schwarz_operators_t *schwarz_ops,d4est_solver_schwarz_apply_lhs_t *apply_lhs,const char *input_file,const char *input_section);


#endif
