#ifndef CG_H
#define CG_H 

#include <problem_data.h>
#include <d4est_elliptic_eqns.h>

typedef struct {

  int max_iter;
  int mpi_rank;
  int monitor;
  double rtol;
  double atol;

  int final_iter;
  double final_fnrm;
  
} cg_solver_params_t;

void
cg_nr_solve
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_operators_t* d4est_ops,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double eta_max,
 double fnrm,
 int max_iter,
 int mpi_rank,
 int* final_iter,
 double* final_fnrm
);

void
cg_solve
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_operators_t* d4est_ops,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 cg_solver_params_t* params
);

void curved_cg_solve
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 cg_solver_params_t* params
);

#endif
