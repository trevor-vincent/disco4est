#ifndef NEWTON_PETSC_H
#define NEWTON_PETSC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/d4est_elliptic_eqns.h"
#include <d4est_petsc.h>
#include <krylov_petsc.h>

typedef struct {

  int snes_monitor;
  int snes_linesearch_monitor;
  int snes_view;
  int snes_converged_reason;
  
  char snes_atol [25];
  char snes_type [25];
  char snes_rtol [25];
  char snes_stol [25];
  char snes_max_funcs [25];
  char snes_max_it [25];
  char snes_trtol [25];
  char snes_linesearch_order [25];  
  
} newton_petsc_params_t;


void newton_petsc_solve
(
 p4est_t* p4est,
 d4est_elliptic_problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 krylov_petsc_params_t* krylov_options,
 newton_petsc_params_t* newton_options,
 krylov_pc_t* krylov_pc
);

void
newton_petsc_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* printf_prefix,
 newton_petsc_params_t* input
);

#endif
