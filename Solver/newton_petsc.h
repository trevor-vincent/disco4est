#ifndef NEWTON_PETSC_H
#define NEWTON_PETSC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include <d4est_petsc.h>
#include <krylov_petsc.h>

typedef struct {

  int snes_monitor;
  int snes_linesearch_monitor;
  int snes_linesearch_order;
  double snes_atol;
  double snes_rtol;
  int snes_max_funcs;
  PetscOptions petsc_options;
  
  
} newton_petsc_params_t;

void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 PetscOptions krylov_options,
 PetscOptions newton_options,
 krylov_pc_t* pc
);

#endif
