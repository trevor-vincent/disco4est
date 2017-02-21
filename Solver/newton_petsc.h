#ifndef NEWTON_PETSC_H
#define NEWTON_PETSC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include <krylov_petsc.h>

typedef struct {

  p4est_t* p4est;
  problem_data_t* vecs;
  void* fcns;
  p4est_ghost_t* ghost;
  void* ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom;
  
} newton_petsc_ctx_t;

/* PetscErrorCode get_residual(SNES,Vec,Vec,void*); */

void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 const char* input_file,
 krylov_pc_create_fcn_t pc_create,
 krylov_pc_destroy_fcn_t pc_destroy,
 void* pc_data
);

#endif
