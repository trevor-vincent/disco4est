#ifndef NEWTON_PETSC_H
#define NEWTON_PETSC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include <krylov_petsc.h>

typedef struct {

  p4est_t* p4est;
  problem_data_t* vecs;
  weakeqn_ptrs_t* fcns;
  p4est_ghost_t* ghost;
  element_data_t* ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  
} newton_petsc_ctx_t;

/* PetscErrorCode get_residual(SNES,Vec,Vec,void*); */

void newton_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t** ghost,
 element_data_t** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 krylov_petsc_params_t* krylov_params
);

#endif
