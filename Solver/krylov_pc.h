#ifndef KRYLOV_PC_H
#define KRYLOV_PC_H 

#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"

typedef struct {

  p4est_t* p4est;
  problem_data_t* vecs;
  void* fcns;
  p4est_ghost_t** ghost;
  void** ghost_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom;
  void* pc_data;
  
} krylov_pc_ctx_t;

typedef
void
(*krylov_pc_apply_fcn_t)
(
 krylov_pc_ctx_t*,
 double*,
 double*
);

typedef
void
(*krylov_pc_setup_fcn_t)
(
 krylov_pc_ctx_t*
);

typedef struct {

  krylov_pc_ctx_t* pc_ctx;
  krylov_pc_apply_fcn_t pc_apply;
  krylov_pc_setup_fcn_t pc_setup;
  
} krylov_pc_t;

typedef
krylov_pc_t*
(*krylov_pc_create_fcn_t)
(
 krylov_pc_ctx_t*
);

typedef
void
(*krylov_pc_destroy_fcn_t)
(
 krylov_pc_t*
);


#endif
