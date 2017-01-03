#ifndef KSP_PETSC_SOLVE_H
#define KSP_PETSC_SOLVE_H 

#include "petscsnes.h"
#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include <krylov_pc.h>

typedef struct {

  int iterations;
  double residual_norm;
  
} krylov_petsc_info_t;

typedef struct {

  KSPType ksp_type;
  int ksp_monitor;
  int ksp_view;
  double ksp_atol;
  double ksp_rtol;
  int ksp_maxit;
  int count;
  
  int user_defined_pc;
  krylov_pc_create_fcn_t pc_create;
  krylov_pc_destroy_fcn_t pc_destroy;
  void* pc_data;
  
} krylov_petsc_params_t;

krylov_petsc_info_t
krylov_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 krylov_pc_create_fcn_t pc_create,
 krylov_pc_destroy_fcn_t pc_destroy,
 void* pc_data
 /* krylov_petsc_params_t* krylov_params */
);

#endif
