#ifndef KSP_PETSC_SOLVE_H
#define KSP_PETSC_SOLVE_H 

#include "petscsnes.h"
#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"

#include <d4est_petsc.h>
#include <krylov_pc.h>

typedef struct {

  int iterations;
  double residual_norm;
  
} krylov_petsc_info_t;

typedef struct {

  /* the biggest type is "richardson" = 10 */
  char ksp_type [15]; 
  int ksp_monitor;
  int ksp_view;
  double ksp_atol;
  double ksp_rtol;
  int ksp_maxit;
  krylov_pc_t* pc;

  int count;
  
} krylov_petsc_params_t;

/* This file was automatically generated.  Do not edit! */
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
 const char* input_file,
 krylov_pc_t* pc
);

krylov_petsc_params_t
krylov_petsc_input
(
 p4est_t* p4est,
 const char* input_file
);

#endif
