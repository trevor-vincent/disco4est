#ifndef KSP_PETSC_SOLVE_H
#define KSP_PETSC_SOLVE_H 

#include "petscsnes.h"
#include "../pXest/pXest.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"

#include <d4est_petsc.h>
#include <krylov_pc.h>
#include <krylov_info.h>

typedef struct {

  char input_section [50]; /* useful when using KSP in multiple contexts */
  
  char ksp_type [25];   /* the biggest type is "richardson" = 10 */
  char ksp_atol [25];
  char ksp_rtol [25];
  char ksp_max_it [25];

  int ksp_monitor;
  int ksp_view;
  int ksp_converged_reason;
  int ksp_initial_guess_nonzero;

  char ksp_chebyshev_esteig_steps [25];
  char ksp_chebyshev_esteig [25];
  int ksp_chebyshev_esteig_random;
  
    
} krylov_petsc_params_t;

/* This file was automatically generated.  Do not edit! */
krylov_info_t
krylov_petsc_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 void* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 krylov_petsc_params_t* params,
 krylov_pc_t* krylov_pc
);

void
krylov_petsc_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 krylov_petsc_params_t* params
);

void
krylov_petsc_set_options_database_from_params
(
 krylov_petsc_params_t* input
);

#endif
