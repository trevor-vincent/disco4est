#ifndef KRYLOV_PC_MULTIGRID_H
#define KRYLOV_PC_MULTIGRID_H 

#include <d4est_solver_petsc.h>
#include <d4est_solver_multigrid.h>
#include <d4est_krylov_pc.h>


typedef struct {

   d4est_solver_multigrid_t* mg_data;
  void(*user_setup_fcn)(d4est_krylov_pc_t* kpc);

  int ratio_is_getting_bad_counts;
  
} d4est_krylov_pc_multigrid_data_t;


/* This file was automatically generated.  Do not edit! */
void d4est_krylov_pc_multigrid_destroy(d4est_krylov_pc_t *kpc);
void d4est_krylov_pc_multigrid_apply(d4est_krylov_pc_t *kpc,double *xp,double *yp);
d4est_krylov_pc_t*
d4est_krylov_pc_multigrid_create
(
  d4est_solver_multigrid_t* mg_data,
 void(*d4est_krylov_pc_multigrid_user_setup_fcn)(d4est_krylov_pc_t* kpc)
);


#endif
