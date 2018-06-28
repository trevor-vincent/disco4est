#ifndef KRYLOV_PC_MULTIGRID_H
#define KRYLOV_PC_MULTIGRID_H 

#include <d4est_petsc.h>
#include <multigrid.h>
#include <krylov_pc.h>


typedef struct {

  multigrid_data_t* mg_data;
  void(*user_setup_fcn)(krylov_pc_t* kpc);

  int ratio_is_getting_bad_counts;
  
} krylov_pc_multigrid_data_t;


/* This file was automatically generated.  Do not edit! */
void krylov_pc_multigrid_destroy(krylov_pc_t *kpc);
void krylov_pc_multigrid_apply(krylov_pc_t *kpc,double *xp,double *yp);
krylov_pc_t*
krylov_pc_multigrid_create
(
 multigrid_data_t* mg_data,
 void(*krylov_pc_multigrid_user_setup_fcn)(krylov_pc_t* kpc)
);


#endif
