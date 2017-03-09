#ifndef KRYLOV_PC_MULTIGRID_H
#define KRYLOV_PC_MULTIGRID_H 

#include <d4est_petsc.h>
#include <multigrid.h>
#include <krylov_pc.h>

/* This file was automatically generated.  Do not edit! */
void krylov_pc_multigrid_destroy(krylov_pc_t *pc);
void krylov_pc_multigrid_apply(void *kpc_in,double *xp,double *yp);
krylov_pc_t *krylov_pc_multigrid_create(multigrid_data_t *mg_data);


#endif
