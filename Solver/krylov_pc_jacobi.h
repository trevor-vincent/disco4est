#ifndef KRYLOV_PC_JACOBI_H
#define KRYLOV_PC_JACOBI_H 

#include <krylov_pc.h>

typedef struct {

  int local_nodes;
  double* inv_aii;

} krylov_pc_jacobi_data_t;

krylov_pc_t*
krylov_pc_jacobi_create
(krylov_pc_ctx_t* ctx);

void
krylov_pc_jacobi_destroy
(krylov_pc_t* pc);

void
krylov_pc_jacobi_setup
(krylov_pc_ctx_t* kct);

void
krylov_pc_jacobi_apply(krylov_pc_ctx_t* kct, double* xp, double* yp);

#endif
