#ifndef KRYLOV_PC_MULTIGRID_H
#define KRYLOV_PC_MULTIGRID_H 

#include <krylov_pc.h>

krylov_pc_t*
krylov_pc_multigrid_create(krylov_pc_ctx_t* kct);

void
krylov_pc_multigrid_setup(krylov_pc_ctx_t* pc_ctx);

void
krylov_pc_multigrid_destroy(krylov_pc_t* pc);

void
krylov_pc_multigrid_apply(krylov_pc_ctx_t* kct, double* xp, double* yp);


#endif
