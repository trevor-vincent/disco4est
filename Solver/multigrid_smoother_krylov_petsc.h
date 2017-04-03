#ifndef MULTIGRID_SMOOTHER_KRYLOV_PETSC_H
#define MULTIGRID_SMOOTHER_KRYLOV_PETSC_H 

#include <multigrid.h>

/* This file was automatically generated.  Do not edit! */
void multigrid_smoother_krylov_petsc_destroy(multigrid_smoother_t *solver);
multigrid_smoother_t *multigrid_smoother_krylov_petsc_init(p4est_t *p4est,const char *input_file);

#endif
