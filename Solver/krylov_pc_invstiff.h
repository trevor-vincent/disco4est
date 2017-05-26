#ifndef KRYLOV_PC_INVSTIFF_H
#define KRYLOV_PC_INVSTIFF_H 

#include <pXest.h>
#include <d4est_operators.h>
#include <krylov_pc.h>


typedef struct {

  int local_matrix_nodes;
  double* inv_stiff;
  p4est_t* p4est;
  d4est_operators_t* d4est_ops;

} krylov_pc_invstiff_data_t;

krylov_pc_t*
krylov_pc_invstiff_create(p4est_t* p4est, d4est_operators_t* d4est_ops);

void
krylov_pc_invstiff_destroy
(krylov_pc_t* pc);

void
krylov_pc_invstiff_setup
(krylov_pc_t* pc);

void
krylov_pc_invstiff_apply(krylov_pc_t* kpc, double* xp, double* yp);


#endif
