#ifndef KRYLOV_PC_INVSTIFF_H
#define KRYLOV_PC_INVSTIFF_H 

#include <pXest.h>
#include <dgmath.h>
#include <krylov_pc.h>


typedef struct {

  int local_matrix_nodes;
  double* inv_stiff;
  p4est_t* p4est;
  dgmath_jit_dbase_t* dgmath_jit_dbase;

} krylov_pc_invstiff_data_t;

krylov_pc_t*
krylov_pc_invstiff_create(p4est_t* p4est, dgmath_jit_dbase_t* dgmath_jit_dbase);

void
krylov_pc_invstiff_destroy
(krylov_pc_t* pc);

void
krylov_pc_invstiff_setup
(krylov_pc_t* pc);

void
krylov_pc_invstiff_apply(krylov_pc_t* kpc, double* xp, double* yp);


#endif
