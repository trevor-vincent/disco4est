#ifndef GRAD_H
#define GRAD_H 

#include "../pXest/pXest.h"
#include "../dGMath/dgmath.h"

void
grad
(
 double*u,
 double* gradu [(P4EST_DIM)],
 double h,
 int deg,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
grad_euc_norm
(
 double* u,
 double* gradu_norm,
 double h,
 int deg,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

#endif
