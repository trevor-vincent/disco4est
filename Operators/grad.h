#ifndef GRAD_H
#define GRAD_H 

#include "../pXest/pXest.h"
#include "../dGMath/d4est_operators.h"

void
grad
(
 double*u,
 double* gradu [(P4EST_DIM)],
 double h,
 int deg,
 d4est_operators_t* d4est_ops
);

void
grad_euc_norm
(
 double* u,
 double* gradu_norm,
 double h,
 int deg,
 d4est_operators_t* d4est_ops
);

#endif
