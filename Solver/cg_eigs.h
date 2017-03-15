#ifndef CG_EIGS_H
#define CG_EIGS_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../dGMath/dgmath.h"

void
cg_eigs
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int imax,
 double* eig_max
);



#endif
