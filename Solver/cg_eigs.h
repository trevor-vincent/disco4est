#ifndef CG_EIGS_H
#define CG_EIGS_H 

#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../dGMath/d4est_operators.h"

void
cg_eigs
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 int imax,
 double* eig_max
);



#endif
