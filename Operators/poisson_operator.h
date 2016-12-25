#ifndef POISSON_CALLBACKS_H
#define POISSON_CALLBACKS_H 

#include <poisson_debug_vecs.h>

void
poisson_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
poisson_build_rhs_with_strongbc
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 problem_data_t* problem_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

poisson_debug_vecs_t*
poisson_apply_aij_debug(
                         p4est_t* p4est,
                         p4est_ghost_t* ghost,
                         element_data_t* ghost_data,
                         problem_data_t* prob_vecs,
			 dgmath_jit_dbase_t* dgmath_jit_dbase,
                         int local_element_id
);
#endif
