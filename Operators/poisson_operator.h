#ifndef POISSON_CALLBACKS_H
#define POISSON_CALLBACKS_H 

#include <poisson_debug_vecs.h>

void poisson_apply_aij(
                         p4est_t* p4est,
                         p4est_ghost_t* ghost,
                         void* ghost_data,
                         problem_data_t* prob_vecs,
			 d4est_operators_t* d4est_ops,
                         d4est_geometry_t* d4est_geom
);
poisson_debug_vecs_t*
poisson_apply_aij_debug(
                         p4est_t* p4est,
                         p4est_ghost_t* ghost,
                         element_data_t* ghost_data,
                         problem_data_t* prob_vecs,
			 d4est_operators_t* d4est_ops,
                         int local_element_id
);
#endif
