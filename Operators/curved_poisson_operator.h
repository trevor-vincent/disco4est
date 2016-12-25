#ifndef CURVED_GAUSS_POISSON_CALLBACKS_H
#define CURVED_GAUSS_POISSON_CALLBACKS_H

#include "curved_poisson_debug_vecs.h"

void
curved_Gauss_poisson_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom
);

curved_poisson_debug_vecs_t*
curved_Gauss_poisson_apply_aij_debug
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 int local_element_id
);

void
curved_Gauss_poisson_build_rhs_with_strongbc
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void curved_Gauss_poisson_init_vecs
(
 p4est_iter_volume_info_t * info,
 void *user_data
);

void curved_Gauss_poisson_destroy_vecs
(
 p4est_iter_volume_info_t * info,
 void *user_data
);

void curved_Gauss_poisson_compute_q_elem
(
 p4est_iter_volume_info_t * info,
 void *user_data
);

void curved_Gauss_poisson_compute_Au_elem
(
 p4est_iter_volume_info_t* info,
 void* user_data
);

#endif
