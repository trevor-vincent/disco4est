#ifndef POISSON_DEBUG_VECS_H
#define POISSON_DEBUG_VECS_H 

#include <pXest.h>
#include <d4est_operators.h>

typedef struct {

  int elem_id;
  int deg;

  double* Mdu [(P4EST_DIM)];
  double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)];
  double* Mdivq;
  double* liftedqflux [(P4EST_FACES)];
  double* Au;
  double* q [(P4EST_DIM)];
  double* Mq [(P4EST_DIM)];
  
} poisson_debug_vecs_t;


poisson_debug_vecs_t*
poisson_debug_vecs_init
(
 p4est_t* p4est,
 int local_element_id
);

void
poisson_debug_vecs_destroy
(
 poisson_debug_vecs_t* debug_vecs
);

void
poisson_debug_vecs_set_Mdivq
(
 double* Mdivq,
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);


void
poisson_debug_vecs_set_Mdu
(
 double* Mdu [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);

void
poisson_debug_vecs_set_Au
(
 double* Au,
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
); 

void
poisson_debug_vecs_set_liftedqflux
(
 double* liftedqflux [(P4EST_FACES)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);

void
poisson_debug_vecs_set_lifteduflux
(
 double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);

void
poisson_debug_vecs_set_q
(
 double* q [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);

void
poisson_debug_vecs_set_Mq
(
 double* Mq [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
);

#endif
