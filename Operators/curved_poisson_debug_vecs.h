#ifndef CURVED_POISSON_DEBUG_VECS_H
#define CURVED_POISSON_DEBUG_VECS_H 

#include <pXest.h>
#include <dgmath.h>

typedef struct {

  int elem_id;
  int deg;
  int deg_integ;

  double* Mdu [(P4EST_DIM)];

  double* Mq [(P4EST_DIM)];
  double* q [(P4EST_DIM)];  
  double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)];

  double* Mdivq;

  double* liftedqflux [(P4EST_FACES)];
  
  double* Au;
  double* u;
  
} curved_poisson_debug_vecs_t;

void
curved_poisson_debug_vecs_set_lifteduflux
(
 double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
curved_poisson_debug_vecs_set_liftedqflux
(
 double* liftedqflux [(P4EST_FACES)],
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
curved_poisson_debug_vecs_set_Mdu
(
 double* Mdu [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
curved_poisson_debug_vecs_set_Au
(
 double* Au,
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);


void
curved_poisson_debug_vecs_set_u
(
 double* u,
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);


void
curved_poisson_debug_vecs_set_Mdivq
(
 double* Mdivq,
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

curved_poisson_debug_vecs_t*
curved_poisson_debug_vecs_init
(
 p4est_t* p4est,
 int local_element_id
);

void
curved_poisson_debug_vecs_destroy
(
 curved_poisson_debug_vecs_t* debug_vecs
);

void
curved_poisson_debug_vecs_set_q
(
 double* q [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
curved_poisson_debug_vecs_set_Mq
(
 double* Mq [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
curved_poisson_debug_vecs_print
(
 curved_poisson_debug_vecs_t* debug_vecs
);

void
curved_poisson_debug_vecs_2print
(
 curved_poisson_debug_vecs_t* debug_vecs1,
 curved_poisson_debug_vecs_t* debug_vecs2
);

#endif
