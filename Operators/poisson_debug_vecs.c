#include "poisson_debug_vecs.h"
#include <element_data.h>
#include <linalg.h>


void
poisson_debug_vecs_set_lifteduflux
(
 double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int f = 0; f < (P4EST_FACES); f++){
      linalg_copy_1st_to_2nd(lifteduflux[f][d],
                             debug_vecs->lifteduflux[f][d],
                             volume_nodes);
    }
  }
}

void
poisson_debug_vecs_set_liftedqflux
(
 double* liftedqflux [(P4EST_FACES)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int f = 0; f < (P4EST_FACES); f++){
    linalg_copy_1st_to_2nd(liftedqflux[f], debug_vecs->liftedqflux[f], volume_nodes);
  }
}

void
poisson_debug_vecs_set_Mdu
(
 double* Mdu [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_copy_1st_to_2nd(Mdu[d], debug_vecs->Mdu[d], volume_nodes);
  }
}

void
poisson_debug_vecs_set_Mq
(
 double* Mq [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_copy_1st_to_2nd(Mq[d], debug_vecs->Mq[d], volume_nodes);
  }
}

void
poisson_debug_vecs_set_q
(
 double* q [(P4EST_DIM)],
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_copy_1st_to_2nd(q[d], debug_vecs->q[d], volume_nodes);
  }
}

void
poisson_debug_vecs_set_Au
(
 double* Au,
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  linalg_copy_1st_to_2nd(Au, debug_vecs->Au, volume_nodes);
}

void
poisson_debug_vecs_set_Mdivq
(
 double* Mdivq,
 poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  linalg_copy_1st_to_2nd(Mdivq, debug_vecs->Mdivq, volume_nodes);
}

poisson_debug_vecs_t*
poisson_debug_vecs_init
(
 p4est_t* p4est,
 int local_element_id
)
{
  poisson_debug_vecs_t* debug_vecs
    = P4EST_ALLOC(poisson_debug_vecs_t, 1);
  element_data_t* elem_data
    = element_data_get_element_data
    (
     p4est,
     local_element_id
    );
  
  int id = elem_data->id;
  int deg = elem_data->deg;
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), deg);

  debug_vecs->deg = deg;
  debug_vecs->elem_id = id;
  debug_vecs->deg = deg;
  debug_vecs->Mdivq = P4EST_ALLOC(double, volume_nodes);
  debug_vecs->Au = P4EST_ALLOC(double, volume_nodes);
  for (int i = 0; i < (P4EST_FACES); i++){
    debug_vecs->liftedqflux[i] = P4EST_ALLOC(double, volume_nodes);

    for(int j = 0; j < (P4EST_DIM); j++){
      if (i == 0){
        debug_vecs->q[j] = P4EST_ALLOC(double, volume_nodes);
        debug_vecs->Mq[j] = P4EST_ALLOC(double, volume_nodes);
        debug_vecs->Mdu[j] = P4EST_ALLOC(double, volume_nodes);
      }
      debug_vecs->lifteduflux[i][j] = P4EST_ALLOC(double, volume_nodes);
    }
  }
  return debug_vecs;
}

void
poisson_debug_vecs_destroy
(
 poisson_debug_vecs_t* debug_vecs
)
{
  P4EST_FREE(debug_vecs->Mdivq);
  P4EST_FREE(debug_vecs->Au);
  for (int i = 0; i < (P4EST_FACES); i++){
    P4EST_FREE(debug_vecs->liftedqflux[i]);
    for(int j = 0; j < (P4EST_DIM); j++){
      if (i == 0){
        P4EST_FREE(debug_vecs->q[j]);
        P4EST_FREE(debug_vecs->Mq[j]);
        P4EST_FREE(debug_vecs->Mdu[j]);
      }
      P4EST_FREE(debug_vecs->lifteduflux[i][j]);
    }
  }
  P4EST_FREE(debug_vecs);
}
