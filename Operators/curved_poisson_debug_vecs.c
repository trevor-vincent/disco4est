#include <util.h>
#include "curved_poisson_debug_vecs.h"
#include <d4est_element_data.h>
#include <d4est_linalg.h>


void
curved_poisson_debug_vecs_set_lifteduflux
(
 double* lifteduflux [(P4EST_FACES)][(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int f = 0; f < (P4EST_FACES); f++){
      d4est_linalg_copy_1st_to_2nd
        (
         lifteduflux[f][d],
         debug_vecs->lifteduflux[f][d],
         volume_nodes
        );
    }
  }
}

void
curved_poisson_debug_vecs_set_liftedqflux
(
 double* liftedqflux [(P4EST_FACES)],
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int f = 0; f < (P4EST_FACES); f++){
    d4est_linalg_copy_1st_to_2nd(liftedqflux[f], debug_vecs->liftedqflux[f], volume_nodes);
  }
}

void
curved_poisson_debug_vecs_set_Mdu
(
 double* Mdu [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_copy_1st_to_2nd(Mdu[d], debug_vecs->Mdu[d], volume_nodes);
  }
}



void
curved_poisson_debug_vecs_set_Mq
(
 double* Mq [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_copy_1st_to_2nd(Mq[d], debug_vecs->Mq[d], volume_nodes);
  }
}

void
curved_poisson_debug_vecs_set_q
(
 double* q [(P4EST_DIM)],
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  /* printf("[WARNING] set_q only valid for constant jacobian meshes\n"); */
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_copy_1st_to_2nd(q[d], debug_vecs->q[d], volume_nodes);    
  }
}


void
curved_poisson_debug_vecs_set_Au
(
 double* Au,
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  d4est_linalg_copy_1st_to_2nd(Au, debug_vecs->Au, volume_nodes);
}


void
curved_poisson_debug_vecs_set_u
(
 double* u,
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg);
  d4est_linalg_copy_1st_to_2nd(u, debug_vecs->u, volume_nodes);
}

void
curved_poisson_debug_vecs_set_Mdivq
(
 double* Mdivq,
 curved_poisson_debug_vecs_t* debug_vecs,
 d4est_operators_t* d4est_ops
)
{
  int volume_nodes_quad = d4est_operators_get_nodes((P4EST_DIM), debug_vecs->deg_quad);
  d4est_linalg_copy_1st_to_2nd(Mdivq, debug_vecs->Mdivq, volume_nodes_quad);
}

curved_poisson_debug_vecs_t*
curved_poisson_debug_vecs_init
(
 p4est_t* p4est,
 int local_element_id
)
{
  curved_poisson_debug_vecs_t* debug_vecs
    = P4EST_ALLOC(curved_poisson_debug_vecs_t, 1);
  d4est_element_data_t* elem_data
    = d4est_element_data_get_element_data
    (
     p4est,
     local_element_id
    );
  
  int id = elem_data->id;
  int deg = elem_data->deg;
  int deg_quad = elem_data->deg_quad;
  int volume_nodes_quad = d4est_operators_get_nodes((P4EST_DIM), deg_quad);
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), deg);

  debug_vecs->deg = deg;
  debug_vecs->elem_id = id;
  debug_vecs->deg_quad = deg_quad;
  
  debug_vecs->Mdivq = P4EST_ALLOC(double, volume_nodes_quad);
  debug_vecs->Au = P4EST_ALLOC(double, volume_nodes);
  debug_vecs->u = P4EST_ALLOC(double, volume_nodes);
  for (int i = 0; i < (P4EST_FACES); i++){
    debug_vecs->liftedqflux[i] = P4EST_ALLOC(double, volume_nodes_quad);
    for(int j = 0; j < (P4EST_DIM); j++){
      if (i == 0){
        debug_vecs->q[j] = P4EST_ALLOC(double, volume_nodes);
        debug_vecs->Mq[j] = P4EST_ALLOC(double, volume_nodes);
        debug_vecs->Mdu[j] = P4EST_ALLOC(double, volume_nodes);
      }
      debug_vecs->lifteduflux[i][j] = P4EST_ALLOC(double, volume_nodes_quad);
    }
  }
  return debug_vecs;
}

void
curved_poisson_debug_vecs_destroy
(
 curved_poisson_debug_vecs_t* debug_vecs
)
{
  P4EST_FREE(debug_vecs->Mdivq);
  P4EST_FREE(debug_vecs->Au);
  P4EST_FREE(debug_vecs->u);
  for (int i = 0; i < (P4EST_FACES); i++){
    P4EST_FREE(debug_vecs->liftedqflux[i]);
    for(int j = 0; j < (P4EST_DIM); j++){
      if (i == 0){
        P4EST_FREE(debug_vecs->Mq[j]);
        P4EST_FREE(debug_vecs->q[j]);
        P4EST_FREE(debug_vecs->Mdu[j]);
      }
      P4EST_FREE(debug_vecs->lifteduflux[i][j]);
    }
  }
  P4EST_FREE(debug_vecs);
}


void
curved_poisson_debug_vecs_print
(
 curved_poisson_debug_vecs_t* debug_vecs
)
{
  printf("** Debug Info for Element %d **\n", debug_vecs->elem_id);
  
  DEBUG_PRINT_2ARR_DBL(debug_vecs->Mdu[0],
                       debug_vecs->Mdu[1],
                       d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg)
                      );

  DEBUG_PRINT_4ARR_DBL(debug_vecs->lifteduflux[0][0],
                       debug_vecs->lifteduflux[1][0],
                       debug_vecs->lifteduflux[2][0],
                       debug_vecs->lifteduflux[3][0],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg))
                      );


  DEBUG_PRINT_4ARR_DBL(debug_vecs->lifteduflux[0][1],
                       debug_vecs->lifteduflux[1][1],
                       debug_vecs->lifteduflux[2][1],
                       debug_vecs->lifteduflux[3][1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg))
                      );


  DEBUG_PRINT_2ARR_DBL(debug_vecs->Mq[0],
                       debug_vecs->Mq[1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg))
                      );


  DEBUG_PRINT_2ARR_DBL(debug_vecs->q[0],
                       debug_vecs->q[1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg))
                      );
  

  DEBUG_PRINT_4ARR_DBL(debug_vecs->liftedqflux[0],
                       debug_vecs->liftedqflux[1],
                       debug_vecs->liftedqflux[2],
                       debug_vecs->liftedqflux[3],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg))
                      );
   
  DEBUG_PRINT_ARR_DBL(debug_vecs->Mdivq, (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg)));

  DEBUG_PRINT_ARR_DBL(debug_vecs->Au, (d4est_operators_get_nodes((P4EST_DIM),debug_vecs->deg)));
}


void
curved_poisson_debug_vecs_2print
(
 curved_poisson_debug_vecs_t* debug_vecs1,
 curved_poisson_debug_vecs_t* debug_vecs2
)
{
  mpi_assert(debug_vecs1->deg == debug_vecs2->deg);
  printf("** Debug Info for Element %d and Element %d**\n", debug_vecs1->elem_id, debug_vecs2->elem_id);

  DEBUG_PRINT_2ARR_DBL(debug_vecs1->u, debug_vecs2->u, (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg)));

  
  DEBUG_PRINT_4ARR_DBL(debug_vecs1->Mdu[0],
                       debug_vecs1->Mdu[1],
                       debug_vecs2->Mdu[0],
                       debug_vecs2->Mdu[1],
                       d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg)
                      );

  DEBUG_PRINT_8ARR_DBL(
                       debug_vecs1->lifteduflux[0][0],
                       debug_vecs1->lifteduflux[1][0],
                       debug_vecs1->lifteduflux[2][0],
                       debug_vecs1->lifteduflux[3][0],
                       debug_vecs2->lifteduflux[0][0],
                       debug_vecs2->lifteduflux[1][0],
                       debug_vecs2->lifteduflux[2][0],
                       debug_vecs2->lifteduflux[3][0],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg))
                      );


  DEBUG_PRINT_8ARR_DBL(
                       debug_vecs1->lifteduflux[0][1],
                       debug_vecs1->lifteduflux[1][1],
                       debug_vecs1->lifteduflux[2][1],
                       debug_vecs1->lifteduflux[3][1],
                       debug_vecs2->lifteduflux[0][1],
                       debug_vecs2->lifteduflux[1][1],
                       debug_vecs2->lifteduflux[2][1],
                       debug_vecs2->lifteduflux[3][1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg))
                      );


  DEBUG_PRINT_4ARR_DBL(
                       debug_vecs1->Mq[0],
                       debug_vecs1->Mq[1],
                       debug_vecs2->Mq[0],
                       debug_vecs2->Mq[1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg))
                      );


  DEBUG_PRINT_4ARR_DBL(
                       debug_vecs1->q[0],
                       debug_vecs1->q[1],
                       debug_vecs2->q[0],
                       debug_vecs2->q[1],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg))
                      );
  

  DEBUG_PRINT_8ARR_DBL(
                       debug_vecs1->liftedqflux[0],
                       debug_vecs1->liftedqflux[1],
                       debug_vecs1->liftedqflux[2],
                       debug_vecs1->liftedqflux[3],
                       debug_vecs2->liftedqflux[0],
                       debug_vecs2->liftedqflux[1],
                       debug_vecs2->liftedqflux[2],
                       debug_vecs2->liftedqflux[3],
                       (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg))
                      );
   
  DEBUG_PRINT_2ARR_DBL(debug_vecs1->Mdivq, debug_vecs2->Mdivq, (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg)));
  DEBUG_PRINT_2ARR_DBL(debug_vecs1->Au, debug_vecs2->Au, (d4est_operators_get_nodes((P4EST_DIM),debug_vecs1->deg)));
}


void
curved_poisson_debug_vecs_2print_1node
(
 curved_poisson_debug_vecs_t* debug_vecs1,
 int n1,
 curved_poisson_debug_vecs_t* debug_vecs2,
 int n2
)
{
  mpi_assert(debug_vecs1->deg == debug_vecs2->deg);
  printf("** Debug Info for Element %d, Node %d and Element %d, Node %d **\n", debug_vecs1->elem_id, n1, debug_vecs2->elem_id, n2);

  DEBUG_PRINT_2DBL
    (
     debug_vecs1->u[n1],
     debug_vecs2->u[n2]
    );

  
  DEBUG_PRINT_4DBL(
                   debug_vecs1->Mdu[0][n1],
                   debug_vecs1->Mdu[1][n1],
                   debug_vecs2->Mdu[0][n2],
                   debug_vecs2->Mdu[1][n2]
                  );

  DEBUG_PRINT_8DBL(
                       debug_vecs1->lifteduflux[0][0][n1],
                       debug_vecs1->lifteduflux[1][0][n1],
                       debug_vecs1->lifteduflux[2][0][n1],
                       debug_vecs1->lifteduflux[3][0][n1],
                       debug_vecs2->lifteduflux[0][0][n2],
                       debug_vecs2->lifteduflux[1][0][n2],
                       debug_vecs2->lifteduflux[2][0][n2],
                       debug_vecs2->lifteduflux[3][0][n2]
                      );


  DEBUG_PRINT_8DBL(
                       debug_vecs1->lifteduflux[0][1][n1],
                       debug_vecs1->lifteduflux[1][1][n1],
                       debug_vecs1->lifteduflux[2][1][n1],
                       debug_vecs1->lifteduflux[3][1][n1],
                       debug_vecs2->lifteduflux[0][1][n2],
                       debug_vecs2->lifteduflux[1][1][n2],
                       debug_vecs2->lifteduflux[2][1][n2],
                       debug_vecs2->lifteduflux[3][1][n2]
                      );


  DEBUG_PRINT_4DBL(
                       debug_vecs1->Mq[0][n1],
                       debug_vecs1->Mq[1][n1],
                       debug_vecs2->Mq[0][n2],
                       debug_vecs2->Mq[1][n2]
                      );


  DEBUG_PRINT_4DBL
    (
     debug_vecs1->q[0][n1],
     debug_vecs1->q[1][n1],
     debug_vecs2->q[0][n2],
     debug_vecs2->q[1][n2]
    );
  

  DEBUG_PRINT_8DBL
    (
     debug_vecs1->liftedqflux[0][n1],
     debug_vecs1->liftedqflux[1][n1],
     debug_vecs1->liftedqflux[2][n1],
     debug_vecs1->liftedqflux[3][n1],
     debug_vecs2->liftedqflux[0][n2],
     debug_vecs2->liftedqflux[1][n2],
     debug_vecs2->liftedqflux[2][n2],
     debug_vecs2->liftedqflux[3][n2]
    );
   
  DEBUG_PRINT_2DBL(debug_vecs1->Mdivq[n1], debug_vecs2->Mdivq[n2]);
  DEBUG_PRINT_2DBL(debug_vecs1->Au[n1], debug_vecs2->Au[n2]);
}

