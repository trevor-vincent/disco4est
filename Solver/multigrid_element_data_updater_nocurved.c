#include <element_data.h>
#include <multigrid_element_data_updater_nocurved.h>

multigrid_element_data_updater_t*
multigrid_element_data_updater_nocurved_init
(
 p4est_ghost_t** ghost,
 element_data_t** ghost_data
)
{
  multigrid_element_data_updater_t* updater = P4EST_ALLOC(multigrid_element_data_updater_t,1);
  updater->get_local_nodes = element_data_get_local_nodes;
  updater->update = multigrid_element_data_updater_nocurved_update;
  updater->ghost = ghost;
  updater->ghost_data = (void**)ghost_data;
  
  /* for curved infrastructure*/
  updater->element_data_init_user_fcn = NULL;
  updater->d4est_geom = NULL;
  updater->geometric_factors = NULL;
  updater->user = NULL;
  return updater;
}

void
multigrid_element_data_updater_nocurved_destroy
(
 multigrid_element_data_updater_t* updater
)
{
  P4EST_FREE(updater);
  
}

void
multigrid_element_data_updater_nocurved_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  
  if (mg_data->mg_state == DOWNV_POST_COARSEN){
    element_data_init(p4est, -1);
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
    element_data_init(p4est, -1);

    /* update ghost data */
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count); 
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){    
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count); 

    element_data_init(p4est, -1);
  }
  else {
    return;
  }
}
