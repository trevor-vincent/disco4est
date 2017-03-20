#include <curved_element_data.h>
#include <util.h>
#include <multigrid_element_data_updater_curved.h>


multigrid_element_data_updater_t*
multigrid_element_data_updater_curved_init
(
 int num_of_levels,
 p4est_ghost_t** ghost,
 curved_element_data_t** ghost_data,
 geometric_factors_t* geometric_factors_toplevel,
 d4est_geometry_t* d4est_geom,
 void(*element_data_init_user_fcn)(void*,void*),
 void* user
)
{
  multigrid_element_data_updater_t* updater = P4EST_ALLOC(multigrid_element_data_updater_t,1);
  updater->get_local_nodes = curved_element_data_get_local_nodes;
  updater->update = multigrid_element_data_updater_curved_update;
  updater->ghost = ghost;
  updater->ghost_data = (void**)ghost_data;
  
  /* for curved infrastructure*/
  updater->element_data_init_user_fcn = element_data_init_user_fcn;
  updater->d4est_geom = d4est_geom;
  updater->geometric_factors = P4EST_ALLOC(geometric_factors_t*, num_of_levels);

  for (int level = 0; level < num_of_levels-1; level++){
    updater->geometric_factors[level] = NULL;
  }
  
  updater->geometric_factors[num_of_levels-1] = geometric_factors_toplevel;

  if (element_data_init_user_fcn == NULL){
    mpi_abort("You must set the element data init user fcn for the element data updater\n");
  }

  if (d4est_geom == NULL){
    mpi_abort("You must set the d4est geometry for the element data updater\n");
  }


  if (geometric_factors_toplevel == NULL){
    mpi_abort("You must set the geometric_factors for the element data updater\n");
  }
  
  updater->user = user;
  return updater;
}

void
multigrid_element_data_updater_curved_destroy
(
 multigrid_element_data_updater_t* updater,
 int num_of_levels
)
{
  /* Do not free the toplevel geometric factors because these were sent in 
   * and may be used afterwards */
  
  for (int level = 0; level < num_of_levels-1; level++){
    geometric_factors_destroy(updater->geometric_factors[level]);
  }
  
  P4EST_FREE(updater->geometric_factors);
  P4EST_FREE(updater);
}

void
multigrid_element_data_updater_curved_update
(
 p4est_t* p4est,
 int level,
 problem_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;


    
  
  if (mg_data->mg_state == DOWNV_POST_COARSEN){
    curved_element_data_init_new(p4est,
                                 NULL,
                                 mg_data->dgmath_jit_dbase,
                                 updater->d4est_geom,
                                 updater->element_data_init_user_fcn,
                                 updater->user,
                                 0,
                                 0
                                );
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_geometric_factors = (updater->geometric_factors[level - 1] == NULL);
    if (compute_geometric_factors) {
      updater->geometric_factors[level-1] = geometric_factors_init(p4est);
    }
    curved_element_data_init_new(p4est,
                                 updater->geometric_factors[level - 1],
                                 mg_data->dgmath_jit_dbase,
                                 updater->d4est_geom,
                                 updater->element_data_init_user_fcn,
                                 updater->user,
                                 compute_geometric_factors,
                                 1
                                );
    

    /* update ghost data */
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (curved_element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count); 
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){    
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (curved_element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count); 
    curved_element_data_init_new(p4est,
                                 updater->geometric_factors[level + 1],
                                 mg_data->dgmath_jit_dbase,
                                 updater->d4est_geom,
                                 updater->element_data_init_user_fcn,
                                 updater->user,
                                 0, /* do not computer geometric factors */
                                 1 /* set aliases for geometric factors because they are stored */
                                );
    
  }
  else {
    return;
  }
}
