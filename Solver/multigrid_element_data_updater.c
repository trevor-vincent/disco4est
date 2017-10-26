#include <d4est_element_data.h>
#include <d4est_util.h>
#include <d4est_mesh.h>
#include <multigrid_element_data_updater.h>


multigrid_element_data_updater_t*
multigrid_element_data_updater_init
(
 int num_of_levels,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_mesh_geometry_storage_t* d4est_mesh_geometry_storage_toplevel,
 void(*element_data_init_user_fcn)(d4est_element_data_t*,void*),
 void* user
)
{
  multigrid_element_data_updater_t* updater = P4EST_ALLOC(multigrid_element_data_updater_t,1);
  updater->update = multigrid_element_data_updater_update;
  updater->ghost = ghost;
  updater->ghost_data = ghost_data;
  
  /* for curved infrastructure*/
  updater->element_data_init_user_fcn = element_data_init_user_fcn;
  updater->user = user;
  updater->geometric_factors = P4EST_ALLOC(d4est_mesh_geometry_storage_t*, num_of_levels);

  for (int level = 0; level < num_of_levels-1; level++){
    updater->geometric_factors[level] = NULL;
  }
  
  updater->geometric_factors[num_of_levels-1] = d4est_mesh_geometry_storage_toplevel;

  updater->current_geometric_factors = updater->geometric_factors[num_of_levels-1];
  
  if (element_data_init_user_fcn == NULL){
    D4EST_ABORT("You must set the element data init user fcn for the element data updater\n");
  }

  if (d4est_mesh_geometry_storage_toplevel == NULL){
    D4EST_ABORT("You must set the geometric_factors for the element data updater\n");
  }
  
  return updater;
}

void
multigrid_element_data_updater_destroy
(
 multigrid_element_data_updater_t* updater,
 int num_of_levels
)
{
  /* Do not free the toplevel geometric factors because these were sent in 
   * and may be used afterwards */
  
  for (int level = 0; level < num_of_levels-1; level++){
    d4est_mesh_geometry_storage_destroy(updater->geometric_factors[level]);
  }

  updater->current_geometric_factors = NULL;
  
  P4EST_FREE(updater->geometric_factors);
  P4EST_FREE(updater);
}

void
multigrid_element_data_updater_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;

  if (mg_data->mg_state == DOWNV_POST_COARSEN){
    d4est_mesh_update
      (
       p4est,
       NULL,
       NULL,
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       NULL,
       DO_NOT_INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
       updater->element_data_init_user_fcn,
       updater->user    
      );
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_geometric_factors = (updater->geometric_factors[level - 1] == NULL);
    if (compute_geometric_factors) {
      updater->geometric_factors[level-1] = d4est_mesh_geometry_storage_init();
    }

   /* update ghost data */
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (d4est_element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count); 

    d4est_mesh_update
      (
       p4est,
       *(updater->ghost),
       *(updater->ghost_data),
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->geometric_factors[level - 1],
       INITIALIZE_QUADRATURE_DATA,
       INITIALIZE_GEOMETRY_DATA,
       INITIALIZE_GEOMETRY_ALIASES,
       updater->element_data_init_user_fcn,
       updater->user    
      );

    updater->current_geometric_factors = updater->geometric_factors[level-1];
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){    
    P4EST_FREE(*(updater->ghost_data));
    p4est_ghost_destroy (*(updater->ghost));
    *(updater->ghost) = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
    *(updater->ghost_data) = P4EST_ALLOC (d4est_element_data_t,
                                          (*(updater->ghost))->ghosts.elem_count);



    d4est_mesh_update
      (
       p4est,
       *(updater->ghost),
       *(updater->ghost_data),
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->geometric_factors[level + 1],
       INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       INITIALIZE_GEOMETRY_ALIASES,
       updater->element_data_init_user_fcn,
       updater->user    
      );

    updater->current_geometric_factors = updater->geometric_factors[level+1];
       
  }
  else {
    return;
  }
}
