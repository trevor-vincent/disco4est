#include <d4est_element_data.h>
#include <d4est_util.h>
#include <d4est_mesh.h>
#include <multigrid_element_data_updater.h>


multigrid_element_data_updater_t*
multigrid_element_data_updater_init
(
 int num_of_levels,
 d4est_ghost_t* d4est_ghost_toplevel,
 d4est_ghost_data_t* d4est_ghost_data_toplevel,
 d4est_mesh_data_t* d4est_mesh_data_toplevel,
 void(*element_data_init_user_fcn)(d4est_element_data_t*,void*),
 d4est_mesh_initial_extents_t* initial_extents
)
{
  multigrid_element_data_updater_t* updater = P4EST_ALLOC(multigrid_element_data_updater_t,1);
  updater->update = multigrid_element_data_updater_update;
  
  /* for curved infrastructure*/
  updater->element_data_init_user_fcn = element_data_init_user_fcn;
  updater->user = initial_extents;
  updater->d4est_factors_on_level = P4EST_ALLOC(d4est_mesh_data_t*, num_of_levels);
  updater->d4est_ghost_on_level = P4EST_ALLOC(d4est_ghost_t*, num_of_levels);
  updater->d4est_ghost_data_on_level = P4EST_ALLOC(d4est_ghost_data_t*, num_of_levels);

  for (int level = 0; level < num_of_levels-1; level++){
    updater->d4est_factors_on_level[level] = NULL;
    updater->d4est_ghost_on_level[level] = NULL;
    updater->d4est_ghost_data_on_level[level] = NULL;
  }
  
  updater->d4est_factors_on_level[num_of_levels-1] = d4est_mesh_data_toplevel;
  updater->d4est_ghost_on_level[num_of_levels-1] = d4est_ghost_toplevel;
  updater->d4est_ghost_data_on_level[num_of_levels-1] = d4est_ghost_data_toplevel;

  updater->current_d4est_factors = updater->d4est_factors_on_level[num_of_levels-1];
  updater->current_d4est_ghost = updater->d4est_ghost_on_level[num_of_levels-1];
  updater->current_d4est_ghost_data = updater->d4est_ghost_data_on_level[num_of_levels-1];
  updater->initial_extents = initial_extents;
  
  if (element_data_init_user_fcn == NULL){
    D4EST_ABORT("You must set the element data init user fcn for the element data updater\n");
  }

  if (d4est_mesh_data_toplevel == NULL){
    D4EST_ABORT("You must set the d4est_mesh_data_toplevel for the element data updater\n");
  }

  if (d4est_ghost_toplevel == NULL){
    D4EST_ABORT("You must set the d4est_ghost_toplevel for the element data updater\n");
  }
  
  if (d4est_ghost_data_toplevel == NULL){
    D4EST_ABORT("You must set the d4est_ghost_data_toplevel for the element data updater\n");
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
    d4est_mesh_data_destroy(updater->d4est_factors_on_level[level]);
    d4est_ghost_data_destroy(updater->d4est_ghost_data_on_level[level]);
    d4est_ghost_destroy(updater->d4est_ghost_on_level[level]);
  }

  updater->current_d4est_factors = NULL;
  updater->current_d4est_ghost = NULL;
  updater->current_d4est_ghost_data = NULL;
  
  P4EST_FREE(updater->d4est_factors_on_level);
  P4EST_FREE(updater->d4est_ghost_on_level);
  P4EST_FREE(updater->d4est_ghost_data_on_level);
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
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       NULL,
       updater->initial_extents,
       DO_NOT_INITIALIZE_GHOST,
       DO_NOT_INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
       updater->element_data_init_user_fcn,
       updater->user    
      );
  }
  else if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_geometric_factors = (updater->d4est_factors_on_level[level - 1] == NULL);
    if (compute_geometric_factors) {
      updater->d4est_factors_on_level[level-1] = d4est_mesh_data_init();

      d4est_mesh_update
        (
         p4est,
         &updater->d4est_ghost_on_level[level-1],
         mg_data->d4est_ops,
         mg_data->d4est_geom,
         mg_data->d4est_quad,
         updater->d4est_factors_on_level[level - 1],
         updater->initial_extents,
         INITIALIZE_GHOST,
         INITIALIZE_QUADRATURE_DATA,
         INITIALIZE_GEOMETRY_DATA,
         INITIALIZE_GEOMETRY_ALIASES,
         updater->element_data_init_user_fcn,
         updater->user    
        );

   
      updater->d4est_ghost_data_on_level[level-1] = d4est_ghost_data_init(p4est,
                                                                          updater->d4est_ghost_on_level[level-1],
                                                                          vecs->field_types,
                                                                          vecs->num_of_fields);
    }
    else {

      D4EST_ASSERT(updater->d4est_ghost_on_level[level-1] != NULL);
      d4est_mesh_update
        (
         p4est,
         &updater->d4est_ghost_on_level[level-1],
         mg_data->d4est_ops,
         mg_data->d4est_geom,
         mg_data->d4est_quad,
         updater->d4est_factors_on_level[level - 1],
         updater->initial_extents,
         DO_NOT_INITIALIZE_GHOST,
         DO_NOT_INITIALIZE_QUADRATURE_DATA,
         DO_NOT_INITIALIZE_GEOMETRY_DATA,
         DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
         /* INITIALIZE_GHOST, */
         /* INITIALIZE_QUADRATURE_DATA, */
         /* INITIALIZE_GEOMETRY_DATA, */
         /* INITIALIZE_GEOMETRY_ALIASES, */
         updater->element_data_init_user_fcn,
         updater->user    
        );
    }
    
    updater->current_d4est_ghost = updater->d4est_ghost_on_level[level-1];
    updater->current_d4est_ghost_data = updater->d4est_ghost_data_on_level[level-1];
    updater->current_d4est_factors = updater->d4est_factors_on_level[level-1];
  }
  else if (mg_data->mg_state == UPV_POST_REFINE){    

    d4est_mesh_update
      (
       p4est,
       &updater->d4est_ghost_on_level[level + 1],
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->d4est_factors_on_level[level + 1],
       updater->initial_extents,
       /* DO_NOT_INITIALIZE_GHOST, */
       /* DO_NOT_INITIALIZE_QUADRATURE_DATA, */
       /* DO_NOT_INITIALIZE_GEOMETRY_DATA, */
       /* DO_NOT_INITIALIZE_GEOMETRY_ALIASES, */
       DO_NOT_INITIALIZE_GHOST,
       DO_NOT_INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
       updater->element_data_init_user_fcn,
       updater->user    
      );

    updater->current_d4est_ghost = updater->d4est_ghost_on_level[level+1];
    updater->current_d4est_ghost_data = updater->d4est_ghost_data_on_level[level+1];
    updater->current_d4est_factors = updater->d4est_factors_on_level[level+1];       
  }
  else {
    return;
  }
}
