#include <pXest.h>
#include <d4est_solver_multigrid_smoother_schwarz.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_linalg.h>
#include <d4est_xyz_functions.h>

static void
d4est_solver_multigrid_smoother_schwarz
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int level
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = mg_data->smoother->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;

  for (int i = 0; i < smoother_data->schwarz_metadata_on_level[level]->schwarz_iter; i++){
    double* temp_Au = vecs->Au;
    vecs->Au = r;

    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       updater->current_d4est_ghost,
       updater->current_d4est_ghost_data,
       fcns,
       vecs,
       mg_data->d4est_ops,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->current_d4est_factors
      );

    d4est_linalg_vec_xpby(vecs->rhs, -1., r, vecs->local_nodes);    
    vecs->Au = temp_Au;

    if (p4est->mpisize != 1){
      D4EST_ABORT("We need to transfer the residual, to get the ghost portions");
    }
    
    d4est_solver_schwarz_compute_and_add_correction_version2
      (
       p4est,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->current_d4est_factors,
       updater->current_d4est_ghost,
       updater->current_d4est_ghost_data, /* should contain ghost data for residual only */
       smoother_data->schwarz_metadata_on_level[level],
       smoother_data->schwarz_ops,
       smoother_data->flux_data_for_lhs,
       smoother_data->laplacian_mortar_data_on_level[level],
       vecs->u,
       r
      );

    if (smoother_data->schwarz_metadata_on_level[level]->print_info){
      double r2 = d4est_linalg_vec_dot(r, r, vecs->local_nodes);  
      printf("[SCHWARZ_SMOOTHER_INFO] mg_level = %d, schwarz iter = %d, r = %.15f\n", level, i, sqrt(r2));
    }
  }

  double* temp_Au = vecs->Au;
  vecs->Au = r;

  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     updater->current_d4est_ghost,
     updater->current_d4est_ghost_data,
     fcns,
     vecs,
     mg_data->d4est_ops,
     mg_data->d4est_geom,
     mg_data->d4est_quad,
     updater->current_d4est_factors
    );

  d4est_linalg_vec_xpby(vecs->rhs, -1., r, vecs->local_nodes);    
  vecs->Au = temp_Au;  
}
  

d4est_solver_multigrid_smoother_t*
d4est_solver_multigrid_smoother_schwarz_init
(
 p4est_t* p4est,
 int num_of_levels,
 d4est_operators_t* d4est_ops,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 const char* input_file
)
{
  d4est_solver_multigrid_smoother_t* smoother = P4EST_ALLOC(d4est_solver_multigrid_smoother_t, 1); 
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = P4EST_ALLOC(d4est_solver_multigrid_smoother_schwarz_t,1);

  smoother_data->schwarz_metadata_on_level = P4EST_ALLOC(d4est_solver_schwarz_metadata_t*, num_of_levels); 
  smoother_data->laplacian_mortar_data_on_level = P4EST_ALLOC(d4est_solver_schwarz_laplacian_mortar_data_t*, num_of_levels); 
  
  for (int level = 0; level < num_of_levels; level++){
    smoother_data->schwarz_metadata_on_level[level] = NULL;
    smoother_data->laplacian_mortar_data_on_level[level] = NULL;
  }
  
  smoother_data->schwarz_metadata_on_level[num_of_levels-1]
    = d4est_solver_schwarz_metadata_init(p4est, d4est_ghost, input_file);

  smoother_data->laplacian_mortar_data_on_level[num_of_levels-1]
    = d4est_solver_schwarz_laplacian_mortar_data_init
    (
     p4est,
     smoother_data->schwarz_metadata_on_level[num_of_levels-1],
     d4est_factors
    );

  smoother_data->schwarz_ops = d4est_solver_schwarz_operators_init
    (d4est_ops);

  smoother_data->num_of_levels = num_of_levels;
  
  int str_len = strlen(input_file);
  int str_padding = 10;
  
  smoother_data->input_file = P4EST_ALLOC(char, str_len + str_padding);
  sprintf(smoother_data->input_file, "%s", input_file);

  dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO;  
  smoother_data->bc_data_dirichlet_for_lhs.dirichlet_fcn = zero_fcn;
  smoother_data->bc_data_dirichlet_for_lhs.eval_method = eval_method;
  smoother_data->bc_data_dirichlet_for_lhs.user = NULL;
  smoother_data->flux_data_for_lhs = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &smoother_data->bc_data_dirichlet_for_lhs);
  
  smoother->smooth = d4est_solver_multigrid_smoother_schwarz;
  smoother->update = d4est_solver_multigrid_smoother_schwarz_update;
  smoother->user = smoother_data;

 
  
  return smoother;
}


void
d4est_solver_multigrid_smoother_schwarz_destroy
(
 d4est_solver_multigrid_smoother_t* smoother
)
{
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = smoother->user;

  for (int level = 0;
       level < smoother_data->num_of_levels;
       level++){

    if(smoother_data->schwarz_metadata_on_level[level] != NULL){
      d4est_solver_schwarz_metadata_destroy
        (
         smoother_data->schwarz_metadata_on_level[level]
        );
    }
    if(smoother_data->laplacian_mortar_data_on_level[level] != NULL){
      d4est_solver_schwarz_laplacian_mortar_data_destroy
        (
         smoother_data->laplacian_mortar_data_on_level[level]
        );
    }
    
  }

  d4est_solver_schwarz_operators_destroy
    (
     smoother_data->schwarz_ops
    );


  
  d4est_laplacian_flux_destroy(smoother_data->flux_data_for_lhs);
  
  P4EST_FREE(smoother_data->input_file);
  P4EST_FREE(smoother_data);
  P4EST_FREE(smoother);
}

void
d4est_solver_multigrid_smoother_schwarz_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = mg_data->smoother->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  
  if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_solver_schwarz_data = (smoother_data->schwarz_metadata_on_level[level - 1] == NULL);
    if (compute_solver_schwarz_data) {
      smoother_data->schwarz_metadata_on_level[level-1] = d4est_solver_schwarz_metadata_init
                                                          (
                                                           p4est,
                                                           updater->current_d4est_ghost,
                                                           smoother_data->input_file
                                                          );
    }
    compute_solver_schwarz_data = (smoother_data->laplacian_mortar_data_on_level[level-1] == NULL);
    if (compute_solver_schwarz_data) {
      smoother_data->laplacian_mortar_data_on_level[level-1]
        = d4est_solver_schwarz_laplacian_mortar_data_init
        (
         p4est,
         smoother_data->schwarz_metadata_on_level[level-1],
         updater->current_d4est_factors
        );
    }    
  }
  
}
