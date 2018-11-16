static
int
d4est_solver_multigrid_smoother_cheby_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{  
  d4est_solver_multigrid_smoother_cheby_t* pconfig = ((d4est_solver_multigrid_smoother_cheby_t*)user);

  if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_imax")) {
    D4EST_ASSERT(pconfig->cheby_imax == -1);
    pconfig->cheby_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_cg_imax")) {
    D4EST_ASSERT(pconfig->cheby_eigs_cg_imax == -1);
    pconfig->cheby_eigs_cg_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_lmax_lmin_ratio")) {
    D4EST_ASSERT(pconfig->cheby_eigs_lmax_lmin_ratio == -1);
    pconfig->cheby_eigs_lmax_lmin_ratio = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_max_multiplier")) {
    D4EST_ASSERT(pconfig->cheby_eigs_max_multiplier == -1);
    pconfig->cheby_eigs_max_multiplier = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_reuse_fromdownvcycle")) {
    D4EST_ASSERT(pconfig->cheby_eigs_reuse_fromdownvcycle == -1);
    pconfig->cheby_eigs_reuse_fromdownvcycle = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_eigs_reuse_fromlastvcycle")) {
    D4EST_ASSERT(pconfig->cheby_eigs_reuse_fromlastvcycle == -1);
    pconfig->cheby_eigs_reuse_fromlastvcycle = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_residual_norm")) {
    D4EST_ASSERT(pconfig->cheby_print_residual_norm == 0);
    pconfig->cheby_print_residual_norm = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_spectral_bound")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound == 0);
    pconfig->cheby_print_spectral_bound = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_print_spectral_bound_iterations")) {
    D4EST_ASSERT(pconfig->cheby_print_spectral_bound_iterations == 0);
    pconfig->cheby_print_spectral_bound_iterations = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_smoother_cheby",name,"cheby_use_new_cg_eigs")) {
    D4EST_ASSERT(pconfig->cheby_use_new_cg_eigs == 0);
    pconfig->cheby_use_new_cg_eigs = atoi(value);
  }
  
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


d4est_solver_multigrid_schwarz_smoother_t*
d4est_solver_multigrid_schwarz_smoother_init
(
 int num_of_levels,
 d4est_operators_t* d4est_ops,
 d4est_ghost_t* d4est_ghost,
 const char* input_file
)
{

  d4est_solver_multigrid_smoother_t* smoother = P4EST_ALLOC(d4est_solver_multigrid_smoother_t, 1);
  
  d4est_solver_multigrid_schwarz_smoother_t* smoother_data = P4EST_ALLOC(d4est_solver_multigrid_element_data_updater_t,1);

  smoother->update = d4est_solver_multigrid_element_data_updater_update;

  smoother_data->schwarz_metadata_on_level = P4EST_ALLOC(d4est_solver_schwarz_metadata_t*, num_of_levels); 
  
  for (int level = 0; level < num_of_levels; level++){
    smoother->schwarz_metadata_on_level[level] = NULL;
  }
  
  updater->solver_schwarz_data_on_level[num_of_levels-1] = d4est_solver_schwarz_metadata_toplevel;
  

  if (d4est_solver_schwarz_metadata_toplevel == NULL){
    D4EST_ABORT("You must set the d4est_solver_schwarz_metadata_toplevel for the element data updater\n");
  }

  if (d4est_ghost_toplevel == NULL){
    D4EST_ABORT("You must set the d4est_ghost_toplevel for the element data updater\n");
  }

  return smoother;
}


void
d4est_solver_multigrid_smoother_schwarz_destroy(d4est_solver_multigrid_smoother_t* smoother)
{
 for (int level = 0; level < num_of_levels; level++){ d4est_solver_schwarz_data_destroy(updater->d4est_factors_on_level[level]);
  }

  //d4est_destroy schwarz ops
  
  d4est_solver_multigrid_smoother_cheby_t* cheby = smoother->user;
  P4EST_FREE(cheby->eigs);
  P4EST_FREE(cheby);
  P4EST_FREE(smoother);
}

void
d4est_solver_multigrid_element_data_updater_update
(
 p4est_t* p4est,
 int level,
 d4est_elliptic_data_t* vecs
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;

  if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_solver_schwarz_data = (updater->d4est_solver_schwarz_data_on_level[level - 1] == NULL);
    if (compute_geometric_factors) {
      updater->d4est_solver_schwarz_data_on_level[level-1] = d4est_mesh_data_init();
    }
    else {
      D4EST_ASSERT(updater->d4est_ghost_on_level[level-1] != NULL);
    }    
  }
}
