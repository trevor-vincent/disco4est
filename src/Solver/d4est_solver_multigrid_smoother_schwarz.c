#include <pXest.h>
#include <d4est_solver_multigrid_smoother_schwarz.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_schwarz_helpers.h>
#include <d4est_solver_schwarz.h>
#include <d4est_linalg.h>
#include <d4est_xyz_functions.h>
#include <ini.h>



static
int d4est_solver_multigrid_smoother_schwarz_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_multigrid_smoother_schwarz_t* pconfig = (d4est_solver_multigrid_smoother_schwarz_t*)user;
  const char* input_section = "mg_smoother_schwarz";
  
  if (
      d4est_util_match_couple(section,input_section,name,"smoother_iterations")) {
    D4EST_ASSERT(pconfig->iterations == -1);
    pconfig->iterations = atoi(value);
  }
  else if (
 d4est_util_match_couple(section,input_section,name,"smoother_verbose")) {
    D4EST_ASSERT(pconfig->verbose == -1);
    pconfig->verbose = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_multigrid_smoother_schwarz_set_apply_lhs
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_apply_lhs_t* apply_lhs,
 d4est_solver_multigrid_t* mg_data,
 const char* input_file
){
  
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = mg_data->smoother->user;

  smoother_data->apply_lhs = apply_lhs;
  
  if (smoother_data->schwarz_on_level[mg_data->num_of_levels - 1] != NULL){
    D4EST_ABORT("smoother_data->schwarz[mg_data->num_of_levels - 1] needs to be NULL before you call set_apply_lhs");
  }
  
  smoother_data->schwarz_on_level[mg_data->num_of_levels - 1] =
    d4est_solver_schwarz_init
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_ghost,
     d4est_factors,
     NULL,
     apply_lhs,
     input_file,
     "mg_smoother_schwarz"
    );

  smoother_data->apply_lhs_is_set = 1;
}


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
  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = mg_data->smoother->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;

  if(smoother_data->apply_lhs_is_set == 0){
    D4EST_ABORT("Please call d4est_solver_multigrid_smoother_schwarz_set_apply_lhs before solving");
  }
  
  for (int i = 0; i < smoother_data->iterations; i++){

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


    d4est_solver_schwarz_iterate
      (
       p4est,
       mg_data->d4est_geom,
       mg_data->d4est_quad,
       updater->current_d4est_factors,
       updater->current_d4est_ghost,
       smoother_data->schwarz_on_level[level],
       vecs->u,
       r
      );

    
    if (smoother_data->verbose){
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

  smoother_data->schwarz_on_level = P4EST_ALLOC(d4est_solver_schwarz_t*, num_of_levels); 

  smoother_data->iterations = -1;
  smoother_data->verbose = -1;


  if(
     ini_parse(input_file,
               d4est_solver_multigrid_smoother_schwarz_input_handler,
               smoother_data) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  
  
  for (int level = 0; level < num_of_levels; level++){
    smoother_data->schwarz_on_level[level] = NULL;
  }

  smoother_data->apply_lhs_is_set = 0;

  smoother_data->schwarz_ops = d4est_solver_schwarz_operators_init
    (d4est_ops);

  smoother_data->num_of_levels = num_of_levels;
  
  int str_len = strlen(input_file);
  int str_padding = 10;
  
  smoother_data->input_file = P4EST_ALLOC(char, str_len + str_padding);
  sprintf(smoother_data->input_file, "%s", input_file);
  
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

    if(smoother_data->schwarz_on_level[level] != NULL){
      d4est_solver_schwarz_destroy
        (
         smoother_data->schwarz_on_level[level]
        );
    }    
  }

  d4est_solver_schwarz_operators_destroy
    (
     smoother_data->schwarz_ops
    );

  /* d4est_solver_schwarz_apply_lhsdestroy(smoother_data->apply_lhs); */
  
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
  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_smoother_schwarz_t* smoother_data = mg_data->smoother->user;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  
  if (mg_data->mg_state == DOWNV_POST_BALANCE){
    int compute_solver_schwarz_data = (smoother_data->schwarz_on_level[level - 1] == NULL);
    if (compute_solver_schwarz_data) {

      if(smoother_data->apply_lhs_is_set == 0){
        D4EST_ABORT("apply lhs is not set");
      }
      
      if (p4est->mpirank == 0){
        zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_smoother_schwarz");    
        zlog_debug(c_default, "Computing schwarz data for smoother on level %d", level);
      }
      smoother_data->schwarz_on_level[level - 1] =
        d4est_solver_schwarz_init
        (
         p4est,
         mg_data->d4est_ops,
         mg_data->d4est_geom,
         mg_data->d4est_quad,
         updater->current_d4est_ghost,
         updater->current_d4est_factors,
         NULL,
         smoother_data->apply_lhs,
         smoother_data->input_file,
         "mg_smoother_schwarz"
        );
    }    
  }
}
