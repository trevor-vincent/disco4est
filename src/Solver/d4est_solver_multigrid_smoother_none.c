#include <pXest.h>
#include <d4est_solver_multigrid_smoother_none.h>
#include <d4est_linalg.h>
#include <ini.h>
#include <d4est_util.h>

void
d4est_solver_multigrid_smoother_none_destroy(d4est_solver_multigrid_smoother_t* smoother)
{
  P4EST_FREE(smoother);
}

static void
d4est_solver_multigrid_smoother_none
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r,
 int level
)
{  
  d4est_solver_multigrid_t* mg_data = p4est->user_pointer;
  /* d4est_solver_multigrid_smoother_none_t* none = mg_data->smoother->user; */
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  d4est_ghost_t* d4est_ghost = updater->current_d4est_ghost;
  d4est_ghost_data_t* d4est_ghost_data = updater->current_d4est_ghost_data;
  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     fcns,
     vecs,
     mg_data->d4est_ops,
     mg_data->d4est_geom,
     mg_data->d4est_quad,
     updater->current_d4est_factors
    );
  
  d4est_util_copy_1st_to_2nd(vecs->Au, r, vecs->local_nodes);
  d4est_linalg_vec_xpby(vecs->rhs, -1., r, vecs->local_nodes);

}



d4est_solver_multigrid_smoother_t*
d4est_solver_multigrid_smoother_none_init
(
 p4est_t* p4est,
 int num_of_levels,
 const char* input_file
)
{
  d4est_solver_multigrid_smoother_t* smoother = P4EST_ALLOC(d4est_solver_multigrid_smoother_t, 1);

  smoother->user = NULL;
  smoother->smooth = d4est_solver_multigrid_smoother_none;
  smoother->update = NULL;

  return smoother;
}


