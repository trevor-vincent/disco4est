#include <d4est_solver_multigrid_bottom_solver_cg.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_linalg.h>
#include <sc_reduce.h>

static int
d4est_solver_multigrid_bottom_solver_cg_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_multigrid_bottom_solver_cg_t* pconfig = ((d4est_solver_multigrid_bottom_solver_cg_t*)user);
  
  if (d4est_util_match_couple(section,"mg_bottom_solver_cg",name,"bottom_iter")) {
    D4EST_ASSERT(pconfig->bottom_imax == -1);
    pconfig->bottom_imax = atoi(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cg",name,"bottom_rtol")) {
    D4EST_ASSERT(pconfig->bottom_rtol == -1);
    pconfig->bottom_rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cg",name,"bottom_atol")) {
    D4EST_ASSERT(pconfig->bottom_atol == -1);
    pconfig->bottom_atol = atof(value);
  }
  else if (d4est_util_match_couple(section,"mg_bottom_solver_cg",name,"bottom_print_residual_norm")) {
    D4EST_ASSERT(pconfig->bottom_print_residual_norm == -1);
    pconfig->bottom_print_residual_norm = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
d4est_solver_multigrid_bottom_solver_cg_destroy(d4est_solver_multigrid_bottom_solver_t* bottom)
{
  P4EST_FREE(bottom->user);
  P4EST_FREE(bottom);
}

static void 
d4est_solver_multigrid_bottom_solver_cg
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 double* r
)
{
  d4est_solver_multigrid_data_t* mg_data = p4est->user_pointer;
  d4est_solver_multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  zlog_category_t *c_default = zlog_get_category("d4est_solver_multigrid_bottom_solver_cg");   
  d4est_ghost_t* d4est_ghost = updater->current_d4est_ghost;
  d4est_ghost_data_t* d4est_ghost_data = updater->current_d4est_ghost_data;
  d4est_geometry_t* d4est_geom = mg_data->d4est_geom;
  
  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  d4est_solver_multigrid_bottom_solver_cg_t* cg_params = mg_data->bottom_solver->user;
  d4est_operators_t* d4est_ops = mg_data->d4est_ops;

  double* Au; 
  double* u;
  double* rhs;

  double* d;
  
  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  u = vecs->u;
  rhs = vecs->rhs;

  const int imax = cg_params->bottom_imax;
  const double rtol = cg_params->bottom_rtol;
  const double atol = cg_params->bottom_atol;
  double d_dot_Au;
  
  d = P4EST_ALLOC(double, local_nodes);
  
  /* first iteration data, store Au in r */ 
  /* fcns->apply_lhs(p4est, ghost, ghost_data, vecs, d4est_ops, d4est_geom); */

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
  
  d4est_util_copy_1st_to_2nd(Au, r, local_nodes);
  
  /* r = f - Au ; Au is stored in r so r = rhs - r */
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);

  d4est_util_copy_1st_to_2nd(r, d, local_nodes);
 
  delta_new = d4est_linalg_vec_dot(r,r,local_nodes);
  /* delta_new = (element_data_compute_l2_norm(p4est, r)); */
  double delta_new_global;
  sc_allreduce
    (
     &delta_new,
     &delta_new_global,
     1,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     sc_MPI_COMM_WORLD
    );

  delta_new = delta_new_global;
  delta_0 = delta_new;

  /* start working on d */
  vecs->u = d;
  int i;
  for (i = 0; i < imax && (delta_new > atol*atol + delta_0 * rtol * rtol); i++){
  
    /* Au = A*d; */
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
       mg_data->elem_data_updater->current_d4est_factors
      );
  
    
    /* sc_MPI_Barrier(sc_MPI_COMM_WORLD); */

    d_dot_Au = d4est_linalg_vec_dot(d,Au,local_nodes);
    double d_dot_Au_global;
    
    sc_allreduce
      (
       &d_dot_Au,
       &d_dot_Au_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
    );
    d_dot_Au = d_dot_Au_global;
    alpha = delta_new/d_dot_Au;
    d4est_linalg_vec_axpy(alpha, d, u, local_nodes);

    /* r = r - Au*alpha */
    d4est_linalg_vec_axpy(-alpha, Au, r, local_nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, local_nodes);
    
    sc_allreduce
      (
       &delta_new,
       &delta_new_global,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
    );

    delta_new = delta_new_global;
    
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, local_nodes);

    if (p4est->mpirank == 0 && cg_params->bottom_print_residual_norm == 1){
      zlog_info(c_default, "iter %03d rnrmsqr %.25f", i, delta_new);
    }
  }
  
  vecs->u = u;
  P4EST_FREE(d);
}

d4est_solver_multigrid_bottom_solver_t*
d4est_solver_multigrid_bottom_solver_cg_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  d4est_solver_multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_t, 1);
  d4est_solver_multigrid_bottom_solver_cg_t* bottom_data = P4EST_ALLOC(d4est_solver_multigrid_bottom_solver_cg_t, 1);
  
  bottom_data->bottom_atol = -1;
  bottom_data->bottom_rtol = -1;
  bottom_data->bottom_imax = -1;
  bottom_data->bottom_print_residual_norm = -1;
  
  if (ini_parse(input_file, d4est_solver_multigrid_bottom_solver_cg_input_handler, bottom_data) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("mg_bottom_solver_cg", bottom_data->bottom_atol, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cg", bottom_data->bottom_rtol, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cg", bottom_data->bottom_imax, -1);
  D4EST_CHECK_INPUT("mg_bottom_solver_cg", bottom_data->bottom_print_residual_norm, -1);
  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: d4est_bottom_solver_cg parameters\n");
    printf("[D4EST_INFO]: bottom imax = %d\n", bottom_data->bottom_imax);
    printf("[D4EST_INFO]: bottom rtol = %.25f\n", bottom_data->bottom_rtol);
    printf("[D4EST_INFO]: bottom atol = %.25f\n", bottom_data->bottom_atol);
    printf("[D4EST_INFO]: bottom atol = %d\n", bottom_data->bottom_print_residual_norm);
  }

  bottom_solver->user = bottom_data;
  bottom_solver->solve = d4est_solver_multigrid_bottom_solver_cg;
  bottom_solver->update = NULL;

  return bottom_solver;
}
