#include <multigrid_bottom_solver_cg.h>
#include <util.h>
#include <ini.h>
#include <linalg.h>
#include <sc_reduce.h>

static int
multigrid_bottom_solver_cg_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  multigrid_bottom_solver_cg_t* pconfig = ((multigrid_bottom_solver_cg_t*)user);
  
  if (util_match_couple(section,"multigrid",name,"bottom_iter")) {
    mpi_assert(pconfig->bottom_imax == -1);
    pconfig->bottom_imax = atoi(value);
  }
  else if (util_match_couple(section,"multigrid",name,"bottom_rtol")) {
    mpi_assert(pconfig->bottom_rtol == -1);
    pconfig->bottom_rtol = atof(value);
  }
  else if (util_match_couple(section,"multigrid",name,"bottom_atol")) {
    mpi_assert(pconfig->bottom_atol == -1);
    pconfig->bottom_atol = atof(value);
  }
  else if (util_match_couple(section,"multigrid",name,"bottom_print_residual_norm")) {
    mpi_assert(pconfig->bottom_print_residual_norm == -1);
    pconfig->bottom_print_residual_norm = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


void
multigrid_bottom_solver_cg_destroy(multigrid_bottom_solver_t* bottom)
{
  P4EST_FREE(bottom->user);
  P4EST_FREE(bottom);
}

static void 
multigrid_bottom_solver_cg
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 double* r
)
{
  multigrid_data_t* mg_data = p4est->user_pointer;
  multigrid_element_data_updater_t* updater = mg_data->elem_data_updater;
  p4est_ghost_t* ghost = *(updater->ghost);
  void* ghost_data = *(updater->ghost_data);
  d4est_geometry_t* d4est_geom = updater->d4est_geom;
  
  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  multigrid_bottom_solver_cg_t* cg_params = mg_data->bottom_solver->user;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;

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
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase, d4est_geom);

  linalg_copy_1st_to_2nd(Au, r, local_nodes);
  
  /* r = f - Au ; Au is stored in r so r = rhs - r */
  linalg_vec_xpby(rhs, -1., r, local_nodes);

  linalg_copy_1st_to_2nd(r, d, local_nodes);
 
  delta_new = linalg_vec_dot(r,r,local_nodes);
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
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase, d4est_geom);
    /* sc_MPI_Barrier(sc_MPI_COMM_WORLD); */

    d_dot_Au = linalg_vec_dot(d,Au,local_nodes);
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
    linalg_vec_axpy(alpha, d, u, local_nodes);

    /* r = r - Au*alpha */
    linalg_vec_axpy(-alpha, Au, r, local_nodes);

    delta_old = delta_new;
    delta_new = linalg_vec_dot(r, r, local_nodes);
    
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
    linalg_vec_xpby(r, beta, d, local_nodes);

    if (p4est->mpirank == 0 && cg_params->bottom_print_residual_norm == 1){
      printf("[CG_BOTTOM_SOLVER]: CG ITER %03d RNRMSQR %.10f\n", i, delta_new);
    }
  }
  
  vecs->u = u;
  P4EST_FREE(d);
}

multigrid_bottom_solver_t*
multigrid_bottom_solver_cg_init
(
 p4est_t* p4est,
 const char* input_file
)
{
  multigrid_bottom_solver_t* bottom_solver = P4EST_ALLOC(multigrid_bottom_solver_t, 1);
  multigrid_bottom_solver_cg_t* bottom_data = P4EST_ALLOC(multigrid_bottom_solver_cg_t, 1);
  
  bottom_data->bottom_atol = -1;
  bottom_data->bottom_rtol = -1;
  bottom_data->bottom_imax = -1;
  bottom_data->bottom_print_residual_norm = -1;
  
  if (ini_parse(input_file, multigrid_bottom_solver_cg_input_handler, bottom_data) < 0) {
    mpi_abort("Can't load input file");
  }
  if(bottom_data->bottom_atol == -1){
    mpi_abort("[D4EST_ERROR]: bottom_atol not set in multigrid input");
  }
  if(bottom_data->bottom_rtol == -1){
    mpi_abort("[D4EST_ERROR]: bottom_rtol not set in multigrid input");
  }
  if(bottom_data->bottom_imax == -1){
    mpi_abort("[D4EST_ERROR]: bottom_imax not set in multigrid input");
  }  
  if(p4est->mpirank == 0){
    printf("[D4EST_INFO]: Multigrid_Bottom_Solver_CG Parameters\n");
    printf("[D4EST_INFO]: bottom imax = %d\n", bottom_data->bottom_imax);
    printf("[D4EST_INFO]: bottom rtol = %.25f\n", bottom_data->bottom_rtol);
    printf("[D4EST_INFO]: bottom atol = %.25f\n", bottom_data->bottom_atol);
  }

  bottom_solver->user = bottom_data;
  bottom_solver->solve = multigrid_bottom_solver_cg;
  bottom_solver->update = NULL;

  return bottom_solver;
}
