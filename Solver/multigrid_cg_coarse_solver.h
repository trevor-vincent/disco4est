#ifndef MULTIGRID_CG_COARSE_SOLVER_H
#define MULTIGRID_CG_COARSE_SOLVER_H 

#define NDEBUG
/* #define VERBOSE */

static void 
multigrid_cg_coarse_solver
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double* r,
 multigrid_cg_params_t* cg_params
)
{
  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;

  double* Au; 
  double* u;
  double* rhs;

  double* d;
  
  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  u = vecs->u;
  rhs = vecs->rhs;

  const int imax = cg_params->iter;
  const double rtol = cg_params->rtol;
  double d_dot_Au;
  
  d = P4EST_ALLOC(double, local_nodes);
  
  /* first iteration data, store Au in r */ 
  /* debug("First Au calculation in CG solve starts"); */
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
  /* debug("First Au calculation in CG solve ends"); */

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
  for (i = 0; i < imax && (delta_new > delta_0 * rtol * rtol); i++){
  
    /* Au = A*d; */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgmath_jit_dbase);
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

#ifdef VERBOSE
    printf("[CG_BOTTOM_SOLVER]: %d CG ITER %03d RNRMSQR %.10f\n",
                   mg_data->mpi_rank, i, delta_new);
#endif
  }
  
  vecs->u = u;
  P4EST_FREE(d);
}

#endif
