#define CG_VERBOSE

#include "../Solver/cg.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../pXest/pXest.h"
#include "sc_reduce.h"

void
cg_nr_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 double eta_max,
 double fnrm,
 int max_iter,
 int mpi_rank,
 int* final_iter,
 double* final_fnrm
)
{
  cg_solver_params_t cg_params;
  cg_params.rtol = 0.;
  cg_params.atol = eta_max*fnrm;
  cg_params.max_iter = max_iter;
  cg_params.mpi_rank = mpi_rank;

  cg_solve
    (
     p4est,
     vecs,
     fcns,
     dgmath_jit_dbase,
     ghost,
     ghost_data,
     &cg_params
    );

  *final_iter = cg_params.final_iter;
  *final_fnrm = cg_params.final_fnrm;
}

void cg_solve(p4est_t* p4est, problem_data_t* vecs, weakeqn_ptrs_t* fcns,
              dgmath_jit_dbase_t* dgbase, p4est_ghost_t* ghost,
              element_data_t* ghost_data, cg_solver_params_t* params) {

  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  double* Au;
  double* u;
  double* rhs;

  double* d;
  double* r;

  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  u = vecs->u;
  rhs = vecs->rhs;

  const int imax = params->max_iter;
  const double atol = params->atol;
  const double rtol = params->rtol;
  double d_dot_Au;

  d = P4EST_ALLOC(double, local_nodes);
  r = P4EST_ALLOC(double, local_nodes);

  /* first iteration data, store Au in r */
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgbase, NULL);


  
  linalg_copy_1st_to_2nd(Au, r, local_nodes);

  /* r = f - Au ; Au is stored in r so r = rhs - r */
  linalg_vec_xpby(rhs, -1., r, local_nodes);
  linalg_copy_1st_to_2nd(r, d, local_nodes);
  delta_new = linalg_vec_dot(r, r, local_nodes);

  double delta_new_global;
  double d_dot_Au_global;

  sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
               sc_MPI_COMM_WORLD);

  delta_new = delta_new_global;
  delta_0 = delta_new;

  /* start working on d */
  vecs->u = d;

  /* DEBUG_PRINT_2ARR_DBL(d, rhs, vecs->local_nodes); */
  /* printf("d sum = %.25f\n", linalg_vec_sum(d, vecs->local_nodes)); */
  
  int i;

  for (i = 0; i < imax && (delta_new > atol*atol + delta_0 * rtol * rtol); i++) {
    /* Au = A*d; */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgbase, NULL);
    d_dot_Au = linalg_vec_dot(d, Au, local_nodes);

    sc_allreduce(&d_dot_Au, &d_dot_Au_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);

    d_dot_Au = d_dot_Au_global;
    alpha = delta_new / d_dot_Au;

    linalg_vec_axpy(alpha, d, u, local_nodes);

    /* r = r - Au*alpha */
    linalg_vec_axpy(-alpha, Au, r, local_nodes);

    delta_old = delta_new;
    delta_new = linalg_vec_dot(r, r, local_nodes);

    sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);
    delta_new = delta_new_global;

    beta = delta_new / delta_old;
    linalg_vec_xpby(r, beta, d, local_nodes);

    if (params->mpi_rank == 0 && params->monitor == 1)
      printf("[CG_SOLVER]: ITER %03d FNRMSQR %.30f\n", i, delta_new);
    
  }

  params->final_iter = i;
  params->final_fnrm = sqrt(delta_new);

  /* set back from d */
  vecs->u = u;

  P4EST_FREE(d);
  P4EST_FREE(r);
}

void curved_cg_solve(
                     p4est_t* p4est,
                     problem_data_t* vecs,
                     weakeqn_ptrs_t* fcns,
                     dgmath_jit_dbase_t* dgbase,
                     d4est_geometry_t* geom,
                     p4est_ghost_t* ghost,
                     curved_element_data_t* ghost_data,
                     cg_solver_params_t* params) {

  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  double* Au;
  double* u;
  double* rhs;

  double* d;
  double* r;

  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  u = vecs->u;
  rhs = vecs->rhs;

  const int imax = params->max_iter;
  const double atol = params->atol;
  const double rtol = params->rtol;
  double d_dot_Au;

  d = P4EST_ALLOC(double, local_nodes);
  r = P4EST_ALLOC(double, local_nodes);

  /* first iteration data, store Au in r */
  fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgbase, geom);

  /* DEBUG_PRINT_2ARR_DBL(vecs->u, vecs->Au, vecs->local_nodes); */
  
  linalg_copy_1st_to_2nd(Au, r, local_nodes);

  /* r = f - Au ; Au is stored in r so r = rhs - r */
  linalg_vec_xpby(rhs, -1., r, local_nodes);
  linalg_copy_1st_to_2nd(r, d, local_nodes);
  delta_new = linalg_vec_dot(r, r, local_nodes);

  double delta_new_global;
  double d_dot_Au_global;

  sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
               sc_MPI_COMM_WORLD);

  delta_new = delta_new_global;
  delta_0 = delta_new;

  /* start working on d */
  vecs->u = d;
  /* DEBUG_PRINT_2ARR_DBL(d, rhs, vecs->local_nodes); */
  /* printf("d sum = %.25f\n", linalg_vec_sum(d, vecs->local_nodes)); */
  
  int i;

  for (i = 0; i < imax && (delta_new > atol*atol + delta_0 * rtol * rtol); i++) {
    /* Au = A*d; */
    fcns->apply_lhs(p4est, ghost, ghost_data, vecs, dgbase, geom);
    d_dot_Au = linalg_vec_dot(d, Au, local_nodes);

    sc_allreduce(&d_dot_Au, &d_dot_Au_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);

    d_dot_Au = d_dot_Au_global;
    alpha = delta_new / d_dot_Au;

    linalg_vec_axpy(alpha, d, u, local_nodes);

    /* r = r - Au*alpha */
    linalg_vec_axpy(-alpha, Au, r, local_nodes);

    delta_old = delta_new;
    delta_new = linalg_vec_dot(r, r, local_nodes);

    sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);
    delta_new = delta_new_global;

    beta = delta_new / delta_old;
    linalg_vec_xpby(r, beta, d, local_nodes);

    if (params->mpi_rank == 0 && params->monitor == 1)
      printf("[CG_SOLVER]: ITER %03d FNRMSQR %.30f\n", i, delta_new);
    
  }

  params->final_iter = i;
  params->final_fnrm = sqrt(delta_new);

  /* set back from d */
  vecs->u = u;

  P4EST_FREE(d);
  P4EST_FREE(r);
}
