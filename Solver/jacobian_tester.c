#include <pXest.h>
#include <jacobian_tester.h>
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"

void
jacobian_tester
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 curved_weakeqn_ptrs_t* prob_fcns,
 problem_data_t* prob_vecs
)
{
  int local_nodes = prob_vecs->local_nodes;

  problem_data_t prob_vecs_jac;
  problem_data_t prob_vecs_F_u0;
  problem_data_t prob_vecs_F_u0_p_u;
  problem_data_copy_ptrs(prob_vecs,&prob_vecs_jac);
  problem_data_copy_ptrs(prob_vecs,&prob_vecs_F_u0);
  problem_data_copy_ptrs(prob_vecs,&prob_vecs_F_u0_p_u);

  
  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_p_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* J_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0 = P4EST_ALLOC(double, local_nodes);
  double* F_u0_p_u = P4EST_ALLOC(double, local_nodes);
  double* J_u_from_res = P4EST_ALLOC(double, local_nodes);

  prob_vecs_jac.u = u;
  prob_vecs_jac.u0 = u0;
  prob_vecs_jac.Au = J_u;

  prob_vecs_F_u0.u = u0;
  prob_vecs_F_u0.Au = F_u0; /* doesn't matter, not used */

  prob_vecs_F_u0_p_u.u = u0_p_u;
  prob_vecs_F_u0_p_u.Au = F_u0_p_u; /* doesn't matter, not used */
  
  int num_vecs_to_try = 5;
  double eps = .0001;
  linalg_fill_vec(u0, 0., local_nodes);
  mpi_assert(num_vecs_to_try <= local_nodes);

  double max_err = -1;
  int max_err_i = -1;
  int max_err_j = -1;
  
  for (int i = 0; i < num_vecs_to_try; i++){
    u[i] = 1.;
    prob_fcns->apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       &prob_vecs_jac,
       dgmath_jit_dbase,
       d4est_geom
      );

    /* double sum1 = linalg_vec_sum(J_u,local_nodes); */
    linalg_vec_axpyeqz(eps, u, u0, u0_p_u, local_nodes);
    
    prob_fcns->build_residual
      (
       p4est,
       ghost,
       ghost_data,
       &prob_vecs_F_u0,
       dgmath_jit_dbase,
       d4est_geom
      );

    prob_fcns->build_residual
      (
       p4est,
       ghost,
       ghost_data,
       &prob_vecs_F_u0_p_u,
       dgmath_jit_dbase,
       d4est_geom
      );

    for (int j = 0; j < local_nodes; j++){
      J_u_from_res[j] = (F_u0_p_u[j] - F_u0[j])/eps;
      double err = fabs(J_u_from_res[j] - J_u[j]);
      max_err_i = (err > max_err) ? i : max_err_i;
      max_err_j = (err > max_err) ? j : max_err_j;   
      max_err = (err > max_err) ? err : max_err;
    }

    /* double sum2 = linalg_vec_sum(J_u_from_res,local_nodes);     */

    /* printf("J_u sum = %.25f\n", sum1); */
    /* printf("J_u_from_res sum = %.25f\n", sum2); */

    
    u[i] = 0.;
  }

  printf("max_err, i, j = %.25f, %d, %d\n", max_err, max_err_i, max_err_j);
  
  P4EST_FREE(u);
  P4EST_FREE(u0);
  P4EST_FREE(u0_p_u);
  P4EST_FREE(J_u);
  P4EST_FREE(F_u0);
  P4EST_FREE(F_u0_p_u);
  P4EST_FREE(J_u_from_res);
  
}

