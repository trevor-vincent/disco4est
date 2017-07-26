#include <pXest.h>
#include <d4est_solver_analytic_jacobian_tester.h>
#include <d4est_linalg.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_elliptic_data.h>

/* @article{an2011finite, */
/*   title={On finite difference approximation of a matrix-vector product in the Jacobian-free Newton--Krylov method}, */
/*   author={An, Heng-Bin and Wen, Ju and Feng, Tao}, */
/*   journal={Journal of Computational and Applied Mathematics}, */
/*   volume={236}, */
/*   number={6}, */
/*   pages={1399--1409}, */
/*   year={2011}, */
/*   publisher={Elsevier} */
/* } */

void
d4est_solver_analytic_jacobian_tester_forward_difference
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_vecs
)
{
  int local_nodes = elliptic_vecs->local_nodes;

  d4est_elliptic_data_t elliptic_vecs_jac;
  d4est_elliptic_data_t elliptic_vecs_F_u0;
  d4est_elliptic_data_t elliptic_vecs_F_u0_p_u;
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_jac);
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_F_u0);
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_F_u0_p_u);

  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_p_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* J_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0 = P4EST_ALLOC(double, local_nodes);
  double* F_u0_p_u = P4EST_ALLOC(double, local_nodes);
  double* J_u_from_res = P4EST_ALLOC(double, local_nodes);

  elliptic_vecs_jac.u = u;
  elliptic_vecs_jac.u0 = u0;
  elliptic_vecs_jac.Au = J_u;

  elliptic_vecs_F_u0.u = u0;
  elliptic_vecs_F_u0.Au = F_u0; /* doesn't matter, not used */

  elliptic_vecs_F_u0_p_u.u = u0_p_u;
  elliptic_vecs_F_u0_p_u.Au = F_u0_p_u; /* doesn't matter, not used */
  
  int num_vecs_to_try = 5;
  double eps = .0001;
  d4est_linalg_fill_vec(u0, 0., local_nodes);
  D4EST_ASSERT(num_vecs_to_try <= local_nodes);

  double max_err = -1;
  int max_err_i = -1;
  int max_err_j = -1;
  
  for (int i = 0; i < num_vecs_to_try; i++){
    u[i] = 1.;
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_jac,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    /* double sum1 = d4est_linalg_vec_sum(J_u,local_nodes); */
    d4est_linalg_vec_axpyeqz(eps, u, u0, u0_p_u, local_nodes);
    
    d4est_elliptic_eqns_build_residual
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_F_u0,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    d4est_elliptic_eqns_build_residual
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_F_u0_p_u,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    for (int j = 0; j < local_nodes; j++){
      J_u_from_res[j] = (F_u0_p_u[j] - F_u0[j])/eps;
      double err = fabs(J_u_from_res[j] - J_u[j]);
      max_err_i = (err > max_err) ? i : max_err_i;
      max_err_j = (err > max_err) ? j : max_err_j;   
      max_err = (err > max_err) ? err : max_err;
    }

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



void
d4est_solver_analytic_jacobian_tester_central_difference
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_elliptic_eqns_t* elliptic_eqns,
 d4est_elliptic_data_t* elliptic_vecs
)
{
  int local_nodes = elliptic_vecs->local_nodes;

  d4est_elliptic_data_t elliptic_vecs_jac;
  d4est_elliptic_data_t elliptic_vecs_F_u0_m_u;
  d4est_elliptic_data_t elliptic_vecs_F_u0_p_u;
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_jac);
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_F_u0_m_u);
  d4est_elliptic_data_copy_ptrs(elliptic_vecs,&elliptic_vecs_F_u0_p_u);

  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_m_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_p_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* J_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0_m_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0_p_u = P4EST_ALLOC(double, local_nodes);
  double* J_u_from_res = P4EST_ALLOC(double, local_nodes);

  elliptic_vecs_jac.u = u;
  elliptic_vecs_jac.u0 = u0;
  elliptic_vecs_jac.Au = J_u;

  elliptic_vecs_F_u0_m_u.u = u0_m_u;
  elliptic_vecs_F_u0_m_u.Au = F_u0_m_u; /* doesn't matter, not used */

  elliptic_vecs_F_u0_p_u.u = u0_p_u;
  elliptic_vecs_F_u0_p_u.Au = F_u0_p_u; /* doesn't matter, not used */
  
  int num_vecs_to_try = 5;
  double eps = .0001;
  d4est_linalg_fill_vec(u0, 0., local_nodes);
  D4EST_ASSERT(num_vecs_to_try <= local_nodes);

  double max_err = -1;
  int max_err_i = -1;
  int max_err_j = -1;
  
  for (int i = 0; i < num_vecs_to_try; i++){
    u[i] = 1.;
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_jac,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    /* double sum1 = d4est_linalg_vec_sum(J_u,local_nodes); */
    d4est_linalg_vec_axpyeqz(eps, u, u0, u0_p_u, local_nodes);
    d4est_linalg_vec_axpyeqz(-eps, u, u0, u0_m_u, local_nodes);
    
    d4est_elliptic_eqns_build_residual
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_F_u0,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    d4est_elliptic_eqns_build_residual
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_F_u0_p_u,
       d4est_ops,
       d4est_geom,
       d4est_quad
      );

    for (int j = 0; j < local_nodes; j++){
      J_u_from_res[j] = (F_u0_p_u[j] - F_u0[j])/(2*eps);
      double err = fabs(J_u_from_res[j] - J_u[j]);
      max_err_i = (err > max_err) ? i : max_err_i;
      max_err_j = (err > max_err) ? j : max_err_j;   
      max_err = (err > max_err) ? err : max_err;
    }

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

