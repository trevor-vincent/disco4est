#include <pXest.h>
#include <d4est_solver_jacobian_tester.h>
#include <d4est_linalg.h>
#include <d4est_elliptic_eqns.h>
#include <d4est_elliptic_data.h>
#include <d4est_mesh.h>

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
d4est_solver_jacobian_tester
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_eqns_t* elliptic_eqns,
 int local_nodes,
 d4est_xyz_fcn_t initial_guess,
 void* initial_guess_ctx,
 double eps,
 d4est_analytic_jacobian_tester_type_t type,
 int num_vecs_to_try
)
{
  D4EST_ASSERT(num_vecs_to_try <= local_nodes);

  printf("[D4EST_SOLVER_JACOBIAN_TESTER]: eps = %.25f\n", eps);
  printf("[D4EST_SOLVER_JACOBIAN_TESTER]: num_vecs_to_try = %d\n", num_vecs_to_try);

  if (type == JAC_TEST_FORWARD_DIFFERENCE){
    printf("[D4EST_SOLVER_JACOBIAN_TESTER]: using forward difference\n");
  }

  if (type == JAC_TEST_CENTRAL_DIFFERENCE){
    printf("[D4EST_SOLVER_JACOBIAN_TESTER]: using central difference\n");
  }
  
  
  d4est_elliptic_data_t elliptic_vecs_jac;
  d4est_elliptic_data_t elliptic_vecs_F_u0_m_u;
  d4est_elliptic_data_t elliptic_vecs_F_u0_p_u;

  double* u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_m_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* u0_p_u = P4EST_ALLOC_ZERO(double, local_nodes);
  double* J_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0_m_u = P4EST_ALLOC(double, local_nodes);
  double* F_u0_p_u = P4EST_ALLOC(double, local_nodes);
  double* J_u_from_res = P4EST_ALLOC(double, local_nodes);

  d4est_mesh_init_field
    (
     p4est,
     u0,
     (initial_guess != NULL) ? initial_guess : zero_fcn,
     d4est_ops,
     d4est_geom,
     (initial_guess != NULL) ? initial_guess_ctx : NULL
    );
  
  elliptic_vecs_jac.u = u;
  elliptic_vecs_jac.u0 = u0;
  elliptic_vecs_jac.Au = J_u;
  elliptic_vecs_jac.local_nodes = local_nodes;

  elliptic_vecs_F_u0_m_u.u = u0_m_u;
  elliptic_vecs_F_u0_m_u.Au = F_u0_m_u; /* doesn't matter, not used */
  elliptic_vecs_F_u0_m_u.local_nodes = local_nodes; /* doesn't matter, not used */

  elliptic_vecs_F_u0_p_u.u = u0_p_u;
  elliptic_vecs_F_u0_p_u.Au = F_u0_p_u; /* doesn't matter, not used */
  elliptic_vecs_F_u0_p_u.local_nodes = local_nodes; /* doesn't matter, not used */
  
  double max_err = -1;
  int max_err_i = -1;
  int max_err_j = -1;

  double max_err_avg = 0.;
  double l2_err_avg = 0.;
  
  for (int i = 0; i < num_vecs_to_try; i++){
    u[i] = (type == JAC_TEST_CENTRAL_DIFFERENCE) ? 2*eps : eps;
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

    if (type == JAC_TEST_CENTRAL_DIFFERENCE){
      d4est_linalg_vec_axpyeqz(1.0, u, u0, u0_p_u, local_nodes);
      d4est_linalg_vec_axpyeqz(-1.0, u, u0, u0_m_u, local_nodes);
    }
    else if (type == JAC_TEST_FORWARD_DIFFERENCE){
      d4est_linalg_vec_axpyeqz(1.0, u, u0, u0_p_u, local_nodes);
      d4est_linalg_vec_axpyeqz(0., u, u0, u0_m_u, local_nodes);
    }
    else{
      D4EST_ABORT("Not a supported jacobian tester type");
    }
    
    d4est_elliptic_eqns_build_residual
      (
       p4est,
       ghost,
       ghost_data,
       elliptic_eqns,
       &elliptic_vecs_F_u0_m_u,
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

    double* error_vec = P4EST_ALLOC(double, local_nodes);
    for (int j = 0; j < local_nodes; j++){
      J_u_from_res[j] = (F_u0_p_u[j] - F_u0_m_u[j]);
      /* printf("J_u_from_res[j], Ju[j] =  %.25f, %.25f\n",J_u_from_res[j], J_u[j]); */
      double err = fabs(J_u_from_res[j] - J_u[j]);
      max_err_i = (err > max_err) ? i : max_err_i;
      max_err_j = (err > max_err) ? j : max_err_j;   
      max_err = (err > max_err) ? err : max_err;
    }

    d4est_linalg_vec_axpyeqz(-1., J_u_from_res, J_u, error_vec, local_nodes);
    double l2_err = d4est_mesh_compute_l2_norm_sqr(p4est, d4est_ops, d4est_geom, d4est_quad, error_vec, local_nodes, DO_NOT_STORE_LOCALLY);
    l2_err_avg += l2_err;
    max_err_avg += max_err;

    P4EST_FREE(error_vec);
    u[i] = 0.;
  }

  l2_err_avg /= num_vecs_to_try;
  max_err_avg /= num_vecs_to_try;
  printf("[D4EST_SOLVER_JACOBIAN_TESTER]: l2_err_avg, max_err_avg = %.25f, %.25f\n", sqrt(l2_err_avg), max_err_avg);
  
  P4EST_FREE(u);
  P4EST_FREE(u0);
  P4EST_FREE(u0_p_u);
  P4EST_FREE(u0_m_u);
  P4EST_FREE(J_u);
  P4EST_FREE(F_u0_m_u);
  P4EST_FREE(F_u0_p_u);
  P4EST_FREE(J_u_from_res);
  
}
 
