double
d4est_solver_schwarz_metadata.cg_single_core
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_ghost_t** ghost,
 d4est_ghost_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* d4est_schwarz_data,
 d4est_solver_schwarz_operators_t* d4est_schwarz_ops,
 double* du_subdomain,
 double* rhs_subdomain,
 int subdomain,
 int nodes,
 int iter,
 double tol
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
  double delta_new, delta_old, temp_max, temp_min, d_dot_Ad;
  double alpha = -1.;
  double beta = -1.;
  double alpha_old;
  double beta_old;

  double* Ax = D4EST_TEST_ALLOC(double, nodes); 
  double* d = D4EST_TEST_ALLOC(double, nodes); 
  double* r = D4EST_TEST_ALLOC(double, nodes);

  apply_mat(A, x0, Ax, nodes);
  d4est_util_copy_1st_to_2nd(Au, r, nodes);
  d4est_linalg_vec_xpby(b, -1., r, nodes);
  d4est_util_copy_1st_to_2nd(r, d, nodes);
  delta_new = d4est_linalg_vec_dot(r,r,nodes);

  for (int i = 0; i < iter; i++){

    apply_mat(A, d, Ax, nodes);
    d_dot_Ad = d4est_linalg_vec_dot(d, Ax, nodes);
    alpha_old = alpha;
    alpha = delta_new/d_dot_Ad;
        
    d4est_linalg_vec_axpy(alpha, d, u, nodes);
    d4est_linalg_vec_axpy(-alpha, Au, r, nodes);

    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, nodes);
  }
  
  free(Ax);
  free(d);
  free(r);
  return spectral_bound;
}
