#include "../pXest/pXest.h"
#include "../Solver/cg_eigs.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Utilities/d4est_util.h"
#include "sc_reduce.h"
#include <signal.h>


static void
tridiag_gershgorin
(
 int i,
 int local_nodes,
 double a0,
 double b0,
 double a1,
 double b1,
 double* max,
 double* min
)
{
  double diag, offdiag_sum;
  if (i != 0 && i < local_nodes) {
    diag = (1./a1 + b0/a0);
    offdiag_sum = fabs(sqrt(b1)/a1) + fabs(sqrt(b0)/a0);
  }
  else if (i == 0){
    diag = 1./a1;
    offdiag_sum = sqrt(b1)/a1;
  }
  else {
    diag = 1./a1 + b0/a0;
    offdiag_sum = fabs(sqrt(b0)/a0);
  }
  *max = diag + offdiag_sum;
  *min = diag - offdiag_sum;
}


void
cg_eigs
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 int imax,
 int print_spectral_bound_iterations,
 double* spectral_bound
)
{

  zlog_category_t *c_default = zlog_get_category("cg_eigs");
  int local_nodes;
  double delta_new, delta_old;
  double alpha = -1.;
  double beta = -1.;
  
  double* Au; 
  double* rhs;
  double* d;
  double* r;
  
  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  double* tmp = vecs->u;
  double* u;
  u = vecs->u;
  rhs = vecs->rhs; 

  double d_dot_Au;
  
  d = P4EST_ALLOC(double, local_nodes);
  r = P4EST_ALLOC(double, local_nodes);

  /* first iteration data, store Au in r */  
  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     ghost,
     ghost_data,
     fcns,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );

  d4est_util_copy_1st_to_2nd(Au, r, local_nodes);
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);
  d4est_util_copy_1st_to_2nd(r, d, local_nodes);
  delta_new = d4est_linalg_vec_dot(r,r,local_nodes);

  double delta_new_global = -1;
  double d_dot_Au_global = -1;
  
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
  
  /* start working on d */
  vecs->u = d;

  int i;
  double alpha_old, beta_old;
  double temp_max, temp_min;

  for (i = 0; i < imax; i++){
    
    d4est_elliptic_eqns_apply_lhs
      (
       p4est,
       ghost,
       ghost_data,
       fcns,
       vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors
      );

    d_dot_Au = d4est_linalg_vec_dot(d,Au,local_nodes);

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
    alpha_old = alpha;
    alpha = delta_new/d_dot_Au;
    
    d4est_linalg_vec_axpy(alpha, d, u, local_nodes);
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

    beta_old = beta;
    beta = delta_new/delta_old;
    d4est_linalg_vec_xpby(r, beta, d, local_nodes);

    tridiag_gershgorin(i, local_nodes, alpha_old, beta_old, alpha, beta, &temp_max, &temp_min);

    if (i > 0){
      *spectral_bound = d4est_util_max( *spectral_bound, temp_max );  
    }
    else{
      *spectral_bound = temp_max;
    }

    if (print_spectral_bound_iterations){
      zlog_info(c_default, "iter %d, CG eig spectral bound = %f", i, *spectral_bound);
    }
  }
  
  /* set back from d */
  vecs->u = u;

  P4EST_FREE(d);
  P4EST_FREE(r);

}
