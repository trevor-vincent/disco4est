#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_linalg.h>
#include <util.h>

#define D4EST_REAL_EPS 1e-15

void
test_d4est_operators_mass_1d
(
 d4est_operators_t* d4est_ops
)
{
  int deg = 2;
  int nodes = deg + 1;
  int N = nodes - 1;
  
  /* matrix of size nodes_H x nodes_h */
  double* mij = P4EST_ALLOC(double, nodes*nodes);

  double* w = d4est_operators_fetch_lobatto_weights_1d(d4est_ops, deg);
  double* x = d4est_operators_fetch_lobatto_nodes_1d(d4est_ops, deg);

  /* Compute mass from Teukolsky mass matrix paper */
  for (int i = 0; i < nodes; i++){
    for (int j = 0; j < nodes; j++){
      double hN = 2./(2.*(double)N + 1.);
      double gN = 2./(double)N;
      double alpha = (hN - gN)/(gN*gN);
      double M = 0.;
      double pNxi = d4est_lgl_jacobi(x[i], 0., 0., N);
      double pNxj = d4est_lgl_jacobi(x[j], 0., 0., N);
      if (i == j)
        M += w[i];
      M += alpha*w[i]*w[j]*pNxi*pNxj*hN;
      mij[i*nodes + j] = M;
    }
  }
  double* mij_2 = d4est_operators_fetch_mij_1d(d4est_ops, deg);
  if (util_compare_vecs(mij, mij_2, nodes*nodes, D4EST_REAL_EPS) == 0){
    util_print_matrices(w,x, nodes, 1, "w,x = ");
    util_print_matrices(mij, mij_2, nodes*nodes, 1, "mij, mij_2 = ");
    P4EST_FREE(mij);
    mpi_abort("[D4EST_ERROR]:mij and mij_2 not equal\n");
  }
  P4EST_FREE(mij);
}

void
test_d4est_operators_mass_nd
(
 d4est_operators_t* d4est_ops,
 int dim
)
{
  int deg = 2;
  int volume_nodes = d4est_lgl_get_nodes(dim, deg);

  double* mij = P4EST_ALLOC_ZERO(double, volume_nodes*volume_nodes);
  double* u = P4EST_ALLOC_ZERO(double, volume_nodes);
  double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes);

  for (int i = 0; i < volume_nodes; i++){
    u[i] = 1.;
    d4est_operators_apply_mij(d4est_ops, u, dim, deg, Mu);
    d4est_linalg_set_column(mij, Mu, i, volume_nodes, volume_nodes);
    u[i] = 0.;
  }

  for (int i = 0; i < volume_nodes; i++){
    for (int j = 0; j < volume_nodes; j++){
      printf("M[%d][%d] = %.15f\n", i,j,mij[i*volume_nodes + j]);
    }
  }

  P4EST_FREE(Mu);
  P4EST_FREE(mij);
  P4EST_FREE(u);
}


int main(int argc, char *argv[])
{

  d4est_operators_t* d4est_ops = d4est_ops_init(20);  
  test_d4est_operators_mass_1d(d4est_ops);
  test_d4est_operators_mass_nd(d4est_ops,2);
  d4est_ops_destroy(d4est_ops);
  return 0;
}
