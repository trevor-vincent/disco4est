#include <pXest.h>
#include <d4est_lgl.h>
#include <math.h>
#include <util.h>

int d4est_lgl_get_nodes(int dim, int deg) {
  D4EST_ASSERT(dim > 0 && dim < 4 && deg > 0);
  int volume_nodes = 1;
  for (int i = 0; i < dim; i++) volume_nodes *= deg + 1;
  return volume_nodes;
}

double d4est_lgl_jacobi(double r, double alpha, double beta, int N) {
  int i;
  double J0, J1, J2 = 0;
  double h1, anew, bnew;
  double gamma0 = (pow(2.0, alpha + beta + 1) / (alpha + beta + 1.0)) *
                  tgamma(alpha + 1.0) * tgamma(beta + 1.0) /
                  tgamma(alpha + beta + 1.0);
  double gamma1 =
      (alpha + 1.0) * ((beta + 1.0) / (alpha + beta + 3.0)) * gamma0;
  double aold = (2.0 / (2.0 + alpha + beta)) *
                sqrt((alpha + 1.0) * (beta + 1.0) / (alpha + beta + 3.0));

  J0 = 1.0 / sqrt(gamma0);
  if (N == 0) return J0;

  J1 = ((alpha + beta + 2.0) * r / 2.0 + (alpha - beta) / 2) / sqrt(gamma1);
  if (N == 1) return J1;

  // Calculate new polynomials from Jacobi polynomial recurrence relation
  for (i = 0; i < N - 1; i++) {
    h1 = 2.0 * (i + 1) + alpha + beta;
    anew = (2.0 / (h1 + 2.0)) *
           sqrt(((i + 1) + 1.0) * ((i + 1) + 1.0 + alpha + beta) *
                ((i + 1) + 1.0 + alpha) *
                (((i + 1) + 1.0 + beta) / (h1 + 1.0)) / (h1 + 3.0));
    bnew = -((alpha * alpha - beta * beta) / h1) / (h1 + 2.0);
    J2 = (1.0 / anew) * (-aold * J0 + (r - bnew) * J1);

    /* cycle */
    aold = anew;
    J0 = J1;
    J1 = J2;
  }
  return J2;
}

double d4est_lgl_gradjacobi(double r, double alpha, double beta, int N) {
  if (N == 0)
    return 0;
  else
    return sqrt(N * (N + alpha + beta + 1.0)) *
           d4est_lgl_jacobi(r, alpha + 1.0, beta + 1.0, N - 1);
}


