#include "../Operators/grad.h"
#include "../pXest/pXest.h"
#include "../dGMath/d4est_operators.h"
#include "../LinearAlgebra/d4est_linalg.h"


void
grad(double*u, double* gradu [(P4EST_DIM)], double h, int deg, d4est_operators_t* d4est_ops)
{
  int vol_nodes = d4est_operators_get_nodes((P4EST_DIM),deg);
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_Dij(
                     d4est_ops,
                     u,
                     (P4EST_DIM),
                     deg,
                     i,
                     gradu[i]
                    );
    d4est_linalg_vec_scale(2./h, gradu[i], vol_nodes);
  }
}

void
grad_euc_norm(double* u, double* gradu_norm,  double h, int deg, d4est_operators_t* d4est_ops)
{
  int vol_nodes = d4est_operators_get_nodes((P4EST_DIM),deg);

  double* gradu [(P4EST_DIM)];
  int i,j;
  for (i = 0; i < (P4EST_DIM); i++)
    gradu[i] = P4EST_ALLOC(double, vol_nodes);

  grad(u, gradu, h, deg, d4est_ops);
  
  for (i = 0; i < vol_nodes; i++){
    gradu_norm[i] = 0.;
    for (j = 0; j < (P4EST_DIM); j++){
      gradu_norm[i] += gradu[j][i]*gradu[j][i];
    }
    gradu_norm[i] = sqrt(gradu_norm[i]);
  }

  for (i = 0; i < (P4EST_DIM); i++)
    P4EST_FREE(gradu[i]);
}
