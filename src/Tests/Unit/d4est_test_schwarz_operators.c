#undef NDEBUG

#include <pXest.h>
#include <d4est_operators.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_linalg.h>
#include <d4est_kron.h>
#include <d4est_util.h>

#define D4EST_REAL_EPS 100*1e-15

static void
d4est_test_operators_schwarz_nd
(
 d4est_solver_schwarz_operators_t* schwarz_ops,
 int dim
)
{
  int res_deg = 1;
  int deg = 1;
  
  double* rst_lobatto [3];
  for (int i = 0; i < dim; i++){
    rst_lobatto[i] = d4est_operators_fetch_lobatto_rst_nd(schwarz_ops->d4est_ops, dim, deg, i);  
  }

  int volume_nodes = d4est_lgl_get_nodes((dim), deg);
  double* frst = P4EST_ALLOC(double, volume_nodes);
  double* res_frst = P4EST_ALLOC(double, volume_nodes);
  double* res_frst_trans = P4EST_ALLOC(double, volume_nodes);
  double* res_rst_lobatto [3];

  for (int i = 0; i < volume_nodes; i++){
    frst[i] = 1.;
    for (int d = 0; d < (dim); d++){
      frst[i] *= exp(rst_lobatto[d][i]);
    }
  }

  int res_volume_nodes = 1;
  int faces [3];
  if (dim == 2){
    faces [0] = 1;
    faces [1] = 2;
    faces [2] = -1;
  }
  else {
    faces [0] = 0;
    faces [1] = 3;
    faces [2] = 4;
  }
  for (int i = 0; i < dim; i++){
    if (faces[i] != -1){
      res_volume_nodes *= (res_deg + 1);
    }
    else {
      res_volume_nodes *= (deg + 1);
    }
    res_rst_lobatto[i] = P4EST_ALLOC(double, volume_nodes);

    d4est_solver_schwarz_operators_apply_schwarz_restrictor
      (
       schwarz_ops,
       rst_lobatto[i],
       dim,
       &faces[0],
       deg,
       res_deg+1,
       0,
       res_rst_lobatto[i]
      );
  }
  
  d4est_solver_schwarz_operators_apply_schwarz_restrictor
    (
     schwarz_ops,
     frst,
     dim,
     &faces[0],
     deg,
     res_deg+1,
     D4OPS_NO_TRANSPOSE,
     res_frst
    );

  for (int i = 0; i < res_volume_nodes; i++){
    double exprst = 1.;
    for (int d = 0; d < (dim); d++){
      exprst *= exp(res_rst_lobatto[d][i]);
    }
    if (fabs(res_frst[i] - exprst) > 100*1e-14){
      D4EST_ABORT("fabs(res_frst[i] - exp(r)*exp(s)*exp(t)) > 100*1e-14");
    }
  }

  d4est_solver_schwarz_operators_apply_schwarz_restrictor
    (
     schwarz_ops,
     res_frst,
     dim,
     &faces[0],
     deg,
     res_deg+1,
     D4OPS_TRANSPOSE,
     res_frst_trans
    );
  
  DEBUG_PRINT_ARR_DBL(frst, volume_nodes);
  DEBUG_PRINT_ARR_DBL(res_frst, res_volume_nodes);
  DEBUG_PRINT_ARR_DBL(res_frst_trans, volume_nodes);
  if (dim == 3){
    DEBUG_PRINT_3ARR_DBL(res_rst_lobatto[0],res_rst_lobatto[1],res_rst_lobatto[2], res_volume_nodes);
  }
  else {
    DEBUG_PRINT_2ARR_DBL(res_rst_lobatto[0],res_rst_lobatto[1], res_volume_nodes);
  }
  P4EST_FREE(frst);
  P4EST_FREE(res_frst);
  P4EST_FREE(res_frst_trans);

  for (int i = 0; i < dim; i++){
    P4EST_FREE(res_rst_lobatto[i]);
  }
}

static void
d4est_test_operators_schwarz_weights
(
 d4est_solver_schwarz_operators_t* schwarz_ops
)
{
  int deg = 2;
  int restricted_size = 2;
  
  double* schwarz_weights_1d =
    d4est_solver_schwarz_operators_fetch_schwarz_weights_1d
    (
     schwarz_ops,
     deg,
     restricted_size
    );

  int size = (deg + 1) + 2*restricted_size;
  DEBUG_PRINT_ARR_DBL(schwarz_weights_1d, size);
}


static void
d4est_test_operators_schwarz
(
 d4est_solver_schwarz_operators_t* schwarz_ops
){


  int deg = 3;
  int res_deg = 3;

  int nodes_res = res_deg + 1;
  int nodes_deg = deg + 1;
  
  double* x = d4est_operators_fetch_lobatto_nodes_1d(schwarz_ops->d4est_ops, deg);
  double* fx = P4EST_ALLOC(double, deg + 1);

  double* res_x_side0 = P4EST_ALLOC(double, res_deg + 1);
  double* res_fx_side0 = P4EST_ALLOC(double, res_deg + 1);

  double* res_x_side1 = P4EST_ALLOC(double, res_deg + 1);
  double* res_fx_side1 = P4EST_ALLOC(double, res_deg + 1);

  double* res_x_side0_trans = P4EST_ALLOC(double, deg + 1);
  double* res_fx_side0_trans = P4EST_ALLOC(double, deg + 1);
  double* res_x_side1_trans = P4EST_ALLOC(double, deg + 1);
  double* res_fx_side1_trans = P4EST_ALLOC(double, deg + 1);

  for (int i = 0; i < deg + 1; i++){
    fx[i] = exp(x[i]);
  }
  
  double* schwarz_1d = d4est_solver_schwarz_operators_fetch_schwarz_restrictor_1d
                       (
                        schwarz_ops,
                        deg,
                        res_deg + 1
                       );
  
  double* schwarz_1d_trans = d4est_solver_schwarz_operators_fetch_schwarz_restrictor_transpose_1d
                             (
                              schwarz_ops,
                              deg,
                              res_deg + 1
                             );


  
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d, x, 0., res_x_side0, nodes_res, nodes_deg);       
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d, fx, 0., res_fx_side0, nodes_res, nodes_deg);
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d[nodes_deg*nodes_res], x, 0., res_x_side1, nodes_res, nodes_deg);       
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d[nodes_deg*nodes_res], fx, 0., res_fx_side1, nodes_res, nodes_deg);

  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d_trans, res_x_side0, 0., res_x_side0_trans, nodes_deg, nodes_res);       
  d4est_linalg_matvec_plus_vec(1.0, schwarz_1d_trans, res_fx_side0, 0., res_fx_side0_trans, nodes_deg, nodes_res);

  
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d_trans[nodes_deg*nodes_res], res_x_side1, 0., res_x_side1_trans, nodes_deg, nodes_res);       
  d4est_linalg_matvec_plus_vec(1.0, &schwarz_1d_trans[nodes_deg*nodes_res], res_fx_side1, 0., res_fx_side1_trans, nodes_deg, nodes_res);

  for (int i = 0; i < res_deg + 1; i++){
    D4EST_ASSERT(fabs(exp(res_x_side0[i]) - res_fx_side0[i]) < 1e-14);
    D4EST_ASSERT(fabs(exp(res_x_side1[i]) - res_fx_side1[i]) < 1e-14);
    D4EST_ASSERT(res_fx_side0_trans[i] == 0. || fabs(exp(res_x_side0_trans[i]) - res_fx_side0_trans[i]) < 1e-14);
    D4EST_ASSERT(res_fx_side1_trans[i] == 0. || fabs(exp(res_x_side1_trans[i]) - res_fx_side1_trans[i]) < 1e-14);
  }
  P4EST_FREE(fx);
  P4EST_FREE(res_x_side0);
  P4EST_FREE(res_fx_side0);
  P4EST_FREE(res_x_side1);
  P4EST_FREE(res_fx_side1);

  P4EST_FREE(res_x_side0_trans);
  P4EST_FREE(res_fx_side0_trans);
  P4EST_FREE(res_x_side1_trans);
  P4EST_FREE(res_fx_side1_trans);
  
}

int main(int argc, char *argv[])
{
  d4est_operators_t* d4est_ops = d4est_ops_init(20);  
  d4est_solver_schwarz_operators_t* d4est_schwarz_ops =
    d4est_solver_schwarz_operators_init(d4est_ops);
  
  d4est_test_operators_schwarz_nd
    (
     d4est_schwarz_ops,
     3
    );

  d4est_test_operators_schwarz
    (
     d4est_schwarz_ops
    );

  d4est_test_operators_schwarz_weights
    (
     d4est_schwarz_ops
    );
  
  d4est_solver_schwarz_operators_destroy(d4est_schwarz_ops);
  d4est_ops_destroy(d4est_ops);
  return 0;
}

