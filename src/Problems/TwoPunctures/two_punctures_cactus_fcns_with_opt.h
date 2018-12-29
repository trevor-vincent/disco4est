#ifndef TWO_PUNCTURES_CACTUS_FCNS_WITH_OPT_H
#define TWO_PUNCTURES_CACTUS_FCNS_WITH_OPT_H 

#include <d4est_util.h>
#include <d4est_element_data.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_laplacian_with_opt.h>
#include <d4est_laplacian_with_opt_flux.h>
#include <d4est_solver_multigrid.h>
#include <d4est_solver_multigrid_matrix_operator.h>
#include <d4est_mesh.h>
#include <d4est_hessian.h>

typedef struct {

  double TP_epsilon;
  double TP_Tiny;
  double par_b;
  double par_S_plus [3];
  double par_P_plus [3];
  double par_S_minus [3];
  double par_P_minus [3];
  double par_m_plus;
  double par_m_minus;

  d4est_mesh_data_t* d4est_factors;
  
} two_punctures_params_t;


typedef struct {

  int use_matrix_operator;
  d4est_solver_multigrid_t* mg_data;
  two_punctures_params_t* two_punctures_params;
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_jac;
  d4est_laplacian_with_opt_flux_data_t* flux_data_for_res;
  d4est_amr_smooth_pred_params_t* smooth_pred_params;
  
} problem_ctx_t;


static
void
init_two_punctures_data
(
 two_punctures_params_t* params
)
{
  double M = 1.;

  params->TP_Tiny = 0.;
  params->TP_epsilon = 0.;
  
  params->par_m_minus = .5*M;
  params->par_m_plus = .5*M;
  params->par_b = 3*M;
  
  params->par_P_minus[0] = 0.;
  params->par_P_minus[1] = -0.2*M;
  params->par_P_minus[2] = 0.;

  params->par_P_plus[0] = 0.;
  params->par_P_plus[1] = 0.2*M;
  params->par_P_plus[2] = 0.;

  params->par_S_minus[0] = 0.;
  params->par_S_minus[1] = 0.;
  params->par_S_minus[2] = 0.;

  params->par_S_plus[0] = 0.;
  params->par_S_plus[1] = 0.;
  params->par_S_plus[2] = 0.;
}


void
BY_Aijofxyz (double x, double y, double z, double Aij[3][3], void* user)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;

  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - params->par_b) * (x - params->par_b) + y * y + z * z;
  r2_minus = (x + params->par_b) * (x + params->par_b) + y * y + z * z;
  
  /* r2_plus = sqrt (pow (r2_plus, 2) + pow (params->TP_epsilon, 4)); */
  /* r2_minus = sqrt (pow (r2_minus, 2) + pow (params->TP_epsilon, 4)); */
  /* r2_plus = sqrt(r2_plus*r2_plus + params->TP_epsilon*params->TP_epsilon*params->TP_epsilon*params->TP_epsilon); */
  /* r2_minus = sqrt(r2_minus*r2_minus + params->TP_epsilon*params->TP_epsilon*params->TP_epsilon*params->TP_epsilon); */  
  /* if (r2_plus < params->TP_Tiny*params->TP_Tiny) */
    /* r2_plus = params->TP_Tiny*params->TP_Tiny; */
  /* if (r2_minus < pow(params->TP_Tiny,2)) */
    /* r2_minus = pow(params->TP_Tiny,2); */
  
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - params->par_b) / r_plus;
  n_minus[0] = (x + params->par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * params->par_P_plus[i];
    nm_Pm += n_minus[i] * params->par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * params->par_S_plus[2] - n_plus[2] * params->par_S_plus[1];
  np_Sp[1] = n_plus[2] * params->par_S_plus[0] - n_plus[0] * params->par_S_plus[2];
  np_Sp[2] = n_plus[0] * params->par_S_plus[1] - n_plus[1] * params->par_S_plus[0];
  nm_Sm[0] = n_minus[1] * params->par_S_minus[2] - n_minus[2] * params->par_S_minus[1];
  nm_Sm[1] = n_minus[2] * params->par_S_minus[0] - n_minus[0] * params->par_S_minus[2];
  nm_Sm[2] = n_minus[0] * params->par_S_minus[1] - n_minus[1] * params->par_S_minus[0];
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (params->par_P_plus[i] * n_plus[j] + params->par_P_plus[j] * n_plus[i]
		 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (params->par_P_minus[i] * n_minus[j] + params->par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }
}

double
BY_KKofxyz (double x, double y, double z, void* user)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;

  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];
  
  r2_plus = (x - params->par_b) * (x - params->par_b) + y * y + z * z;
  r2_minus = (x + params->par_b) * (x + params->par_b) + y * y + z * z;
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - params->par_b) / r_plus;
  n_minus[0] = (x + params->par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * params->par_P_plus[i];
    nm_Pm += n_minus[i] * params->par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * params->par_S_plus[2] - n_plus[2] * params->par_S_plus[1];
  np_Sp[1] = n_plus[2] * params->par_S_plus[0] - n_plus[0] * params->par_S_plus[2];
  np_Sp[2] = n_plus[0] * params->par_S_plus[1] - n_plus[1] * params->par_S_plus[0];
  nm_Sm[0] = n_minus[1] * params->par_S_minus[2] - n_minus[2] * params->par_S_minus[1];
  nm_Sm[1] = n_minus[2] * params->par_S_minus[0] - n_minus[0] * params->par_S_minus[2];
  nm_Sm[2] = n_minus[0] * params->par_S_minus[1] - n_minus[1] * params->par_S_minus[0];
  AijAij = 0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij =
	+ 1.5 * (params->par_P_plus[i] * n_plus[j] + params->par_P_plus[j] * n_plus[i]
                 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (params->par_P_minus[i] * n_minus[j] + params->par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);

      
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}


static
double two_punctures_plus_7o8_K2_psi_neg8
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;


  double r_plus = sqrt ((x - params->par_b) * (x - params->par_b) + y * y + z * z);
  double r_minus = sqrt ((x + params->par_b) * (x + params->par_b) + y * y + z * z);

  double psi =
    1. + 0.5 * params->par_m_plus / r_plus + 0.5 * params->par_m_minus / r_minus + u;
  double psi2 = psi * psi;
  double psi4 = psi2 * psi2;
  double psi8 = psi4 * psi4;

  return 0.875 * BY_KKofxyz (x, y, z,user) / psi8;
}


static
double two_punctures_neg_1o8_K2_psi_neg7
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;

  double r_plus, r_minus, psi, psi2, psi4, psi7;
  r_plus = sqrt ((x - params->par_b) * (x - params->par_b) + y * y + z * z);
  r_minus = sqrt ((x + params->par_b) * (x + params->par_b) + y * y + z * z);

  psi =
    1. + 0.5 * params->par_m_plus / r_plus + 0.5 * params->par_m_minus / r_minus + u;
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;

  /* printf("cactus r_plus r_minus m_plus m_minus = %.6f %.6f %.6f %.6f\n", r_plus, r_minus, params->par_m_plus, params->par_m_minus); */
  /* printf("cactus BY_KKofxyz(x,y,z,user) = %.25f\n", BY_KKofxyz (x, y, z, user)); */
  /* printf("cactus psi, psi7 = %.25f, %.25f\n", psi, psi7); */
  return -0.125 * BY_KKofxyz (x, y, z, user) / psi7;
}

static
void two_punctures_apply_jac_add_nonlinear_term_using_matrix
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  double* M_plus_7o8_K2_psi_neg8_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);
  problem_ctx_t* ctx = user;
  d4est_solver_multigrid_t* mg_data = ctx->mg_data;
  d4est_solver_multigrid_matrix_op_t* matrix_op = mg_data->user_callbacks->user;

  int matrix_stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;

        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        d4est_linalg_matvec_plus_vec(1.,&(matrix_op->matrix[matrix_stride]), &prob_vecs->u[ed->nodal_stride], 0., &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->nodal_stride], volume_nodes, volume_nodes);
        
        matrix_stride += volume_nodes*volume_nodes;
        
      }
    }

  d4est_linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);
}



static
void two_punctures_apply_jac_add_nonlinear_term
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  double* M_plus_7o8_K2_psi_neg8_of_u0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;

        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
        mesh_object.q[2] = ed->q[2];
    
        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);
                               
          d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                                 (
                                                  d4est_factors,
                                                  ed
                                                 );
        
        d4est_quadrature_apply_fofufofvlilj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed->deg,
           md_on_e.xyz_quad,
           J_quad,
           ed->deg_quad,
           &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->nodal_stride],
           two_punctures_plus_7o8_K2_psi_neg8,
           user,
           NULL,
           NULL,
           QUAD_APPLY_MATRIX
          );

      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);
}


static
void two_punctures_apply_jac
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_t* d4est_ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;
  d4est_laplacian_with_opt_flux_data_t* flux_data = ctx->flux_data_for_jac;
  d4est_laplacian_with_opt_apply_aij(p4est,
                          d4est_ghost,
                          d4est_ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad,
                          d4est_factors,
                          0
                         );

  if (ctx->use_matrix_operator == 0)
    two_punctures_apply_jac_add_nonlinear_term
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       user
      );
  else {
    
    d4est_solver_multigrid_t* mg_data = ctx->mg_data;
    d4est_solver_multigrid_matrix_op_t* matrix_op = mg_data->user_callbacks->user;

    if (matrix_op->matrix != matrix_op->matrix_at0){
    two_punctures_apply_jac_add_nonlinear_term_using_matrix
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       user
      );
    }
    else {
    two_punctures_apply_jac_add_nonlinear_term
      (
       p4est,
       d4est_ghost,
       d4est_ghost_data,
       prob_vecs,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       user
      );
    }
  }
}

static void
two_punctures_build_residual_add_nonlinear_term
(
 p4est_t* p4est,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  double* M_neg_1o8_K2_psi_neg7_vec= P4EST_ALLOC(double, prob_vecs->local_nodes);
 
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;


        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
        mesh_object.q[2] = ed->q[2];

    
        double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,
                                                                      ed);

          d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                                 (
                                                  d4est_factors,
                                                  ed
                                                 );

        
        
        d4est_quadrature_apply_fofufofvlj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed->deg,
           J_quad,
           md_on_e.xyz_quad,
           ed->deg_quad,
           &M_neg_1o8_K2_psi_neg7_vec[ed->nodal_stride],
           two_punctures_neg_1o8_K2_psi_neg7,
           user,
           NULL,
           NULL
          );
      }
    }

  d4est_linalg_vec_axpy(1.0,
                  M_neg_1o8_K2_psi_neg7_vec,
                  prob_vecs->Au,
                  prob_vecs->local_nodes);

  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec);

}


static void
two_punctures_pointwise_residual
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double* u,
 double* residual_pointwise_quad,
 int interpolation_option,
 void* user
){
  //interpolation option = 0 -> interpolate u to quad and compute
  //f(u) on quad

  //interpolation option = 1 -> interpolate f(u) to quad  
  d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     HESSIAN_ANALYTICAL,
     u,
     residual_pointwise_quad
    );

  if (interpolation_option == 0){
    double* u_quad = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt){
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;

      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif

        d4est_mesh_data_on_element_t md_on_e
          = d4est_mesh_data_on_element
          (
           d4est_factors,
           ed
          );



        d4est_quadrature_interpolate
          (
           d4est_ops,
           d4est_quad,
           d4est_geom,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &u[ed->nodal_stride],
           ed->deg,
           &u_quad[ed->quad_stride],
           ed->deg_quad
          );

        int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);


        for (int i = 0; i < volume_nodes_quad; i++){
          double x = md_on_e.xyz_quad[0][i];
          double y = md_on_e.xyz_quad[1][i];
          double z = md_on_e.xyz_quad[2][i];
          double u = u_quad[ed->quad_stride + i];
           
          residual_pointwise_quad[ed->quad_stride +  i] =
            two_punctures_neg_1o8_K2_psi_neg7(x,y,z,u,user)
            - residual_pointwise_quad[ed->quad_stride +  i];
        }
      }
    }
    P4EST_FREE(u_quad);
  }
  else {
    double* f_quad = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
    double* f_lob = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes);
    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt){
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;

      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
        
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif

        d4est_mesh_data_on_element_t md_on_e
          = d4est_mesh_data_on_element
          (
           d4est_factors,
           ed
          );
        
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);

        for (int i = 0; i < volume_nodes; i++){
          double x = md_on_e.xyz[0][i];
          double y = md_on_e.xyz[1][i];
          double z = md_on_e.xyz[2][i];
          f_lob[ed->nodal_stride + i] =
            two_punctures_neg_1o8_K2_psi_neg7(x,y,z,u[ed->nodal_stride + i],user);
        }


        d4est_quadrature_interpolate
          (
           d4est_ops,
           d4est_quad,
           d4est_geom,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &f_lob[ed->nodal_stride],
           ed->deg,
           &f_quad[ed->quad_stride],
           ed->deg_quad
          );
        
        for (int i = 0; i < volume_nodes_quad; i++){
          residual_pointwise_quad[ed->quad_stride +  i] =
            f_quad[ed->quad_stride + i]
            - residual_pointwise_quad[ed->quad_stride +  i];
        }
    
      }
    }
    P4EST_FREE(f_quad);
    P4EST_FREE(f_lob);
  }


    
}


static void
two_punctures_build_residual
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_ghost_data_t* d4est_ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* user
)
{
  problem_ctx_t* ctx = user;
  two_punctures_params_t* params = ctx->two_punctures_params;
  d4est_laplacian_with_opt_flux_data_t* flux_data = ctx->flux_data_for_res;
  
  d4est_laplacian_with_opt_apply_aij(p4est,
                          d4est_ghost,
                          d4est_ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad,
                          d4est_factors,
                          0
                         );

  two_punctures_build_residual_add_nonlinear_term
    (
     p4est,
     d4est_ghost,
     d4est_ghost_data,
     prob_vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     user
    );
}

double
two_punctures_robin_coeff_brick_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
 int mortar_node
)
{
  
  double r2 = x*x + y*y + z*z;
  
  if (fabs(boundary_data->n_on_f_m_quad[0][mortar_node]) > 0.)
    return boundary_data->n_on_f_m_quad[0][mortar_node]*x/r2;
  else if (fabs(boundary_data->n_on_f_m_quad[1][mortar_node]) > 0.)
    return boundary_data->n_on_f_m_quad[1][mortar_node]*y/r2;
  else if (fabs(boundary_data->n_on_f_m_quad[2][mortar_node]) > 0.)
    return boundary_data->n_on_f_m_quad[2][mortar_node]*z/r2;
  else{
    D4EST_ABORT("Should not occur");
    return NAN;
  }
}

double
two_punctures_robin_coeff_sphere_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
 int mortar_node
)
{
  double r2 = x*x + y*y + z*z;
  return 1/sqrt(r2);
}



double
two_punctures_robin_bc_rhs_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user,
 d4est_laplacian_with_opt_flux_boundary_data_t* boundary_data,
 int mortar_node
)
{
  return 0.;
}

static double
two_punctures_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 1.;
}


static
void two_punctures_pc_setup_fcn
(
 d4est_krylov_pc_t* d4est_krylov_pc
)
{

  zlog_category_t* c_default = zlog_get_category("two_punctures_pc_setup_fcn");
  d4est_krylov_pc_multigrid_data_t* d4est_krylov_pcmgdata = d4est_krylov_pc->pc_data;
  d4est_solver_multigrid_t* mg_data = d4est_krylov_pcmgdata->mg_data;
  krylov_ctx_t* ctx = d4est_krylov_pc->pc_ctx;

  if (ctx->p4est->mpirank == 0)
    zlog_info(c_default, "Initializing Matrix Operator");

  d4est_solver_multigrid_matrix_setup_fofufofvlilj_operator
      (
       ctx->p4est,
       ctx->d4est_ops,
       ctx->d4est_geom,
       ctx->d4est_quad,
       ctx->d4est_factors,
       ctx->vecs->u0,
       NULL,
       two_punctures_plus_7o8_K2_psi_neg8,
       ctx->fcns->user,
       NULL,
       NULL,
       mg_data->user_callbacks->user
      ); 
}

#endif
