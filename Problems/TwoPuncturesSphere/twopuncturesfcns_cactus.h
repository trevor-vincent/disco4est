#ifndef TWOPUNCTURESFCNS_CACTUS_H
#define TWOPUNCTURESFCNS_CACTUS_H 

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
  int deg_offset_for_puncture_nonlinearity_quad;
  
} twopunctures_cactus_params_t;

static
void
init_cactus_puncture_data
(
 twopunctures_cactus_params_t* params,
 int deg_offset_for_puncture_nonlinearity_quad
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
  
  params->deg_offset_for_puncture_nonlinearity_quad = deg_offset_for_puncture_nonlinearity_quad;
}


void
BY_Aijofxyz (double x, double y, double z, double Aij[3][3], void* user)
{
  twopunctures_cactus_params_t* params = user;
  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - params->par_b) * (x - params->par_b) + y * y + z * z;
  r2_minus = (x + params->par_b) * (x + params->par_b) + y * y + z * z;
  r2_plus = sqrt (pow (r2_plus, 2) + pow (params->TP_epsilon, 4));
  r2_minus = sqrt (pow (r2_minus, 2) + pow (params->TP_epsilon, 4));
  if (r2_plus < pow(params->TP_Tiny,2))
    r2_plus = pow(params->TP_Tiny,2);
  if (r2_minus < pow(params->TP_Tiny,2))
    r2_minus = pow(params->TP_Tiny,2);
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
  twopunctures_cactus_params_t* params = user;
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
double twopunctures_cactus_plus_7o8_K2_psi_neg8
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_cactus_params_t* params = user;

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
double twopunctures_cactus_neg_1o8_K2_psi_neg7
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_cactus_params_t* params = user;
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
void twopunctures_cactus_apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  twopunctures_cactus_params_t* params = prob_vecs->user;
  curved_poisson_operator_primal_apply_aij(p4est,
                                           ghost,
                                           ghost_data,
                                           prob_vecs,
                                           d4est_ops,
                                           d4est_geom);
  
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
        curved_element_data_t* ed = quad->p.user_data;        
        curved_element_data_apply_fofufofvlilj_Gaussnodes
          (
           d4est_ops,
           d4est_geom,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed,
           ed->deg_quad + params->deg_offset_for_puncture_nonlinearity_quad,
           (P4EST_DIM),
           &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->nodal_stride],
           twopunctures_cactus_plus_7o8_K2_psi_neg8,
           prob_vecs->user,
           NULL,
           NULL
          );
      }
    }
  
  linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);
}

static
void
twopunctures_cactus_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  twopunctures_cactus_params_t* params = prob_vecs->user;
  curved_poisson_operator_primal_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);

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
        curved_element_data_t* ed = quad->p.user_data;        
        curved_element_data_apply_fofufofvlj_Gaussnodes
          (
           d4est_ops,
           d4est_geom,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed,
           ed->deg_quad + params->deg_offset_for_puncture_nonlinearity_quad,
           (P4EST_DIM),
           &M_neg_1o8_K2_psi_neg7_vec[ed->nodal_stride],
           twopunctures_cactus_neg_1o8_K2_psi_neg7,
           prob_vecs->user,
           NULL,
           NULL
          );
      }
    }

  linalg_vec_axpy(1.0,
                  M_neg_1o8_K2_psi_neg7_vec,
                  prob_vecs->Au,
                  prob_vecs->local_nodes);

  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec); 
}

#endif
