#ifndef TWOPUNCTURESFCNS_NOMG_H
#define TWOPUNCTURESFCNS_NOMG_H 

#include <d4est_element_data.h>
#include <d4est_poisson.h>

#define MAX_PUNCTURES 10
 
typedef struct {
  double C_bh [MAX_PUNCTURES][3];
  double P_bh [MAX_PUNCTURES][3];
  double S_bh [MAX_PUNCTURES][3];
  double M_bh [MAX_PUNCTURES];
  int num_punctures;
  double puncture_eps;
  
} twopunctures_params_t;

static
void
init_onepuncture_data
(
 twopunctures_params_t* params
)
{
  double M = 1.;
  params->num_punctures = 1;
  params->puncture_eps = 1e-10;
  mpi_assert(params->num_punctures < (MAX_PUNCTURES));
  
  params->M_bh[0] = M;

  params->C_bh[0][0] = 0.;
  params->C_bh[0][1] = 0.;
  params->C_bh[0][2] = 0.;

  params->P_bh[0][0] = 0.;
  params->P_bh[0][1] = 0.;
  params->P_bh[0][2] = 0.;

  params->S_bh[0][0] = 0.;
  params->S_bh[0][1] = 0.;
  params->S_bh[0][2] = 0.2*M;
   
}


static
void
init_twopunctures_data
(
 twopunctures_params_t* params
)
{
  double M = 1.;
  params->num_punctures = 2;
  params->puncture_eps = 1e-10;
  mpi_assert(params->num_punctures < (MAX_PUNCTURES));
  
  params->M_bh[0] = .5*M;
  params->M_bh[1] = .5*M;

  params->C_bh[0][0] = -3*M;
  params->C_bh[0][1] = 0.;
  params->C_bh[0][2] = 0.;

  params->C_bh[1][0] = 3*M;
  params->C_bh[1][1] = 0;
  params->C_bh[1][2] = 0;

  params->P_bh[0][0] = 0.;
  params->P_bh[0][1] = -0.2*M;
  params->P_bh[0][2] = 0.;

  params->P_bh[1][0] = 0.;
  params->P_bh[1][1] = 0.2*M;
  params->P_bh[1][2] = 0.;

  params->S_bh[0][0] = 0.;
  params->S_bh[0][1] = 0.;
  params->S_bh[0][2] = 0.;
  
  params->S_bh[1][0] = 0.;
  params->S_bh[1][1] = 0.;
  params->S_bh[1][2] = 0.;

}




static
void
init_random_puncture_data
(
 p4est_t* p4est,
 twopunctures_params_t* params,
 int num_punctures
)
{
  mpi_assert(num_punctures < (MAX_PUNCTURES));
  params->puncture_eps = 1e-10;
  params->num_punctures = num_punctures;
  double M = 1.;

  double rand_x [MAX_PUNCTURES];
  double rand_y [MAX_PUNCTURES];
  double rand_px [MAX_PUNCTURES];
  double rand_py [MAX_PUNCTURES];

  util_gen_rand_vec(&rand_x[0], num_punctures, 1532413243, -5., 5.);
  util_gen_rand_vec(&rand_y[0], num_punctures, 1532413243, -5., 5.);
  util_gen_rand_vec(&rand_px[0], num_punctures, 13232413243, -.2, .2);
  util_gen_rand_vec(&rand_py[0], num_punctures, 14432413243, -.2, .2);
  
  for (int i = 0; i < num_punctures; i++){
    params->M_bh[i] = M/(double)(num_punctures);
    params->C_bh[i][0] = rand_x[i];
    params->C_bh[i][1] = rand_y[i];
    params->C_bh[i][2] = 0.;
    params->P_bh[i][0] = rand_px[i];
    params->P_bh[i][1] = rand_py[i];
    params->P_bh[i][2] = 0.;
    params->S_bh[i][0] = 0.;
    params->S_bh[i][1] = 0.;
    params->S_bh[i][2] = 0.;
    if (p4est->mpirank == 0){
    printf("Puncture %d: M_bh, x, y, z, px, py, pz, sx, sy ,sz \n", i);
    printf(" %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f ,%.6f \n",
           params->M_bh[i],
           params->C_bh[i][0],
           params->C_bh[i][1],
           params->C_bh[i][2],
           params->P_bh[i][0],
           params->P_bh[i][1],
           params->P_bh[i][2],
           params->S_bh[i][0],
           params->S_bh[i][1],
           params->S_bh[i][2]);
    }
    
  }

}

static
double levi_civita(int a, int b, int c)
{
  double eps;
  if( ( ( a == 0 )&&( b == 1 )&&( c == 2 ) ) ||
      ( ( a == 1 )&&( b == 2 )&&( c == 0 ) ) ||
      ( ( a == 2 )&&( b == 0 )&&( c == 1 ) ) ) {
    eps = 1.;
  } else
    if( ( ( a == 1 )&&( b == 0 )&&( c == 2 ) ) ||
        ( ( a == 0 )&&( b == 2 )&&( c == 1 ) ) ||
        ( ( a == 2 )&&( b == 1 )&&( c == 0 ) ) ) {
      eps = -1.;
    } else {
      eps = 0.;
    }
  return eps;
}

static
double kronecker(int a, int b)
{
  if (a != b)
    return 0.0;
  else
    return 1.0;
}

static
double compute_confAij
(
 int i,
 int j,
 double x,
 double y,
 double z,
 void* user
)
{
  twopunctures_params_t* params = user;
  
  double nvec[MAX_PUNCTURES][3];
  double confAij = 0.;
  
  /* initialize and check if (x,y,z) is a puncture point*/
  for (int n = 0; n < params->num_punctures; n++){
    double dxn = (x - params->C_bh[n][0]);
    double dyn = (y - params->C_bh[n][1]);
    double dzn = (z - params->C_bh[n][2]);
    double rn = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    if (rn == 0.)
      rn += params->puncture_eps;
    nvec[n][0] = (x - params->C_bh[n][0])/rn;
    nvec[n][1] = (y - params->C_bh[n][1])/rn;
    nvec[n][2] = (z - params->C_bh[n][2])/rn;

    double lcikl_dot_Sk_dot_nln_times_njn = 0.;
    double lcjkl_dot_Sk_dot_nln_times_nin = 0.;
    double Pn_dot_nvec = 0.;
    
    for (int k = 0; k < 3; k++){
      Pn_dot_nvec += params->P_bh[n][k]*nvec[n][k];
      for (int l = 0; l < 3; l++){
        lcikl_dot_Sk_dot_nln_times_njn += (levi_civita(i,k,l))*params->S_bh[n][k]*nvec[n][l]*nvec[n][j];
        lcjkl_dot_Sk_dot_nln_times_nin += (levi_civita(j,k,l))*params->S_bh[n][k]*nvec[n][l]*nvec[n][i];
      }
    }

    double a = nvec[n][i]*params->P_bh[n][j] + nvec[n][j]*params->P_bh[n][i] - (kronecker(i,j) - nvec[n][i]*nvec[n][j])*Pn_dot_nvec;
    double b = lcikl_dot_Sk_dot_nln_times_njn + lcjkl_dot_Sk_dot_nln_times_nin;

    confAij += (1./(rn*rn))*(a + (2./rn)*b);
  }

  confAij *= 3./2.;
  return confAij;
}

static
double compute_confAij_sqr
(
 double x,
 double y,
 double z,
 void *user
)
{
  double confAij;
  double confAij_sqr = 0.;
  int i,j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      confAij = compute_confAij(i,j,x,y,z,user);
      confAij_sqr += confAij*confAij;
    }
  return confAij_sqr;
}


static
double twopunctures_neg_1o8_K2_psi_neg7
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_params_t* params = user;
  double sumn_mn_o_2rn = 0.;
  double dxn, dyn, dzn, r;
  for (int n = 0; n < params->num_punctures; n++){
    dxn = x - params->C_bh[n][0];
    dyn = y - params->C_bh[n][1];
    dzn = z - params->C_bh[n][2];
    r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    sumn_mn_o_2rn += params->M_bh[n]/(2.*r);
  }


  double psi_0 = 1. + u + sumn_mn_o_2rn;
  double confAij_sqr = compute_confAij_sqr(x,y,z,user);

  if (r > params->puncture_eps)
    return (-1./8.)*confAij_sqr/(psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0);
  else{
    return 0.;
  }
}


static
double twopunctures_plus_7o8_K2_psi_neg8
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_params_t* params = user;
  double sumn_mn_o_2rn = 0.;
  double dxn, dyn, dzn, r;
  for (int n = 0; n < params->num_punctures; n++){
    dxn = x - params->C_bh[n][0];
    dyn = y - params->C_bh[n][1];
    dzn = z - params->C_bh[n][2];
    r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    sumn_mn_o_2rn += params->M_bh[n]/(2.*r);
  }

  double psi_0 = 1. + u + sumn_mn_o_2rn;
  double confAij_sqr = compute_confAij_sqr(x,y,z,user);

  if (r > params->puncture_eps)
    return (7./8.)*confAij_sqr/(psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0*psi_0);
  else{
    return 0.;
  }
}


static
void twopunctures_apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_mortar_fcn_ptrs_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  twopunctures_params_t* params = prob_vecs->user;
  mpi_assert(params->num_punctures > 0);
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_fcn_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad);
  
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
           ed->xyz_quad,
           ed->J_quad,
           ed->deg_quad,
           &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->nodal_stride],
           twopunctures_plus_7o8_K2_psi_neg8,
           prob_vecs->user,
           NULL,
           NULL
          );
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_plus_7o8_K2_psi_neg8_of_u0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_plus_7o8_K2_psi_neg8_of_u0_u_vec);
}

void
twopunctures_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_problem_data_t* prob_vecs,
 d4est_mortar_fcn_ptrs_t* flux_fcn_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  twopunctures_params_t* params = prob_vecs->user;
  mpi_assert(params->num_punctures > 0);
  
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_fcn_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad
                         );

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
           ed->J_quad,
           ed->xyz_quad,
           ed->deg_quad,
           &M_neg_1o8_K2_psi_neg7_vec[ed->nodal_stride],
           twopunctures_neg_1o8_K2_psi_neg7,
           prob_vecs->user,
           NULL,
           NULL
          );

        /* DEBUG_PRINT_ARR_DBL */
        
        /* int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), ed->deg); */
        /* int volume_nodes_gauss = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad); */
        /* double* u_tmp = &prob_vecs->u[ed->nodal_stride]; */
        /* DEBUG_PRINT_ARR_DBL(u_tmp, volume_nodes_lobatto); */
        /* DEBUG_PRINT_ARR_DBL(ed->J_quad, volume_nodes_gauss); */
        
      }
    }

  d4est_linalg_vec_axpy(1.0,
                  M_neg_1o8_K2_psi_neg7_vec,
                  prob_vecs->Au,
                  prob_vecs->local_nodes);

  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec);
}




#endif
