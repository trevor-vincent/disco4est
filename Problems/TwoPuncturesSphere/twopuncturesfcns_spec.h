#ifndef TWOPUNCTURESFCNS_SPEC_H
#define TWOPUNCTURESFCNS_SPEC_H 

static
double twopunctures_spec_plus_7o8_K2_psi_neg8
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_params_t* params = user;
  double oneoveralpha = twopunctures_spec_compute_oneoveralpha(
                                                               x,
                                                               y,
                                                               z,
                                                               params->C_bh,
                                                               params->M_bh,
                                                               params->num_punctures
                                                              );
  
  double psi_0 = 1. + u + oneoveralpha;
  double confAij_sqr = twopunctures_spec_compute_tildeAijsqr
                       (
                        x,
                        y,
                        z,
                        params->P_bh,
                        params->S_bh,
                        params->C_bh,
                        params->num_punctures
                       );

  return (.875)*confAij_sqr*pow(psi_0, -8);
}

static
double twopunctures_spec_neg_1o8_K2_psi_neg7
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  twopunctures_params_t* params = user;
  double oneoveralpha = twopunctures_spec_compute_oneoveralpha(
                                                               x,
                                                               y,
                                                               z,
                                                               params->C_bh,
                                                               params->M_bh,
                                                               params->num_punctures
                                                              );
  
  double psi_0 = 1. + u + oneoveralpha;
  double confAij_sqr = twopunctures_spec_compute_tildeAijsqr
                       (
                        x,
                        y,
                        z,
                        params->P_bh,
                        params->S_bh,
                        params->C_bh,
                        params->num_punctures
                       );
  
  return (-.125)*confAij_sqr*pow(psi_0, -7);
}

double
twopunctures_spec_compute_oneoveralpha
(
 double global_x,
 double global_y,
 double global_z,
 double C [MAX_PUNCTURES][3],
 double M [MAX_PUNCTURES][3],
 int num_punctures
)
{
  double oneoveralpha = 0.;
  for (int n = 0; n < num_punctures; n++){
    double dxn = x - C[n][0];
    double dyn = y - C[n][1];
    double dzn = z - C[n][2];
    double r = sqrt(dxn*dxn + dyn*dyn + dzn*dzn);
    oneoveralpha += M[n]/(2.0*r);
  }
  return oneoveralpha;
}


double
twopunctures_spec_compute_tildeAijsqr
(
double global_x,
double global_y,
double global_z,
double P [MAX_PUNCTURES][3],
double S [MAX_PUNCTURES][3],
double C [MAX_PUNCTURES][3],
int num_punctures
)
{ 
  double tildeAijsqr = 0.;
  double tildeAij_n [3][3];
  double tildeAij [3][3];


  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      tildeAij[i][j] = 0.;
    }
  }

  
  for (int n = 0; n < num_punctures; n++){

    double P_n [] = {P[n][0], P[n][1], P[n][2]};
    double S_n [] = {S[n][0], S[n][1], S[n][2]};
    double C_n [] = {C[n][0], C[n][1], C[n][2]};
    
    twopunctures_spec_compute_tildeAij
      (
       global_x,
       global_y,
       global_z,
       P_n,
       S_n,
       C_n,
       tildeAij_n
      );

    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
        tildeAij[i][j] += tildeAij_n[i][j];
      }
    }
    
  }
  
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      tildeAijsqr += tildeAij[i][j]*tildeAij[i][j]; 
    }
  }

  return tildeAijsqr;
}

void
twopunctures_spec_compute_oneBH_tildeAij
(
 double global_x,
 double global_y,
 double global_z,
 double P [3],
 double S [3],
 double C [3],
 double tildeAij[3][3]
)
{
  const int dim=3;
  double x [] =
    {
      global_x - C[0],
      global_y - C[1],
      global_z - C[2]
    };
        
    /* for(int d=0; d<dim; ++d){ */
    /*   x[d] -= mC1[d]; */
    /* } */

    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r = sqrt(r2);
    /* const DataMesh r2=sqr(x[0])+sqr(x[1])+sqr(x[2]); */
    /* const DataMesh r=sqrt(r2); */
    double normal [] = {x[0]/r, x[1]/r, x[2]/r};
    
    /* const MyVector<DataMesh> n(MV::fill, x[0]/r, x[1]/r, x[2]/r); */
    double Pdotn = P[0]*n[0] + P[1]*n[1] + P[2]*n[2];
    /* const DataMesh Pdotn = mP1[0]*n[0]+mP1[1]*n[1]+mP1[2]*n[2]; */

    double temp [] = {S[1]*n[2] - S[2]*n[1],
                      S[2]*n[0] - S[0]*n[2],
                      S[0]*n[1] - S[1]*n[0]};
    
    // temp[i] = (S\times n)_i
    /* MyVector<DataMesh> temp(MV::fill, */
			    /* mS1[1]*n[2] - mS1[2]*n[1], */
			    /* mS1[2]*n[0] - mS1[0]*n[2], */
			    /* mS1[0]*n[1] - mS1[1]*n[0]); */

    for(int i=0; i<dim; ++i) {
      for(int j=0; j<=i; ++j) {
	tildeAij[i][j]=mP1[i]*n[j]+mP1[j]*n[i] + n[i]*n[j]*Pdotn;
      }
      tildeAij(i,i) -= Pdotn;
    }
    for(int i=0; i<3; ++i) {
      for(int j=0; j<=i; ++j) {
	tildeAij[i][j] *= 1.5/r2;
      }
    }
    for(int i=0; i<3; ++i) {
      for(int j=0; j<=i; ++j) {
	tildeAij[i][j] += 3/(r*r2)*(temp[i]*n[j]+temp[j]*n[i]);
      }
    }
}

static
void twopunctures_spec_apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom
)
{
  twopunctures_params_t* params = prob_vecs->user;
  curved_poisson_operator_primal_apply_aij(p4est,
                                           ghost,
                                           ghost_data,
                                           prob_vecs,
                                           dgmath_jit_dbase,
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
           dgmath_jit_dbase,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed,
           ed->deg_integ + params->deg_offset_for_puncture_nonlinearity_integ,
           (P4EST_DIM),
           &M_plus_7o8_K2_psi_neg8_of_u0_u_vec[ed->nodal_stride],
           twopunctures_spec_plus_7o8_K2_psi_neg8,
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
twopunctures_spec_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom
)
{
  twopunctures_params_t* params = prob_vecs->user;
  curved_poisson_operator_primal_apply_aij(p4est, ghost, ghost_data, prob_vecs, dgmath_jit_dbase, d4est_geom);

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
           dgmath_jit_dbase,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed,
           ed->deg_integ + params->deg_offset_for_puncture_nonlinearity_integ,
           (P4EST_DIM),
           &M_neg_1o8_K2_psi_neg7_vec[ed->nodal_stride],
           twopunctures_spec_neg_1o8_K2_psi_neg7,
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
