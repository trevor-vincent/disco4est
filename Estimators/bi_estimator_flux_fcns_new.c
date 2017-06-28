#include "../Estimators/bi_estimator_flux_fcns.h"
#include "../LinearAlgebra/d4est_linalg.h"

double bi_est_sipg_flux_penalty_prefactor;
penalty_calc_t bi_est_u_prefactor_calculate_fcn;
penalty_calc_t bi_est_u_dirichlet_prefactor_calculate_fcn;
penalty_calc_t bi_est_gradu_prefactor_calculate_fcn;

static void
bi_est_dirichlet
(
 element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 void* params
)
{
  grid_fcn_t u_at_bndry = bndry_fcn;
  double prefactor = bi_est_u_dirichlet_prefactor_calculate_fcn
                 (
                  e_m->deg,
                  e_m->h,
                  e_m->deg,
                  e_m->h,
                  bi_est_sipg_flux_penalty_prefactor
                 );

  /* printf("prefactor = %f\n", prefactor); */
  
  int face_nodes_m = d4est_lgl_get_nodes ( (P4EST_DIM) - 1, e_m->deg );
  double* tmp = P4EST_ALLOC(double, face_nodes_m);
  double* xyz_on_f_m [(P4EST_DIM)];
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);
  
  /* get solution variable on this face */
  d4est_operators_apply_slicer(
                      d4est_ops,
                      e_m->u_elem,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);

  /* get xyz on this face */
  int dir;
  for (dir = 0; dir < (P4EST_DIM); dir++){
    xyz_on_f_m[dir] = P4EST_ALLOC(double, face_nodes_m);

    double* rst = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), e_m->deg, dir);
    d4est_operators_apply_slicer(d4est_ops, rst, (P4EST_DIM), f_m, e_m->deg, tmp);
    
    d4est_reference_rtox_array(tmp, e_m->xyz_corner[dir], e_m->h, xyz_on_f_m[dir], face_nodes_m);
  }

  /* get boundary values on this face */
  int i;
  for (i = 0; i < face_nodes_m; i++){
    tmp[i] = u_at_bndry(xyz_on_f_m[0][i],
                             xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                             ,
                             xyz_on_f_m[2][i]
#endif                          
                            );
  }

  double n [3];
  d4est_reference_get_normal(f_m, (P4EST_DIM), &n[0]);

  double* Je2 = P4EST_ALLOC(double, face_nodes_m);
  double* MJe2 = P4EST_ALLOC(double, face_nodes_m);

  
  for (dir = 0; dir < (P4EST_DIM); dir++){
    /* calculate qstar - q(-) */
    for(i = 0; i < face_nodes_m; i++){
      /* e_m->ustar_min_u[f_m*face_nodes_m + i] = 0.; */

      
      /* printf("tmp[%d] = %f", tmp[i]); */
      /* the ||mu \del u|| term does not contribute at the boundary to eta2*/
      /* e_m->qstar_min_q[dir][f_m*face_nodes_m + i] = n[dir]*(u_m_on_f_m[i] - tmp[i]);// */
      Je2[i] = n[dir]*(u_m_on_f_m[i] - tmp[i]);//
      /* e_m->qstar_min_q[dir][f_m*face_nodes_m + i] *= prefactor; */
      Je2[i] *= prefactor;
    }


    d4est_operators_apply_mij(d4est_ops, Je2, (P4EST_DIM)-1, e_m->deg, MJe2);
    double Je2MJe2 = d4est_linalg_vec_dot(Je2, MJe2, face_nodes_m);
    Je2MJe2 *= e_m->surface_jacobian;
    e_m->local_estimator += Je2MJe2;
  }

  for (dir = 0; dir < (P4EST_DIM); dir++){
    P4EST_FREE(xyz_on_f_m[dir]);
  }

  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(Je2);
  P4EST_FREE(MJe2);
  P4EST_FREE(tmp);
}

static void
bi_est_interface
(
 element_data_t** e_m,
 int faces_m,
 int f_m,
 element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 d4est_operators_t* d4est_ops,
 void* params
)
{
  int deg_p [(P4EST_HALF)];
  int max_deg_p = -1;
  int face_nodes_p [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int face_nodes_m [(P4EST_HALF)];
  int max_deg_m = -1;
  
  int deg_mortar [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  int nodes_mortar [(P4EST_HALF)];
  
  double gradu_prefactor_mortar [(P4EST_HALF)];
  double u_prefactor_mortar [(P4EST_HALF)];
  
  double n [(P4EST_DIM)];
  d4est_reference_get_normal(f_m, (P4EST_DIM), &n[0]);
  
  int i,j;
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m = 0;
  for (i = 0; i < faces_m; i++){
    deg_m[i] = e_m[i]->deg;
    if (e_m[i]->deg > max_deg_m) max_deg_m = e_m[i]->deg;
    face_nodes_m[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg );
    total_side_nodes_m += face_nodes_m[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p = 0;
  for (i = 0; i < faces_p; i++){
    deg_p[i] = e_p[i]->deg;
    if (e_p[i]->deg > max_deg_p) max_deg_p = e_p[i]->deg;
    face_nodes_p[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p[i]->deg );
    total_side_nodes_p += face_nodes_p[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar = 0;
  for (i = 0; i < faces_m; i++)
    for (j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar[i+j] = util_max_int( e_m[i]->deg, e_p[j]->deg );

      /* find IP penalty parameter for each face pair of the two sides*/
      gradu_prefactor_mortar[i+j] = bi_est_gradu_prefactor_calculate_fcn
                            (
                             e_m[i]->deg,
                             e_m[i]->h,
                             e_p[j]->deg,
                             e_p[j]->h,
                             bi_est_sipg_flux_penalty_prefactor
                            );

      u_prefactor_mortar[i+j] = bi_est_u_prefactor_calculate_fcn
                            (
                             e_m[i]->deg,
                             e_m[i]->h,
                             e_p[j]->deg,
                             e_p[j]->h,
                             bi_est_sipg_flux_penalty_prefactor
                            );    

      /* printf("Je1_prefactor_mortar, Je2_prefactor_mortar = %f,%f\n", gradu_prefactor_mortar[i], u_prefactor_mortar[i]); */
      
      nodes_mortar[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar[i+j] );
      
      total_nodes_mortar += nodes_mortar[i+j];
    }

  /* scalar and vector fields on each of the (-) and (+) elements */
  double* du_m = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), max_deg_m));
  double* du_p = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), max_deg_p));
  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* du_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m); 
  double* du_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* u_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);

  double* du_m_on_f_m_mortar [(P4EST_DIM)];
  double* du_p_on_f_p_mortar [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    du_m_on_f_m_mortar[d] = P4EST_ALLOC(double, total_nodes_mortar);
    du_p_on_f_p_mortar[d] = P4EST_ALLOC(double, total_nodes_mortar);
  }
  
  double* qstar_min_q_mortar = P4EST_ALLOC(double, total_nodes_mortar);

  /* NOTE THE SUBTLE ALLOC_ZERO here, we want to zero out ustar_min_u_mortar because it's a dot product */
  double* ustar_min_u_mortar = P4EST_ALLOC_ZERO(double, total_nodes_mortar);
  double* qstar_min_q_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* ustar_min_u_m = P4EST_ALLOC(double, total_side_nodes_m);

  int stride = 0;
  for (i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_m[i]->u_elem[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    stride += face_nodes_m[i];
  }
 
  stride = 0;
  for (i = 0; i < faces_p; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_p[i]->u_elem[0]),
       (P4EST_DIM),
       f_p,
       e_p[i]->deg,
       &u_p_on_f_p[stride]
      );
    stride += face_nodes_p[i];
  }
  
  /* project (-)-side u trace vector onto mortar space */ 
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar
    );

  /* project (+)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar
    );

  /* For each component of the vector */
  int dir;
  for (dir = 0; dir < (P4EST_DIM); dir++){

    /* compute the (-)-u-derivative and project on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (i = 0; i < faces_m; i++){
      d4est_operators_apply_dij
        (
         d4est_ops,
         e_m[i]->u_elem,
         (P4EST_DIM),
         e_m[i]->deg,
         dir,
         du_m
        );
      
      d4est_linalg_vec_scale(2./e_m[i]->h, du_m, d4est_lgl_get_nodes((P4EST_DIM), e_m[i]->deg));

      d4est_operators_apply_slicer
        (
         d4est_ops,
         du_m,
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &du_m_on_f_m[stride]
        );
    
       stride += face_nodes_m[i];
    }

    /* compute the (+)-u-derivative and project on the (+)-side faces and project q onto the (+)-side faces */
    stride = 0;
    for (i = 0; i < faces_p; i++){
      d4est_operators_apply_dij(d4est_ops, &(e_p[i]->u_elem[0]), (P4EST_DIM), e_p[i]->deg, dir, du_p);
      d4est_linalg_vec_scale(2./e_p[i]->h, du_p, d4est_lgl_get_nodes((P4EST_DIM), e_p[i]->deg));

      d4est_operators_apply_slicer
        (
         d4est_ops,
         du_p,
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &du_p_on_f_p[stride]
        );

      stride += face_nodes_p[i];
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       du_m_on_f_m,
       faces_m,
       deg_m,
       du_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar
      );

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       du_p_on_f_p,
       faces_p,
       deg_p,
       du_p_on_f_p_mortar[d],
       faces_mortar,
       deg_mortar
      );

  }

    /* printf("e_m[0]->id = %d, e_p[0]->id = %d, f_m = %d, f_p = %d \n", e_m[0]->id, e_p[0]->id, f_m, f_p); */
   /* DEBUG_PRINT_ARR_DBL */
   /*    ( */
   /*     u_m_on_f_m, */
   /*     total_side_nodes_m */
   /*    ); */
   /*  DEBUG_PRINT_ARR_DBL */
   /*    ( */
   /*     du_m_on_f_m, */
   /*     total_side_nodes_m */
   /*    ); */
   /*  DEBUG_PRINT_ARR_DBL */
   /*    ( */
   /*     du_m_on_f_m, */
   /*     total_side_nodes_p */
   /*    ); */
   /* DEBUG_PRINT_ARR_DBL */
   /*    ( */
   /*     u_p_on_f_p, */
   /*     total_side_nodes_p */
   /*    ); */
   /*  DEBUG_PRINT_2ARR_DBL */
   /*    ( */
   /*     du_m_on_f_m_mortar, */
   /*     du_p_on_f_p_mortar, */
   /*     total_nodes_mortar */
   /*    ); */
    

  double* Je1 = P4EST_ALLOC(double, total_nodes_mortar);
  double* MJe1 = P4EST_ALLOC(double, total_nodes_mortar);
  double* Je2 = P4EST_ALLOC(double, total_nodes_mortar);
  double* MJe2 = P4EST_ALLOC(double, total_nodes_mortar);
 
    /* calculate symmetric interior penalty flux */
    int k;
    int f;
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (k = 0; k < nodes_mortar[f]; k++){
        int ks = k + stride;
        Je1[ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
        
          Je2[d][ks] = n[d]*u_m_on_f_m_mortar[ks];
          Je2[d][ks] -= n[d]*u_p_on_f_p_mortar[ks];
          Je2[d][ks] *= u_prefactor_mortar[f];
        
          Je1[ks] += gradu_prefactor_mortar[f]
                     *n[dir]*(du_m_on_f_m_mortar[d][ks] - du_p_on_f_p_mortar[d][ks]);
        }
        /* printf("ustar_min_u_mortar[%d]\n", ks, ustar_min_u_mortar[ks]); */
        /* printf("qstar_min_q_mortar[%d]\n", ks, qstar_min_q_mortar[ks]); */
      }
      stride += nodes_mortar[f];
    }

    double sj_m = e_m[0]->surface_jacobian;
    double sj_p = e_p[0]->surface_jacobian;
    double sj_mortar = (sj_m < sj_p) ? sj_m : sj_p;

       
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){

         d4est_operators_apply_mij(d4est_ops,
                         &Je2[d][stride],
                         (P4EST_DIM)-1,
                         deg_mortar[f],
                         &MJe2[stride]);
        
        double Je2MJe2 = sj_mortar*d4est_linalg_vec_dot(&Je2[d][stride], &MJe2[stride], nodes_mortar[f]);

        if(faces_m == (P4EST_HALF)){
          e_m[f]->local_estimator += Je2MJe2;
        }
        else{
          e_m[0]->local_estimator += Je2MJe2;
        }
      }
      stride += nodes_mortar[f];
    }


    stride = 0;
    for (f = 0; f < faces_mortar; f++){  

      d4est_operators_apply_mij(d4est_ops,
                       &Je1[stride],
                       (P4EST_DIM)-1,
                       deg_mortar[f],
                       &MJe1[stride]);
        
      double Je1MJe1 = sj_mortar*d4est_linalg_vec_dot(&Je1[stride],
                                                &MJe1[stride],
                                                nodes_mortar[f]);



      if(faces_m == (P4EST_HALF)){
        e_m[f]->local_estimator += Je1MJe1;
      }
      else{
        e_m[0]->local_estimator += Je1MJe1;
      }
      stride += nodes_mortar_gauss[f];
    }

    for (int d = 0; d < (P4EST_DIM); d++){
      P4EST_FREE(du_m_on_f_m_mortar[d]);
      P4EST_FREE(du_p_on_f_p_mortar[d]);
    }   
    
    P4EST_FREE(MJe2);
    P4EST_FREE(Je2);
    P4EST_FREE(MJe1);
    P4EST_FREE(Je1);
    P4EST_FREE(u_p_on_f_p_mortar);
    P4EST_FREE(u_m_on_f_m_mortar);
    P4EST_FREE(u_p_on_f_p);
    P4EST_FREE(du_p_on_f_p);
    P4EST_FREE(du_m_on_f_m);
    P4EST_FREE(u_m_on_f_m);
    P4EST_FREE(du_p);
    P4EST_FREE(du_m);
    P4EST_FREE(qstar_min_q_m);
    P4EST_FREE(ustar_min_u_m);
    P4EST_FREE(qstar_min_q_mortar);
    P4EST_FREE(ustar_min_u_mortar);
}

flux_fcn_ptrs_t
bi_est_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 double penalty_prefactor
)
{
  flux_fcn_ptrs_t bi_est_fcns;
  bi_est_fcns.flux_interface_fcn
    = bi_est_interface;

  bi_est_fcns.flux_boundary_fcn
    = bi_est_dirichlet;

  bi_est_fcns.bndry_fcn = bndry_fcn;

  bi_est_sipg_flux_penalty_prefactor = penalty_prefactor;

  bi_est_u_prefactor_calculate_fcn = u_penalty_fcn;
  bi_est_u_dirichlet_prefactor_calculate_fcn = u_dirichlet_penalty_fcn;
  bi_est_gradu_prefactor_calculate_fcn = gradu_penalty_fcn;
  
  return bi_est_fcns;
}
