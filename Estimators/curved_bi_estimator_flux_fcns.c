#include "../Estimators/curved_bi_estimator_flux_fcns.h"
#include "../Estimators/bi_estimator_flux_fcns.h"
#include "../LinearAlgebra/linalg.h"

/* #define D4EST_DEBUG */

double curved_bi_est_sipg_flux_penalty_prefactor;
penalty_calc_t curved_bi_est_u_prefactor_calculate_fcn;
penalty_calc_t curved_bi_est_u_dirichlet_prefactor_calculate_fcn;
penalty_calc_t curved_bi_est_gradu_prefactor_calculate_fcn;

static void
curved_bi_est_dirichlet
(
 curved_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_geometry_t* geom,
 void* params
)
{
  
  grid_fcn_t u_at_bndry = bndry_fcn;
  /* int vol_nodes_m = dgmath_get_nodes ( (P4EST_DIM) , e_m->deg ); */

  int face_nodes_m = dgmath_get_nodes( (P4EST_DIM) - 1, e_m->deg);
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);

  double* Je2 = P4EST_ALLOC(double, face_nodes_m);
  double* sje2 = P4EST_ALLOC(double, face_nodes_m);
  double* MJe2 = P4EST_ALLOC(double, face_nodes_m);

  double* xyz_on_f_m [(P4EST_DIM)];
  double* n_on_f_m [(P4EST_DIM)];
  double* du_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);
  double* sj_on_f_m = P4EST_ALLOC(double, face_nodes_m);

  for (int d = 0; d < (P4EST_DIM); d++){
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m);
    n_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m);
  }

  
  dgmath_apply_slicer(dgmath_jit_dbase,
                      e_m->u_elem,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);


  for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &e_m->xyz[d][0],
       (P4EST_DIM),
       f_m,
       e_m->deg,
       xyz_on_f_m[d]
      );
  }

  
  curved_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     &e_m,
     1,
     1,
     &e_m->deg,
     f_m,
     n_on_f_m,
     sj_on_f_m,
     geom,
     dgmath_jit_dbase
    );  

  
  double h = (e_m->volume/e_m->surface_area[f_m]);
  /* find IP penalty parameter for each face pair of the two sides*/
  double Je2_prefactor = curved_bi_est_u_dirichlet_prefactor_calculate_fcn
                  (
                   e_m->deg,
                   h,
                   e_m->deg,
                   h,
                   curved_bi_est_sipg_flux_penalty_prefactor
                  );
    
  int dim;
  for (dim = 0; dim < (P4EST_DIM); dim++){
    /* calculate qstar - q(-) */
    for(int i = 0; i < face_nodes_m; i++){
      Je2[i] = Je2_prefactor*n_on_f_m[dim][i]*(u_m_on_f_m[i]
                                              -
                                              u_at_bndry
                                              (
                                               xyz_on_f_m[0][i],
                                               xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                                                ,
                                               xyz_on_f_m[2][i]
#endif
                                              )
                                             );
      sje2[i] = sj_on_f_m[i]*Je2[i];
    }

    dgmath_apply_Mij
      (
       dgmath_jit_dbase, 
       sje2,
       (P4EST_DIM) - 1,
       e_m->deg,
       MJe2
      );

    double Je2MJe2 = linalg_vec_dot(Je2, MJe2, face_nodes_m);
    e_m->local_estimator += Je2MJe2;
  }

  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m[d]);
    P4EST_FREE(xyz_on_f_m[d]);
  }

  P4EST_FREE(du_m_on_f_m);
  P4EST_FREE(Je2);
  P4EST_FREE(sje2);
  P4EST_FREE(MJe2);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m);
}

static void
curved_bi_est_interface
(
 curved_element_data_t** e_m,
 int faces_m,
 int f_m,
 curved_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_geometry_t* geom,
 void* params
)
{ 
  int stride;
  int deg_p [(P4EST_HALF)];
  int max_deg_p = -1;
  int face_nodes_p [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int face_nodes_m [(P4EST_HALF)];
  int max_deg_m = -1;
  int nodes_mortar [(P4EST_HALF)];
  int deg_mortar [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double Je1_prefactor_mortar [(P4EST_HALF)];
  double Je2_prefactor_mortar [(P4EST_HALF)];

  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m[i] = e_m[i]->deg;
    if (e_m[i]->deg > max_deg_m) max_deg_m = e_m[i]->deg;
    face_nodes_m[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg );
    total_side_nodes_m += face_nodes_m[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p[i] = e_p_oriented[i]->deg;
    if (e_p_oriented[i]->deg > max_deg_p) max_deg_p = e_p_oriented[i]->deg;
    face_nodes_p[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    total_side_nodes_p += face_nodes_p[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar[i+j] = util_max_int( e_m[i]->deg, e_p_oriented[j]->deg );
      nodes_mortar[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar[i+j] );     
      total_nodes_mortar += nodes_mortar[i+j];
      double h = util_min((e_m[i]->volume/e_m[i]->surface_area[f_m]), (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]));

      Je1_prefactor_mortar[i+j] =  curved_bi_est_gradu_prefactor_calculate_fcn
                                   (
                                    deg_mortar[i+j],
                                    h,
                                    deg_mortar[i+j],
                                    h,
                                    curved_bi_est_sipg_flux_penalty_prefactor
                                   );    
    
      Je2_prefactor_mortar[i+j] = curved_bi_est_u_prefactor_calculate_fcn
                                  (
                                   deg_mortar[i+j],
                                   h,
                                   deg_mortar[i+j],
                                   h,
                                   curved_bi_est_sipg_flux_penalty_prefactor
                                  );    
      
    }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  double* du_m_on_f_m [(P4EST_DIM)]; 
  double* du_p_on_f_p [(P4EST_DIM)]; 

  for (int d = 0; d < (P4EST_DIM); d++){
    du_m_on_f_m[d] = P4EST_ALLOC(double, total_side_nodes_m);
    du_p_on_f_p[d] = P4EST_ALLOC(double, total_side_nodes_p);
  }
  
  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* du_m_on_f_m_mortar [(P4EST_DIM)];
  double* du_p_on_f_p_mortar [(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    du_m_on_f_m_mortar[d] = P4EST_ALLOC(double, total_nodes_mortar);
    du_p_on_f_p_mortar[d] = P4EST_ALLOC(double, total_nodes_mortar);
  }
  
  double* Je1 = P4EST_ALLOC_ZERO(double, total_nodes_mortar);
  double* sje1 = P4EST_ALLOC_ZERO(double, total_nodes_mortar);
  double* MJe1 = P4EST_ALLOC(double, total_nodes_mortar);

  double* Je2 [(P4EST_DIM)]; /* = P4EST_ALLOC_ZERO(double, total_nodes_mortar); */
  double* sje2 [(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++) {
    Je2[d] = P4EST_ALLOC(double, total_nodes_mortar);
    sje2[d] = P4EST_ALLOC(double, total_nodes_mortar);
  }

  double* MJe2 = P4EST_ALLOC(double, total_nodes_mortar);
  
  double* tmp = P4EST_ALLOC(double, total_side_nodes_p);
  double* sj_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* n_on_f_m_mortar [(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    n_on_f_m_mortar[d] = P4EST_ALLOC(double, total_nodes_mortar);
  }

  stride = 0;
  for (int i = 0; i < faces_m; i++){   
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_m[i]->u_elem[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    stride += face_nodes_m[i];
  }
 
  stride = 0;
  for (int i = 0; i < faces_p; i++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_p_oriented[i]->u_elem[0]),
       (P4EST_DIM),
       f_p,
       e_p_oriented[i]->deg,
       tmp
      );
    
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       tmp,
       ((P4EST_DIM) - 1),
       e_p_oriented[i]->deg,
       orientation,
       f_m,
       f_p,
       &u_p_on_f_p[stride]
      );
    stride += face_nodes_p[i];
  }


  /* project (-)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     faces_m,
     deg_m,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar
    );

  /* project (+)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_p_on_f_p,
     faces_p,
     deg_p,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar
    );
  
  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (int i = 0; i < faces_m; i++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_m[i]->du_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &du_m_on_f_m[d][stride]
        );
    
      stride += face_nodes_m[i];
    }
    
    for (int i = 0; i < faces_p; i++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_p_oriented[i]->du_elem[d][0],
         (P4EST_DIM),
         f_p,
         e_p_oriented[i]->deg,
         tmp
        );

      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         tmp,
         ((P4EST_DIM) - 1),
         e_p_oriented[i]->deg,
         orientation,
         f_m,
         f_p,
         &du_p_on_f_p[d][stride]
        );

      stride += face_nodes_p[i];
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       du_m_on_f_m[d],
       faces_m,
       deg_m,
       du_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar
      );

    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       du_p_on_f_p[d],
       faces_p,
       deg_p,
       du_p_on_f_p_mortar[d],
       faces_mortar,
       deg_mortar
      );

  }
  
  curved_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     e_m,
     faces_m,
     faces_mortar,
     &deg_mortar[0],
     f_m,
     n_on_f_m_mortar,
     sj_on_f_m_mortar,
     geom,
     dgmath_jit_dbase
    );

    /* calculate symmetric interior penalty flux */
    int k;
    int f;
    int ks;
    double n_ks;
    double sj_ks;
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (k = 0; k < nodes_mortar[f]; k++){
        ks = k + stride;
        sj_ks = sj_on_f_m_mortar[ks];
        Je1[ks] = 0.;
        sje1[ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          n_ks = n_on_f_m_mortar[d][ks];        
          Je1[ks] += Je1_prefactor_mortar[f]*n_ks*
                     (du_m_on_f_m_mortar[d][ks] - du_p_on_f_p_mortar[d][ks]);
          sje1[ks] += sj_ks*Je1_prefactor_mortar[f]*n_ks*
                      (du_m_on_f_m_mortar[d][ks] - du_p_on_f_p_mortar[d][ks]);

          Je2[d][ks] = n_ks*u_m_on_f_m_mortar[ks];
          Je2[d][ks] -= n_ks*u_p_on_f_p_mortar[ks];
          Je2[d][ks] *= Je2_prefactor_mortar[f];
          sje2[d][ks] = sj_ks*Je2[d][ks];
        }
      }
      stride += nodes_mortar[f];
    }

    /* the contribution in every direction must be added up due to it being a vector norm */
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_apply_Mij
          (
           dgmath_jit_dbase,
           &sje2[d][stride],
           (P4EST_DIM) - 1,
           deg_mortar[f],
           &MJe2[stride]
          );

        /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */
        double Je2MJe2 = linalg_vec_dot(&Je2[d][stride], &MJe2[stride], nodes_mortar[f]);
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
    dgmath_apply_Mij
      (
       dgmath_jit_dbase,
       &sje1[stride],
       (P4EST_DIM) - 1,
       deg_mortar[f],
       &MJe1[stride]
      );

    double Je1MJe1 = linalg_vec_dot(&Je1[stride], &MJe1[stride], nodes_mortar[f]);

    /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */
      if(faces_m == (P4EST_HALF)){
        e_m[f]->local_estimator += Je1MJe1;   
      }
      else{
        e_m[0]->local_estimator += Je1MJe1;   
      }
    stride += nodes_mortar[f];
  }

  P4EST_FREE(tmp);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar[i]);
    P4EST_FREE(du_m_on_f_m_mortar[i]);
    P4EST_FREE(du_m_on_f_m[i]);
    P4EST_FREE(du_p_on_f_p[i]);
    P4EST_FREE(du_p_on_f_p_mortar[i]);
    P4EST_FREE(Je2[i]);
    P4EST_FREE(sje2[i]);
  }
  
  P4EST_FREE(sje1);
  P4EST_FREE(Je1);
  P4EST_FREE(MJe1);
  P4EST_FREE(MJe2);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(sj_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m);
}

curved_flux_fcn_ptrs_t
curved_bi_est_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 double penalty_prefactor
)
{
  curved_flux_fcn_ptrs_t curved_bi_est_fcns;
  curved_bi_est_fcns.flux_interface_fcn
    = curved_bi_est_interface;

  curved_bi_est_fcns.flux_boundary_fcn
    = curved_bi_est_dirichlet;

  curved_bi_est_fcns.bndry_fcn = bndry_fcn;

  curved_bi_est_sipg_flux_penalty_prefactor = penalty_prefactor;

  curved_bi_est_u_prefactor_calculate_fcn = u_penalty_fcn;
  curved_bi_est_u_dirichlet_prefactor_calculate_fcn = u_dirichlet_penalty_fcn;
  curved_bi_est_gradu_prefactor_calculate_fcn = gradu_penalty_fcn;
  
  return curved_bi_est_fcns;
}
