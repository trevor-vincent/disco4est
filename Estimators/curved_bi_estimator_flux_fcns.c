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
 d4est_geometry_t* geom,
 void* params
)
{
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_on_f_m_min_u_at_bndry_Lobatto = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_on_f_m_min_u_at_bndry_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  

  double* MJe2 = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* Je2 = P4EST_ALLOC(double, face_nodes_m_Gauss);

  double* xyz_on_f_m [(P4EST_DIM)];
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  /* double* sj_n_on_f_m_Gauss [(P4EST_DIM)]; */

  
  for (int d = 0; d < (P4EST_DIM); d++) {
    /* xyz_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss); */
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_Lobatto);
    n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    /* sj_n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss); */


    dgmath_apply_slicer(dgmath_jit_dbase,
                        e_m->xyz[d],
                        (P4EST_DIM),
                        f_m,
                        e_m->deg,
                        xyz_on_f_m[d]);
    
  }

  
  
  dgmath_apply_slicer(dgmath_jit_dbase,
                      e_m->u_storage,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);

  /* dgmath_interp_GLL_to_GL */
  /*   ( */
  /*    dgmath_jit_dbase, */
  /*    u_m_on_f_m, */
  /*    e_m->deg, */
  /*    e_m->deg_integ, */
  /*    u_m_on_f_m_Gauss, */
  /*    (P4EST_DIM)-1 */
  /*   ); */
 
  /* curved_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    &e_m, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_Gauss, */
  /*    sj_on_f_m_Gauss, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    xyz_on_f_m_Gauss */
  /*   ); */

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     NULL,
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     NULL,
     NULL,
     GAUSS,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
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
    

  for (int i = 0; i < face_nodes_m_Lobatto; i++){
    u_on_f_m_min_u_at_bndry_Lobatto[i] = u_m_on_f_m[i]
                                         - u_at_bndry
                                         (
                                          xyz_on_f_m[0][i],
                                          xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                                          ,
                                          xyz_on_f_m[2][i]
#endif
                                         );
  }
  dgmath_interp_GLL_to_GL
    (
     dgmath_jit_dbase,
     u_on_f_m_min_u_at_bndry_Lobatto,
     e_m->deg,
     e_m->deg_integ,
     u_on_f_m_min_u_at_bndry_Gauss,
     (P4EST_DIM)-1
    );
  
  for (int dim = 0; dim < (P4EST_DIM); dim++){
    /* calculate qstar - q(-) */

    for(int i = 0; i < face_nodes_m_Gauss; i++){
      Je2[i] = n_on_f_m_Gauss[dim][i]*Je2_prefactor*(u_on_f_m_min_u_at_bndry_Gauss[i]);
     
    }
    
    
    /* dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase, */
    /*                                             Je2, */
    /*                                             e_m->deg_integ, */
    /*                                             sj_n_on_f_m_Gauss[dim], */
    /*                                             e_m->deg_integ, */
    /*                                             (P4EST_DIM)-1, */
    /*                                             MJe2); */
    

    /* double Je2MJe2 = linalg_vec_dot(Je2, MJe2, face_nodes_m_Gauss); */
    

    double Je2MJe2 = dgmath_Gauss_quadrature(
                                             dgmath_jit_dbase,
                                             Je2,
                                             Je2,
                                             sj_on_f_m_Gauss,
                                             e_m->deg_integ,
                                             (P4EST_DIM)-1);
    
    e_m->local_estimator += Je2MJe2;
    /* if (e_m->id == 1){ */
      /* printf("element id = %d, f = %d, Je2MJe2 = %.25f\n", e_m->id, f_m, Je2MJe2); */
      /* printf("element id = %d, f = %d, Je1MJe1 = %.25f\n", e_m->id, f_m, 0.); */
    /* } */
  }


  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(Je2);
  P4EST_FREE(MJe2);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_Lobatto);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_Gauss);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(xyz_on_f_m[d]);
    P4EST_FREE(n_on_f_m_Gauss[d]);
    /* P4EST_FREE(sj_n_on_f_m_Gauss[d]); */
  }
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
 d4est_geometry_t* geom,
 void* params
)
{
  int stride;
  int deg_p_Lobatto [(P4EST_HALF)];
  int face_nodes_p_Lobatto [(P4EST_HALF)];
  int face_nodes_p_Gauss [(P4EST_HALF)];
  int deg_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Gauss [(P4EST_HALF)];
  int nodes_mortar_Gauss [(P4EST_HALF)];
  int nodes_mortar_Lobatto [(P4EST_HALF)];
  int deg_mortar_Gauss [(P4EST_HALF)];
  int deg_mortar_Lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double Je1_prefactor_mortar [(P4EST_HALF)];
  double Je2_prefactor_mortar [(P4EST_HALF)];
  int deg_p_Lobatto_porder [(P4EST_HALF)];
  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  int total_side_nodes_m_Lobatto = 0;
  int total_side_nodes_m_Gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_Lobatto[i] = e_m[i]->deg;
    /* deg_m_Gauss[i] = e_m[i]->deg_integ; */
    
    face_nodes_m_Lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_Gauss[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_integ);
    
    total_side_nodes_m_Lobatto += face_nodes_m_Lobatto[i];
    total_side_nodes_m_Gauss += face_nodes_m_Gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_Lobatto = 0;
  int total_side_nodes_p_Gauss = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_Lobatto[i] = e_p_oriented[i]->deg;
    deg_p_Lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_Gauss[i] = e_p_oriented[i]->deg_integ; */

    face_nodes_p_Lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_Gauss[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_integ);
    
    total_side_nodes_p_Lobatto += face_nodes_p_Lobatto[i];
    total_side_nodes_p_Gauss += face_nodes_p_Gauss[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_Gauss = 0;
  int total_nodes_mortar_Lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_Gauss[i+j] = util_max_int( e_m[i]->deg_integ,
                                            e_p_oriented[j]->deg_integ);
      deg_mortar_Lobatto[i+j] = util_max_int( e_m[i]->deg,
                                            e_p_oriented[j]->deg );      
      nodes_mortar_Gauss[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Gauss[i+j] );     
      nodes_mortar_Lobatto[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Lobatto[i+j] );     
      total_nodes_mortar_Gauss += nodes_mortar_Gauss[i+j];
      total_nodes_mortar_Lobatto += nodes_mortar_Lobatto[i+j];
      Je1_prefactor_mortar[i+j] =  curved_bi_est_gradu_prefactor_calculate_fcn
                                   (
                                    e_m[i]->deg,
                                    (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                                    e_p[j]->deg,
                                    (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                                    curved_bi_est_sipg_flux_penalty_prefactor
                                   );    
    
      Je2_prefactor_mortar[i+j] = curved_bi_est_u_prefactor_calculate_fcn
                                  (
                                   e_m[i]->deg,
                                   (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                                   e_p[j]->deg,
                                   (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                                   curved_bi_est_sipg_flux_penalty_prefactor
                                  );    
      
    }


  int deg_mortar_Gauss_porder [(P4EST_HALF)];
  int nodes_mortar_Gauss_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_Gauss_porder[inew] = deg_mortar_Gauss[i];
    nodes_mortar_Gauss_porder[inew] = nodes_mortar_Gauss[i];
  }

  
  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_m_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss [(P4EST_DIM)];
  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  /* double* sj_n_on_f_m_mortar_Gauss [(P4EST_DIM)];   */


  double* drst_dxyz_m_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_Gauss_porder [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder, total_nodes_mortar_Gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dudr_p_on_f_p_porder[i] = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
    dudr_p_on_f_p_mortar_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudr_p_on_f_p_mortar_Gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_p_on_f_p_mortar_Gauss_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_p_on_f_p_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    
    dudr_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
    dudr_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudr_m_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    dudx_m_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    /* sj_n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss); */
  }
 
  double* Je1 = P4EST_ALLOC_ZERO(double, total_nodes_mortar_Gauss);
  double* MJe1 = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* Je2 [(P4EST_DIM)]; 
  for (int d = 0; d < (P4EST_DIM); d++) {
    Je2[d] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }
  double* MJe2 = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_Gauss);

  stride = 0;
  for (int i = 0; i < faces_m; i++){   
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_m[i]->u_storage[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    stride += face_nodes_m_Lobatto[i];
  }
 
  stride = 0;
  for (int i = 0; i < faces_p; i++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_p_oriented[i]->u_storage[0]),
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
    stride += face_nodes_p_Lobatto[i];
  }

  P4EST_FREE(tmp);

  /* project (-)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     faces_m,
     deg_m_Lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );

  /* project (+)-side u trace vector onto mortar space */
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     u_p_on_f_p,
     faces_p,
     deg_p_Lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_Gauss
    );


  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_m_on_f_m_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_m_on_f_m_mortar_Gauss[stride], (P4EST_DIM)-1);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_p_on_f_p_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_p_on_f_p_mortar_Gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_Gauss[f];
  }
  
  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (int i = 0; i < faces_m; i++){
    dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_m[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &dudr_m_on_f_m[d][stride]
        );
      

      stride += face_nodes_m_Lobatto[i];
    }
    
    /* compute the (+)-u-derivative 
   * and project on the (+)-side faces and project q onto the (+)-side faces */
    stride = 0;
    for (int i = 0; i < faces_p; i++){
      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_apply_slicer
          (
           dgmath_jit_dbase,
           &e_p[i]->dudr_elem[d][0],
           (P4EST_DIM),
           f_p,
           e_p[i]->deg,
           &dudr_p_on_f_p_porder[d][stride]
          );
      }
      stride += dgmath_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_Lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_Gauss_porder
      );
  
    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_Lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_Gauss
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_m_on_f_m_mortar[d][stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &dudr_m_on_f_m_mortar_Gauss[d][stride], (P4EST_DIM)-1);
      stride += nodes_mortar_Gauss[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_p_on_f_p_mortar_porder[d][stride], deg_mortar_Gauss_porder[f], deg_mortar_Gauss_porder[f], &dudr_p_on_f_p_mortar_Gauss_porder[d][stride], (P4EST_DIM)-1);
      }
      stride += nodes_mortar_Gauss_porder[f];
    }
    

  }
  


 /* curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data */
 /*    ( */
 /*     e_m, */
 /*     faces_m, */
 /*     faces_mortar, */
 /*     &deg_mortar_Gauss[0], */
 /*     f_m, */
 /*     drst_dxyz_m_on_mortar_Gauss, */
 /*     sj_on_f_m_mortar_Gauss, */
 /*     n_on_f_m_mortar_Gauss, */
 /*     geom->p4est_geom, */
 /*     dgmath_jit_dbase, */
 /*     NULL */
 /*    ); */

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     drst_dxyz_m_on_mortar_Gauss,
     sj_on_f_m_mortar_Gauss,
     n_on_f_m_mortar_Gauss,
     NULL,
     NULL,
     GAUSS,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     faces_p,
     faces_mortar,
     &deg_mortar_Gauss_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_Gauss_porder,
     NULL,
     NULL,
     NULL,
     NULL,
     GAUSS,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  

  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_Gauss[d],
       0.0,
       total_nodes_mortar_Gauss
      );

    linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_Gauss_porder[d],
       0.0,
       total_nodes_mortar_Gauss
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_Gauss; k++){
        dudx_m_on_f_m_mortar_Gauss[j][k] +=
          drst_dxyz_m_on_mortar_Gauss[i][j][k]*dudr_m_on_f_m_mortar_Gauss[i][k];
        dudx_p_on_f_p_mortar_Gauss_porder[j][k] += drst_dxyz_p_on_mortar_Gauss_porder[i][j][k]*dudr_p_on_f_p_mortar_Gauss_porder[i][k];
      }
    }    
  }


  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
      face_p = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &dudx_p_on_f_p_mortar_Gauss_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_Gauss[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_Gauss[d][face_mortar_stride]
        );
    }
    
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]);
  }


  
  
  /* DEBUG_PRINT_ARR_DBL */
  /*   ( */
  /*    u_m_on_f_m, */
  /*    total_side_nodes_m_Lobatto */
  /*   ); */
  /* DEBUG_PRINT_ARR_DBL */
  /*   ( */
  /*    u_p_on_f_p, */
  /*    total_side_nodes_p_Lobatto */
  /*   ); */
  /* DEBUG_PRINT_3ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m[0], */
  /*    dudr_m_on_f_m[1], */
  /*    dudr_m_on_f_m[2], */
  /*    total_side_nodes_m_Lobatto */
  /*   ); */
  /* DEBUG_PRINT_3ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m[0], */
  /*    dudr_m_on_f_m[1], */
  /*    dudr_m_on_f_m[2], */
  /*    total_side_nodes_p_Lobatto */
  /*   ); */

  /* DEBUG_PRINT_6ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m_mortar[0], */
  /*    dudr_m_on_f_m_mortar[1], */
  /*    dudr_m_on_f_m_mortar[2], */
  /*    dudr_p_on_f_p_mortar[0], */
  /*    dudr_p_on_f_p_mortar[1], */
  /*    dudr_p_on_f_p_mortar[2], */
  /*    total_nodes_mortar_Gauss */
  /*   ); */


  /* for (int i = 0; i < total_nodes_mortar_Gauss; i++){ */
  /*   printf("dudx_mortar_Lobatto = %.25f %.25f %.25f\n", 4*dudr_m_on_f_m_mortar[0][i],4*dudr_m_on_f_m_mortar[1][i], 4*dudr_m_on_f_m_mortar[2][i]);     */
  /* } */

  /* double dudx_div_dudr = dudx_m_on_f_m_mortar_Gauss[0][0]/dudr_m_on_f_m_mortar_Gauss[0][0]; */
  /* double nx = n_on_f_m_mortar_Gauss[0][0]; */
  /* double ny = n_on_f_m_mortar_Gauss[1][0]; */
  /* double nz = n_on_f_m_mortar_Gauss[2][0]; */
  
  
  /* double* tmpxyz [(P4EST_DIM)]; */
  /* D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_Gauss); */
  /* curved_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    e_m, */
  /*    faces_m, */
  /*    faces_mortar, */
  /*    &deg_mortar_Gauss[0], */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_mortar_Gauss, */
  /*    sj_on_f_m_mortar_Gauss, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    tmpxyz */
  /*   ); */
  /* D4EST_FREE_DIM_VEC(tmpxyz); */
  

  /* for(int d = 0; d < (P4EST_DIM); d++){ */
  /*   stride = 0; */
  /*   for (int f = 0; f < faces_mortar; f++){ */
  /*     dgmath_reorient_face_data */
  /*       ( */
  /*        dgmath_jit_dbase, */
  /*        &dudx_p_on_f_p_mortar_Gauss[d][stride], */
  /*        ((P4EST_DIM) - 1), */
  /*        deg_mortar_Gauss[f], */
  /*        orientation, */
  /*        f_m, */
  /*        f_p, */
  /*        &dudx_p_on_f_p_mortar_Gauss_reoriented[d][stride] */
  /*       );       */
  /*     stride += nodes_mortar_Gauss[f]; */
  /*   }    */
  /* } */



    /* calculate symmetric interior penalty flux */
    int k;
    int f;
    int ks;
    double n_ks;
    double sj_ks;
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (k = 0; k < nodes_mortar_Gauss[f]; k++){
        ks = k + stride;
        sj_ks = sj_on_f_m_mortar_Gauss[ks];
        Je1[ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          n_ks = n_on_f_m_mortar_Gauss[d][ks];        
          Je1[ks] += Je1_prefactor_mortar[f]*n_ks*
                     (dudx_m_on_f_m_mortar_Gauss[d][ks] - dudx_p_on_f_p_mortar_Gauss[d][ks]);
          /* sje1[ks] += sj_ks*Je1_prefactor_mortar[f]*n_ks* */
                      /* (du_m_on_f_m_mortar[d][ks] - du_p_on_f_p_mortar[d][ks]); */

          Je2[d][ks] = n_ks*u_m_on_f_m_mortar_Gauss[ks];
          Je2[d][ks] -= n_ks*u_p_on_f_p_mortar_Gauss[ks];
          Je2[d][ks] *= Je2_prefactor_mortar[f];
          
        }
        /* printf("Je1_test = %.25f\n", Je1_test); */
      }
      stride += nodes_mortar_Gauss[f];
    }


    
    /* the contribution in every direction must be added up due to it being a vector norm */
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        


        double Je2MJe2 = dgmath_Gauss_quadrature(
                                                 dgmath_jit_dbase,
                                                 &Je2[d][stride],
                                                 &Je2[d][stride],
                                                 &sj_on_f_m_mortar_Gauss[stride],
                                                 deg_mortar_Gauss[f],
                                                 (P4EST_DIM)-1);
        

        /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */

        if(faces_m == (P4EST_HALF)){
          e_m[f]->local_estimator += Je2MJe2;
        }
        else{
          e_m[0]->local_estimator += Je2MJe2;
        }
      }
      stride += nodes_mortar_Gauss[f];
    }


    
  stride = 0;
  for (f = 0; f < faces_mortar; f++){  

    /* dgmath_apply_curvedGaussMass_onGaussNodeVec */
    /*   ( */
    /*    dgmath_jit_dbase, */
    /*    &Je1[stride], */
    /*    deg_mortar_Gauss[f], */
    /*    &sj_on_f_m_mortar_Gauss[stride], */
    /*    deg_mortar_Gauss[f], */
    /*    (P4EST_DIM)-1, */
    /*    &MJe1[stride] */
    /*   ); */

    /* double Je1MJe1 = linalg_vec_dot(&Je1[stride], &MJe1[stride], nodes_mortar_Gauss[f]); */

    double Je1MJe1 = dgmath_Gauss_quadrature(
                                             dgmath_jit_dbase,
                                             &Je1[stride],
                                             &Je1[stride],
                                             &sj_on_f_m_mortar_Gauss[stride],
                                             deg_mortar_Gauss[f],
                                             (P4EST_DIM)-1);


    /* dgmath_apply_Mij(dgmath_jit_dbase, &Je1_test_mortar[stride], (P4EST_DIM)-1, deg_mortar_Gauss[f], &MJe1_test_mortar[stride]); */
    /* linalg_vec_scale(sj_on_f_m_mortar_Gauss[0], &MJe1_test_mortar[stride], nodes_mortar_Gauss[f]); */
    /* double Je1MJe1_test = linalg_vec_dot(&Je1_test_mortar[stride], &MJe1_test_mortar[stride], nodes_mortar_Gauss[f]); */


    /* if (fabs(Je1MJe1 - Je1MJe1_test) > .00001){ */
    /* } */
    /* printf("id, f, Je1MJe1, Je1MJe1_test = %d %d %.25f, %.25f\n", e_m[0]->id, f_m, Je1MJe1, Je1MJe1_test); */

    
    /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */
      if(faces_m == (P4EST_HALF)){
        e_m[f]->local_estimator += Je1MJe1;
      }
      else{
        e_m[0]->local_estimator += Je1MJe1;
      }
      stride += nodes_mortar_Gauss[f];
  }


  /* for (int f = 0; f < faces_m; f++){ */
    /* printf("eta2 on face %d = %.25f\n", f, e_m[f]->local_estimator); */
  /* } */
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder);
  
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_Gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(Je2[i]);
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_Gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_Gauss_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_Gauss[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_Gauss[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_Gauss[i]);
    P4EST_FREE(n_on_f_m_mortar_Gauss[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  P4EST_FREE(Je1);
  P4EST_FREE(MJe1);
  P4EST_FREE(MJe2);

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
