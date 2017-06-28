#include "../Estimators/curved_bi_estimator_flux_fcns.h"
#include "../Estimators/bi_estimator_flux_fcns.h"
#include "../LinearAlgebra/d4est_linalg.h"

/* #define D4EST_DEBUG */

double curved_bi_est_sipg_flux_penalty_prefactor;
penalty_calc_t curved_bi_est_u_prefactor_calculate_fcn;
penalty_calc_t curved_bi_est_u_dirichlet_prefactor_calculate_fcn;
penalty_calc_t curved_bi_est_gradu_prefactor_calculate_fcn;

static void
curved_bi_est_dirichlet
(
 d4est_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_on_f_m_min_u_at_bndry_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  

  double* MJe2 = P4EST_ALLOC(double, face_nodes_m_quad);
  double* sj_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* Je2 = P4EST_ALLOC(double, face_nodes_m_quad);

  double* xyz_on_f_m [(P4EST_DIM)];
  double* n_on_f_m_quad [(P4EST_DIM)];
  /* double* sj_n_on_f_m_quad [(P4EST_DIM)]; */

  
  for (int d = 0; d < (P4EST_DIM); d++) {
    /* xyz_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad); */
    xyz_on_f_m[d] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    n_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad);
    /* sj_n_on_f_m_quad[d] = P4EST_ALLOC(double, face_nodes_m_quad); */


    d4est_operators_apply_slicer(d4est_ops,
                        e_m->xyz[d],
                        (P4EST_DIM),
                        f_m,
                        e_m->deg,
                        xyz_on_f_m[d]);
    
  }

  
  
  d4est_operators_apply_slicer(d4est_ops,
                      e_m->u_elem,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);

  /* d4est_operators_interp_lobatto_to_GL */
  /*   ( */
  /*    d4est_ops, */
  /*    u_m_on_f_m, */
  /*    e_m->deg, */
  /*    e_m->deg_quad, */
  /*    u_m_on_f_m_quad, */
  /*    (P4EST_DIM)-1 */
  /*   ); */
 
  /* d4est_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    &e_m, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_quad, */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_quad, */
  /*    sj_on_f_m_quad, */
  /*    geom, */
  /*    d4est_ops, */
  /*    xyz_on_f_m_quad */
  /*   ); */

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     NULL,
     sj_on_f_m_quad,
     n_on_f_m_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     d4est_ops,
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
    

  for (int i = 0; i < face_nodes_m_lobatto; i++){
    u_on_f_m_min_u_at_bndry_lobatto[i] = u_m_on_f_m[i]
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


  
    d4est_operators_interp(d4est_ops,
                  u_on_f_m_min_u_at_bndry_lobatto,
                  QUAD_LOBATTO,
                  e_m->deg,
                  u_on_f_m_min_u_at_bndry_quad,
                  geom->geom_quad_type,
                  e_m->deg_quad,
                  (P4EST_DIM)-1);

  
  for (int dim = 0; dim < (P4EST_DIM); dim++){
    /* calculate qstar - q(-) */

    for(int i = 0; i < face_nodes_m_quad; i++){
      Je2[i] = n_on_f_m_quad[dim][i]*Je2_prefactor*(u_on_f_m_min_u_at_bndry_quad[i]);
     
    }
    
    
    /* d4est_operators_apply_curvedquadMass_onquadNodeVec(d4est_ops, */
    /*                                             Je2, */
    /*                                             e_m->deg_quad, */
    /*                                             sj_n_on_f_m_quad[dim], */
    /*                                             e_m->deg_quad, */
    /*                                             (P4EST_DIM)-1, */
    /*                                             MJe2); */
    

    /* double Je2MJe2 = d4est_linalg_vec_dot(Je2, MJe2, face_nodes_m_quad); */
    

    double Je2MJe2 = d4est_operators_quadrature(
                                       d4est_ops,
                                       Je2,
                                       Je2,
                                       sj_on_f_m_quad,
                                       e_m->deg_quad,
                                       geom->geom_quad_type,
                                       (P4EST_DIM)-1);
    
    e_m->local_estimator += Je2MJe2;
    /* if (e_m->id == 1){ */
      /* printf("element id = %d, f = %d, Je2MJe2 = %.25f\n", e_m->id, f_m, Je2MJe2); */
      /* printf("element id = %d, f = %d, Je1MJe1 = %.25f\n", e_m->id, f_m, 0.); */
    /* } */
  }


  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m_quad);
  P4EST_FREE(Je2);
  P4EST_FREE(MJe2);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_lobatto);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_quad);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(xyz_on_f_m[d]);
    P4EST_FREE(n_on_f_m_quad[d]);
    /* P4EST_FREE(sj_n_on_f_m_quad[d]); */
  }
}

static void
curved_bi_est_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{
  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_quad [(P4EST_HALF)];
  int nodes_mortar_quad [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double Je1_prefactor_mortar [(P4EST_HALF)];
  double Je2_prefactor_mortar [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    /* deg_m_quad[i] = e_m[i]->deg_quad; */
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_quad[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_quad = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                            e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
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


  int deg_mortar_quad_porder [(P4EST_HALF)];
  int nodes_mortar_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    nodes_mortar_quad_porder[inew] = nodes_mortar_quad[i];
  }

  
  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_m_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* u_p_on_f_p_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];
  /* double* sj_n_on_f_m_mortar_quad [(P4EST_DIM)];   */


  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_quad, total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder, total_nodes_mortar_quad);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dudr_p_on_f_p_porder[i] = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
    dudr_p_on_f_p_mortar_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudr_p_on_f_p_mortar_quad_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_p_on_f_p_mortar_quad_porder[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_p_on_f_p_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    
    dudr_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
    dudr_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudr_m_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    dudx_m_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    n_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad);
    /* sj_n_on_f_m_mortar_quad[i] = P4EST_ALLOC(double, total_nodes_mortar_quad); */
  }
 
  double* Je1 = P4EST_ALLOC_ZERO(double, total_nodes_mortar_quad);
  double* MJe1 = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* Je2 [(P4EST_DIM)]; 
  for (int d = 0; d < (P4EST_DIM); d++) {
    Je2[d] = P4EST_ALLOC(double, total_nodes_mortar_quad);
  }
  double* MJe2 = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_quad);

  stride = 0;
  for (int i = 0; i < faces_m; i++){   
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_m[i]->u_elem[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    stride += face_nodes_m_lobatto[i];
  }
 
  stride = 0;
  for (int i = 0; i < faces_p; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_p_oriented[i]->u_elem[0]),
       (P4EST_DIM),
       f_p,
       e_p_oriented[i]->deg,
       tmp
      );
    
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       tmp,
       ((P4EST_DIM) - 1),
       e_p_oriented[i]->deg,
       orientation,
       f_m,
       f_p,
       &u_p_on_f_p[stride]
      );
    stride += face_nodes_p_lobatto[i];
  }

  P4EST_FREE(tmp);

  /* project (-)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m_lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_quad
    );

  /* project (+)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p_lobatto,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar_quad
    );


  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    d4est_operators_interp(d4est_ops,
                  &u_m_on_f_m_mortar[stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &u_m_on_f_m_mortar_quad[stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);
    
    d4est_operators_interp(d4est_ops,
                  &u_p_on_f_p_mortar[stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &u_p_on_f_p_mortar_quad[stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);
    
    stride += nodes_mortar_quad[f];
  }
  
  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (int i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_m[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &dudr_m_on_f_m[d][stride]
        );
      

      stride += face_nodes_m_lobatto[i];
    }
    
    /* compute the (+)-u-derivative 
   * and project on the (+)-side faces and project q onto the (+)-side faces */
    stride = 0;
    for (int i = 0; i < faces_p; i++){
        d4est_operators_apply_slicer
          (
           d4est_ops,
           &e_p[i]->dudr_elem[d][0],
           (P4EST_DIM),
           f_p,
           e_p[i]->deg,
           &dudr_p_on_f_p_porder[d][stride]
          );
      stride += d4est_lgl_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_p_on_f_p_porder[d],
       faces_p,
       deg_p_lobatto_porder,
       dudr_p_on_f_p_mortar_porder[d],
       faces_mortar,
       deg_mortar_quad_porder
      );
  
    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       dudr_m_on_f_m[d],
       faces_m,
       deg_m_lobatto,
       dudr_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_quad
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

      d4est_operators_interp(d4est_ops,
                    &dudr_m_on_f_m_mortar[d][stride],
                    QUAD_LOBATTO,
                    deg_mortar_quad[f],
                    &dudr_m_on_f_m_mortar_quad[d][stride],
                    geom->geom_quad_type,
                    deg_mortar_quad[f],
                    (P4EST_DIM)-1);
      
      stride += nodes_mortar_quad[f];
    }

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){

      d4est_operators_interp(d4est_ops,
                    &dudr_p_on_f_p_mortar_porder[d][stride],
                    QUAD_LOBATTO,
                    deg_mortar_quad_porder[f],
                    &dudr_p_on_f_p_mortar_quad_porder[d][stride],
                    geom->geom_quad_type,
                    deg_mortar_quad_porder[f],
                    (P4EST_DIM)-1);
              
      stride += nodes_mortar_quad_porder[f];
    }
    

  }
  


 /* curved_data_compute_drst_dxyz_quad_on_mortar_using_volume_data */
 /*    ( */
 /*     e_m, */
 /*     faces_m, */
 /*     faces_mortar, */
 /*     &deg_mortar_quad[0], */
 /*     f_m, */
 /*     drst_dxyz_m_on_mortar_quad, */
 /*     sj_on_f_m_mortar_quad, */
 /*     n_on_f_m_mortar_quad, */
 /*     geom->p4est_geom, */
 /*     d4est_ops, */
 /*     NULL */
 /*    ); */

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     drst_dxyz_m_on_mortar_quad,
     sj_on_f_m_mortar_quad,
     n_on_f_m_mortar_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     faces_p,
     faces_mortar,
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder,
     NULL,
     NULL,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     d4est_ops,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad[d],
       0.0,
       total_nodes_mortar_quad
      );

    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_quad_porder[d],
       0.0,
       total_nodes_mortar_quad
      );
  }
  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_quad; k++){
        dudx_m_on_f_m_mortar_quad[j][k] +=
          drst_dxyz_m_on_mortar_quad[i][j][k]*dudr_m_on_f_m_mortar_quad[i][k];
        dudx_p_on_f_p_mortar_quad_porder[j][k] += drst_dxyz_p_on_mortar_quad_porder[i][j][k]*dudr_p_on_f_p_mortar_quad_porder[i][k];
      }
    }    
  }


  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
      face_p = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
    }


    for (int d = 0; d < (P4EST_DIM); d++){
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &dudx_p_on_f_p_mortar_quad_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_quad[d][face_mortar_stride]
        );
    }
    
    face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
  }


  
  
  /* DEBUG_PRINT_ARR_DBL */
  /*   ( */
  /*    u_m_on_f_m, */
  /*    total_side_nodes_m_lobatto */
  /*   ); */
  /* DEBUG_PRINT_ARR_DBL */
  /*   ( */
  /*    u_p_on_f_p, */
  /*    total_side_nodes_p_lobatto */
  /*   ); */
  /* DEBUG_PRINT_3ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m[0], */
  /*    dudr_m_on_f_m[1], */
  /*    dudr_m_on_f_m[2], */
  /*    total_side_nodes_m_lobatto */
  /*   ); */
  /* DEBUG_PRINT_3ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m[0], */
  /*    dudr_m_on_f_m[1], */
  /*    dudr_m_on_f_m[2], */
  /*    total_side_nodes_p_lobatto */
  /*   ); */

  /* DEBUG_PRINT_6ARR_DBL */
  /*   ( */
  /*    dudr_m_on_f_m_mortar[0], */
  /*    dudr_m_on_f_m_mortar[1], */
  /*    dudr_m_on_f_m_mortar[2], */
  /*    dudr_p_on_f_p_mortar[0], */
  /*    dudr_p_on_f_p_mortar[1], */
  /*    dudr_p_on_f_p_mortar[2], */
  /*    total_nodes_mortar_quad */
  /*   ); */


  /* for (int i = 0; i < total_nodes_mortar_quad; i++){ */
  /*   printf("dudx_mortar_lobatto = %.25f %.25f %.25f\n", 4*dudr_m_on_f_m_mortar[0][i],4*dudr_m_on_f_m_mortar[1][i], 4*dudr_m_on_f_m_mortar[2][i]);     */
  /* } */

  /* double dudx_div_dudr = dudx_m_on_f_m_mortar_quad[0][0]/dudr_m_on_f_m_mortar_quad[0][0]; */
  /* double nx = n_on_f_m_mortar_quad[0][0]; */
  /* double ny = n_on_f_m_mortar_quad[1][0]; */
  /* double nz = n_on_f_m_mortar_quad[2][0]; */
  
  
  /* double* tmpxyz [(P4EST_DIM)]; */
  /* D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_quad); */
  /* d4est_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    e_m, */
  /*    faces_m, */
  /*    faces_mortar, */
  /*    &deg_mortar_quad[0], */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_mortar_quad, */
  /*    sj_on_f_m_mortar_quad, */
  /*    geom, */
  /*    d4est_ops, */
  /*    tmpxyz */
  /*   ); */
  /* D4EST_FREE_DIM_VEC(tmpxyz); */
  

  /* for(int d = 0; d < (P4EST_DIM); d++){ */
  /*   stride = 0; */
  /*   for (int f = 0; f < faces_mortar; f++){ */
  /*     d4est_operators_reorient_face_data */
  /*       ( */
  /*        d4est_ops, */
  /*        &dudx_p_on_f_p_mortar_quad[d][stride], */
  /*        ((P4EST_DIM) - 1), */
  /*        deg_mortar_quad[f], */
  /*        orientation, */
  /*        f_m, */
  /*        f_p, */
  /*        &dudx_p_on_f_p_mortar_quad_reoriented[d][stride] */
  /*       );       */
  /*     stride += nodes_mortar_quad[f]; */
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
      for (k = 0; k < nodes_mortar_quad[f]; k++){
        ks = k + stride;
        sj_ks = sj_on_f_m_mortar_quad[ks];
        Je1[ks] = 0.;
        for (int d = 0; d < (P4EST_DIM); d++){
          n_ks = n_on_f_m_mortar_quad[d][ks];        
          Je1[ks] += Je1_prefactor_mortar[f]*n_ks*
                     (dudx_m_on_f_m_mortar_quad[d][ks] - dudx_p_on_f_p_mortar_quad[d][ks]);
          /* sje1[ks] += sj_ks*Je1_prefactor_mortar[f]*n_ks* */
                      /* (du_m_on_f_m_mortar[d][ks] - du_p_on_f_p_mortar[d][ks]); */

          Je2[d][ks] = n_ks*u_m_on_f_m_mortar_quad[ks];
          Je2[d][ks] -= n_ks*u_p_on_f_p_mortar_quad[ks];
          Je2[d][ks] *= Je2_prefactor_mortar[f];
          
        }
        /* printf("Je1_test = %.25f\n", Je1_test); */
      }
      stride += nodes_mortar_quad[f];
    }


    
    /* the contribution in every direction must be added up due to it being a vector norm */
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (int d = 0; d < (P4EST_DIM); d++){
        


        double Je2MJe2 = d4est_operators_quadrature(
                                                 d4est_ops,
                                                 &Je2[d][stride],
                                                 &Je2[d][stride],
                                                 &sj_on_f_m_mortar_quad[stride],
                                                 deg_mortar_quad[f],
                                                 geom->geom_quad_type,
                                                 (P4EST_DIM)-1);
        

        /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */

        if(faces_m == (P4EST_HALF)){
          e_m[f]->local_estimator += Je2MJe2;
        }
        else{
          e_m[0]->local_estimator += Je2MJe2;
        }
      }
      stride += nodes_mortar_quad[f];
    }


    
  stride = 0;
  for (f = 0; f < faces_mortar; f++){  

    /* d4est_operators_apply_curvedquadMass_onquadNodeVec */
    /*   ( */
    /*    d4est_ops, */
    /*    &Je1[stride], */
    /*    deg_mortar_quad[f], */
    /*    &sj_on_f_m_mortar_quad[stride], */
    /*    deg_mortar_quad[f], */
    /*    (P4EST_DIM)-1, */
    /*    &MJe1[stride] */
    /*   ); */

    /* double Je1MJe1 = d4est_linalg_vec_dot(&Je1[stride], &MJe1[stride], nodes_mortar_quad[f]); */

    double Je1MJe1 = d4est_operators_quadrature(
                                             d4est_ops,
                                             &Je1[stride],
                                             &Je1[stride],
                                             &sj_on_f_m_mortar_quad[stride],
                                             deg_mortar_quad[f],
                                             geom->geom_quad_type,
                                             (P4EST_DIM)-1);


    /* d4est_operators_apply_mij(d4est_ops, &Je1_test_mortar[stride], (P4EST_DIM)-1, deg_mortar_quad[f], &MJe1_test_mortar[stride]); */
    /* d4est_linalg_vec_scale(sj_on_f_m_mortar_quad[0], &MJe1_test_mortar[stride], nodes_mortar_quad[f]); */
    /* double Je1MJe1_test = d4est_linalg_vec_dot(&Je1_test_mortar[stride], &MJe1_test_mortar[stride], nodes_mortar_quad[f]); */


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
      stride += nodes_mortar_quad[f];
  }


  /* for (int f = 0; f < faces_m; f++){ */
    /* printf("eta2 on face %d = %.25f\n", f, e_m[f]->local_estimator); */
  /* } */
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder);
  
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_quad);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_quad);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(Je2[i]);
    P4EST_FREE(dudr_p_on_f_p_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_porder[i]);
    P4EST_FREE(dudr_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad_porder[i]);
    P4EST_FREE(dudx_p_on_f_p_mortar_quad[i]);
    P4EST_FREE(dudr_m_on_f_m[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar[i]);
    P4EST_FREE(dudr_m_on_f_m_mortar_quad[i]);
    P4EST_FREE(dudx_m_on_f_m_mortar_quad[i]);
    P4EST_FREE(n_on_f_m_mortar_quad[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_quad);
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
