#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../ElementData/d4est_element_data.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Flux/curved_gauss_central_flux_vector_fcns.h"

/* #define DEALIASING */

static void
curved_gauss_central_flux_vector_dirichlet
(
 d4est_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* geom,
 void* params
)
{

  /* int deg_face_gauss_quad = e_m->deg_quad; */
  
  central_flux_params_t* central_params = (central_flux_params_t*) params;
  double penalty = central_params->central_flux_penalty_prefactor;
  
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m_lobatto = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_gauss = d4est_lgl_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* u_m_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* q_m_on_f_m [(P4EST_DIM)];
  double* q_m_on_f_m_gauss [(P4EST_DIM)];
  /* double* xyz_on_f_m [(P4EST_DIM)]; */
  /* double* xyz_on_f_m_gauss [(P4EST_DIM)]; */
  double* xyz_on_f_m [(P4EST_DIM)];
  for (int i = 0; i < (P4EST_DIM); i++) {
    q_m_on_f_m[i] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    q_m_on_f_m_gauss[i] = P4EST_ALLOC(double, face_nodes_m_gauss);
    /* xyz_on_f_m[i] = P4EST_ALLOC(double, face_nodes_m_lobatto); */
    xyz_on_f_m[i] = P4EST_ALLOC(double, face_nodes_m_lobatto);
    /* xyz_on_f_m_gauss[i] = P4EST_ALLOC(double, face_nodes_m_gauss); */
  }

  double* M_qstar_min_q_n = P4EST_ALLOC_ZERO(double, face_nodes_m_lobatto);
  double* u_on_f_m_min_u_at_bndry_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* u_on_f_m_min_u_at_bndry_lobatto = P4EST_ALLOC(double, face_nodes_m_lobatto);
  double* sj_on_f_m_gauss = P4EST_ALLOC(double, face_nodes_m_gauss);
  double* qstar_min_q_gauss [(P4EST_DIM)];
  double* n_on_f_m_gauss [(P4EST_DIM)];
  double* sj_n_on_f_m_gauss [(P4EST_DIM)];

  for (int d = 0; d < (P4EST_DIM); d++){
    n_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    sj_n_on_f_m_gauss[d] = P4EST_ALLOC(double, face_nodes_m_gauss);
    qstar_min_q_gauss[d] = P4EST_ALLOC_ZERO(double, face_nodes_m_gauss);
  }

  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_slicer(d4est_ops,
                        e_m->xyz[d],
                        (P4EST_DIM),
                        f_m,
                        e_m->deg,
                        xyz_on_f_m[d]);
    
    /* d4est_operators_interp_lobatto_to_GL */
    /*   ( */
    /*    d4est_ops, */
    /*    xyz_on_f_m[d], */
    /*    /\* e_m->deg, *\/ */
    /*    e_m->deg_quad, */
    /*    e_m->deg_quad, */
    /*    xyz_on_f_m_gauss[d], */
    /*    (P4EST_DIM)-1 */
    /*   ); */
    
    d4est_operators_apply_slicer
      (
       d4est_ops,
       e_m->q_elem[d],
       (P4EST_DIM),
       f_m,
       e_m->deg,
       q_m_on_f_m[d]
      );

   d4est_operators_interp_lobatto_to_GL
     (
      d4est_ops,
      q_m_on_f_m[d],
      e_m->deg,
      e_m->deg_quad,
      q_m_on_f_m_gauss[d],
      (P4EST_DIM)-1
     );
    
  }

  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  d4est_operators_interp_lobatto_to_GL
    (
     d4est_ops,
     u_m_on_f_m,
     e_m->deg,
     e_m->deg_quad,
     u_m_on_f_m_gauss,
     (P4EST_DIM)-1
    );

  d4est_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     &e_m,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     geom->dxdr_method,
     1,
     n_on_f_m_gauss,
     sj_on_f_m_gauss,
     geom,
     d4est_ops,
     (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
         , NULL
#endif
         }
    );
  /* get boundary values on this face */
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

    /* printf("u_at_bndry = %.25f\n", u_at_bndry(xyz_on_f_m[0][i], xyz_on_f_m[1][i] */
/* #if (P4EST_DIM)==3 */
                                              /* , xyz_on_f_m[2][i] */
/* #endif */
                                             /* )); */
  }
  
  d4est_operators_interp_lobatto_to_GL
    (
     d4est_ops,
     u_on_f_m_min_u_at_bndry_lobatto,
     e_m->deg,
     e_m->deg_quad,
     u_on_f_m_min_u_at_bndry_gauss,
     (P4EST_DIM)-1
    );
  
  
  
  for(int i = 0; i < face_nodes_m_gauss; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      qstar_min_q_gauss[d][i] = q_m_on_f_m_gauss[d][i]
                       - penalty*n_on_f_m_gauss[d][i]*(u_on_f_m_min_u_at_bndry_gauss[i])
                       - q_m_on_f_m_gauss[d][i];
      sj_n_on_f_m_gauss[d][i] = sj_on_f_m_gauss[i]*n_on_f_m_gauss[d][i];
    }
  }

  for (int i = 0; i < face_nodes_m_lobatto; i++){
    e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto + i] = 0.;
  }
  
  for (int d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_curvedgaussMass_ongaussNodeVec(d4est_ops,
                                                qstar_min_q_gauss[d],
                                                e_m->deg,
                                                sj_n_on_f_m_gauss[d],
                                                e_m->deg_quad,
                                                (P4EST_DIM)-1,
                                                M_qstar_min_q_n);
    
    for (int i = 0; i < face_nodes_m_lobatto; i++){
      e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto + i] += M_qstar_min_q_n[i];
      /* printf("e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto + i] = %.25f\n", e_m->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto + i]); */
    }
  }

  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_gauss);
  for (int i = 0; i < (P4EST_DIM); i++) {
    P4EST_FREE(q_m_on_f_m[i]);
    P4EST_FREE(q_m_on_f_m_gauss[i]);
    P4EST_FREE(xyz_on_f_m[i]);
    /* P4EST_FREE(xyz_on_f_m_gauss[i]); */
    P4EST_FREE(qstar_min_q_gauss[i]);
  }

  P4EST_FREE(sj_on_f_m_gauss);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_gauss);
  P4EST_FREE(u_on_f_m_min_u_at_bndry_lobatto);
  P4EST_FREE(M_qstar_min_q_n);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m_gauss[d]);
    P4EST_FREE(sj_n_on_f_m_gauss[d]);
  }
}

static void
curved_gauss_central_flux_vector_interface
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
  central_flux_params_t* central_params = (central_flux_params_t*) params;
  double curved_gauss_central_flux_penalty_prefactor = central_params->central_flux_penalty_prefactor;
  
  int stride;

  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  /* int deg_p_gauss [(P4EST_HALF)]; */
  int face_nodes_p_gauss [(P4EST_HALF)];

  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  /* int deg_m_gauss [(P4EST_HALF)]; */
  int face_nodes_m_gauss [(P4EST_HALF)];
  
  int nodes_mortar_gauss [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_gauss [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    /* deg_m_gauss[i] = e_m[i]->deg_quad; */
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_gauss[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_gauss += face_nodes_m_gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_gauss = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    /* deg_p_gauss[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_gauss[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_gauss += face_nodes_p_gauss[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_gauss = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_gauss[i+j] = util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                            e_p_oriented[j]->deg );      
      nodes_mortar_gauss[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_gauss[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_gauss += nodes_mortar_gauss[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
    }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  double* q_p_on_f_p [(P4EST_DIM)];
  double* q_m_on_f_m [(P4EST_DIM)];

  for (int i = 0; i < (P4EST_DIM); i++) {
    q_m_on_f_m[i] = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
    q_p_on_f_p[i] = P4EST_ALLOC(double, total_side_nodes_p_lobatto);
  }
  
  /* projections of f_m/f_p on to mortar space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_m_on_f_m_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* u_p_on_f_p_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  double* q_p_on_f_p_mortar [(P4EST_DIM)];
  double* q_p_on_f_p_mortar_gauss [(P4EST_DIM)];
  double* q_m_on_f_m_mortar [(P4EST_DIM)];
  double* q_m_on_f_m_mortar_gauss [(P4EST_DIM)];

  for (int i = 0; i < (P4EST_DIM); i++){
    q_p_on_f_p_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    q_p_on_f_p_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    q_m_on_f_m_mortar[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    q_m_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  }
  
  double* sj_on_f_m_mortar_gauss = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  /* double* sj_on_f_m_mortar_gauss_m = P4EST_ALLOC(double, total_nodes_mortar_gauss); */
  double* n_on_f_m_mortar_gauss [(P4EST_DIM)];
  /* double* n_on_f_m_mortar_gauss_m [(P4EST_DIM)]; */
  double* sj_n_on_f_m_mortar_gauss [(P4EST_DIM)];

  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_gauss);
  double* qstar_min_q_mortar_gauss [(P4EST_DIM)];
  double* proj_M_qstar_min_q_dot_n_mortar = P4EST_ALLOC_ZERO(double, total_side_nodes_m_lobatto);
  double* M_qstar_min_q_dot_n_mortar = P4EST_ALLOC_ZERO(double, total_nodes_mortar_gauss);
  double* M_qstar_min_q_n_mortar = P4EST_ALLOC_ZERO(double, total_nodes_mortar_gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++) {
    n_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    /* n_on_f_m_mortar_gauss_m[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss); */
    sj_n_on_f_m_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
    qstar_min_q_mortar_gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_gauss);
  }  
  
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


  /* project (-)-side u trace vector onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m_lobatto,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar_gauss
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
     deg_mortar_gauss
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    d4est_operators_interp_lobatto_to_GL(d4est_ops, &u_m_on_f_m_mortar[stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &u_m_on_f_m_mortar_gauss[stride], (P4EST_DIM)-1);
    d4est_operators_interp_lobatto_to_GL(d4est_ops, &u_p_on_f_p_mortar[stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &u_p_on_f_p_mortar_gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_gauss[f];
  }

  
  /* For each component of the vector */
  for (int d = 0; d < (P4EST_DIM); d++){

    /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (int i = 0; i < faces_m; i++){    
      d4est_operators_apply_slicer
        (
         d4est_ops,
         e_m[i]->q_elem[d],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &q_m_on_f_m[d][stride]
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
         &e_p_oriented[i]->q_elem[d][0],
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
         &q_p_on_f_p[d][stride]
        );

      stride += face_nodes_p_lobatto[i];
    }

    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       q_p_on_f_p[d],
       faces_p,
       deg_p_lobatto,
       q_p_on_f_p_mortar[d],
       faces_mortar,
       deg_mortar_gauss
      );

    /* project q from the (-) side onto the mortar space */
    d4est_mortars_project_side_onto_mortar_space
      (
       d4est_ops,
       q_m_on_f_m[d],
       faces_m,
       deg_m_lobatto,
       q_m_on_f_m_mortar[d],
       faces_mortar,
       deg_mortar_gauss
      );

    stride = 0;
    for (int f = 0; f < faces_mortar; f++){
      d4est_operators_interp_lobatto_to_GL(d4est_ops, &q_m_on_f_m_mortar[d][stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &q_m_on_f_m_mortar_gauss[d][stride], (P4EST_DIM)-1);
      d4est_operators_interp_lobatto_to_GL(d4est_ops, &q_p_on_f_p_mortar[d][stride], deg_mortar_gauss[f], deg_mortar_gauss[f], &q_p_on_f_p_mortar_gauss[d][stride], (P4EST_DIM)-1);
      stride += nodes_mortar_gauss[f];
    }
  }

  double* tmpxyz [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_gauss);
  d4est_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     e_m,
     faces_m,
     faces_mortar,
     &deg_mortar_gauss[0],
     f_m,
     geom->dxdr_method,
     1,
     n_on_f_m_mortar_gauss,
     sj_on_f_m_mortar_gauss,
     geom,
     d4est_ops,
     tmpxyz
    );
  D4EST_FREE_DIM_VEC(tmpxyz);

  
    /* calculate symmetric interior penalty flux */
    stride = 0;
    int stride_lobatto = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_gauss[f]; k++){
        int ks = k + stride;
        double sigma = curved_gauss_central_flux_penalty_prefactor;
        
        for (int d = 0; d < (P4EST_DIM); d++){
          qstar_min_q_mortar_gauss[d][ks] = .5*(q_p_on_f_p_mortar_gauss[d][ks] + q_m_on_f_m_mortar_gauss[d][ks]);
          qstar_min_q_mortar_gauss[d][ks] -= -sigma*n_on_f_m_mortar_gauss[d][ks]*u_p_on_f_p_mortar_gauss[ks];
          qstar_min_q_mortar_gauss[d][ks] -= sigma*n_on_f_m_mortar_gauss[d][ks]*u_m_on_f_m_mortar_gauss[ks];
          qstar_min_q_mortar_gauss[d][ks] -= q_m_on_f_m_mortar_gauss[d][ks];
          sj_n_on_f_m_mortar_gauss[d][ks] = sj_on_f_m_mortar_gauss[ks]*n_on_f_m_mortar_gauss[d][ks];
        }        
      }
  
      for (int d = 0; d < (P4EST_DIM); d++){
        d4est_operators_apply_curvedgaussMass_ongaussNodeVec(
                                     d4est_ops,
                                     &qstar_min_q_mortar_gauss[d][stride],
                                     deg_mortar_lobatto[f],
                                     &sj_n_on_f_m_mortar_gauss[d][stride],
                                     deg_mortar_gauss[f],
                                     (P4EST_DIM)-1,
                                     &M_qstar_min_q_n_mortar[stride_lobatto]
                                    );
    
        for (int i = 0; i < nodes_mortar_lobatto[f]; i++){
          M_qstar_min_q_dot_n_mortar[stride_lobatto + i] += M_qstar_min_q_n_mortar[stride_lobatto + i];
        }
      }
      stride += nodes_mortar_gauss[f];
      stride_lobatto += nodes_mortar_lobatto[f];
      /* stride_lobatto += d4est_lgl_get_nodes((P4EST_DIM)-1, e_m[f]->deg); */
    }

    /* if (faces_mortar == P4EST_HALF) */
      /* mpi_abort("Only testing p-nonconforming atm"); */
    
    d4est_mortars_project_mass_mortar_onto_side
      (
       d4est_ops,
       M_qstar_min_q_dot_n_mortar,
       faces_mortar,
       deg_mortar_lobatto,
       proj_M_qstar_min_q_dot_n_mortar,
       faces_m,
       deg_m_lobatto
      );
    
  /* copy result back to element */
  stride = 0;
  for (int i = 0; i < faces_m; i++){
    if(e_m_is_ghost[i] == 0)
      d4est_linalg_copy_1st_to_2nd
        (
         &proj_M_qstar_min_q_dot_n_mortar[stride],
         &(e_m[i]->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto[i]]),
         face_nodes_m_lobatto[i]
        );
    stride += face_nodes_m_lobatto[i];
  }

  /* stride = 0; */
  /* for (int i = 0; i < faces_m; i++){ */
  /*   if(e_m_is_ghost[i] == 0) */
  /*     d4est_linalg_copy_1st_to_2nd */
  /*       ( */
  /*        &M_qstar_min_q_dot_n_mortar[stride], */
  /*        &(e_m[i]->M_qstar_min_q_dot_n[f_m*face_nodes_m_lobatto[i]]), */
  /*        face_nodes_m_lobatto[i] */
  /*       ); */
  /*   stride += face_nodes_m_lobatto[i]; */
  /* } */
    
 
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(q_p_on_f_p[i]);
    P4EST_FREE(q_m_on_f_m[i]);
  }

  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_gauss);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(q_p_on_f_p_mortar[i]);
    P4EST_FREE(q_p_on_f_p_mortar_gauss[i]);
    P4EST_FREE(q_m_on_f_m_mortar[i]);
    P4EST_FREE(q_m_on_f_m_mortar_gauss[i]);
  }

  P4EST_FREE(sj_on_f_m_mortar_gauss);
  /* P4EST_FREE(sj_on_f_m_mortar_gauss_m); */
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(n_on_f_m_mortar_gauss[i]);
    /* P4EST_FREE(n_on_f_m_mortar_gauss_m[i]); */
    P4EST_FREE(sj_n_on_f_m_mortar_gauss[i]);
    P4EST_FREE(qstar_min_q_mortar_gauss[i]);
  }    
  P4EST_FREE(M_qstar_min_q_n_mortar);
  P4EST_FREE(M_qstar_min_q_dot_n_mortar);
  P4EST_FREE(proj_M_qstar_min_q_dot_n_mortar);
  P4EST_FREE(tmp);  
  
}

curved_flux_fcn_ptrs_t
curved_gauss_central_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 central_flux_params_t* curved_gauss_central_params
)
{  
  curved_flux_fcn_ptrs_t curved_gauss_central_flux_vector_fcns;
  curved_gauss_central_flux_vector_fcns.flux_interface_fcn = curved_gauss_central_flux_vector_interface;
  curved_gauss_central_flux_vector_fcns.flux_boundary_fcn = curved_gauss_central_flux_vector_dirichlet;
  curved_gauss_central_flux_vector_fcns.bndry_fcn = bndry_fcn;
  curved_gauss_central_flux_vector_fcns.params = (void*)curved_gauss_central_params;
  
  return curved_gauss_central_flux_vector_fcns;
}
