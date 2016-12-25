#include "../Flux/curved_Gauss_sipg_flux_scalar_fcns.h"
#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../LinearAlgebra/linalg.h"

/* TODO: Add in dealiasing to scalar_dirichlet */
/* #define DEALIASING */

static void
curved_Gauss_sipg_flux_scalar_dirichlet
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
  /* int vol_nodes_m = dgmath_get_nodes( (P4EST_DIM), e_m->deg); */
  int face_nodes_m_Lobatto = dgmath_get_nodes( (P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = dgmath_get_nodes( (P4EST_DIM) - 1, e_m->deg_integ);


  double* xyz_on_f_m_Gauss [(P4EST_DIM)];
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* sj_n_on_f_m_Gauss [(P4EST_DIM)];
  
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* ustar_min_u_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);

  for (int i = 0; i < (P4EST_DIM); i++) {
    /* xyz_on_f_m[i] = P4EST_ALLOC(double, face_nodes_m_Lobatto);     */
    /* xyz_on_f_m[i] = P4EST_ALLOC(double, face_nodes_m_Gauss);     */
    xyz_on_f_m_Gauss[i] = P4EST_ALLOC(double, face_nodes_m_Gauss);    
    n_on_f_m_Gauss[i] = P4EST_ALLOC(double, face_nodes_m_Gauss);    
    sj_n_on_f_m_Gauss[i] = P4EST_ALLOC(double, face_nodes_m_Gauss);    
  }
  

  curved_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     &e_m,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     geom->dxdr_method,
     1,
     n_on_f_m_Gauss,
     sj_on_f_m_Gauss,
     geom,
     dgmath_jit_dbase,
     xyz_on_f_m_Gauss
    );
  
  dgmath_apply_slicer(dgmath_jit_dbase, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  dgmath_interp_GLL_to_GL
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     e_m->deg,
     e_m->deg_integ,
     u_m_on_f_m_Gauss,
     (P4EST_DIM)-1
    );
  
  for (int i = 0; i < face_nodes_m_Gauss; i++){
    ustar_min_u_on_f_m_Gauss[i] = u_at_bndry(xyz_on_f_m_Gauss[0][i],
                                             xyz_on_f_m_Gauss[1][i]
#if (P4EST_DIM)==3
                                             ,
                                             xyz_on_f_m_Gauss[2][i]
#endif                          
                                            ) - u_m_on_f_m_Gauss[i]; 
  }


  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int i = 0; i < face_nodes_m_Gauss; i++){
      sj_n_on_f_m_Gauss[d][i] = sj_on_f_m_Gauss[i]
                                *n_on_f_m_Gauss[d][i];
    }
    
    dgmath_apply_curvedGaussMass_onGaussNodeVec
      (
       dgmath_jit_dbase,
       ustar_min_u_on_f_m_Gauss,
       e_m->deg,
       sj_n_on_f_m_Gauss[d],
       e_m->deg_integ,
       (P4EST_DIM)-1,
       &e_m->M_ustar_min_u_n[d][f_m*face_nodes_m_Lobatto]
      );
  }
  
  for (int dir = 0; dir < (P4EST_DIM); dir++){
    /* P4EST_FREE(xyz_on_f_m[dir]); */
    P4EST_FREE(xyz_on_f_m_Gauss[dir]);
    P4EST_FREE(n_on_f_m_Gauss[dir]);
    P4EST_FREE(sj_n_on_f_m_Gauss[dir]);
  }

  P4EST_FREE(ustar_min_u_on_f_m_Gauss);
  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(u_m_on_f_m_Gauss);
  P4EST_FREE(u_m_on_f_m);
}

static void
curved_Gauss_sipg_flux_scalar_interface
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
    
  /* int deg_m_Gauss [(P4EST_HALF)]; */
  int deg_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Gauss [(P4EST_HALF)];
  
  int deg_mortar_Gauss [(P4EST_HALF)];
  int deg_mortar_Lobatto [(P4EST_HALF)];
  int nodes_mortar_Gauss [(P4EST_HALF)];
  int nodes_mortar_Lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  
  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(
                                                         e_p,
                                                         (P4EST_DIM)-1,
                                                         f_m,
                                                         f_p,
                                                         orientation,
                                                         faces_p,
                                                         e_p_oriented
                                                        );
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_Lobatto = 0;
  int total_side_nodes_m_Gauss = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_Lobatto[i] = e_m[i]->deg;
    /* deg_m_Gauss[i] = e_m[i]->deg_integ; */

    face_nodes_m_Lobatto[i] = dgmath_get_nodes
                              (
                               (P4EST_DIM) - 1,
                               e_m[i]->deg
                              );
    face_nodes_m_Gauss[i] = dgmath_get_nodes
                            (
                             (P4EST_DIM) - 1,
                             e_m[i]->deg_integ
                            );
    
    total_side_nodes_m_Lobatto += face_nodes_m_Lobatto[i];
    total_side_nodes_m_Gauss += face_nodes_m_Gauss[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_Lobatto = 0;
  /* int total_side_nodes_p_Gauss = 0; */
  for (int i = 0; i < faces_p; i++){
    deg_p_Lobatto[i] = e_p_oriented[i]->deg;
    face_nodes_p_Lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    total_side_nodes_p_Lobatto += face_nodes_p_Lobatto[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_Gauss = 0;
  int total_nodes_mortar_Lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_Gauss[i+j] = util_max_int( e_m[i]->deg_integ, e_p_oriented[j]->deg_integ);
      deg_mortar_Lobatto[i+j] = util_max_int( e_m[i]->deg, e_p_oriented[j]->deg);
      nodes_mortar_Gauss[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Gauss[i+j] );     
      nodes_mortar_Lobatto[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_Lobatto[i+j] );     
      total_nodes_mortar_Gauss += nodes_mortar_Gauss[i+j];
      total_nodes_mortar_Lobatto += nodes_mortar_Lobatto[i+j];
    }

  /* scalar and scalar fields on each of the (-) and (+) elements */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);

  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_m_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* u_p_on_f_p_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);

  /* printf("total_nodes_mortar_Gauss = %d\n",total_nodes_mortar_Gauss); */
  double* ustar_min_u_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* M_sj_n_ustar_min_u_mortar = P4EST_ALLOC(double, total_nodes_mortar_Lobatto);
  double* proj_M_sj_n_ustar_min_u_mortar = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);

  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* sj_n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  
  for (int i = 0; i < (P4EST_DIM); i++){
    sj_n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
    n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
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
    
     stride += face_nodes_m_Lobatto[i];
  }

  stride = 0;
  double* tmp = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);
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
    stride += face_nodes_p_Lobatto[i];
  }
  P4EST_FREE(tmp);
  
  /* project (-)-side u trace scalar onto mortar space */ 
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

  /* project (+)-side u trace scalar onto mortar space */
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
 
  double* tmpxyz [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_Gauss);
  curved_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     e_m,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     geom->dxdr_method,
     1,
     n_on_f_m_mortar_Gauss,
     sj_on_f_m_mortar_Gauss,
     geom,
     dgmath_jit_dbase,
     tmpxyz
    );
  D4EST_FREE_DIM_VEC(tmpxyz);  
  /* calculate symmetric interior penalty flux */
  int k;
  int f;
  stride = 0;
  for (f = 0; f < faces_mortar; f++){
    for (k = 0; k < nodes_mortar_Gauss[f]; k++){
      int ks = k + stride;
      ustar_min_u_mortar_Gauss[ks] = .5*(u_p_on_f_p_mortar_Gauss[ks] + u_m_on_f_m_mortar_Gauss[ks]);
      ustar_min_u_mortar_Gauss[ks] -= u_m_on_f_m_mortar_Gauss[ks];
    }
    stride += nodes_mortar_Gauss[f];
  }
    
  for (int d = 0; d < (P4EST_DIM); d++){

    stride = 0;
    int stride_Lobatto = 0;
    for (int f = 0; f < faces_mortar; f++){
      for (int k = 0; k < nodes_mortar_Gauss[f]; k++){
        int ks = k + stride;
        sj_n_on_f_m_mortar_Gauss[d][ks] = sj_on_f_m_mortar_Gauss[ks]*n_on_f_m_mortar_Gauss[d][ks];
      }
      dgmath_apply_curvedGaussMass_onGaussNodeVec(
                                                  dgmath_jit_dbase,
                                                  &ustar_min_u_mortar_Gauss[stride],
                                                  deg_mortar_Lobatto[f],
                                                  &sj_n_on_f_m_mortar_Gauss[d][stride],
                                                  deg_mortar_Gauss[f],
                                                  (P4EST_DIM)-1,
                                                  &M_sj_n_ustar_min_u_mortar[stride_Lobatto]
                                                 );

      /* dgmath_apply_curvedGaussMass_onGaussNodeVec( */
      /*                                             dgmath_jit_dbase, */
      /*                                             &ustar_min_u_mortar_Gauss[stride], */
      /*                                             deg_m_Lobatto[f], */
      /*                                             &sj_n_on_f_m_mortar_Gauss[d][stride], */
      /*                                             deg_mortar_Gauss[f], */
      /*                                             (P4EST_DIM)-1, */
      /*                                             &M_sj_n_ustar_min_u_mortar[stride_Lobatto] */
      /*                                            ); */
      
      stride_Lobatto += nodes_mortar_Lobatto[f];
      /* stride_Lobatto += dgmath_get_nodes((P4EST_DIM)-1, e_m[f]->deg); */
      stride += nodes_mortar_Gauss[f];
    }
      
    /* DEBUG_PRINT_ARR_INT(deg_mortar_Gauss, faces_mortar); */
    /* DEBUG_PRINT_ARR_INT(deg_m_Gauss, faces_m); */


    /* if (faces_mortar == P4EST_HALF) */
      /* mpi_abort("Only testing p-nonconforming atm"); */
    

    /* project mortar data back onto the (-) side */
    dgmath_project_mass_mortar_onto_side
      (
       dgmath_jit_dbase,
       M_sj_n_ustar_min_u_mortar,
       faces_mortar,
       deg_mortar_Lobatto,
       proj_M_sj_n_ustar_min_u_mortar,
       faces_m,
       deg_m_Lobatto
      );
  
    /* copy result back to element */
    stride = 0;
    for (int i = 0; i < faces_m; i++){
      if(e_m_is_ghost[i] == 0)
        linalg_copy_1st_to_2nd
          (
           &proj_M_sj_n_ustar_min_u_mortar[stride],
           &(e_m[i]->M_ustar_min_u_n[d][f_m*face_nodes_m_Lobatto[i]]),
           face_nodes_m_Lobatto[i]
          );
      stride += face_nodes_m_Lobatto[i];
    }


    /* for (int i = 0; i < faces_m; i++){ */
    /*   if(e_m_is_ghost[i] == 0) */
    /*     linalg_copy_1st_to_2nd */
    /*       ( */
    /*        &M_sj_n_ustar_min_u_mortar[stride], */
    /*        &(e_m[i]->M_ustar_min_u_n[d][f_m*face_nodes_m_Lobatto[i]]), */
    /*        face_nodes_m_Lobatto[i] */
    /*       ); */
    /*   stride += face_nodes_m_Lobatto[i]; */
    /* }     */
    
  }


/* #if D4EST_DEBUG */
/*   if( */
/*       (e_m[0]->id == 2 && e_p[0]->id == 3) || */
/*       (e_m[0]->id == 3 && e_p[0]->id == 2) */
/*      ) */
/*     { */
/*       printf("e_m = %d, e_p = %d\n", e_m[0]->id, e_p[0]->id); */
/*       DEBUG_PRINT_2ARR_DBL(n_on_f_m_mortar_Gauss[0], n_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
/*       DEBUG_PRINT_ARR_DBL(ustar_min_u_mortar_Gauss, total_nodes_mortar_Gauss); */
/*       double* em0 = &(e_m[0]->M_ustar_min_u_n[0][f_m*face_nodes_m_Lobatto[0]]); */
/*       double* ep0 = &(e_p[0]->M_ustar_min_u_n[0][f_m*face_nodes_p_Lobatto[0]]); */
      
/*       DEBUG_PRINT_ARR_DBL(em0, face_nodes_m_Lobatto[0]); */
/*       DEBUG_PRINT_ARR_DBL(ep0, face_nodes_p_Lobatto[0]); */
/*     } */
/* #endif */


  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(n_on_f_m_mortar_Gauss[d]);
    P4EST_FREE(sj_n_on_f_m_mortar_Gauss[d]);
  }
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  P4EST_FREE(proj_M_sj_n_ustar_min_u_mortar);  
  P4EST_FREE(M_sj_n_ustar_min_u_mortar);
  P4EST_FREE(ustar_min_u_mortar_Gauss);
  P4EST_FREE(u_p_on_f_p_mortar_Gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_Gauss);
  P4EST_FREE(u_m_on_f_m_mortar);
}

curved_flux_fcn_ptrs_t
curved_Gauss_sipg_flux_scalar_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn
)
{
  curved_flux_fcn_ptrs_t curved_Gauss_sipg_flux_scalar_fcns;
  curved_Gauss_sipg_flux_scalar_fcns.flux_interface_fcn
    = curved_Gauss_sipg_flux_scalar_interface;
  curved_Gauss_sipg_flux_scalar_fcns.flux_boundary_fcn
    = curved_Gauss_sipg_flux_scalar_dirichlet;
  curved_Gauss_sipg_flux_scalar_fcns.bndry_fcn = bndry_fcn;
  
  return curved_Gauss_sipg_flux_scalar_fcns;
}
