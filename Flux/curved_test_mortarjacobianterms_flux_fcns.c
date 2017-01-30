#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../ElementData/curved_element_data.h"
#include "../LinearAlgebra/linalg.h"
#include "../Flux/curved_test_mortarjacobianterms_flux_fcns.h"


#define TEST_ERROR_EPS .00000001

/* #define DEALIASING */

static void
curved_test_mortarjacobianterms_dirichlet
(
 curved_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 void* params
)
{
  test_mortarjacobianterms_data_t* data = params;
  int face_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* sjvol_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* nvol_on_f_m_Gauss [(P4EST_DIM)];
  double* xyz_on_f_m_Gauss [(P4EST_DIM)];
  double* drst_dxyz_on_f_m_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(n_on_f_m_Gauss,face_nodes_m_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_on_f_m_Gauss,face_nodes_m_Gauss);
  D4EST_ALLOC_DIM_VEC(xyz_on_f_m_Gauss,face_nodes_m_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_on_f_m_Gauss,face_nodes_m_Gauss);
  
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

  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     &e_m,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     drst_dxyz_on_f_m_Gauss,
     sjvol_on_f_m_Gauss,
     nvol_on_f_m_Gauss,
     geom->p4est_geom,
     dgmath_jit_dbase
    );

  DEBUG_PRINT_2ARR_DBL(sj_on_f_m_Gauss,
                       sjvol_on_f_m_Gauss,
                       face_nodes_m_Gauss);

  double maxerror1 = util_max_error(sj_on_f_m_Gauss, sjvol_on_f_m_Gauss, face_nodes_m_Gauss);
  if (maxerror1 > TEST_ERROR_EPS){
    printf("Holy shit batman, LOTS OF ERROR HERE\n");
  }
  
  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(sjvol_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(n_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(nvol_on_f_m_Gauss);
  D4EST_FREE_DIM_VEC(xyz_on_f_m_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_on_f_m_Gauss);
}

static void
curved_test_mortarjacobianterms_interface
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
  test_mortarjacobianterms_data_t*  data = params;

  int stride;
  int deg_p_Lobatto [(P4EST_HALF)];
  int face_nodes_p_Lobatto [(P4EST_HALF)];
  /* int deg_p_Gauss [(P4EST_HALF)]; */
  int face_nodes_p_Gauss [(P4EST_HALF)];

  int deg_m_Lobatto [(P4EST_HALF)];
  int face_nodes_m_Lobatto [(P4EST_HALF)];
  /* int deg_m_Gauss [(P4EST_HALF)]; */
  int face_nodes_m_Gauss [(P4EST_HALF)];
  
  int nodes_mortar_Gauss [(P4EST_HALF)];
  int nodes_mortar_Lobatto [(P4EST_HALF)];
  int deg_mortar_Gauss [(P4EST_HALF)];
  int deg_mortar_Lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];


  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
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
      
    }
  
  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_m_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_p_on_f_p_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_p_on_f_p_mortar_Gauss_reoriented = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* nvol_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_Gauss [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* tmpxyz [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(n_on_f_m_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_m_on_f_m_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss,total_nodes_mortar_Gauss);
  
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

  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e_m,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     drst_dxyz_m_on_mortar_Gauss,
     sjvol_m_on_f_m_mortar_Gauss,
     nvol_m_on_f_m_mortar_Gauss,
     geom->p4est_geom,
     dgmath_jit_dbase
    );



  int deg_mortar_Gauss_p [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF))
      inew = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    deg_mortar_Gauss_p[inew] = deg_mortar_Gauss[i];
  }


  
  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e_p,
     faces_p,
     faces_mortar,
     &deg_mortar_Gauss_p[0], /* CAREFUL THIS NEEDS REORIENTING */
     f_p,
     drst_dxyz_p_on_mortar_Gauss,
     sjvol_p_on_f_p_mortar_Gauss,
     nvol_p_on_f_p_mortar_Gauss,
     geom->p4est_geom,
     dgmath_jit_dbase
    );


  /* int oriented_stride = 0; */
  /* int oriented_strides [P4EST_HALF]; */
  /* for (int i = 0; i < faces_p; i++){ */
  /*   for (int j = 0; j < faces_p; j++){ */
  /*     if (e_p[j]->tree_quadid == e_p_oriented[i]->tree_quadid){ */
  /*       oriented_strides[j] = oriented_stride; */
  /*       break; */
  /*     } */
  /*   } */
  /*   oriented_stride += nodes_mortar_Gauss[i]; */
  /* } */

  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
        face_p = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss_p[b]);
    }

    printf("face, face_p, deg_mortar_Gauss[face], deg_mortar_Gauss_p[face_p] = %d,%d,%d,%d\n", face, face_p, deg_mortar_Gauss[face], deg_mortar_Gauss_p[face_p]);
    
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       &sjvol_p_on_f_p_mortar_Gauss[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_Gauss[face],
       orientation,
       f_m,
       f_p,
       &sjvol_p_on_f_p_mortar_Gauss_reoriented[face_mortar_stride]
      );
  
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]);
  }
  
  
  /* printf("faces_m, faces_p = %d,%d\n", faces_m, faces_p); */

  
  /* stride = 0; */
  /* for (int f = 0; f < faces_mortar; f++){ */
  /*   dgmath_reorient_face_data */
  /*     ( */
  /*      dgmath_jit_dbase, */
  /*      &sjvol_p_on_f_p_mortar_Gauss[stride], */
  /*      ((P4EST_DIM) - 1), */
  /*      deg_mortar_Gauss[f], */
  /*      orientation, */
  /*      f_m, */
  /*      f_p, */
  /*      &sjvol_p_on_f_p_mortar_Gauss_reoriented[stride] */
  /*     );       */
  /*   stride += nodes_mortar_Gauss[f]; */
  /* }    */

  

  DEBUG_PRINT_4ARR_DBL(sj_on_f_m_mortar_Gauss,
                       sjvol_m_on_f_m_mortar_Gauss,
                       sjvol_p_on_f_p_mortar_Gauss,
                       sjvol_p_on_f_p_mortar_Gauss_reoriented,
                       total_nodes_mortar_Gauss);


  double maxerror1 = util_max_error(sj_on_f_m_mortar_Gauss, sjvol_m_on_f_m_mortar_Gauss, total_nodes_mortar_Gauss);
  double maxerror2 = util_max_error(sj_on_f_m_mortar_Gauss, sjvol_p_on_f_p_mortar_Gauss_reoriented, total_nodes_mortar_Gauss);
  if (maxerror1 > TEST_ERROR_EPS || maxerror2 > TEST_ERROR_EPS){
    printf("Holy shit batman, LOTS OF ERROR HERE\n");
  }

  
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  P4EST_FREE(sjvol_m_on_f_m_mortar_Gauss);
  P4EST_FREE(sjvol_p_on_f_p_mortar_Gauss);
  P4EST_FREE(sjvol_p_on_f_p_mortar_Gauss_reoriented);
  D4EST_FREE_DIM_VEC(tmpxyz);
  D4EST_FREE_DIM_VEC(n_on_f_m_mortar_Gauss);
  D4EST_FREE_DIM_VEC(nvol_m_on_f_m_mortar_Gauss);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss);
}

curved_flux_fcn_ptrs_t
curved_test_mortarjacobianterms_fetch_fcns
(
 test_mortarjacobianterms_data_t* data
)
{
  
  curved_flux_fcn_ptrs_t curved_test_mortarjacobianterms_fcns;
  curved_test_mortarjacobianterms_fcns.flux_interface_fcn = curved_test_mortarjacobianterms_interface;
  curved_test_mortarjacobianterms_fcns.flux_boundary_fcn = curved_test_mortarjacobianterms_dirichlet;
  curved_test_mortarjacobianterms_fcns.params = (void*)data;

  return curved_test_mortarjacobianterms_fcns;
}
