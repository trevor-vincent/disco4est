#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../ElementData/curved_element_data.h"
#include "../LinearAlgebra/linalg.h"
#include "../Flux/curved_test_mortarjacobianterms_flux_fcns.h"

void curved_test_mortarjacobianterms_init_vecs
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;
  test_mortarjacobianterms_data_t* test_data = (test_mortarjacobianterms_data_t*) user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = test_data->dgmath_jit_dbase;
  
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg);
  int face_nodes_Lobatto = dgmath_get_nodes(dim-1,deg);
  int volume_nodes_Gauss = dgmath_get_nodes(dim, elem_data->deg_integ);
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_Dij(dgmath_jit_dbase,
                     &(test_data->u[elem_data->nodal_stride]),
                     dim,
                     elem_data->deg,
                     i,
                     &elem_data->dudr_elem[i][0]);  
  }



  linalg_copy_1st_to_2nd
    (
     &(test_data->u[elem_data->nodal_stride]),
     &(elem_data->u_storage)[0],
     volume_nodes_Lobatto
    );
  

}

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
     dgmath_jit_dbase,
     NULL
    );

  /* DEBUG_PRINT_2ARR_DBL(sj_on_f_m_Gauss, */
  /*                      sjvol_on_f_m_Gauss, */
  /*                      face_nodes_m_Gauss); */

  double maxerror = util_max_error(sj_on_f_m_Gauss, sjvol_on_f_m_Gauss, face_nodes_m_Gauss);
  maxerror += util_max_error(n_on_f_m_Gauss[0], nvol_on_f_m_Gauss[0], face_nodes_m_Gauss);
  maxerror += util_max_error(n_on_f_m_Gauss[1], nvol_on_f_m_Gauss[1], face_nodes_m_Gauss);
#if (P4EST_DIM)==3
  maxerror += util_max_error(n_on_f_m_Gauss[2], nvol_on_f_m_Gauss[2], face_nodes_m_Gauss);
#endif
  if (maxerror > data->local_eps){
    printf("Holy shit batman, LOTS OF ERROR HERE\n");
  }

  data->global_err += maxerror;
  
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
  int deg_p_Lobatto_porder [(P4EST_HALF)];
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


  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_Gauss [(P4EST_DIM)];
  double* mortar_flux [(P4EST_DIM)];

  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_Lobatto);
  double* u_p_on_f_p_porder = P4EST_ALLOC(double, total_side_nodes_p_Lobatto);

  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_porder, total_side_nodes_p_Lobatto); 
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_mortar_porder, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_mortar_Gauss_porder, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_Gauss_porder, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_Gauss, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m, total_side_nodes_m_Lobatto);
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_mortar, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(mortar_flux, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_mortar_Gauss, total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(dudx_m_on_f_m_mortar_Gauss, total_nodes_mortar_Gauss);

  stride = 0;
  for (int i = 0; i < faces_m; i++){    
    for (int d = 0; d < (P4EST_DIM); d++){

      /* project (-)-u-derivative on the (-)-side faces and project q onto the (-)-side faces */   
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_m[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &dudr_m_on_f_m[d][stride]
        );
      
    }
    stride += face_nodes_m_Lobatto[i];
  }


  stride = 0;
  for (int i = 0; i < faces_m; i++){    
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_m[i]->u_storage[0],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &u_m_on_f_m[stride]
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

  stride = 0;
  for (int i = 0; i < faces_p; i++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         &e_p[i]->u_storage[0],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &u_p_on_f_p_porder[stride]
        );
    stride += dgmath_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
  }
  

  /* if (e_m[0]->id == 0 && e_p[0]->id == 4){ */
  /*   DEBUG_PRINT_2ARR_DBL(u_p_on_f_p_porder, u_m_on_f_m, total_side_nodes_m_Lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[0], dudr_m_on_f_m[0], total_side_nodes_m_Lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[1], dudr_m_on_f_m[1], total_side_nodes_m_Lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[2], dudr_m_on_f_m[2], total_side_nodes_m_Lobatto); */

  /* } */
  

  for (int d = 0; d < (P4EST_DIM); d++){
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
  }
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_m_on_f_m_mortar[d][stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &dudr_m_on_f_m_mortar_Gauss[d][stride], (P4EST_DIM)-1);
    }
    stride += nodes_mortar_Gauss[f];
  }

  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dudr_p_on_f_p_mortar_porder[d][stride], deg_mortar_Gauss_porder[f], deg_mortar_Gauss_porder[f], &dudr_p_on_f_p_mortar_Gauss_porder[d][stride], (P4EST_DIM)-1);
    }
    stride += nodes_mortar_Gauss_porder[f];
  }

  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_m_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_p_on_f_p_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* sjvol_p_on_f_p_mortar_Gauss_reoriented = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* nvol_m_on_f_m_mortar_Gauss [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_Gauss [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_Gauss_reoriented [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_Gauss_porder [(P4EST_DIM)][(P4EST_DIM)];
  double* tmpxyz [(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(n_on_f_m_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_m_on_f_m_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_Gauss_reoriented,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss,total_nodes_mortar_Gauss);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder,total_nodes_mortar_Gauss);
  
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
     dgmath_jit_dbase,
     NULL
    );





  
  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e_p,
     faces_p,
     faces_mortar,
     &deg_mortar_Gauss_porder[0], 
     f_p,
     drst_dxyz_p_on_mortar_Gauss_porder,
     sjvol_p_on_f_p_mortar_Gauss,
     nvol_p_on_f_p_mortar_Gauss,
     geom->p4est_geom,
     dgmath_jit_dbase,
     NULL
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
        dudx_m_on_f_m_mortar_Gauss[j][k] += drst_dxyz_m_on_mortar_Gauss[i][j][k]*dudr_m_on_f_m_mortar_Gauss[i][k];
        dudx_p_on_f_p_mortar_Gauss_porder[j][k] += drst_dxyz_p_on_mortar_Gauss_porder[i][j][k]*dudr_p_on_f_p_mortar_Gauss_porder[i][k];

        /* if (e_m[0]->id == 0 && e_p[0]->id == 4){ */
        /*   printf("drst_dxyz_p_on_mortar_Gauss_porder[%d][%d][%d], drst_dxyz_m_on_mortar_Gauss[%d][%d][%d] = %.25f, %.25f\n",i,j,k,i,j,k, drst_dxyz_p_on_mortar_Gauss_porder[i][j][k], drst_dxyz_m_on_mortar_Gauss[i][j][k]); */
        /* } */
        
      }
    }    
  }

  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_Gauss_porder[0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_Gauss_porder[1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_Gauss_porder[2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_Gauss_porder[0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_Gauss_porder[1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_Gauss_porder[2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[0][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[0][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[0][2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[1][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[1][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[1][2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[2][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[2][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_Gauss_porder[2][2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[2], total_side_nodes_m_Lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[1], total_side_nodes_m_Lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[0], total_side_nodes_m_Lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[0][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[0][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[0][2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[1][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[1][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[1][2], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[2][0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[2][1], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_Gauss[2][2], total_nodes_mortar_Gauss); */

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
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss_porder[b]);
    }

    /* printf("face, face_p, deg_mortar_Gauss[face], deg_mortar_Gauss_p[face_p] = %d,%d,%d,%d\n", face, face_p, deg_mortar_Gauss[face], deg_mortar_Gauss_p[face_p]); */
    
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
      
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &nvol_p_on_f_p_mortar_Gauss[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_Gauss[face],
         orientation,
         f_m,
         f_p,
         &nvol_p_on_f_p_mortar_Gauss_reoriented[d][face_mortar_stride]
        );

      

      linalg_vec_scale(-1., &nvol_p_on_f_p_mortar_Gauss_reoriented[d][face_mortar_stride], dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]));
    }
    
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face]);
  }
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int i = 0; i < total_nodes_mortar_Gauss; i++){
      mortar_flux[d][i] = dudx_p_on_f_p_mortar_Gauss[d][i] + dudx_m_on_f_m_mortar_Gauss[d][i];
    }
  }
  printf("\n***NEW FACE FLUX CALCULATION***\n");
  printf("Elements on m-face %d\n", f_m);
  for (int f = 0; f < faces_m; f++){
    printf("e_m[%d]->id = %d, ",f, e_m[f]->id);
  }
  printf("\nElements on p-face %d\n", f_p);
  for (int f = 0; f < faces_p; f++){
    printf("e_p[%d]->id = %d, ",f, e_p[f]->id);
  }
  printf("\n");
  for (int d = 0; d < (P4EST_DIM); d++){
    DEBUG_PRINT_ARR_DBL_SUM(mortar_flux[d], total_nodes_mortar_Gauss);
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

  

  /* DEBUG_PRINT_4ARR_DBL(sj_on_f_m_mortar_Gauss, */
  /*                      sjvol_m_on_f_m_mortar_Gauss, */
  /*                      sjvol_p_on_f_p_mortar_Gauss, */
  /*                      sjvol_p_on_f_p_mortar_Gauss_reoriented, */
  /*                      total_nodes_mortar_Gauss); */

  double maxerror = util_max_error(sj_on_f_m_mortar_Gauss, sjvol_m_on_f_m_mortar_Gauss, total_nodes_mortar_Gauss);
  maxerror += util_max_error(sj_on_f_m_mortar_Gauss, sjvol_p_on_f_p_mortar_Gauss_reoriented, total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[0], nvol_p_on_f_p_mortar_Gauss_reoriented[0], total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[0], nvol_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[1], nvol_p_on_f_p_mortar_Gauss_reoriented[1], total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[1], nvol_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss);
  maxerror += util_max_error(dudx_p_on_f_p_mortar_Gauss[0], dudx_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss);
  maxerror += util_max_error(dudx_p_on_f_p_mortar_Gauss[1], dudx_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss);
#if (P4EST_DIM)==3
  maxerror += util_max_error(dudx_p_on_f_p_mortar_Gauss[2], dudx_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[2], nvol_p_on_f_p_mortar_Gauss_reoriented[2], total_nodes_mortar_Gauss);
  maxerror += util_max_error(n_on_f_m_mortar_Gauss[2], nvol_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss);
#endif
  
  /* DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[0], dudx_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss); */
  /* DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[1], dudx_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
/* #if (P4EST_DIM)==3 */
/*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[2], dudx_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss); */
/* #endif */
  // if (e_m[0]->id == 0 && e_p[0]->id == 4){
    /* if (maxerror > data->local_eps){ */
      /* printf("Holy shit batman, LOTS OF ERROR HERE, Error = %.25f\n", maxerror); */
      /* printf("faces_m, faces_p = %d,%d\n", faces_m, faces_p); */
      /* printf("e_m[0]->tree, e_p[0]->tree = %d,%d\n", e_m[0]->tree, e_p[0]->tree); */
      /* printf("e_m[0]->id, e_p[0]->id = %d,%d\n", e_m[0]->id, e_p[0]->id); */
      /* printf("f_m, f_p = %d,%d\n", f_m, f_p); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[0], dudx_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[1], dudx_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_Gauss[2], dudx_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss); */
    /* DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_Gauss_porder[0], dudr_m_on_f_m_mortar_Gauss[0], total_nodes_mortar_Gauss); */
    /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_Gauss_porder[1], dudr_m_on_f_m_mortar_Gauss[1], total_nodes_mortar_Gauss); */
    /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_Gauss_porder[2], dudr_m_on_f_m_mortar_Gauss[2], total_nodes_mortar_Gauss); */
    /* } */
    /* else { */
      /* printf("Holy fucks duckbutt, NO ERROR HERE, Error = %.25f\n", maxerror); */
      /* printf("faces_m, faces_p = %d,%d\n", faces_m, faces_p); */
      /* printf("e_m[0]->tree, e_p[0]->tree = %d,%d\n", e_m[0]->tree, e_p[0]->tree); */
      /* printf("e_m[0]->id, e_p[0]->id = %d,%d\n", e_m[0]->id, e_p[0]->id); */
    /* } */
    // }
  data->global_err += maxerror;

  P4EST_FREE(u_p_on_f_p_porder);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  P4EST_FREE(sjvol_m_on_f_m_mortar_Gauss);
  P4EST_FREE(sjvol_p_on_f_p_mortar_Gauss);
  P4EST_FREE(sjvol_p_on_f_p_mortar_Gauss_reoriented);
  D4EST_FREE_DIM_VEC(tmpxyz);
  D4EST_FREE_DIM_VEC(n_on_f_m_mortar_Gauss);
  D4EST_FREE_DIM_VEC(mortar_flux);
  D4EST_FREE_DIM_VEC(nvol_m_on_f_m_mortar_Gauss);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_Gauss);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_Gauss_reoriented);
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_Gauss);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_Gauss_porder);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_porder); 
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_mortar_porder);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_mortar_Gauss_porder);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_Gauss_porder);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_Gauss);
  D4EST_FREE_DIM_VEC(dudr_m_on_f_m);
  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_mortar);
  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_mortar_Gauss);
  D4EST_FREE_DIM_VEC(dudx_m_on_f_m_mortar_Gauss);
  
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
