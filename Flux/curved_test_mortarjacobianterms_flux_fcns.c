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
  int volume_nodes_lobatto = dgmath_get_nodes(dim,deg);
  int face_nodes_lobatto = dgmath_get_nodes(dim-1,deg);
  int volume_nodes_quad = dgmath_get_nodes(dim, elem_data->deg_integ);
  
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
     volume_nodes_lobatto
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
  int face_nodes_m_lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

  double* sj_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* sjvol_on_f_m_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  double* sjvol_on_f_m_quad_analytic = P4EST_ALLOC(double, face_nodes_m_quad);
  double* n_on_f_m_quad [(P4EST_DIM)];
  double* nvol_on_f_m_quad [(P4EST_DIM)];
  double* nvol_on_f_m_quad_analytic [(P4EST_DIM)];
  double* xyz_on_f_m_quad [(P4EST_DIM)];
  double* drst_dxyz_on_f_m_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_on_f_m_quad_analytic [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DIM_VEC(n_on_f_m_quad,face_nodes_m_quad);
  D4EST_ALLOC_DIM_VEC(nvol_on_f_m_quad,face_nodes_m_quad);
  D4EST_ALLOC_DIM_VEC(nvol_on_f_m_quad_analytic,face_nodes_m_quad);
  D4EST_ALLOC_DIM_VEC(xyz_on_f_m_quad,face_nodes_m_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_on_f_m_quad,face_nodes_m_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_on_f_m_quad_analytic,face_nodes_m_quad);

  /* mpi_assert(geom->X_mapping_type == MAP_ISOPARAMETRIC); */
  /* curved_element_data_compute_mortar_normal_and_sj_using_face_data */
  /*   ( */
  /*    &e_m, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    geom->dxdr_method, */
  /*    1, */
  /*    n_on_f_m_quad, */
  /*    sj_on_f_m_quad, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    xyz_on_f_m_quad */
  /*   ); */


  geometric_quantity_compute_method_t mapping_orig  = geom->DX_compute_method;


  geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     drst_dxyz_on_f_m_quad ,
     sjvol_on_f_m_quad,
     nvol_on_f_m_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );


  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m->tree,
     e_m->q,
     e_m->dq,
     1,
     1,
     &e_m->deg_integ,
     f_m,
     drst_dxyz_on_f_m_quad_analytic,
     sjvol_on_f_m_quad_analytic,
     nvol_on_f_m_quad_analytic,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  
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
     sj_on_f_m_quad,
     n_on_f_m_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_CROSS_PRODUCT
    );


  /* geom->X_mapping_type = MAP_ISOPARAMETRIC; */
  /* d4est_geometry_compute_geometric_data_on_mortar */
  /*   ( */
  /*    e_m->tree, */
  /*    e_m->q, */
  /*    e_m->dq, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    drst_dxyz_on_f_m_quad , */
  /*    sjvol_on_f_m_quad, */
  /*    nvol_on_f_m_quad, */
  /*    NULL, */
  /*    LOBATTO, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    COMPUTE_NORMAL_USING_JACOBIAN */
  /*   ); */


  /* geom->X_mapping_type = MAP_ANALYTIC; */
  /* d4est_geometry_compute_geometric_data_on_mortar */
  /*   ( */
  /*    e_m->tree, */
  /*    e_m->q, */
  /*    e_m->dq, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    drst_dxyz_on_f_m_quad_analytic, */
  /*    sjvol_on_f_m_quad_analytic, */
  /*    nvol_on_f_m_quad_analytic, */
  /*    NULL, */
  /*    LOBATTO, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    COMPUTE_NORMAL_USING_JACOBIAN */
  /*   ); */

  
  /* geom->X_mapping_type = MAP_ISOPARAMETRIC; */
  /* d4est_geometry_compute_geometric_data_on_mortar */
  /*   ( */
  /*    e_m->tree, */
  /*    e_m->q, */
  /*    e_m->dq, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    NULL, */
  /*    sj_on_f_m_quad, */
  /*    n_on_f_m_quad, */
  /*    NULL, */
  /*    LOBATTO, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    COMPUTE_NORMAL_USING_CROSS_PRODUCT */
  /*   ); */


  /* d4est_geometry_compute_geometric_data_on_mortar_TESTINGONLY */
  /*   ( */
  /*    e_m->tree, */
  /*    e_m->q, */
  /*    e_m->dq, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    LOBATTO, */
  /*    n_on_f_m_quad, */
  /*    sj_on_f_m_quad, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    xyz_on_f_m_quad */
  /*   ); */

  geom->DX_compute_method = mapping_orig;

  /* DEBUG_PRINT_2ARR_DBL(sjvol_on_f_m_quad_analytic, sjvol_on_f_m_quad, face_nodes_m_quad); */
  
  /* curved_data_compute_drst_dxyz_quad_on_mortar_using_volume_data */
  /*   ( */
  /*    &e_m, */
  /*    1, */
  /*    1, */
  /*    &e_m->deg_integ, */
  /*    f_m, */
  /*    drst_dxyz_on_f_m_quad, */
  /*    sjvol_on_f_m_quad, */
  /*    nvol_on_f_m_quad, */
  /*    geom->p4est_geom, */
  /*    dgmath_jit_dbase, */
  /*    NULL */
  /*   ); */

  /* DEBUG_PRINT_2ARR_DBL(sj_on_f_m_quad, */
  /*                      sjvol_on_f_m_quad, */
  /*                      face_nodes_m_quad); */

  double* n_error [(P4EST_DIM)];
  double* sj_error = P4EST_ALLOC(double, face_nodes_m_quad);
  D4EST_ALLOC_DIM_VEC(n_error, face_nodes_m_quad);
  
  double maxerror = util_max_error(sj_on_f_m_quad, sjvol_on_f_m_quad, face_nodes_m_quad);
  maxerror += util_max_error(n_on_f_m_quad[0], nvol_on_f_m_quad[0], face_nodes_m_quad);

  maxerror += util_max_error(n_on_f_m_quad[1], nvol_on_f_m_quad[1], face_nodes_m_quad);
  util_compute_error_array(n_on_f_m_quad[0], nvol_on_f_m_quad[0], n_error[0], face_nodes_m_quad);
  util_compute_error_array(n_on_f_m_quad[1], nvol_on_f_m_quad[1], n_error[1], face_nodes_m_quad);
  util_compute_error_array(sj_on_f_m_quad, sjvol_on_f_m_quad, sj_error, face_nodes_m_quad);
  
#if (P4EST_DIM)==3
  util_compute_error_array(n_on_f_m_quad[2], nvol_on_f_m_quad[2], n_error[2], face_nodes_m_quad);
  maxerror += util_max_error(n_on_f_m_quad[2], nvol_on_f_m_quad[2],face_nodes_m_quad);
#endif
  
  if (maxerror > data->local_eps){
    printf("Holy shit batman, LOTS OF N_ERROR ON THE BOUNDARY HERE, NOT SURPRISING BECAUSE THIS IS ISOPARAMETRIC, OR IS IT?\n");
    printf("face = %d\n", f_m);
    
    DEBUG_PRINT_4ARR_DBL(n_on_f_m_quad[0], nvol_on_f_m_quad[0], n_error[0],nvol_on_f_m_quad_analytic[0],face_nodes_m_quad);
    DEBUG_PRINT_4ARR_DBL(n_on_f_m_quad[1], nvol_on_f_m_quad[1], n_error[1],nvol_on_f_m_quad_analytic[1],face_nodes_m_quad);
    DEBUG_PRINT_4ARR_DBL(n_on_f_m_quad[2], nvol_on_f_m_quad[2], n_error[2],nvol_on_f_m_quad_analytic[2],face_nodes_m_quad);
    DEBUG_PRINT_4ARR_DBL(sj_on_f_m_quad, sjvol_on_f_m_quad, sj_error, sjvol_on_f_m_quad_analytic, face_nodes_m_quad);
  }

  D4EST_FREE_DIM_VEC(n_error);
  P4EST_FREE(sj_error);  
  
  data->global_err += maxerror;
  
  P4EST_FREE(sj_on_f_m_quad);
  P4EST_FREE(sjvol_on_f_m_quad);
  P4EST_FREE(sjvol_on_f_m_quad_analytic);
  D4EST_FREE_DIM_VEC(n_on_f_m_quad);
  D4EST_FREE_DIM_VEC(nvol_on_f_m_quad);
  D4EST_FREE_DIM_VEC(nvol_on_f_m_quad_analytic);
  D4EST_FREE_DIM_VEC(xyz_on_f_m_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_on_f_m_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_on_f_m_quad_analytic);
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
  int deg_p_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  /* int deg_p_quad [(P4EST_HALF)]; */
  int face_nodes_p_quad [(P4EST_HALF)];

  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  /* int deg_m_quad [(P4EST_HALF)]; */
  int face_nodes_m_quad [(P4EST_HALF)];
  
  int nodes_mortar_quad [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];


  curved_element_data_t* e_p_oriented [(P4EST_HALF)];
  curved_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    /* deg_m_quad[i] = e_m[i]->deg_integ; */
    
    face_nodes_m_lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_integ);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_quad[i] = e_p_oriented[i]->deg_integ; */

    face_nodes_p_lobatto[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_integ);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_quad = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = util_max_int( e_m[i]->deg_integ,
                                            e_p_oriented[j]->deg_integ);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      
    }


  
  int deg_mortar_quad_porder [(P4EST_HALF)];
  int nodes_mortar_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    nodes_mortar_quad_porder[inew] = nodes_mortar_quad[i];
  }


  double* dudr_p_on_f_p_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar [(P4EST_DIM)];
  double* dudr_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudr_m_on_f_m_mortar_quad [(P4EST_DIM)];

  double* dudx_p_on_f_p_mortar_quad_porder [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];

  double* dudx_p_on_f_p_mortar_quad_porder_analytic [(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad_analytic [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad_analytic [(P4EST_DIM)];

  double* mortar_flux [(P4EST_DIM)];

  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
  double* u_p_on_f_p_porder = P4EST_ALLOC(double, total_side_nodes_p_lobatto);

  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_porder, total_side_nodes_p_lobatto); 
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_mortar_porder, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudr_p_on_f_p_mortar_quad_porder, total_nodes_mortar_quad);

  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m, total_side_nodes_m_lobatto);
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_mortar, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(mortar_flux, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudr_m_on_f_m_mortar_quad, total_nodes_mortar_quad);

  D4EST_ALLOC_DIM_VEC(dudx_m_on_f_m_mortar_quad, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_quad_porder, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_quad, total_nodes_mortar_quad);

  D4EST_ALLOC_DIM_VEC(dudx_m_on_f_m_mortar_quad_analytic, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_quad_porder_analytic, total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(dudx_p_on_f_p_mortar_quad_analytic, total_nodes_mortar_quad);
  
  
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
    stride += face_nodes_m_lobatto[i];
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
      
    stride += face_nodes_m_lobatto[i];
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
  /*   DEBUG_PRINT_2ARR_DBL(u_p_on_f_p_porder, u_m_on_f_m, total_side_nodes_m_lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[0], dudr_m_on_f_m[0], total_side_nodes_m_lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[1], dudr_m_on_f_m[1], total_side_nodes_m_lobatto); */
  /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_porder[2], dudr_m_on_f_m[2], total_side_nodes_m_lobatto); */

  /* } */
  

  for (int d = 0; d < (P4EST_DIM); d++){
  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     dudr_p_on_f_p_porder[d],
     faces_p,
     deg_p_lobatto_porder,
     dudr_p_on_f_p_mortar_porder[d],
     faces_mortar,
     deg_mortar_quad_porder
    );
  

  dgmath_project_side_onto_mortar_space
    (
     dgmath_jit_dbase,
     dudr_m_on_f_m[d],
     faces_m,
     deg_m_lobatto,
     dudr_m_on_f_m_mortar[d],
     faces_mortar,
     deg_mortar_quad
    );
  }
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){

    dgmath_interp(dgmath_jit_dbase,
                  &dudr_m_on_f_m_mortar[d][stride],
                  QUAD_LOBATTO,
                  deg_mortar_quad[f],
                  &dudr_m_on_f_m_mortar_quad[d][stride],
                  geom->geom_quad_type,
                  deg_mortar_quad[f],
                  (P4EST_DIM)-1);
    }
    stride += nodes_mortar_quad[f];
  }

  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){

      dgmath_interp(dgmath_jit_dbase,
                    &dudr_p_on_f_p_mortar_porder[d][stride],
                    QUAD_LOBATTO,
                    deg_mortar_quad_porder[f],
                    &dudr_p_on_f_p_mortar_quad_porder[d][stride],
                    geom->geom_quad_type,
                    deg_mortar_quad_porder[f],
                    (P4EST_DIM)-1);
    }
    stride += nodes_mortar_quad_porder[f];
  }

  double* sj_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];
  double* tmpxyz [(P4EST_DIM)];

  D4EST_ALLOC_DIM_VEC(tmpxyz,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(n_on_f_m_mortar_quad,total_nodes_mortar_quad);

  double* sjvol_m_on_f_m_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_p_on_f_p_mortar_quad = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_p_on_f_p_mortar_quad_reoriented = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* nvol_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_quad_reoriented [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];

  D4EST_ALLOC_DIM_VEC(nvol_m_on_f_m_mortar_quad,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad_reoriented,total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_quad,total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder,total_nodes_mortar_quad);


  double* sjvol_m_on_f_m_mortar_quad_analytic = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_p_on_f_p_mortar_quad_analytic = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_p_on_f_p_mortar_quad_reoriented_analytic = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* nvol_m_on_f_m_mortar_quad_analytic [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_quad_analytic [(P4EST_DIM)];
  double* nvol_p_on_f_p_mortar_quad_reoriented_analytic [(P4EST_DIM)];
  double* drst_dxyz_m_on_mortar_quad_analytic [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_p_on_mortar_quad_porder_analytic [(P4EST_DIM)][(P4EST_DIM)];

  D4EST_ALLOC_DIM_VEC(nvol_m_on_f_m_mortar_quad_analytic,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad_analytic,total_nodes_mortar_quad);
  D4EST_ALLOC_DIM_VEC(nvol_p_on_f_p_mortar_quad_reoriented_analytic,total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_m_on_mortar_quad_analytic,total_nodes_mortar_quad);
  D4EST_ALLOC_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder_analytic,total_nodes_mortar_quad);


  geometric_quantity_compute_method_t mapping_orig = geom->DX_compute_method;
  geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  /* geom->X_mapping_type = MAP_ISOPARAMETRIC; */
  /* d4est_geometry_compute_geometric_data_on_mortar_TESTINGONLY */
  /*   ( */
  /*    e_m[0]->tree, */
  /*    e_m[0]->q, */
  /*    e_m[0]->dq, */
  /*    faces_m, */
  /*    faces_mortar, */
  /*    &deg_mortar_quad[0], */
  /*    f_m, */
  /*    GAUSS, */
  /*    n_on_f_m_mortar_quad, */
  /*    sj_on_f_m_mortar_quad, */
  /*    geom, */
  /*    dgmath_jit_dbase, */
  /*    tmpxyz */
  /*   ); */


  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     NULL,
     sj_on_f_m_mortar_quad,
     n_on_f_m_mortar_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_CROSS_PRODUCT
    );

  
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
     sjvol_m_on_f_m_mortar_quad,
     nvol_m_on_f_m_mortar_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
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
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder,
     sjvol_p_on_f_p_mortar_quad,
     nvol_p_on_f_p_mortar_quad,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  if(geom->DX != NULL){
  geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  /* geom->X_mapping_type = MAP_ISOPARAMETRIC; */
  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     drst_dxyz_m_on_mortar_quad_analytic,
     sjvol_m_on_f_m_mortar_quad_analytic,
     nvol_m_on_f_m_mortar_quad_analytic,
     NULL,
     NULL,
     geom->geom_quad_type,
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
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder_analytic,
     sjvol_p_on_f_p_mortar_quad_analytic,
     nvol_p_on_f_p_mortar_quad_analytic,
     NULL,
     NULL,
     geom->geom_quad_type,
     geom,
     dgmath_jit_dbase,
     COMPUTE_NORMAL_USING_JACOBIAN
    ); 
  }
  geom->DX_compute_method = mapping_orig;
  
  
  

  for (int d = 0; d < (P4EST_DIM); d++){
    linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad[d],
       0.0,
       total_nodes_mortar_quad
      );


    linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_quad_porder[d],
       0.0,
       total_nodes_mortar_quad
      );


    linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad_analytic[d],
       0.0,
       total_nodes_mortar_quad
      );


    linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_quad_porder_analytic[d],
       0.0,
       total_nodes_mortar_quad
      );
    
  }





  
  for (int j = 0; j < (P4EST_DIM); j++){
    for (int i = 0; i < (P4EST_DIM); i++){   
      for (int k = 0; k < total_nodes_mortar_quad; k++){
        dudx_m_on_f_m_mortar_quad[j][k] += drst_dxyz_m_on_mortar_quad[i][j][k]*dudr_m_on_f_m_mortar_quad[i][k];
        dudx_p_on_f_p_mortar_quad_porder[j][k] += drst_dxyz_p_on_mortar_quad_porder[i][j][k]*dudr_p_on_f_p_mortar_quad_porder[i][k];

          if(geom->DX != NULL){
        dudx_m_on_f_m_mortar_quad_analytic[j][k] += drst_dxyz_m_on_mortar_quad_analytic[i][j][k]*dudr_m_on_f_m_mortar_quad[i][k];
        dudx_p_on_f_p_mortar_quad_porder_analytic[j][k] += drst_dxyz_p_on_mortar_quad_porder_analytic[i][j][k]*dudr_p_on_f_p_mortar_quad_porder[i][k];
          }        
        /* if (e_m[0]->id == 0 && e_p[0]->id == 4){ */
        /*   printf("drst_dxyz_p_on_mortar_quad_porder[%d][%d][%d], drst_dxyz_m_on_mortar_quad[%d][%d][%d] = %.25f, %.25f\n",i,j,k,i,j,k, drst_dxyz_p_on_mortar_quad_porder[i][j][k], drst_dxyz_m_on_mortar_quad[i][j][k]); */
        /* } */
        
      }
    }    
  }

  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_quad_porder[0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_quad_porder[1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_p_on_f_p_mortar_quad_porder[2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_quad_porder[0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_quad_porder[1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_p_on_f_p_mortar_quad_porder[2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[0][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[0][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[0][2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[1][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[1][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[1][2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[2][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[2][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_p_on_mortar_quad_porder[2][2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudx_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[2], total_side_nodes_m_lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[1], total_side_nodes_m_lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(dudr_m_on_f_m[0], total_side_nodes_m_lobatto); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[0][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[0][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[0][2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[1][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[1][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[1][2], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[2][0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[2][1], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_ARR_DBL_SUM(drst_dxyz_m_on_mortar_quad[2][2], total_nodes_mortar_quad); */

  /* int oriented_stride = 0; */
  /* int oriented_strides [P4EST_HALF]; */
  /* for (int i = 0; i < faces_p; i++){ */
  /*   for (int j = 0; j < faces_p; j++){ */
  /*     if (e_p[j]->tree_quadid == e_p_oriented[i]->tree_quadid){ */
  /*       oriented_strides[j] = oriented_stride; */
  /*       break; */
  /*     } */
  /*   } */
  /*   oriented_stride += nodes_mortar_quad[i]; */
  /* } */

  int face_mortar_stride = 0;
  for (int face = 0; face < faces_mortar; face++){
    int face_p = face;
    if (faces_mortar == (P4EST_HALF))
      face_p = dgmath_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

    int oriented_face_mortar_stride = 0;
    for (int b = 0; b < face_p; b++){
      oriented_face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
    }

    /* printf("face, face_p, deg_mortar_quad[face], deg_mortar_quad_p[face_p] = %d,%d,%d,%d\n", face, face_p, deg_mortar_quad[face], deg_mortar_quad_p[face_p]); */
    
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       &sjvol_p_on_f_p_mortar_quad[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_quad[face],
       orientation,
       f_m,
       f_p,
       &sjvol_p_on_f_p_mortar_quad_reoriented[face_mortar_stride]
      );

     if(geom->DX != NULL){ 
    dgmath_reorient_face_data
      (
       dgmath_jit_dbase,
       &sjvol_p_on_f_p_mortar_quad_analytic[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_quad[face],
       orientation,
       f_m,
       f_p,
       &sjvol_p_on_f_p_mortar_quad_reoriented_analytic[face_mortar_stride]
      );
     }    

    for (int d = 0; d < (P4EST_DIM); d++){
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &dudx_p_on_f_p_mortar_quad_porder[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_quad[d][face_mortar_stride]
        );
      
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &nvol_p_on_f_p_mortar_quad[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &nvol_p_on_f_p_mortar_quad_reoriented[d][face_mortar_stride]
        );

  if(geom->DX != NULL){
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &dudx_p_on_f_p_mortar_quad_porder_analytic[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &dudx_p_on_f_p_mortar_quad_analytic[d][face_mortar_stride]
        );
      
      dgmath_reorient_face_data
        (
         dgmath_jit_dbase,
         &nvol_p_on_f_p_mortar_quad_analytic[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &nvol_p_on_f_p_mortar_quad_reoriented_analytic[d][face_mortar_stride]
        );
  }
      

      linalg_vec_scale(-1., &nvol_p_on_f_p_mortar_quad_reoriented[d][face_mortar_stride], dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]));


      linalg_vec_scale(-1., &nvol_p_on_f_p_mortar_quad_reoriented_analytic[d][face_mortar_stride], dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]));
      
    }
    
    face_mortar_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
  }
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int i = 0; i < total_nodes_mortar_quad; i++){
      mortar_flux[d][i] = dudx_p_on_f_p_mortar_quad[d][i] + dudx_m_on_f_m_mortar_quad[d][i];
    }
  }
  /* printf("\n***NEW FACE FLUX CALCULATION***\n"); */
  /* printf("Elements on m-face %d\n", f_m); */
  /* for (int f = 0; f < faces_m; f++){ */
  /*   printf("e_m[%d]->id = %d, ",f, e_m[f]->id); */
  /* } */
  /* printf("\nElements on p-face %d\n", f_p); */
  /* for (int f = 0; f < faces_p; f++){ */
  /*   printf("e_p[%d]->id = %d, ",f, e_p[f]->id); */
  /* } */
  /* printf("\n"); */
  /* for (int d = 0; d < (P4EST_DIM); d++){ */
  /*   DEBUG_PRINT_ARR_DBL_SUM(mortar_flux[d], total_nodes_mortar_quad); */
  /* } */
  /* printf("faces_m, faces_p = %d,%d\n", faces_m, faces_p); */

  
  /* stride = 0; */
  /* for (int f = 0; f < faces_mortar; f++){ */
  /*   dgmath_reorient_face_data */
  /*     ( */
  /*      dgmath_jit_dbase, */
  /*      &sjvol_p_on_f_p_mortar_quad[stride], */
  /*      ((P4EST_DIM) - 1), */
  /*      deg_mortar_quad[f], */
  /*      orientation, */
  /*      f_m, */
  /*      f_p, */
  /*      &sjvol_p_on_f_p_mortar_quad_reoriented[stride] */
  /*     );       */
  /*   stride += nodes_mortar_quad[f]; */
  /* }    */

  

  /* DEBUG_PRINT_4ARR_DBL(sj_on_f_m_mortar_quad, */
  /*                      sjvol_m_on_f_m_mortar_quad, */
  /*                      sjvol_p_on_f_p_mortar_quad, */
  /*                      sjvol_p_on_f_p_mortar_quad_reoriented, */
  /*                      total_nodes_mortar_quad); */

  double maxerror = util_max_error(sj_on_f_m_mortar_quad, sjvol_m_on_f_m_mortar_quad, total_nodes_mortar_quad);
  maxerror += util_max_error(sj_on_f_m_mortar_quad, sjvol_p_on_f_p_mortar_quad_reoriented, total_nodes_mortar_quad);
  maxerror += util_max_error(n_on_f_m_mortar_quad[0], nvol_p_on_f_p_mortar_quad_reoriented[0], total_nodes_mortar_quad);
  maxerror += util_max_error(n_on_f_m_mortar_quad[0], nvol_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad);
  maxerror += util_max_error(n_on_f_m_mortar_quad[1], nvol_p_on_f_p_mortar_quad_reoriented[1], total_nodes_mortar_quad);
  maxerror += util_max_error(n_on_f_m_mortar_quad[1], nvol_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad);
  /* maxerror += util_max_error(dudx_p_on_f_p_mortar_quad[0], dudx_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
  /* maxerror += util_max_error(dudx_p_on_f_p_mortar_quad[1], dudx_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */

  double maxerror_analytic = 0.;
  if(geom->DX != NULL){
  maxerror_analytic += util_max_error(sjvol_m_on_f_m_mortar_quad_analytic, sjvol_p_on_f_p_mortar_quad_reoriented_analytic, total_nodes_mortar_quad);
  maxerror_analytic += util_max_error(nvol_m_on_f_m_mortar_quad_analytic[0], nvol_p_on_f_p_mortar_quad_reoriented_analytic[0], total_nodes_mortar_quad);
  maxerror_analytic += util_max_error(nvol_m_on_f_m_mortar_quad_analytic[1], nvol_p_on_f_p_mortar_quad_reoriented_analytic[1], total_nodes_mortar_quad);
  /* maxerror_analytic += util_max_error(dudx_p_on_f_p_mortar_quad_analytic[0], dudx_m_on_f_m_mortar_quad_analytic[0], total_nodes_mortar_quad); */
  /* maxerror_analytic += util_max_error(dudx_p_on_f_p_mortar_quad_analytic[1], dudx_m_on_f_m_mortar_quad_analytic[1], total_nodes_mortar_quad); */

  /* DEBUG_PRINT_2ARR_DBL(sjvol_m_on_f_m_mortar_quad_analytic, sjvol_p_on_f_p_mortar_quad_reoriented_analytic, total_nodes_mortar_quad); */
  }
  
#if (P4EST_DIM)==3
  /* maxerror += util_max_error(dudx_p_on_f_p_mortar_quad[2], dudx_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
  maxerror += util_max_error(n_on_f_m_mortar_quad[2], nvol_p_on_f_p_mortar_quad_reoriented[2], total_nodes_mortar_quad);
  maxerror += util_max_error(n_on_f_m_mortar_quad[2], nvol_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad);

  if(geom->DX != NULL){
    /* maxerror_analytic += util_max_error(dudx_p_on_f_p_mortar_quad_analytic[2], dudx_m_on_f_m_mortar_quad_analytic[2], total_nodes_mortar_quad); */
    maxerror_analytic += util_max_error(nvol_m_on_f_m_mortar_quad_analytic[2], nvol_p_on_f_p_mortar_quad_reoriented_analytic[2], total_nodes_mortar_quad);
  }
#endif


  maxerror += maxerror_analytic;
  if (maxerror > data->local_eps){
    printf("HOLY FUCKS DONALD TRUMP, LOOK AT THIS ERROR\n");
    /* DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[0], dudr_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
      /* DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[1], dudr_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
      /* DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[2], dudr_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
  }
  
  /* DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[0], dudx_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
  /* DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[1], dudx_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
/* #if (P4EST_DIM)==3 */
/*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[2], dudx_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
/* #endif */
  // if (e_m[0]->id == 0 && e_p[0]->id == 4){
    /* if (maxerror > data->local_eps){ */
      /* printf("Holy shit batman, LOTS OF ERROR HERE, Error = %.25f\n", maxerror); */
      /* printf("faces_m, faces_p = %d,%d\n", faces_m, faces_p); */
      /* printf("e_m[0]->tree, e_p[0]->tree = %d,%d\n", e_m[0]->tree, e_p[0]->tree); */
      /* printf("e_m[0]->id, e_p[0]->id = %d,%d\n", e_m[0]->id, e_p[0]->id); */
      /* printf("f_m, f_p = %d,%d\n", f_m, f_p); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[0], dudx_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[1], dudx_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_2ARR_DBL(dudx_p_on_f_p_mortar_quad[2], dudx_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
    /* DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[0], dudr_m_on_f_m_mortar_quad[0], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[1], dudr_m_on_f_m_mortar_quad[1], total_nodes_mortar_quad); */
    /*   DEBUG_PRINT_2ARR_DBL(dudr_p_on_f_p_mortar_quad_porder[2], dudr_m_on_f_m_mortar_quad[2], total_nodes_mortar_quad); */
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
  P4EST_FREE(sj_on_f_m_mortar_quad);

  D4EST_FREE_DIM_VEC(tmpxyz);
  D4EST_FREE_DIM_VEC(n_on_f_m_mortar_quad);
  D4EST_FREE_DIM_VEC(mortar_flux);
  P4EST_FREE(sjvol_m_on_f_m_mortar_quad);
  P4EST_FREE(sjvol_p_on_f_p_mortar_quad);
  P4EST_FREE(sjvol_p_on_f_p_mortar_quad_reoriented);
  D4EST_FREE_DIM_VEC(nvol_m_on_f_m_mortar_quad);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad_reoriented);
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_quad);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder);


  P4EST_FREE(sjvol_m_on_f_m_mortar_quad_analytic);
  P4EST_FREE(sjvol_p_on_f_p_mortar_quad_analytic);
  P4EST_FREE(sjvol_p_on_f_p_mortar_quad_reoriented_analytic);
  D4EST_FREE_DIM_VEC(nvol_m_on_f_m_mortar_quad_analytic);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad_analytic);
  D4EST_FREE_DIM_VEC(nvol_p_on_f_p_mortar_quad_reoriented_analytic);
  D4EST_FREE_DBYD_MAT(drst_dxyz_m_on_mortar_quad_analytic);
  D4EST_FREE_DBYD_MAT(drst_dxyz_p_on_mortar_quad_porder_analytic);


  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_porder); 
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_mortar_porder);
  D4EST_FREE_DIM_VEC(dudr_p_on_f_p_mortar_quad_porder);

  D4EST_FREE_DIM_VEC(dudr_m_on_f_m);
  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_mortar);
  D4EST_FREE_DIM_VEC(dudr_m_on_f_m_mortar_quad);
  D4EST_FREE_DIM_VEC(dudx_m_on_f_m_mortar_quad);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_quad_porder);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_quad);
 
  D4EST_FREE_DIM_VEC(dudx_m_on_f_m_mortar_quad_analytic);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_quad_porder_analytic);
  D4EST_FREE_DIM_VEC(dudx_p_on_f_p_mortar_quad_analytic);
  
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
