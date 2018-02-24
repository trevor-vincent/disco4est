
#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_mortars_compute_flux.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_util.h>
#include <limits.h>

#define DEG_LOBATTO 2
#define DEG_QUAD 2

void
problem_set_degrees
(
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  elem_data->deg = DEG_LOBATTO;
  elem_data->deg_quad = DEG_QUAD;
}

typedef struct {
  
  double global_err;
  double local_eps;
  d4est_operators_t* d4est_ops;
  
} test_mortarjacobianterms_data_t;

static int
uni_refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}


p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}


static void
curved_test_mortarjacobianterms_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
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


  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
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
      deg_mortar_quad[i+j] = d4est_util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = d4est_util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      
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
      d4est_operators_apply_slicer
        (
         d4est_ops,
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
      d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_m[i]->u_elem[0],
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
      d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_p[i]->dudr_elem[d][0],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &dudr_p_on_f_p_porder[d][stride]
        );
    }
    stride += d4est_lgl_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
  }

  stride = 0;
  for (int i = 0; i < faces_p; i++){
      d4est_operators_apply_slicer
        (
         d4est_ops,
         &e_p[i]->u_elem[0],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &u_p_on_f_p_porder[stride]
        );
    stride += d4est_lgl_get_nodes((P4EST_DIM)-1, e_p[i]->deg);
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

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     mortar_side_id_m,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     drst_dxyz_m_on_mortar_quad,
     sjvol_m_on_f_m_mortar_quad,
     nvol_m_on_f_m_mortar_quad,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_INTEGRAND_UNKNOWN,
     e_p[0]->tree,
     e_p[0]->q,
     e_p[0]->dq,
     mortar_side_id_p,
     faces_p,
     faces_mortar,
     &deg_mortar_quad_porder[0],
     f_p,
     drst_dxyz_p_on_mortar_quad_porder,
     sjvol_p_on_f_p_mortar_quad,
     nvol_p_on_f_p_mortar_quad,
     NULL,
     NULL,
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


    d4est_linalg_fill_vec
      (
       dudx_m_on_f_m_mortar_quad_analytic[d],
       0.0,
       total_nodes_mortar_quad
      );


    d4est_linalg_fill_vec
      (
       dudx_p_on_f_p_mortar_quad_porder_analytic[d],
       0.0,
       total_nodes_mortar_quad
      );
    
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

    /* printf("face, face_p, deg_mortar_quad[face], deg_mortar_quad_p[face_p] = %d,%d,%d,%d\n", face, face_p, deg_mortar_quad[face], deg_mortar_quad_p[face_p]); */
    
    d4est_operators_reorient_face_data
      (
       d4est_ops,
       &sjvol_p_on_f_p_mortar_quad[oriented_face_mortar_stride],
       (P4EST_DIM)-1,
       deg_mortar_quad[face],
       orientation,
       f_m,
       f_p,
       &sjvol_p_on_f_p_mortar_quad_reoriented[face_mortar_stride]
      );

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
      
      d4est_operators_reorient_face_data
        (
         d4est_ops,
         &nvol_p_on_f_p_mortar_quad[d][oriented_face_mortar_stride],
         (P4EST_DIM)-1,
         deg_mortar_quad[face],
         orientation,
         f_m,
         f_p,
         &nvol_p_on_f_p_mortar_quad_reoriented[d][face_mortar_stride]
        );

      

      d4est_linalg_vec_scale(-1., &nvol_p_on_f_p_mortar_quad_reoriented[d][face_mortar_stride], d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]));


      d4est_linalg_vec_scale(-1., &nvol_p_on_f_p_mortar_quad_reoriented_analytic[d][face_mortar_stride], d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]));
      
    }
    
    face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad[face]);
  }
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int i = 0; i < total_nodes_mortar_quad; i++){
      mortar_flux[d][i] = dudx_p_on_f_p_mortar_quad[d][i] + dudx_m_on_f_m_mortar_quad[d][i];
    }
  }
  DEBUG_PRINT_3ARR_DBL(sjvol_m_on_f_m_mortar_quad,
                       sjvol_p_on_f_p_mortar_quad,
                       sjvol_p_on_f_p_mortar_quad_reoriented,
                       total_nodes_mortar_quad);
  

  double max_error = d4est_util_max_error(sjvol_p_on_f_p_mortar_quad_reoriented, sjvol_m_on_f_m_mortar_quad, total_nodes_mortar_quad);
  printf("max error 1 = %f\n", max_error);
  max_error += d4est_util_max_error(nvol_m_on_f_m_mortar_quad[0],  nvol_p_on_f_p_mortar_quad_reoriented[0], total_nodes_mortar_quad);
  max_error += d4est_util_max_error(nvol_m_on_f_m_mortar_quad[1],  nvol_p_on_f_p_mortar_quad_reoriented[1], total_nodes_mortar_quad);
  
#if (P4EST_DIM)==3
    max_error += d4est_util_max_error(nvol_m_on_f_m_mortar_quad[2],  nvol_p_on_f_p_mortar_quad_reoriented[2], total_nodes_mortar_quad);
#endif
  printf("max error 6 = %f\n", max_error);
  
  data->global_err += max_error;

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

d4est_mortars_fcn_ptrs_t
curved_test_mortarjacobianterms_fetch_fcns
(
 test_mortarjacobianterms_data_t* data
)
{
  
  d4est_mortars_fcn_ptrs_t curved_test_mortarjacobianterms_fcns;
  curved_test_mortarjacobianterms_fcns.flux_interface_fcn = curved_test_mortarjacobianterms_interface;
  curved_test_mortarjacobianterms_fcns.flux_boundary_fcn = NULL;
  curved_test_mortarjacobianterms_fcns.params = (void*)data;

  return curved_test_mortarjacobianterms_fcns;
}


int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_CUBED_SPHERE_OUTER_SHELL;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = 10;
  sphere_attrs->R2 = 1000000;
  sphere_attrs->compactify_outer_shell = 1;
  sphere_attrs->compactify_inner_shell = -1;
  d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geom, sphere_attrs);
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init();
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG4;
  d4est_quadrature_compactified_new(p4est, d4est_ops, d4est_geom, d4est_quad, "", "");


  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  int num_unifrefs = 1;
  for (int level = 0; level < num_unifrefs; ++level){

      p4est_refine_ext(p4est,
                       0,
                       -1,
                       uni_refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FACE,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);

  }


  d4est_mesh_update
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     geometric_factors,
     INITIALIZE_QUADRATURE_DATA,
     INITIALIZE_GEOMETRY_DATA,
     INITIALIZE_GEOMETRY_ALIASES,
     problem_set_degrees,
     NULL
    );
  
  test_mortarjacobianterms_data_t test_data;
  test_data.global_err = 0.;
  test_data.local_eps = .00000000001;
  test_data.d4est_ops = d4est_ops;
  d4est_mortars_fcn_ptrs_t ffp = curved_test_mortarjacobianterms_fetch_fcns(&test_data);
  
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &ffp,
     EXCHANGE_GHOST_DATA
    );
    


  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_geometry_destroy(d4est_geom);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
