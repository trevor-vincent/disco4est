#include "../GridFunctions/grid_functions.h"
#include "../ElementData/curved_element_data.h"
#include "../dGMath/dgmath.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Flux/curved_compute_flux.h"
#include <curved_dg_norm.h>
#include <ip_flux_aux.h>

/* See houston2007 for the definition of the dg-norm */

static void
curved_dg_norm_boundary
(
 curved_element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 void* params
)
{
  curved_dg_norm_params_t* curved_dg_norm_params = (curved_dg_norm_params_t*) params;
  ip_flux_params_t* ip_flux_params = curved_dg_norm_params->ip_flux_params;
  double curved_curved_dg_norm_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t curved_dg_norm_u_dirichlet_prefactor_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  
  int face_nodes_m_Lobatto = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_Gauss = dgmath_get_nodes((P4EST_DIM) - 1, e_m->deg_integ);

  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m_Lobatto);
  double* u_m_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  

  double* Mfaceterm = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* sj_on_f_m_Gauss = P4EST_ALLOC(double, face_nodes_m_Gauss);
  double* faceterm = P4EST_ALLOC(double, face_nodes_m_Gauss);

  double* xyz_on_f_m_Gauss [(P4EST_DIM)];
  double* n_on_f_m_Gauss [(P4EST_DIM)];
  double* sj_n_on_f_m_Gauss [(P4EST_DIM)];

  
  for (int d = 0; d < (P4EST_DIM); d++) {
    xyz_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
    sj_n_on_f_m_Gauss[d] = P4EST_ALLOC(double, face_nodes_m_Gauss);
  }
  
  dgmath_apply_slicer(dgmath_jit_dbase,
                      e_m->u_storage,
                      (P4EST_DIM),
                      f_m,
                      e_m->deg,
                      u_m_on_f_m);

  dgmath_interp_GLL_to_GL
    (
     dgmath_jit_dbase,
     u_m_on_f_m,
     e_m->deg,
     e_m->deg_integ,
     u_m_on_f_m_Gauss,
     (P4EST_DIM)-1
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
     sj_on_f_m_Gauss,
     n_on_f_m_Gauss,
     NULL,
     GAUSS,
     geom,
     dgmath_jit_dbase
    );
  
  double h = (e_m->volume/e_m->surface_area[f_m]);
  /* find IP penalty parameter for each face pair of the two sides*/
  double prefactor = curved_dg_norm_u_dirichlet_prefactor_calculate_fcn
                     (
                      e_m->deg,
                      h,
                      e_m->deg,
                      h,
                      curved_curved_dg_norm_penalty_prefactor
                     );
    


  
  for (int dim = 0; dim < (P4EST_DIM); dim++){
    /* calculate qstar - q(-) */

    for(int i = 0; i < face_nodes_m_Gauss; i++){
      faceterm[i] = n_on_f_m_Gauss[dim][i]*(u_m_on_f_m_Gauss[i]);
    }
    double facetermMfaceterm = dgmath_Gauss_quadrature(
                                                       dgmath_jit_dbase,
                                                       faceterm,
                                                       faceterm,
                                                       sj_on_f_m_Gauss,
                                                       e_m->deg_integ,
                                                       (P4EST_DIM)-1);
    
    curved_dg_norm_params->dg_norm_face_term += prefactor*facetermMfaceterm;



    //printf("dirichlet faceterm_prefactor_mortar[f] = %.25f\n",prefactor);
    //printf("dirichlet facetermMfaceterm = %.25f\n",facetermMfaceterm);
    //printf("dirichlet curved_dg_norm_params->dg_norm_face_term= %.25f\n",curved_dg_norm_params->dg_norm_face_term);
    
  }    

  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_m_on_f_m_Gauss);
  P4EST_FREE(sj_on_f_m_Gauss);
  P4EST_FREE(faceterm);
  P4EST_FREE(Mfaceterm);
  
  for (int d = 0; d < (P4EST_DIM); d++){
    P4EST_FREE(xyz_on_f_m_Gauss[d]);
    P4EST_FREE(n_on_f_m_Gauss[d]);
    P4EST_FREE(sj_n_on_f_m_Gauss[d]);
  }
}

static void
curved_dg_norm_interface
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
  curved_dg_norm_params_t* curved_dg_norm_params = (curved_dg_norm_params_t*) params;
  ip_flux_params_t* ip_flux_params = curved_dg_norm_params->ip_flux_params;
  double curved_dg_norm_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t curved_dg_norm_penalty_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

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
  double faceterm_prefactor_mortar [(P4EST_HALF)];
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
   
    
      faceterm_prefactor_mortar[i+j] = curved_dg_norm_penalty_calculate_fcn
                                  (
                                   e_m[i]->deg,
                                   (e_m[i]->volume/e_m[i]->surface_area[f_m]),
                                   e_p[j]->deg,
                                   (e_p_oriented[j]->volume/e_p_oriented[j]->surface_area[f_p]),
                                   curved_dg_norm_penalty_prefactor
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
  double* sj_on_f_m_mortar_Gauss = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  double* n_on_f_m_mortar_Gauss [(P4EST_DIM)];

  
  for (int i = 0; i < (P4EST_DIM); i++){
    n_on_f_m_mortar_Gauss[i] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }
 
  double* faceterm [(P4EST_DIM)]; 
  for (int d = 0; d < (P4EST_DIM); d++) {
    faceterm[d] = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
  }
  double* Mfaceterm = P4EST_ALLOC(double, total_nodes_mortar_Gauss);
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

  d4est_geometry_compute_geometric_data_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     &deg_mortar_Gauss[0],
     f_m,
     NULL,
     sj_on_f_m_mortar_Gauss,
     n_on_f_m_mortar_Gauss,
     NULL,
     GAUSS,
     geom,
     dgmath_jit_dbase
    );
  
  stride = 0;
  for (int f = 0; f < faces_mortar; f++){
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_m_on_f_m_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_m_on_f_m_mortar_Gauss[stride], (P4EST_DIM)-1);
    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &u_p_on_f_p_mortar[stride], deg_mortar_Gauss[f], deg_mortar_Gauss[f], &u_p_on_f_p_mortar_Gauss[stride], (P4EST_DIM)-1);
    stride += nodes_mortar_Gauss[f];
  }
  
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
      for (int d = 0; d < (P4EST_DIM); d++){
        n_ks = n_on_f_m_mortar_Gauss[d][ks];        
        faceterm[d][ks] = n_ks*u_m_on_f_m_mortar_Gauss[ks];
        faceterm[d][ks] -= n_ks*u_p_on_f_p_mortar_Gauss[ks];
      }
    }
    stride += nodes_mortar_Gauss[f];
  }
    
  /* the contribution in every direction must be added up due to it being a vector norm */
  stride = 0;
  for (f = 0; f < faces_mortar; f++){
    for (int d = 0; d < (P4EST_DIM); d++){


      double facetermMfaceterm = dgmath_Gauss_quadrature(
                                               dgmath_jit_dbase,
                                               &faceterm[d][stride],
                                               &faceterm[d][stride],
                                               &sj_on_f_m_mortar_Gauss[stride],
                                               deg_mortar_Gauss[f],
                                               (P4EST_DIM)-1);
        

      /* even if it's a ghost we can still add it to the ghosts estimator because we will not send it back! */
      /* we need a .5 because this will be counted twice because we visit each interface twice */
      curved_dg_norm_params->dg_norm_face_term += .5*faceterm_prefactor_mortar[f]*facetermMfaceterm;

      //printf("faceterm_prefactor_mortar[f] = %.25f\n",faceterm_prefactor_mortar[f]);
      //printf("facetermMfaceterm = %.25f\n",facetermMfaceterm);
      //printf(" curved_dg_norm_params->dg_norm_face_term= %.25f\n",curved_dg_norm_params->dg_norm_face_term);
      
    }
    stride += nodes_mortar_Gauss[f];
  }




  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar_Gauss);
  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p_mortar_Gauss);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(faceterm[i]);
    P4EST_FREE(n_on_f_m_mortar_Gauss[i]);
  }
  P4EST_FREE(sj_on_f_m_mortar_Gauss);
  P4EST_FREE(Mfaceterm);
}

curved_flux_fcn_ptrs_t
curved_dg_norm_fetch_fcns
(
 curved_dg_norm_params_t* curved_dg_params
)
{
  curved_flux_fcn_ptrs_t curved_dg_norm_fcns;
  curved_dg_norm_fcns.flux_interface_fcn
    = curved_dg_norm_interface;

  curved_dg_norm_fcns.flux_boundary_fcn
    = curved_dg_norm_boundary;

  curved_dg_norm_fcns.bndry_fcn = NULL;
  curved_dg_norm_fcns.params = (void*)curved_dg_params;

  return curved_dg_norm_fcns;
}
