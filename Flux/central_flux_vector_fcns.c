#include "../Flux/central_flux_vector_fcns.h"
#include "../Utilities/util.h"
#include "../dGMath/dgmath.h"
#include "../LinearAlgebra/linalg.h"

static void
central_flux_vector_dirichlet
(
 element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 void* params
)
{
  central_flux_params_t* central_params = (central_flux_params_t*) params;
  double penalty = central_params->central_flux_penalty_prefactor;
  
  grid_fcn_t u_at_bndry = bndry_fcn;
  int face_nodes_m = dgmath_get_nodes ( (P4EST_DIM) - 1, e_m->deg );
  /* int vol_nodes_m = dgmath_get_nodes ( (P4EST_DIM) , e_m->deg ); */
  double* tmp = P4EST_ALLOC(double, face_nodes_m);
  double* xyz_on_f_m [(P4EST_DIM)];
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);
  double* q_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);
  
  /* get solution variable on this face */
  dgmath_apply_slicer(dgmath_jit_dbase, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);

  /* get xyz on this face */
  int dir;
  for (dir = 0; dir < (P4EST_DIM); dir++){
    xyz_on_f_m[dir] = P4EST_ALLOC(double, face_nodes_m);

    double* rst = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), e_m->deg, dir);
    dgmath_apply_slicer(dgmath_jit_dbase, rst, (P4EST_DIM), f_m, e_m->deg, tmp);
    
    dgmath_rtox_array(tmp, e_m->xyz_corner[dir], e_m->h, xyz_on_f_m[dir], face_nodes_m);
  }

  /* get boundary values on this face */
  int i;
  for (i = 0; i < face_nodes_m; i++){
    tmp[i] = u_at_bndry(xyz_on_f_m[0][i],
                             xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                             ,
                             xyz_on_f_m[2][i]
#endif                          
                            );
  }

  double n [3];
  dgmath_get_normal(f_m, (P4EST_DIM), &n[0]);

  for (dir = 0; dir < (P4EST_DIM); dir++){
    dgmath_apply_slicer(dgmath_jit_dbase, e_m->q_elem[dir], (P4EST_DIM), f_m, e_m->deg, q_m_on_f_m);
    
    /* calculate qstar - q(-) */
    for(i = 0; i < face_nodes_m; i++){
      e_m->qstar_min_q[dir][f_m*face_nodes_m + i] = q_m_on_f_m[i] - penalty*n[dir]*(u_m_on_f_m[i] - tmp[i]);
      e_m->qstar_min_q[dir][f_m*face_nodes_m + i] -= q_m_on_f_m[i];
    }
  }

  for (dir = 0; dir < (P4EST_DIM); dir++){
    P4EST_FREE(xyz_on_f_m[dir]);
  }
  P4EST_FREE(q_m_on_f_m);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(tmp);
}

static void
central_flux_vector_interface
(
 element_data_t** e_m,
 int faces_m,
 int f_m,
 element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 void* params
)
{
  central_flux_params_t* central_params = (central_flux_params_t*) params;
  double curved_central_flux_penalty_prefactor = central_params->central_flux_penalty_prefactor;
  
  int stride;
  int deg_p [(P4EST_HALF)];
  int max_deg_p = -1;
  int face_nodes_p [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int face_nodes_m [(P4EST_HALF)];
  int max_deg_m = -1;
  int nodes_mortar [(P4EST_HALF)];

  /* double penalty_mortar [(P4EST_HALF)]; */
  int deg_mortar [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  double n [(P4EST_DIM)];
  dgmath_get_normal(f_m, (P4EST_DIM), &n[0]);
  
  int i,j;
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m = 0;
  for (i = 0; i < faces_m; i++){
    deg_m[i] = e_m[i]->deg;
    if (e_m[i]->deg > max_deg_m) max_deg_m = e_m[i]->deg;
    face_nodes_m[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg );
    total_side_nodes_m += face_nodes_m[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p = 0;
  for (i = 0; i < faces_p; i++){
    deg_p[i] = e_p[i]->deg;
    if (e_p[i]->deg > max_deg_p) max_deg_p = e_p[i]->deg;
    face_nodes_p[i] = dgmath_get_nodes( (P4EST_DIM) - 1, e_p[i]->deg );
    total_side_nodes_p += face_nodes_p[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar = 0;
  for (i = 0; i < faces_m; i++)
    for (j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar[i+j] = util_max_int( e_m[i]->deg, e_p[j]->deg );
      /* printf("penalty_mortar[%d] = %f\n", i+j, penalty_mortar[i+j]); */
      nodes_mortar[i+j] = dgmath_get_nodes( (P4EST_DIM) - 1, deg_mortar[i+j] );
      
      total_nodes_mortar += nodes_mortar[i+j];
    }

  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* q_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  double* q_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* q_m_on_f_m_mortar= P4EST_ALLOC(double, total_nodes_mortar);  
  double* q_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);  
  double* u_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);

  double* qstar_min_q_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* qstar_min_q_m = P4EST_ALLOC(double, total_side_nodes_m);

  stride = 0;
  for (i = 0; i < faces_m; i++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_m[i]->u_storage[0]),
       (P4EST_DIM),
       f_m,
       e_m[i]->deg,
       &u_m_on_f_m[stride]
      );
    stride += face_nodes_m[i];
  }
 
  stride = 0;
  for (i = 0; i < faces_p; i++){
    dgmath_apply_slicer
      (
       dgmath_jit_dbase,
       &(e_p[i]->u_storage[0]),
       (P4EST_DIM),
       f_p,
       e_p[i]->deg,
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
  int dir;
  for (dir = 0; dir < (P4EST_DIM); dir++){

    /* compute the (-)-u-derivative and project on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (i = 0; i < faces_m; i++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         e_m[i]->q_elem[dir],
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &q_m_on_f_m[stride]
        );

      stride += face_nodes_m[i];
    }

    /* compute the (+)-u-derivative and project on the (+)-side faces and project q onto the (+)-side faces */

#ifndef D4EST_DEBUG
    mpi_abort("Must be in debug mode to run this central flux sipg, otherwise change q_elem to q_storage");
#endif
    
    stride = 0;
    for (i = 0; i < faces_p; i++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         e_p[i]->q_elem[dir],
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &q_p_on_f_p[stride]
        );

      stride += face_nodes_p[i];
    }

     /* project q from the (-) side onto the mortar space */
    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       q_m_on_f_m,
       faces_m,
       deg_m,
       q_m_on_f_m_mortar,
       faces_mortar,
       deg_mortar
      );


    dgmath_project_side_onto_mortar_space
      (
       dgmath_jit_dbase,
       q_p_on_f_p,
       faces_p,
       deg_p,
       q_p_on_f_p_mortar,
       faces_mortar,
       deg_mortar
      );


    

    /* calculate symmetric interior penalty flux */
    int k;
    int f;
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      double sigma = curved_central_flux_penalty_prefactor;
      for (k = 0; k < nodes_mortar[f]; k++){
        int ks = k + stride;
        qstar_min_q_mortar[ks] = .5*(q_p_on_f_p_mortar[ks] + q_m_on_f_m_mortar[ks]);
        qstar_min_q_mortar[ks] -= -sigma*n[dir]*u_p_on_f_p_mortar[ks];
        qstar_min_q_mortar[ks] -= sigma*n[dir]*u_m_on_f_m_mortar[ks];
        qstar_min_q_mortar[ks] -= q_m_on_f_m_mortar[ks];
      }
      stride += nodes_mortar[f];
    }

    /* project mortar data back onto the (-) side */
    dgmath_project_mortar_onto_side
      (
       dgmath_jit_dbase,
       qstar_min_q_mortar,
       faces_mortar,
       deg_mortar,
       qstar_min_q_m,
       faces_m,
       deg_m
      );

    /* copy result back to element */
    stride = 0;
    for (i = 0; i < faces_m; i++){
      if(e_m_is_ghost[i] == 0)
        linalg_copy_1st_to_2nd
          (
           &qstar_min_q_m[stride],
           &(e_m[i]->qstar_min_q[dir][f_m*face_nodes_m[i]]),
           face_nodes_m[i]
          );
      stride += face_nodes_m[i];                     
    }
  }

  P4EST_FREE(u_p_on_f_p_mortar);
  P4EST_FREE(q_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(q_p_on_f_p_mortar);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(q_m_on_f_m);
  P4EST_FREE(q_p_on_f_p);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(qstar_min_q_m);
  P4EST_FREE(qstar_min_q_mortar);
}

flux_fcn_ptrs_t
central_flux_vector_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 central_flux_params_t* central_params
)
{  
  flux_fcn_ptrs_t central_flux_vector_fcns;
  central_flux_vector_fcns.flux_interface_fcn = central_flux_vector_interface;
  central_flux_vector_fcns.flux_boundary_fcn = central_flux_vector_dirichlet;
  central_flux_vector_fcns.bndry_fcn = bndry_fcn;
  central_flux_vector_fcns.params = (void*)central_params;
  
  return central_flux_vector_fcns;
}
