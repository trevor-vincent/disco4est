#include "../GridFunctions/grid_functions.h"
#include "../ElementData/element_data.h"
#include "../dGMath/d4est_operators.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../Flux/compute_flux.h"

static void
dg_norm_ip_flux_dirichlet
(
 element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double ip_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t dg_norm_ip_flux_prefactor_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;

  
  /* grid_fcn_t u_at_bndry = bndry_fcn; */
  double penalty = dg_norm_ip_flux_prefactor_calculate_fcn
                             (
                              e_m->deg,
                              e_m->h,
                              e_m->deg,
                              e_m->h,
                              ip_flux_penalty_prefactor
                             );
  
  int face_nodes_m = d4est_operators_get_nodes ( (P4EST_DIM) - 1, e_m->deg );
  double* tmp = P4EST_ALLOC(double, face_nodes_m);
  double* xyz_on_f_m [(P4EST_DIM)];
  double* u_m_on_f_m = P4EST_ALLOC(double, face_nodes_m);
  
  /* get solution variable on this face */
  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);

  /* get xyz on this face */
  int dir;
  for (dir = 0; dir < (P4EST_DIM); dir++){
    xyz_on_f_m[dir] = P4EST_ALLOC(double, face_nodes_m);

    double* rst = d4est_operators_fetch_xyz_nd(d4est_ops, (P4EST_DIM), e_m->deg, dir);
    d4est_operators_apply_slicer(d4est_ops, rst, (P4EST_DIM), f_m, e_m->deg, tmp);
    
    d4est_operators_rtox_array(tmp, e_m->xyz_corner[dir], e_m->h, xyz_on_f_m[dir], face_nodes_m);
  }

  /* get boundary values on this face */
  int i;
/*   for (i = 0; i < face_nodes_m; i++){ */
/*     tmp[i] = u_at_bndry(xyz_on_f_m[0][i], */
/*                              xyz_on_f_m[1][i] */
/* #if (P4EST_DIM)==3 */
/*                              , */
/*                              xyz_on_f_m[2][i] */
/* #endif                           */
/*                             ); */
/*   } */

  double prefactor = sqrt(penalty);
  double n [3];
  d4est_operators_get_normal(f_m, (P4EST_DIM), &n[0]);
  for (dir = 0; dir < (P4EST_DIM); dir++){
    for(i = 0; i < face_nodes_m; i++){
      e_m->qstar_min_q[dir][f_m*face_nodes_m + i] = prefactor*n[dir]*(u_m_on_f_m[i]);
      /* printf("tmp[i] = %f\n", tmp[i]); */
    }
  }

  for (dir = 0; dir < (P4EST_DIM); dir++){
    P4EST_FREE(xyz_on_f_m[dir]);
  }
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(tmp);
}


static void
dg_norm_ip_flux_interface
(
 element_data_t** e_m,
 int faces_m,
 int f_m,
 element_data_t** e_p,
 int faces_p,
 int f_p,
 int* e_m_is_ghost,
 d4est_operators_t* d4est_ops,
 void* params
)
{
  ip_flux_params_t* ip_flux_params = (ip_flux_params_t*) params;
  double ip_flux_penalty_prefactor = ip_flux_params->ip_flux_penalty_prefactor;
  penalty_calc_t dg_norm_ip_flux_prefactor_calculate_fcn = ip_flux_params->ip_flux_penalty_calculate_fcn;


  
  int deg_p [(P4EST_HALF)];
  int max_deg_p = -1;
  int face_nodes_p [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int face_nodes_m [(P4EST_HALF)];
  int max_deg_m = -1;

  int nodes_mortar[(P4EST_HALF)];
  int deg_mortar [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar[(P4EST_HALF)];
  
  double n [(P4EST_DIM)];
  d4est_operators_get_normal(f_m, (P4EST_DIM), &n[0]);
  
  int i,j;
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m = 0;
  for (i = 0; i < faces_m; i++){
    deg_m[i] = e_m[i]->deg;
    if (e_m[i]->deg > max_deg_m) max_deg_m = e_m[i]->deg;
    face_nodes_m[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg );
    total_side_nodes_m += face_nodes_m[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p = 0;
  for (i = 0; i < faces_p; i++){
    deg_p[i] = e_p[i]->deg;
    if (e_p[i]->deg > max_deg_p) max_deg_p = e_p[i]->deg;
    face_nodes_p[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_p[i]->deg );
    total_side_nodes_p += face_nodes_p[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar = 0;
  for (i = 0; i < faces_m; i++)
    for (j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar[i+j] = util_max_int( e_m[i]->deg, e_p[j]->deg );

      /* find IP penalty parameter for each face pair of the two sides*/
      penalty_mortar[i+j] = dg_norm_ip_flux_prefactor_calculate_fcn
                            (
                             e_m[i]->deg,
                             e_m[i]->h,
                             e_p[j]->deg,
                             e_p[j]->h,
                             ip_flux_penalty_prefactor
                            );

      nodes_mortar[i+j] = d4est_operators_get_nodes( (P4EST_DIM) - 1, deg_mortar[i+j] );

      total_nodes_mortar += nodes_mortar[i+j];
    }

  /* scalar and vector fields on each of the (-) and (+) elements */
  double* du_m = P4EST_ALLOC(double, d4est_operators_get_nodes((P4EST_DIM), max_deg_m));
  double* du_p = P4EST_ALLOC(double, d4est_operators_get_nodes((P4EST_DIM), max_deg_p));
  /* slices of scalar/vector fields of (-) onto f_m and (+) onto f_p */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* du_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m); 
  double* du_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* du_m_on_f_m_mortar= P4EST_ALLOC(double, total_nodes_mortar); 
  double* du_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);
  double* u_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);

  double* qstar_min_q_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* qstar_min_q_m = P4EST_ALLOC(double, total_side_nodes_m);

  int stride;
  stride = 0;
  for (i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
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
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_p[i]->u_storage[0]),
       (P4EST_DIM),
       f_p,
       e_p[i]->deg,
       &u_p_on_f_p[stride]
      );
    stride += face_nodes_p[i];
  }
  
  /* project (-)-side u trace vector onto mortar space */ 
  d4est_operators_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar
    );

  /* project (+)-side u trace vector onto mortar space */
  d4est_operators_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar
    );

  int dir;
  /* For each component of the vector */
  for (dir = 0; dir < (P4EST_DIM); dir++){

    /* compute the (-)-u-derivative and project on the (-)-side faces and project q onto the (-)-side faces */
    stride = 0;
    for (i = 0; i < faces_m; i++){
      d4est_operators_apply_Dij(d4est_ops, e_m[i]->u_elem, (P4EST_DIM), e_m[i]->deg, dir, du_m);
      linalg_vec_scale(2./e_m[i]->h, du_m, d4est_operators_get_nodes((P4EST_DIM), e_m[i]->deg));

      d4est_operators_apply_slicer
        (
         d4est_ops,
         du_m,
         (P4EST_DIM),
         f_m,
         e_m[i]->deg,
         &du_m_on_f_m[stride]
        );

      stride += face_nodes_m[i];
    }

    /* compute the (+)-u-derivative and project on the (+)-side faces and project q onto the (+)-side faces */
    stride = 0;
    for (i = 0; i < faces_p; i++){
      d4est_operators_apply_Dij(d4est_ops, &(e_p[i]->u_storage[0]), (P4EST_DIM), e_p[i]->deg, dir, du_p);
      linalg_vec_scale(2./e_p[i]->h, du_p, d4est_operators_get_nodes((P4EST_DIM), e_p[i]->deg));

      d4est_operators_apply_slicer
        (
         d4est_ops,
         du_p,
         (P4EST_DIM),
         f_p,
         e_p[i]->deg,
         &du_p_on_f_p[stride]
        );

      stride += face_nodes_p[i];
    }

    /* project the derivatives from (-) and (+) sides onto the mortar space */
    d4est_operators_project_side_onto_mortar_space
      (
       d4est_ops,
       du_m_on_f_m,
       faces_m,
       deg_m,
       du_m_on_f_m_mortar,
       faces_mortar,
       deg_mortar
      );

    d4est_operators_project_side_onto_mortar_space
      (
       d4est_ops,
       du_p_on_f_p,
       faces_p,
       deg_p,
       du_p_on_f_p_mortar,
       faces_mortar,
       deg_mortar
      );

    /* calculate symmetric interior penalty flux */
    int k;
    int f;
    stride = 0;
    for (f = 0; f < faces_mortar; f++){
      double prefactor = sqrt(penalty_mortar[f]);
      for (k = 0; k < nodes_mortar[f]; k++){
        int ks = k + stride;
        qstar_min_q_mortar[ks] = n[dir]*u_m_on_f_m_mortar[ks];
        qstar_min_q_mortar[ks] -= n[dir]*u_p_on_f_p_mortar[ks];
        qstar_min_q_mortar[ks] *= prefactor;
      }
      stride += nodes_mortar[f];
    }

    /* project mortar data back onto the (-) side */
    d4est_operators_project_mortar_onto_side
      (
       d4est_ops,
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
  P4EST_FREE(du_p_on_f_p_mortar);
  P4EST_FREE(du_m_on_f_m_mortar);
  P4EST_FREE(u_m_on_f_m_mortar);
  P4EST_FREE(u_p_on_f_p);
  P4EST_FREE(du_p_on_f_p);
  P4EST_FREE(du_m_on_f_m);
  P4EST_FREE(u_m_on_f_m);
  P4EST_FREE(du_p);
  P4EST_FREE(du_m);
  P4EST_FREE(qstar_min_q_m);
  P4EST_FREE(qstar_min_q_mortar);
}

flux_fcn_ptrs_t
dg_norm_ip_flux_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* ip_params
)
{
  flux_fcn_ptrs_t dg_norm_ip_flux_fcns;
  dg_norm_ip_flux_fcns.flux_interface_fcn
    = dg_norm_ip_flux_interface;

  dg_norm_ip_flux_fcns.flux_boundary_fcn
    = dg_norm_ip_flux_dirichlet;

  dg_norm_ip_flux_fcns.bndry_fcn = bndry_fcn;
  dg_norm_ip_flux_fcns.params = (void*)ip_params;

  return dg_norm_ip_flux_fcns;
}
