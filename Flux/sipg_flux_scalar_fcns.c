#include "../Flux/sipg_flux_scalar_fcns.h"
#include "../Utilities/util.h"
#include "../dGMath/d4est_operators.h"
#include "../LinearAlgebra/d4est_linalg.h"

static void
sipg_flux_scalar_dirichlet
(
 element_data_t* e_m,
 int f_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 void* params
)
{
  grid_fcn_t u_at_bndry = bndry_fcn;
  /* int vol_nodes_m = d4est_lgl_get_nodes( (P4EST_DIM), e_m->deg); */
  int face_nodes_m = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m->deg);
  double* tmp = P4EST_ALLOC(double, face_nodes_m);
  double* xyz_on_f_m [(P4EST_DIM)];
  int dir, i;
  double* u_m_on_f_m = P4EST_ALLOC(double,
                                   d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg)
                                  );

  d4est_operators_apply_slicer(d4est_ops, e_m->u_elem, (P4EST_DIM), f_m, e_m->deg, u_m_on_f_m);
  
  for (dir = 0; dir < (P4EST_DIM); dir++){
    xyz_on_f_m[dir] = P4EST_ALLOC(double, face_nodes_m);

    double* rst = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), e_m->deg, dir);
    d4est_operators_apply_slicer(d4est_ops, rst, (P4EST_DIM), f_m, e_m->deg, tmp);

    d4est_reference_rtox_array(tmp,
                      e_m->xyz_corner[dir],
                      e_m->h,
                      xyz_on_f_m[dir],
                      face_nodes_m);
  }
  

  for (i = 0; i < face_nodes_m; i++){
    (e_m->ustar_min_u)[f_m*face_nodes_m + i] = u_at_bndry(xyz_on_f_m[0][i],
                             xyz_on_f_m[1][i]
#if (P4EST_DIM)==3
                             ,
                             xyz_on_f_m[2][i]
#endif                          
                            ) - u_m_on_f_m[i];

    /* printf("u_at_bndry(%f,%f) = %f\n", xyz_on_f_m[0][i], xyz_on_f_m[1][i], u_at_bndry(xyz_on_f_m[0][i], */
                                                                                      /* xyz_on_f_m[1][i]) ); */
  }

  for (dir = 0; dir < (P4EST_DIM); dir++){
    P4EST_FREE(xyz_on_f_m[dir]);
  }
  P4EST_FREE(tmp);
  P4EST_FREE(u_m_on_f_m);
  /* P4EST_FREE(rst); */
}


static void
sipg_flux_scalar_interface
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
  int stride;
  int deg_p [(P4EST_HALF)];
  int max_deg_p = -1;
  int face_nodes_p [(P4EST_HALF)];
  int deg_m [(P4EST_HALF)];
  int face_nodes_m [(P4EST_HALF)];
  int max_deg_m = -1;
  
  int deg_mortar [(P4EST_HALF)];
  int nodes_mortar [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  double n [(P4EST_DIM)];
  d4est_reference_get_normal(f_m, (P4EST_DIM), &n[0]);
  
  int i,j;
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m = 0;
  for (i = 0; i < faces_m; i++){
    deg_m[i] = e_m[i]->deg;
    /* if (sum_ghost_array != 0) */
      /* printf(" sum_ghost_array = %d, deg_m[i] = %d\n" , sum_ghost_array, deg_m[i]); */
    if (e_m[i]->deg > max_deg_m) max_deg_m = e_m[i]->deg;
    face_nodes_m[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg );
    total_side_nodes_m += face_nodes_m[i];
  }
  
  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p = 0;
  for (i = 0; i < faces_p; i++){
    deg_p[i] = e_p[i]->deg;
    if (e_p[i]->deg > max_deg_p) max_deg_p = e_p[i]->deg;
    face_nodes_p[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p[i]->deg );
    total_side_nodes_p += face_nodes_p[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar = 0;
  for (i = 0; i < faces_m; i++)
    for (j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar[i+j] = util_max_int( e_m[i]->deg, e_p[j]->deg );
      nodes_mortar[i+j] = d4est_lgl_get_nodes ( (P4EST_DIM) - 1, deg_mortar[i+j] );
      total_nodes_mortar += d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar[i+j] );
    }

  /* scalar and scalar fields on each of the (-) and (+) elements */
  double* u_m_on_f_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* ustar_min_u_m = P4EST_ALLOC(double, total_side_nodes_m);
  double* u_p_on_f_p = P4EST_ALLOC(double, total_side_nodes_p);
  /* projections of f_m slices to max_deg space */
  double* u_m_on_f_m_mortar = P4EST_ALLOC(double, total_nodes_mortar);
  double* u_p_on_f_p_mortar= P4EST_ALLOC(double, total_nodes_mortar);
  double* ustar_min_u_mortar= P4EST_ALLOC(double, total_nodes_mortar);

  stride = 0;
  for (i = 0; i < faces_m; i++){
    d4est_operators_apply_slicer
      (
       d4est_ops,
       &(e_m[i]->u_elem[0]),
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
       &(e_p[i]->u_elem[0]),
       (P4EST_DIM),
       f_p,
       e_p[i]->deg,
       &u_p_on_f_p[stride]
      );
    stride += face_nodes_p[i];
  }
  
  /* project (-)-side u trace scalar onto mortar space */ 
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m,
     u_m_on_f_m_mortar,
     faces_mortar,
     deg_mortar
    );

  /* project (+)-side u trace scalar onto mortar space */
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_p_on_f_p,
     faces_p,
     deg_p,
     u_p_on_f_p_mortar,
     faces_mortar,
     deg_mortar
    );

  /* calculate symmetric interior penalty flux */
  int k;
  int f;
  stride = 0;
    for (f = 0; f < faces_mortar; f++){
      for (k = 0; k < nodes_mortar[f]; k++){
        int ks = k + stride;
        ustar_min_u_mortar[ks] = .5*(u_p_on_f_p_mortar[ks] + u_m_on_f_m_mortar[ks]);
        ustar_min_u_mortar[ks] -= u_m_on_f_m_mortar[ks];
      }
      stride += nodes_mortar[f];
    }

    /* project mortar data back onto the (-) side */
    d4est_mortars_project_mortar_onto_side
      (
       d4est_ops,
       ustar_min_u_mortar,
       faces_mortar,
       deg_mortar,
       ustar_min_u_m,
       faces_m,
       deg_m
      );

    /* copy result back to element */
    stride = 0;
    for (i = 0; i < faces_m; i++){
      if(e_m_is_ghost[i] == 0)
        d4est_linalg_copy_1st_to_2nd
          (
           &ustar_min_u_m[stride],
           &(e_m[i]->ustar_min_u[f_m*face_nodes_m[i]]),
           face_nodes_m[i]
          );
      stride += face_nodes_m[i];                     
    }
    
    P4EST_FREE(ustar_min_u_m);
    P4EST_FREE(ustar_min_u_mortar);
    P4EST_FREE(u_p_on_f_p_mortar);
    P4EST_FREE(u_m_on_f_m_mortar);
    P4EST_FREE(u_p_on_f_p);
    P4EST_FREE(u_m_on_f_m);
}

flux_fcn_ptrs_t
sipg_flux_scalar_dirichlet_fetch_fcns
(
 grid_fcn_t bndry_fcn
)
{
  flux_fcn_ptrs_t sipg_flux_scalar_fcns;
  sipg_flux_scalar_fcns.flux_interface_fcn
    = sipg_flux_scalar_interface;
  sipg_flux_scalar_fcns.flux_boundary_fcn
    = sipg_flux_scalar_dirichlet;
  sipg_flux_scalar_fcns.bndry_fcn = bndry_fcn;
  sipg_flux_scalar_fcns.params = NULL;
  
  return sipg_flux_scalar_fcns;
}
