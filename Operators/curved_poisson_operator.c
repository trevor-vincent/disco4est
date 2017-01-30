#include <pXest.h>
#include <curved_element_data.h>
#include <problem_data.h>
#include <dgmath.h>
#include <linalg.h>
#include <curved_poisson_operator.h>
#include <grid_functions.h>
#include <util.h>

/* #define DEALIASING */


typedef struct {
  
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  problem_data_t* problem_data;
#if D4EST_DEBUG
  curved_poisson_debug_vecs_t* debug_vecs;
#endif
    
} curved_Gauss_poisson_user_data_t;

/* static */
void curved_Gauss_poisson_init_vecs
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;

  curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data;
  problem_data_t* problem_data = (problem_data_t*) curved_Gauss_poisson_user_data->problem_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase;
  
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg);
  int face_nodes_Lobatto = dgmath_get_nodes(dim-1,deg);
  int volume_nodes_Gauss = dgmath_get_nodes(dim, elem_data->deg_integ);
  /* int face_nodes_Gauss = dgmath_get_nodes(dim-1, elem_data->deg_integ); */

  
  int i,j,k;
  for (i = 0; i < (P4EST_DIM); i++){
    elem_data->M_ustar_min_u_n[i] = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*face_nodes_Lobatto);
    //elem_data->q_elem[i] = P4EST_ALLOC_ZERO(double, volume_nodes);
    /* elem_data->du_elem[i] = P4EST_ALLOC_ZERO(double, volume_nodes);    */
  }
  
  elem_data->M_qstar_min_q_dot_n = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*face_nodes_Lobatto);  
  elem_data->Au_elem = &(problem_data->Au[elem_data->nodal_stride]);

  linalg_copy_1st_to_2nd
    (
     &(problem_data->u[elem_data->nodal_stride]),
     &(elem_data->u_elem)[0],
     volume_nodes_Lobatto
    );

  linalg_fill_vec
    (
     &(elem_data->du_elem[0][0]),
     0.0,
     (P4EST_DIM)*(MAX_NODES)
    );

  double* u_array =&(problem_data->u[elem_data->nodal_stride]);
  double* u_elem = &(elem_data->u_elem)[0];
  /* DEBUG_PRINT_2ARR_DBL(u_array, u_elem, volume_nodes_Lobatto); */
  
  /* double* du_di = P4EST_ALLOC(double, volume_nodes_Gauss); */
  double* du_di_prolonged = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* du_di_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  
  for (i = 0; i < (P4EST_DIM); i++){

    /* PROLONG u here first, then take derivative, then interpolate to Gauss, this is most consistent */
    /* dgmath_apply_p_prolong(dgmath_jit_dbase, &(problem_data->u[elem_data->nodal_stride]), deg, (P4EST_DIM), elem_data->deg_integ, u_prolonged); */
    
    /* dgmath_apply_Dij(dgmath_jit_dbase, u_prolonged, dim, elem_data->deg_integ, i, &elem_data->dudr_elem[i][0]); */
    /* dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &elem_data->dudr_elem[i][0], elem_data->deg_integ, elem_data->deg_integ, du_di_Gauss, (P4EST_DIM)); */

    /* dgmath_apply_p_prolong(dgmath_jit_dbase, &(problem_data->u[elem_data->nodal_stride]), deg, (P4EST_DIM), elem_data->deg_integ, u_prolonged); */
    dgmath_apply_Dij(dgmath_jit_dbase, &(problem_data->u[elem_data->nodal_stride]), dim, elem_data->deg, i, &elem_data->dudr_elem[i][0]);
    
    dgmath_apply_p_prolong(dgmath_jit_dbase, &elem_data->dudr_elem[i][0], deg, (P4EST_DIM), elem_data->deg_integ, du_di_prolonged);

    dgmath_interp_GLL_to_GL(dgmath_jit_dbase, du_di_prolonged, elem_data->deg_integ, elem_data->deg_integ, du_di_Gauss, (P4EST_DIM));


    
    for (j = 0; j < (P4EST_DIM); j++){
      for (k = 0; k < volume_nodes_Gauss; k++){
        elem_data->du_elem[j][k] += elem_data->rst_xyz_integ[i][j][k]*du_di_Gauss[k];
      }
    }    
  }

  /* double* tmp_ptr = &elem_data->dudr_elem[0][0]; */
  /* DEBUG_PRINT_ARR_DBL(tmp_ptr, dgmath_get_nodes((P4EST_DIM),elem_data->deg_integ)); */

  
  P4EST_FREE(du_di_prolonged);
  P4EST_FREE(du_di_Gauss);
}

/* static */
/* void curved_Gauss_poisson_compute_q_elem_old */
/* ( */
/*  p4est_iter_volume_info_t * info, */
/*  void *user_data */
/* ) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data; */

/*   curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data; */
/*   dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase; */
/*   int dim = (P4EST_DIM); */
/*   int deg = element_data->deg; */
/*   int faces = 2*dim; */
/*   int face_nodes_Gauss = dgmath_get_nodes(dim-1, elem_data->deg_integ); */
/*   int volume_nodes_Gauss = dgmath_get_nodes(dim, elem_data->deg_integ); */
/*   /\* double* u_elem = element_data->u_elem; *\/ */
/*   double* vol_tmp = P4EST_ALLOC(double, volume_nodes_Gauss); */
/*   double* Si_u [(P4EST_DIM)]; */
/*   int i; */

/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     Si_u[i] = P4EST_ALLOC(double, volume_nodes_Gauss); */
/*   } */
  
/*   for (int i = 0; i < (P4EST_DIM); i++){ */
/*     dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase, */
/*                                                 element_data->du_elem[i], */
/*                                                 elem_data->deg_integ, */
/*                                                 element_data->J_integ, */
/*                                                 elem_data->deg_integ, */
/*                                                 (P4EST_DIM), */
/*                                                 Si_u[i]); */
/*   } */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int f = 0; f < faces; f++){ */
/*       dgmath_apply_LIFT(dgmath_jit_dbase, &element_data->M_ustar_min_u_n[d][f*face_nodes], dim, elem_data->deg_integ, f, vol_tmp); */
/*       linalg_vec_axpy(1.0, vol_tmp, Si_u[d], volume_nodes_Gauss); */
/*     } */
/*   } */

/*   for (int i = 0; i < (P4EST_DIM); i++){ */
/*     dgmath_apply_curvedInverseGaussMass(dgmath_jit_dbase, Si_u[i], elem_data->deg_integ, element_data->J_integ, elem_data->deg_integ, (P4EST_DIM), element_data->q_elem[i]); */
/*   } */
   
/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     P4EST_FREE(Si_u[i]); */
/*   } */
/*   P4EST_FREE(vol_tmp); */
/* } */

void curved_Gauss_poisson_compute_q_elem
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data;

  curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase;
  int dim = (P4EST_DIM);
  /* int deg = element_data->deg; */
  int faces = 2*dim;
  /* int face_nodes_Gauss = dgmath_get_nodes(dim-1,element_data->deg_integ); */
  int face_nodes_Lobatto = dgmath_get_nodes(dim-1,element_data->deg);
  /* int volume_nodes_Gauss = dgmath_get_nodes(dim,element_data->deg_integ); */
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,element_data->deg);
  /* double* u_elem = element_data->u_elem; */

  double* Si_u [(P4EST_DIM)];
  double* Mq [(P4EST_DIM)];
  double* lifted_ustar_min_u [(P4EST_FACES)][(P4EST_DIM)];
  /* double* M_lifted_ustar_min_u [(P4EST_FACES)][(P4EST_DIM)]; */

/* #ifdef D4EST_DEBUG */
  /* curved_poisson_debug_vecs_t* debug_vecs = (curved_poisson_debug_vecs_t*)curved_Gauss_poisson_user_data->debug_vecs; */
/* #endif */


  /* int volume_nodes_integ = dgmath_get_nodes((P4EST_DIM), element_data->deg_integ); */
  /* DEBUG_PRINT_2ARR_DBL(element_data->du_elem[0], element_data->du_elem[1], volume_nodes_integ); */
  
  
  for (int i = 0; i < (P4EST_FACES); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      lifted_ustar_min_u[i][j] = P4EST_ALLOC(double, volume_nodes_Lobatto);
      /* M_lifted_ustar_min_u[i][j] = P4EST_ALLOC(double, volume_nodes_Lobatto); */
    }
  
  for (int i = 0; i < (P4EST_DIM); i++){
    Si_u[i] = P4EST_ALLOC(double, volume_nodes_Lobatto);
    Mq[i] = P4EST_ALLOC(double, volume_nodes_Lobatto);
  }
  
  for (int i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase,
                                                element_data->du_elem[i],
                                                element_data->deg,
                                                element_data->J_integ,
                                                element_data->deg_integ,
                                                (P4EST_DIM),
                                                Si_u[i]);

    /* util_print_matrix(Si_u[i], volume_nodes_Lobatto, 1, "Si_u[i] = ", 0); */
  }


  double* u_ptr = &element_data->u_elem[0];
  /* DEBUG_PRINT_4ARR_DBL(u_ptr, Si_u[0], Si_u[1], Si_u[2],volume_nodes_Lobatto); */
  /* DEBUG_PRINT_3ARR_DBL(element_data->u, element_data->J_integ, Si_u[0], volume_nodes_Lobatto); */

  /* double* Minv_Mustar_min_u_n = P4EST_ALLOC(double, face_nodes_Lobatto); */

  for (int d = 0; d < (P4EST_DIM); d++){
    for (int f = 0; f < faces; f++){

      /* linalg_matvec_plus_vec(1., */
      /*                        element_data->invMface[f], */
      /*                        &element_data->M_ustar_min_u_n[d][f*face_nodes_Lobatto], */
      /*                        0., */
      /*                        Minv_Mustar_min_u_n, */
      /*                        face_nodes_Lobatto, */
      /*                        face_nodes_Lobatto); */

      /* dgmath_apply_LIFT( */
      /*                   dgmath_jit_dbase, */
      /*                   Minv_Mustar_min_u_n, */
      /*                   dim, */
      /*                   element_data->deg, */
      /*                   f, */
      /*                   lifted_ustar_min_u[f][d] */
      /* ); */
      dgmath_apply_LIFT(
                        dgmath_jit_dbase,
                        &element_data->M_ustar_min_u_n[d][f*face_nodes_Lobatto],
                        dim,
                        element_data->deg,
                        f,
                        lifted_ustar_min_u[f][d]
                       );
      /* linalg_vec_axpy(1.0, vol_tmp, Si_u[d], volume_nodes); */
    }
  }

  /* P4EST_FREE(Minv_Mustar_min_u_n); */
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int i = 0; i < volume_nodes_Lobatto; i++) {
      Mq[d][i] = Si_u[d][i]; 
    }
    for (int i = 0; i < volume_nodes_Lobatto; i++) {
      for (int f = 0; f < faces; f++){
        Mq[d][i] += lifted_ustar_min_u[f][d][i];
      }
    }
  }


  
  for (int i = 0; i < (P4EST_DIM); i++){
    if (element_data->deg == element_data->deg_integ){
      dgmath_apply_curvedInverseGaussMass(dgmath_jit_dbase, Mq[i], element_data->deg, element_data->J_integ, element_data->deg_integ, (P4EST_DIM), element_data->q_elem[i]);
    }
    else {
      linalg_matvec_plus_vec(1., element_data->invM, Mq[i], 0., element_data->q_elem[i], volume_nodes_Lobatto, volume_nodes_Lobatto);
    }
  }

  /* DEBUG_PRINT_3ARR_DBL(element_data->q_elem[0],element_data->q_elem[1], element_data->q_elem[2], volume_nodes_Lobatto); */
  
  /* for (int d = 0; d < (P4EST_DIM); d++){ */
  /*   for (int i = 0; i < volume_nodes_Lobatto; i++) { */
  /*     for (int f = 0; f < faces; f++){ */
  /*       element_data->q_elem[d][i] += lifted_ustar_min_u[f][d][i]; */
  /*     } */
  /*   } */
  /* } */

  /* double* testinvMq = P4EST_ALLOC(double, volume_nodes_Lobatto); */
  /* dgmath_apply_invMij(dgmath_jit_dbase, Mq[0], (P4EST_DIM), element_data->deg, testinvMq); */
  /* linalg_vec_scale(1./element_data->J_integ[0], testinvMq, volume_nodes_Lobatto); */


  /* DEBUG_PRINT_3ARR_DBL( */
  /*                      (&element_data->q_elem[0][0]), */
  /*                      (&element_data->q_elem[1][0]), */
  /*                      (&element_data->q_elem[2][0]), */
  /*                      volume_nodes_Lobatto */
  /*                     ); */



  /* DEBUG_PRINT_3ARR_DBL( */
  /*                      (&Mq[0][0]), */
  /*                      (&Mq[1][0]), */
  /*                      (&Mq[2][0]), */
  /*                      volume_nodes_Lobatto */
  /*                     ); */
  
  
  /* P4EST_FREE(testinvMq); */
  
#ifdef D4EST_DEBUG
  curved_poisson_debug_vecs_t* debug_vecs = (curved_poisson_debug_vecs_t*)curved_Gauss_poisson_user_data->debug_vecs;
  if(debug_vecs != NULL && debug_vecs->elem_id == element_data->id){
  curved_poisson_debug_vecs_set_Mdu(Si_u, debug_vecs, dgmath_jit_dbase);
  curved_poisson_debug_vecs_set_u(&element_data->u_elem[0], debug_vecs, dgmath_jit_dbase);


    double* qptrs [(P4EST_DIM)];
    qptrs[0] = &element_data->q_elem[0][0];
    qptrs[1] = &element_data->q_elem[1][0];

#if (P4EST_DIM)==3
    qptrs[2] = &element_data->q_elem[2][0];
#endif
    
    curved_poisson_debug_vecs_set_q(qptrs, debug_vecs, dgmath_jit_dbase);

  curved_poisson_debug_vecs_set_Mq(Mq, debug_vecs, dgmath_jit_dbase);
  curved_poisson_debug_vecs_set_lifteduflux(lifted_ustar_min_u, debug_vecs, dgmath_jit_dbase);
  /* printf("element_id = %d\n", element_data->id); */
  /* DEBUG_PRINT_ARR_DBL(element_data->q_elem[0], volume_nodes_Lobatto); */
  /* DEBUG_PRINT_ARR_DBL(element_data->q_elem[1], volume_nodes_Lobatto); */
  }
#endif
  
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(Si_u[i]);
    P4EST_FREE(Mq[i]);
  }
  for (int i = 0; i < (P4EST_FACES); i++){
    for (int j = 0; j < (P4EST_DIM); j++){
      P4EST_FREE(lifted_ustar_min_u[i][j]);
    }
  }
  
}


/* /\* static *\/ */
/* void curved_Gauss_poisson_compute_Au_elem(p4est_iter_volume_info_t* info, void* user_data) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data; */

/*   curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data; */
/*   dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase; */
  
/*   int dim = (P4EST_DIM); */
/*   int deg = element_data->deg; */
/*   int faces = 2*dim; */
/*   int face_nodes = dgmath_get_nodes(dim-1,deg); */
/*   int volume_nodes = dgmath_get_nodes(dim,deg); */

/*   double* Au = element_data->Au_elem; */
/*   linalg_fill_vec(Au, 0., volume_nodes); */
  
/*   double* dq = P4EST_ALLOC(double, dgmath_get_nodes(dim,deg)); */
/*   double* dq_Gauss = P4EST_ALLOC(double, dgmath_get_nodes(dim,deg)); */
/*   double* vol_tmp = P4EST_ALLOC_ZERO(double, dgmath_get_nodes(dim,deg)); */

/*   /\* compute  *\/ */
/*   int i,j,k; */
/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     for (j = 0; j < (P4EST_DIM); j++){ */
/*       dgmath_apply_Dij(dgmath_jit_dbase, &element_data->q_elem[i][0], dim, deg, j, dq); */
/*       dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dq, deg, deg, dq_Gauss, (P4EST_DIM)); */
/*       for (k = 0; k < volume_nodes; k++){ */
/*         vol_tmp[k] += element_data->rst_xyz_integ[j][i][k]*dq_Gauss[k]; */
/*       } */
/*     }     */
/*   } */
  
/*   /\* linalg_component_mult(vol_tmp2, element_data->J, vol_tmp, volume_nodes); *\/ */
/*   /\* dgmath_apply_Mij(dgmath_jit_dbase, vol_tmp, dim, deg, Au); *\/ */

/*   dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase, vol_tmp, deg, element_data->J_integ, deg, (P4EST_DIM), Au); */
  
/*   /\* compute surface_integral[ n dot M*(q - q*) ] over boundary *\/ */
/*   /\* int f; *\/ */
  
/*   for (int f = 0; f < faces; f++){ */
/*     dgmath_apply_LIFT(dgmath_jit_dbase, &element_data->M_qstar_min_q_dot_n[f*face_nodes], dim, deg, f, vol_tmp); */
/*     /\* printf("element = %d, face = %d, dim = %d, x,y,z = %f,%f\n", element_data->id, f, d, element_data->xyz[0][0], element_data->xyz[1][0]); *\/ */
/*     /\* printf("element = %d, face = %d, x,y,z = %f,%f\n", element_data->id, f, element_data->xyz[0][0], element_data->xyz[1][0]); *\/ */
/*     /\* util_print_matrix(vol_tmp, volume_nodes, 1, "Mqstar_min_q", 0); *\/ */

/*     /\* printf("element = %d, face = %d, dim = %d, x,y,z = %f,%f\n", element_data->id, f, d); *\/ */
/*     /\* printf("element = %d, face = %d\n", element_data->id, f); *\/ */
/*     /\* util_print_matrix(vol_tmp, volume_nodes, 1, "Mqstar_min_q", 0); *\/ */
/*     /\* for (int i = 0; i < volume_nodes; i++){ *\/ */
/*       /\* printf("Mqstar_min_q_dot_n = %.20f, xyz[0] = %.20f, xyz[1] = %.20f\n", vol_tmp[i], element_data->xyz[0][i], element_data->xyz[1][i]); *\/ */
/*     /\* } *\/ */

    
/*     linalg_vec_axpy(1., vol_tmp, Au, volume_nodes); */
/*   } */

/*   /\* Au *= -1 because Au matrix is negative definite before this! *\/ */
/*   linalg_vec_scale(-1., Au, volume_nodes); */
  
/*   P4EST_FREE(vol_tmp); */
/*   P4EST_FREE(dq_Gauss); */
/*   P4EST_FREE(dq); */
/* } */

void curved_Gauss_poisson_compute_Au_elem(p4est_iter_volume_info_t* info, void* user_data)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data;

  curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase;

  
  int dim = (P4EST_DIM);
  /* int deg = element_data->deg; */
  int faces = 2*dim;
  /* int face_nodes_Gauss = dgmath_get_nodes(dim-1,element_data->deg_integ); */
  int face_nodes_Lobatto = dgmath_get_nodes(dim-1,element_data->deg);
  int volume_nodes_Gauss = dgmath_get_nodes(dim,element_data->deg_integ);
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,element_data->deg);

  double* Au = element_data->Au_elem;
  linalg_fill_vec(Au, 0., volume_nodes_Lobatto);
  
  double* dq = P4EST_ALLOC(double, volume_nodes_Lobatto);
  double* dq_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* dxq_Gauss = P4EST_ALLOC_ZERO(double, volume_nodes_Gauss);

  double* Mdq = P4EST_ALLOC(double, volume_nodes_Lobatto);
  double* lifted_qstar_min_q [P4EST_FACES];
  /* double* M_lifted_qstar_min_q [P4EST_FACES]; */

  for (int f = 0; f < (P4EST_FACES); f++){
    lifted_qstar_min_q[f] = P4EST_ALLOC(double, volume_nodes_Lobatto);
    /* M_lifted_qstar_min_q[f] = P4EST_ALLOC(double, volume_nodes_Lobatto); */
  }
  
  /* compute  */
  int i,j,k;
  for (i = 0; i < (P4EST_DIM); i++){
    for (j = 0; j < (P4EST_DIM); j++){
      dgmath_apply_Dij(dgmath_jit_dbase, &element_data->q_elem[i][0], dim, element_data->deg, j, dq);
      dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dq, element_data->deg, element_data->deg_integ, dq_Gauss, (P4EST_DIM));
      for (k = 0; k < volume_nodes_Gauss; k++){
        dxq_Gauss[k] += element_data->rst_xyz_integ[j][i][k]*dq_Gauss[k];
      }
    }    
  }

  /* DEBUG_PRINT_ARR_DBL(dxq_Gauss, volume_nodes_Gauss); */
  
  
  dgmath_apply_curvedGaussMass_onGaussNodeVec(
                                              dgmath_jit_dbase,
                                              dxq_Gauss,
                                              element_data->deg,
                                              element_data->J_integ,
                                              element_data->deg_integ,
                                              (P4EST_DIM),
                                              Mdq
                                             );



  /* DEBUG_PRINT_ARR_DBL(Mdq, volume_nodes_Lobatto); */

  /* double* Minv_Mqstar_min_q_dot_n = P4EST_ALLOC(double, face_nodes_Lobatto); */

  for (int f = 0; f < faces; f++){
    /* linalg_matvec_plus_vec(1., */
    /*                        element_data->invMface[f], */
    /*                        &element_data->M_qstar_min_q_dot_n[f*face_nodes_Lobatto], */
    /*                        0., */
    /*                        Minv_Mqstar_min_q_dot_n, */
    /*                        face_nodes_Lobatto, */
    /*                        face_nodes_Lobatto); */

    
    dgmath_apply_LIFT(
                      dgmath_jit_dbase,
                      &element_data->M_qstar_min_q_dot_n[f*face_nodes_Lobatto],
                      dim,
                      element_data->deg,
                      f,
                      lifted_qstar_min_q[f]);

    /* double* tmp = &element_data->M_qstar_min_q_dot_n[f*face_nodes_Lobatto]; */
    /* DEBUG_PRINT_ARR_DBL(tmp, face_nodes_Lobatto); */
    
    /* dgmath_apply_LIFT( */
    /*                   dgmath_jit_dbase, */
    /*                   Minv_Mqstar_min_q_dot_n, */
    /*                   dim, */
    /*                   element_data->deg, */
    /*                   f, */
    /*                   lifted_qstar_min_q[f]); */



    /* dgmath_apply_curvedGaussMass( */
    /*                              dgmath_jit_dbase, */
    /*                              lifted_qstar_min_q[f], */
    /*                              element_data->deg, */
    /*                              element_data->J_integ, */
    /*                              element_data->deg_integ, */
    /*                              (P4EST_DIM), */
    /*                              M_lifted_qstar_min_q[f] */
    /* );     */
  }

  /* P4EST_FREE(Minv_Mqstar_min_q_dot_n);   */
  
  for (int i = 0; i < volume_nodes_Lobatto; i++){
    Au[i] = Mdq[i];
    for (int f = 0; f < faces; f++){
      Au[i] += lifted_qstar_min_q[f][i];
      /* printf("lifted_qstar_min_q[f][i] = %.25f\n",lifted_qstar_min_q[f][i]); */
    }
  }

  /* for (int i = 0; i < volume_nodes_Lobatto; i++){ */
  /*   Au[i] = Mdq[i]; */
  /*   for (int f = 0; f < faces; f++){ */
  /*     Au[i] += M_lifted_qstar_min_q[f][i]; */
  /*   } */
  /* } */



  /* DEBUG_PRINT_ARR_DBL(Au_integ, volume_nodes_Gauss); */


  
  /* dgmath_apply_p_prolong_transpose(dgmath_jit_dbase, Au_integ, element_data->deg_integ, (P4EST_DIM), deg, Au); */

   
  /* Au *= -1 because Au matrix is negative definite before this! */
  linalg_vec_scale(-1., Au, volume_nodes_Lobatto);

  /* DEBUG_PRINT_ARR_DBL(Au, volume_nodes_Lobatto); */

  
#ifdef D4EST_DEBUG
  curved_poisson_debug_vecs_t* debug_vecs = (curved_poisson_debug_vecs_t*)curved_Gauss_poisson_user_data->debug_vecs;

  if (debug_vecs != NULL && debug_vecs->elem_id == element_data->id){
    curved_poisson_debug_vecs_set_Mdivq(Mdq, debug_vecs, dgmath_jit_dbase);
    curved_poisson_debug_vecs_set_liftedqflux(lifted_qstar_min_q, debug_vecs, dgmath_jit_dbase);
    curved_poisson_debug_vecs_set_Au(Au, debug_vecs, dgmath_jit_dbase);
    
  }
#endif

  
  for (int f = 0; f < (P4EST_FACES); f++){
    P4EST_FREE(lifted_qstar_min_q[f]);
    /* P4EST_FREE(M_lifted_qstar_min_q[f]); */
  }
  
  P4EST_FREE(dq_Gauss);
  P4EST_FREE(dxq_Gauss);
  P4EST_FREE(dq);
  P4EST_FREE(Mdq);
}

/* void curved_Gauss_poisson_compute_Au_elem_old(p4est_iter_volume_info_t* info, void* user_data) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data; */

/*   curved_Gauss_poisson_user_data_t* curved_Gauss_poisson_user_data = (curved_Gauss_poisson_user_data_t*) user_data; */
/*   dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) curved_Gauss_poisson_user_data->dgmath_jit_dbase; */
  
/*   int dim = (P4EST_DIM); */
/*   int deg = element_data->deg; */
/*   int faces = 2*dim; */
/*   int face_nodes = dgmath_get_nodes(dim-1,deg); */
/*   int volume_nodes = dgmath_get_nodes(dim,deg); */

/*   double* Au = element_data->Au_elem; */
/*   linalg_fill_vec(Au, 0., volume_nodes); */
  
/*   double* dq = P4EST_ALLOC(double, dgmath_get_nodes(dim,deg)); */
/*   double* dq_Gauss = P4EST_ALLOC(double, dgmath_get_nodes(dim,deg)); */
/*   double* vol_tmp = P4EST_ALLOC_ZERO(double, dgmath_get_nodes(dim,deg)); */

/*   /\* compute  *\/ */
/*   int i,j,k; */
/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     for (j = 0; j < (P4EST_DIM); j++){ */
/*       dgmath_apply_Dij(dgmath_jit_dbase, &element_data->q_elem[i][0], dim, deg, j, dq); */
/*       dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dq, deg, deg, dq_Gauss, (P4EST_DIM)); */
/*       for (k = 0; k < volume_nodes; k++){ */
/*         vol_tmp[k] += element_data->rst_xyz_integ[j][i][k]*dq_Gauss[k]; */
/*       } */
/*     }     */
/*   } */
  
/*   /\* linalg_component_mult(vol_tmp2, element_data->J, vol_tmp, volume_nodes); *\/ */
/*   /\* dgmath_apply_Mij(dgmath_jit_dbase, vol_tmp, dim, deg, Au); *\/ */

/*   dgmath_apply_curvedGaussMass_onGaussNodeVec(dgmath_jit_dbase, vol_tmp, deg, element_data->J_integ, deg, (P4EST_DIM), Au); */
  
/*   /\* compute surface_integral[ n dot M*(q - q*) ] over boundary *\/ */
/*   /\* int f; *\/ */
  
/*   for (int f = 0; f < faces; f++){ */
/*     dgmath_apply_LIFT(dgmath_jit_dbase, &element_data->M_qstar_min_q_dot_n[f*face_nodes], dim, deg, f, vol_tmp); */
/*     /\* printf("element = %d, face = %d, dim = %d, x,y,z = %f,%f\n", element_data->id, f, d, element_data->xyz[0][0], element_data->xyz[1][0]); *\/ */
/*     /\* printf("element = %d, face = %d, x,y,z = %f,%f\n", element_data->id, f, element_data->xyz[0][0], element_data->xyz[1][0]); *\/ */
/*     /\* util_print_matrix(vol_tmp, volume_nodes, 1, "Mqstar_min_q", 0); *\/ */

/*     /\* printf("element = %d, face = %d, dim = %d, x,y,z = %f,%f\n", element_data->id, f, d); *\/ */
/*     /\* printf("element = %d, face = %d\n", element_data->id, f); *\/ */
/*     /\* util_print_matrix(vol_tmp, volume_nodes, 1, "Mqstar_min_q", 0); *\/ */
/*     /\* for (int i = 0; i < volume_nodes; i++){ *\/ */
/*       /\* printf("Mqstar_min_q_dot_n = %.20f, xyz[0] = %.20f, xyz[1] = %.20f\n", vol_tmp[i], element_data->xyz[0][i], element_data->xyz[1][i]); *\/ */
/*     /\* } *\/ */

    
/*     linalg_vec_axpy(1., vol_tmp, Au, volume_nodes); */
/*   } */

/*   /\* Au *= -1 because Au matrix is negative definite before this! *\/ */
/*   linalg_vec_scale(-1., Au, volume_nodes); */
  
/*   P4EST_FREE(vol_tmp); */
/*   P4EST_FREE(dq_Gauss); */
/*   P4EST_FREE(dq); */
/* } */


void
curved_Gauss_poisson_destroy_vecs(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* element_data = (curved_element_data_t *) q->p.user_data;
  
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    //P4EST_FREE(element_data->q_elem[i]);
    P4EST_FREE(element_data->M_ustar_min_u_n[i]);
    /* P4EST_FREE(element_data->du_elem[i]); */
  }
  P4EST_FREE(element_data->M_qstar_min_q_dot_n);
}



void
curved_Gauss_poisson_apply_aij
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom
)
{
  curved_Gauss_poisson_user_data_t curved_Gauss_poisson_user_data;
  curved_Gauss_poisson_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_Gauss_poisson_user_data.problem_data = prob_vecs;
 #if D4EST_DEBUG
  curved_Gauss_poisson_user_data.debug_vecs = NULL;
#endif
  
  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_compute_flux_user_data.geom = geom;
  
  void* tmpptr = p4est->user_pointer;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &curved_Gauss_poisson_user_data,
		curved_Gauss_poisson_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data(p4est,ghost,ghost_data);
 
  curved_compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->curved_scalar_flux_fcn_data;
  p4est->user_pointer = &curved_compute_flux_user_data;
  
  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  /* DEBUG_PRINT_ARR_DBL(prob_vecs->u, prob_vecs->local_nodes); */
  
  
  p4est_iterate(p4est,
  		NULL,
  		(void*)&curved_Gauss_poisson_user_data,
  		curved_Gauss_poisson_compute_q_elem,
  		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data (p4est, ghost, ghost_data);
  curved_compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->curved_vector_flux_fcn_data;
  p4est->user_pointer = &curved_compute_flux_user_data;

  p4est_iterate (p4est,
  		 ghost,
  		 (void *) ghost_data,
  		 NULL,
                 curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&curved_Gauss_poisson_user_data,
  		 curved_Gauss_poisson_compute_Au_elem,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est_iterate(p4est,
		NULL,
  		(void*)&curved_Gauss_poisson_user_data,
		curved_Gauss_poisson_destroy_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  p4est->user_pointer = tmpptr;
}


curved_poisson_debug_vecs_t*
curved_Gauss_poisson_apply_aij_debug
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 problem_data_t* prob_vecs,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 int local_element_id
)
{
  curved_Gauss_poisson_user_data_t curved_Gauss_poisson_user_data;
  curved_Gauss_poisson_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_Gauss_poisson_user_data.problem_data = prob_vecs;


#ifdef D4EST_DEBUG
  curved_Gauss_poisson_user_data.debug_vecs = curved_poisson_debug_vecs_init
                                            (
                                             p4est,
                                             local_element_id
                                            ); 
  mpi_assert(curved_Gauss_poisson_user_data.debug_vecs != NULL);
#endif
  
  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_compute_flux_user_data.geom = geom;
  
  
  void* tmpptr = p4est->user_pointer;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &curved_Gauss_poisson_user_data,
		curved_Gauss_poisson_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data(p4est,ghost,ghost_data);
 
  curved_compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->curved_scalar_flux_fcn_data;
  p4est->user_pointer = &curved_compute_flux_user_data;
  
  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  /* DEBUG_PRINT_ARR_DBL(prob_vecs->u, prob_vecs->local_nodes); */
  
  
  p4est_iterate(p4est,
  		NULL,
  		(void*)&curved_Gauss_poisson_user_data,
  		curved_Gauss_poisson_compute_q_elem,
  		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data (p4est, ghost, ghost_data);
  curved_compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->curved_vector_flux_fcn_data;
  p4est->user_pointer = &curved_compute_flux_user_data;

  p4est_iterate (p4est,
  		 ghost,
  		 (void *) ghost_data,
  		 NULL,
                 curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&curved_Gauss_poisson_user_data,
  		 curved_Gauss_poisson_compute_Au_elem,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est_iterate(p4est,
		NULL,
  		(void*)&curved_Gauss_poisson_user_data,
		curved_Gauss_poisson_destroy_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  p4est->user_pointer = tmpptr;

#if D4EST_DEBUG
  return curved_Gauss_poisson_user_data.debug_vecs;
#else
  return NULL;
#endif
}
