#include "../pXest/pXest.h"
#include "../ElementData/element_data.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../hpAMR/hp_amr.h"
#include "../Estimators/curved_bi_estimator_flux_fcns.h"
#include "../Estimators/curved_bi_estimator.h"

typedef struct {
  
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  problem_data_t* problem_data;
  /* double* curved_bi_estimator_local_eta2; */
  
} curved_bi_estimator_user_data_t;

static double debug_eta2_after_res_calc;

/* static double curved_bi_estimator_local_eta2 = 0.; */

static void
curved_bi_estimator_init
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;
  curved_bi_estimator_user_data_t* curved_bi_estimator_user_data = (curved_bi_estimator_user_data_t*) user_data;
  problem_data_t* problem_data = curved_bi_estimator_user_data->problem_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = curved_bi_estimator_user_data->dgmath_jit_dbase;
  
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int volume_nodes_Lobatto = dgmath_get_nodes(dim,deg);
  int face_nodes_Lobatto = dgmath_get_nodes(dim-1,deg);
  int volume_nodes_Gauss = dgmath_get_nodes(dim, elem_data->deg_integ);
  /* int face_nodes_Gauss = dgmath_get_nodes(dim-1, elem_data->deg_integ); */
  
  linalg_copy_1st_to_2nd
    (
     &(problem_data->u[elem_data->nodal_stride]),
     &(elem_data->u_storage)[0],
     volume_nodes_Lobatto
    );
  
  double* u_prolonged = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* du_di_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  
  for (int i = 0; i < (P4EST_DIM); i++){

    /* PROLONG u here first, then take derivative, then interpolate to Gauss, this is most consistent */
    dgmath_apply_p_prolong(dgmath_jit_dbase, &(problem_data->u[elem_data->nodal_stride]), deg, (P4EST_DIM), elem_data->deg_integ, u_prolonged);
    
    dgmath_apply_Dij(dgmath_jit_dbase, u_prolonged, dim, elem_data->deg_integ, i, &elem_data->dudr_elem[i][0]);
  /*   dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &elem_data->dudr_elem[i][0], elem_data->deg_integ, elem_data->deg_integ, du_di_Gauss, (P4EST_DIM)); */
    
  /*   for (int j = 0; j < (P4EST_DIM); j++){ */
  /*     for (int k = 0; k < volume_nodes_Gauss; k++){ */
  /*       elem_data->du_elem[j][k] += elem_data->rst_xyz_integ[i][j][k]*du_di_Gauss[k]; */
  /*     } */
  /*   }     */
  }

  /* double* tmp_ptr = &elem_data->dudr_elem[0][0]; */
  /* DEBUG_PRINT_ARR_DBL(tmp_ptr, dgmath_get_nodes((P4EST_DIM),elem_data->deg_integ)); */

  
  P4EST_FREE(u_prolonged);
  P4EST_FREE(du_di_Gauss);
}

static
void curved_bi_estimator_compute_volume_term_callback
(
p4est_iter_volume_info_t* info,
void* user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* element_data = (curved_element_data_t*) q->p.user_data;

  int deg = element_data->deg;
  double* eta2 = &(element_data->local_estimator);

  /* handle ||R||^2 * h^2/p^2 term */
  double h = element_data->diam;
  *eta2 *= h*h/(deg*deg);


  /* DEBUG_PRINT_DBL(*eta2); */
  
  debug_eta2_after_res_calc += *eta2;
}

/* static */
/* void curved_bi_estimator_compute_local_eta_callback */
/* ( */
/*  p4est_iter_volume_info_t * info, */
/*  void *user_data */
/* ) */
/* { */
/*   p4est_quadrant_t *q = info->quad; */
/*   curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data; */
/*   curved_bi_estimator_user_data_t* curved_bi_estimator_user_data = (curved_bi_estimator_user_data_t*) user_data; */

/*   *(curved_bi_estimator_user_data->curved_bi_estimator_local_eta2) += elem_data->local_estimator; */
/* } */

void
curved_bi_estimator_compute
(
 p4est_t* p4est,
 problem_data_t* vecs,
 curved_weakeqn_ptrs_t* fcns,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 grid_fcn_t u_bndry_fcn,
 double penalty_prefactor,
 p4est_ghost_t* ghost,
 curved_element_data_t* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom
)
{
  /* double curved_bi_estimator_local_eta2 = 0.; */
    
  fcns->build_residual
    (
     p4est,
     ghost,
     ghost_data,
     vecs,
     dgmath_jit_dbase,
     geom
    );


  
  /* DEBUG_PRINT_ARR_DBL(vecs->Au, vecs->local_nodes); */
 /* x */
  /* double check = */
  curved_element_data_compute_l2_norm_sqr
    (
     p4est,
     vecs->Au,
     vecs->local_nodes,
     dgmath_jit_dbase,
     STORE_LOCALLY
    );
  curved_element_data_print_local_estimator(p4est);

  curved_bi_estimator_user_data_t curved_bi_estimator_user_data;
  curved_bi_estimator_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_bi_estimator_user_data.problem_data = vecs;
  /* curved_bi_estimator_user_data.curved_bi_estimator_local_eta2 = &curved_bi_estimator_local_eta2; */

  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_compute_flux_user_data.geom = geom;

  debug_eta2_after_res_calc = 0;
  
  p4est_iterate
    (
     p4est,
     NULL,
     (void*) &curved_bi_estimator_user_data,
     curved_bi_estimator_compute_volume_term_callback,
     NULL,
  #if (P4EST_DIM)==3
     NULL,
#endif
     NULL
    );

  /* printf("local eta2 after residual calc = %.25f\n",  debug_eta2_after_res_calc); */
  
  /* subtract  off numerical solution*/
  void* tmpptr = p4est->user_pointer;

  p4est_iterate
    (
     p4est,
     NULL,
     (void *) &curved_bi_estimator_user_data,
     curved_bi_estimator_init,
     NULL,
#if (P4EST_DIM)==3
     NULL,
#endif
     NULL
    );
  
  p4est_ghost_exchange_data
    (
     p4est,
     ghost,
     ghost_data
    );

  curved_flux_fcn_ptrs_t curved_bi_estimator_flux_fcn_ptrs = curved_bi_est_dirichlet_fetch_fcns
                                                      (
                                                       u_bndry_fcn,
                                                       u_penalty_fcn,
                                                       u_dirichlet_penalty_fcn,
                                                       gradu_penalty_fcn,
                                                       penalty_prefactor
                                                      );
  
  curved_compute_flux_user_data.flux_fcn_ptrs = &curved_bi_estimator_flux_fcn_ptrs;
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

  /* printf("curved_bi_estimator_local_eta2 after = %f\n", curved_bi_estimator_local_eta2); */
/*   p4est_iterate(p4est, */
/* 		NULL, */
/* 		(void *) &curved_bi_estimator_user_data, */
/*                 curved_bi_estimator_compute_local_eta_callback, */
/* 		NULL, */
/* #if (P4EST_DIM)==3 */
/*                 NULL, */
/* #endif */
/* 		NULL); */
  
  p4est->user_pointer = tmpptr;

  /* P4EST_FREE(estimator_bi_res);   */
  /* p4est_ghost_destroy(ghost); */
  /* P4EST_FREE(ghost_data); */

  /* printf("local eta2 = %.25f\n",curved_bi_estimator_local_eta2); */
  curved_element_data_print_local_estimator(p4est);

  /* return curved_bi_estimator_local_eta2; */
}
