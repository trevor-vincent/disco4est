#include "../pXest/pXest.h"
#include "../ElementData/element_data.h"
#include "../EllipticSystem/problem_data.h"
#include "../EllipticSystem/problem_weakeqn_ptrs.h"
#include "../LinearAlgebra/linalg.h"
#include "../Utilities/util.h"
#include "../hpAMR/hp_amr.h"
#include "../Estimators/bi_estimator_flux_fcns.h"
#include "../Estimators/bi_estimator.h"

/* #define D4EST_DEBUG */

typedef struct {
  
  dgmath_jit_dbase_t* dgmath_jit_dbase;
  problem_data_t* problem_data;
  double* bi_estimator_local_eta2;
  
} bi_estimator_user_data_t;


/* static double bi_estimator_local_eta2 = 0.; */

static void
bi_estimator_init
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;

  bi_estimator_user_data_t* bi_estimator_user_data = (bi_estimator_user_data_t*) user_data;
  problem_data_t* problem_data = bi_estimator_user_data->problem_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = bi_estimator_user_data->dgmath_jit_dbase;

  
  double h = element_data->h;
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int i;

  element_data->ustar_min_u = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*dgmath_get_nodes(dim-1, deg));
  
  for (i = 0; i < (P4EST_DIM); i++){
    element_data->qstar_min_q[i] = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*dgmath_get_nodes(dim-1, deg));
    element_data->du_elem[i] = P4EST_ALLOC_ZERO(double, dgmath_get_nodes(dim, deg));   
  }

  element_data->u_elem = &(problem_data->u[element_data->stride]);

  int volume_nodes = dgmath_get_nodes(dim,deg);

  linalg_copy_1st_to_2nd(
  element_data->u_elem,
    &(element_data->u_storage)[0],
    volume_nodes
    );

  for (i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_Dij(dgmath_jit_dbase, element_data->u_elem, dim, deg, i, element_data->du_elem[i]);
    linalg_vec_scale(2./h, element_data->du_elem[i], volume_nodes);
  }
}

static
void bi_estimator_compute_callback
(
p4est_iter_volume_info_t* info,
void* user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t*) q->p.user_data;

  bi_estimator_user_data_t* bi_estimator_user_data = (bi_estimator_user_data_t*) user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = bi_estimator_user_data->dgmath_jit_dbase;
  double* bi_estimator_local_eta2 = bi_estimator_user_data->bi_estimator_local_eta2;
  
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int faces = 2*dim;
  int face_nodes = dgmath_get_nodes(dim-1,deg);

  double surface_jacobian = element_data->surface_jacobian;
  double h = element_data->h;

  double* eta = &(element_data->local_estimator);

  /* handle ||R||^2 * h^2/p^2 term */
  *eta *= h*h/(deg*deg);

  double* Je1;
  double* Je2;
  double* MJe1 = P4EST_ALLOC(double,face_nodes);
  double* MJe2 = P4EST_ALLOC(double,face_nodes);

  double Nsqre1 = 0.;
  double Nsqre2 = 0.;

  int f,d;
  for (f = 0; f < faces; f++){
    /* calculate Ne1_sqr term */
    Je1 =  &(element_data->ustar_min_u[f*face_nodes]);
    dgmath_apply_Mij(dgmath_jit_dbase, Je1, dim - 1, deg, MJe1);
    linalg_vec_scale(surface_jacobian, MJe1, face_nodes);
    Nsqre1 += linalg_vec_dot(Je1, MJe1, face_nodes);

#ifndef NDEBUG
  /* printf("id, face, Je1MJe1 = %d,%d,%.25f\n", element_data->id, f, linalg_vec_dot(Je1, MJe1, face_nodes)); */
#endif
    
    /* calculate Ne2_sqr term */
    for (d = 0; d < (P4EST_DIM); d++){
      Je2 = &(element_data->qstar_min_q[d][f*face_nodes]);
      dgmath_apply_Mij(dgmath_jit_dbase, Je2, dim - 1, deg, MJe2);
      linalg_vec_scale(surface_jacobian, MJe2, face_nodes);
      Nsqre2 += linalg_vec_dot(Je2, MJe2, face_nodes);

#ifndef NDEBUG
  /* printf("id, face, Je2MJe2 = %d,%d,%.25f\n", element_data->id, f, linalg_vec_dot(Je2, MJe2, face_nodes)); */
#endif
      
    }

  }

  /* printf("Nsqre1 = %.20f\n", Nsqre1); */
  /* printf("Nsqre2 = %.20f\n", Nsqre2); */

  /* printf("eta ||R|| term = %.25f\n", *eta); */
  
  *eta += Nsqre1;
  *eta += Nsqre2;
  
  /* printf("eta = ||R|| term + Nsqre1 + Nsqre2 = %.25f\n", *eta); */
  
  /* take advantage of opportunity to compute max, sum and max error density */
  *bi_estimator_local_eta2 += (*eta);
  /* printf("eta =%.20f\n", *eta); */

  P4EST_FREE(MJe1);
  P4EST_FREE(MJe2);
}

static
void bi_estimator_destroy(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;

  P4EST_FREE(element_data->ustar_min_u);
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(element_data->qstar_min_q[i]);
    P4EST_FREE(element_data->du_elem[i]);
  }
}


static
void bi_estimator_set_eta_randomly
(
p4est_iter_volume_info_t* info,
void* user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t*) q->p.user_data;  
  double surface_jacobian = element_data->surface_jacobian;
  element_data->local_estimator = .001*info->quadid*surface_jacobian;  
}

double
bi_estimator_test_compute
(
 p4est_t* p4est,
 int level
)
{
  p4est_iterate(
                p4est,
		NULL,
		NULL,
		bi_estimator_set_eta_randomly,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL
               );
  
  return util_dbl_pow_int( .1, level);
}

void
bi_estimator_compute
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 penalty_calc_t u_penalty_fcn,
 penalty_calc_t u_dirichlet_penalty_fcn,
 penalty_calc_t gradu_penalty_fcn,
 grid_fcn_t u_bndry_fcn,
 double penalty_prefactor,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  double bi_estimator_local_eta2 = 0.;
    
  fcns->build_residual
    (
     p4est,
     ghost,
     ghost_data,
     vecs,
     dgmath_jit_dbase
    );
 
  element_data_compute_l2_norm_sqr(p4est,vecs->Au, dgmath_jit_dbase);

  /* /\* element_data_t* ghost *\/_data; */
  bi_estimator_user_data_t bi_estimator_user_data;
  bi_estimator_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  bi_estimator_user_data.problem_data = vecs;
  bi_estimator_user_data.bi_estimator_local_eta2 = &bi_estimator_local_eta2;

  compute_flux_user_data_t compute_flux_user_data;
  compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;

  /* subtract  off numerical solution*/
  void* tmpptr = p4est->user_pointer;

  p4est_iterate(p4est,
		NULL,
		(void *) &bi_estimator_user_data,
		bi_estimator_init,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);
  
  p4est_ghost_exchange_data(p4est,ghost,ghost_data);

  
  
  flux_fcn_ptrs_t bi_estimator_flux_fcn_ptrs = bi_est_dirichlet_fetch_fcns(u_bndry_fcn,
                                                                        u_penalty_fcn,
                                                                        u_dirichlet_penalty_fcn,
                                                                        gradu_penalty_fcn,
                                                                        penalty_prefactor
                                                                        );
  compute_flux_user_data.flux_fcn_ptrs = &bi_estimator_flux_fcn_ptrs;
  p4est->user_pointer = &compute_flux_user_data; 

  
  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est_iterate(p4est,
  		NULL,
  		(void*) &bi_estimator_user_data,
  		bi_estimator_compute_callback,
  		NULL,
  #if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);


  /* printf("bi_estimator_local_eta2 after = %f\n", bi_estimator_local_eta2); */
  p4est_iterate(p4est,
		NULL,
		(void *) &bi_estimator_user_data,
                bi_estimator_destroy,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est->user_pointer = tmpptr;

  /* estimator_stats_t stats; */
  /* estimator_stats_compute(p4est,&stats); */
  /* return stats; */
}
