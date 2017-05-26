#include <pXest.h>
#include <element_data.h>
#include <problem_data.h>
#include <d4est_operators.h>
#include <linalg.h>
#include <poisson_operator.h>
#include <grid_functions.h>
#include <util.h>
#include <poisson_debug_vecs.h>

typedef struct {
  
  d4est_operators_t* d4est_ops;
  problem_data_t* problem_data;

#if D4EST_DEBUG
  poisson_debug_vecs_t* debug_vecs;
#endif
  
} poisson_user_data_t;

static
void poisson_init_vecs(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;

  poisson_user_data_t* poisson_user_data = (poisson_user_data_t*) user_data;
  problem_data_t* problem_data = (problem_data_t*) poisson_user_data->problem_data;
  d4est_operators_t* d4est_ops = (d4est_operators_t*) poisson_user_data->d4est_ops;
  
  double h = element_data->h;
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int i;
  
  for (i = 0; i < (P4EST_DIM); i++){
    element_data->q_elem[i] = P4EST_ALLOC_ZERO(double, d4est_operators_get_nodes(dim, deg));
    element_data->qstar_min_q[i] = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*d4est_operators_get_nodes(dim-1, deg));
    element_data->du_elem[i] = P4EST_ALLOC_ZERO(double, d4est_operators_get_nodes(dim, deg));   
  }

  element_data->ustar_min_u = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*d4est_operators_get_nodes(dim-1,deg));
  element_data->u_elem = &(problem_data->u[element_data->stride]);
  element_data->Au_elem = &(problem_data->Au[element_data->stride]);

  linalg_copy_1st_to_2nd(element_data->u_elem, &(element_data->u_storage)[0], d4est_operators_get_nodes(dim, deg));
  
  for (i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_Dij(d4est_ops, element_data->u_elem, dim, deg, i, element_data->du_elem[i]);
    linalg_vec_scale(2./h, element_data->du_elem[i], d4est_operators_get_nodes(dim, deg));

    /* util_print_matrix */
    /*   ( */
    /*    element_data->du_elem[i], */
    /*    d4est_operators_get_nodes(dim, deg), */
    /*    1, */
    /*    "du_elem[i] = ", */
    /*    0 */
    /*   ); */
  }
}

static
void poisson_compute_q_elem(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t*) q->p.user_data;

  poisson_user_data_t* poisson_user_data = (poisson_user_data_t*) user_data;
  d4est_operators_t* d4est_ops = (d4est_operators_t*) poisson_user_data->d4est_ops;

  
  double jacobian = element_data->jacobian;
  double surface_jacobian = element_data->surface_jacobian;
  double h = element_data->h;
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int faces = 2*dim;
  int face_nodes = d4est_operators_get_nodes(dim-1,deg);
  int vol_nodes = d4est_operators_get_nodes(dim,deg);

  double* u_elem = element_data->u_elem;
  
  double* vol_tmp = P4EST_ALLOC(double, vol_nodes);
  double* Si_u [(P4EST_DIM)];
  double n [(P4EST_DIM)];

  /* double* ustar_minus_u = P4EST_ALLOC(double, faces*d4est_operators_get_nodes(dim-1,deg)); */
  double* M_ustar_minus_u = P4EST_ALLOC(double, faces*d4est_operators_get_nodes(dim-1,deg));

  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    Si_u[i] = P4EST_ALLOC(double, d4est_operators_get_nodes(dim,deg));
  }

  /* compute vol_temp2[i] = S_i * u */
  for (i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_Dij(d4est_ops, u_elem, dim, deg, i, vol_tmp);
    linalg_vec_scale(2./h, vol_tmp, vol_nodes);
    d4est_operators_apply_Mij(d4est_ops, vol_tmp, dim, deg, Si_u[i]);
    linalg_vec_scale(jacobian, Si_u[i], vol_nodes);
/* #ifndef NDEBUG */

    /* util_print_matrix(Si_u[i], vol_nodes, 1, "Si_u[i] = ", 0); */
/* #endif */
  }

  /* double* u_ptr = &element_data->u_elem[0]; */
  /* DEBUG_PRINT_4ARR_DBL(u_elem, Si_u[0], Si_u[1], Si_u[2], vol_nodes); */
  /* printf("jacobian = %.25f\n",jacobian); */

#ifdef D4EST_DEBUG
  double* lifted_ustar_min_u [(P4EST_FACES)][(P4EST_DIM)];
  poisson_debug_vecs_t* debug_vecs
    = (poisson_debug_vecs_t*)poisson_user_data->debug_vecs;
  if(debug_vecs != NULL && element_data->id == debug_vecs->elem_id){
    poisson_debug_vecs_set_Mdu(Si_u, debug_vecs, d4est_ops);
  }
    for (int i = 0; i < (P4EST_FACES); i++)
      for (int j = 0; j < (P4EST_DIM); j++){
        lifted_ustar_min_u[i][j] = P4EST_ALLOC_ZERO(double, vol_nodes);
      }
#endif  

  /* compute M*(u - u*) on boundary */
  int f,d;
  for (d = 0; d < (P4EST_DIM); d++){
  for (f = 0; f < faces; f++){
    /* d4est_operators_apply_slicer(u_elem, dim, f, deg, &ustar_minus_u[f*face_nodes]); */

    /* printf("\nl = %d\n", f); */
    /* util_print_matrix( &u_flux[f*face_nodes], face_nodes, 1,"u_flux = ",0); */
    /* util_print_matrix(&element_data->ustar_min_u[f*face_nodes], face_nodes, 1,"ustar_min_u = ",0); */
    
    /* linalg_vec_xpby(&u_flux[f*face_nodes], -1., &ustar_minus_u[f*face_nodes], face_nodes); */
    d4est_operators_apply_Mij(d4est_ops, &element_data->ustar_min_u[f*face_nodes], dim - 1, deg, &M_ustar_minus_u[f*face_nodes]);
    linalg_vec_scale( surface_jacobian, &M_ustar_minus_u[f*face_nodes], face_nodes);
 
    d4est_operators_apply_LIFT(d4est_ops, &M_ustar_minus_u[f*face_nodes], dim, deg, f, vol_tmp);

/* #ifndef NDEBUG */
    /* printf("f = %d\n",f); */
    /* util_print_matrix(vol_tmp, vol_nodes, 1, "Mustar_min_u lifted = ", 0); */
/* #endif */
    
    /* lift face integral and dot with normal */
    d4est_operators_get_normal(f, dim, &n[0]);

#ifdef D4EST_DEBUG
    linalg_vec_axpy(n[d], vol_tmp, lifted_ustar_min_u[f][d], vol_nodes);
#endif
    
    linalg_vec_axpy(n[d], vol_tmp, Si_u[d], vol_nodes);
  }
  /* util_print_matrix(Si_u[d], vol_nodes, 1, "Sd_u[d] = ", 0); */
  }

/* #ifndef NDEBUG */
  /* util_print_matrix(Si_u[d], nodes, 1, "-Mq = ", 0); */
/* #endif */


  
  for (d = 0; d < (P4EST_DIM); d++){
    d4est_operators_apply_invMij(d4est_ops, Si_u[d], dim, deg, element_data->q_elem[d]);
    linalg_vec_scale(1./jacobian, element_data->q_elem[d], vol_nodes);
/* #ifndef NDEBUG */
/* #endif */

    /* util_print_matrix(element_data->q_elem[d], vol_nodes, 1, "q = ", 0); */
  }


  /* DEBUG_PRINT_3ARR_DBL( */
  /*                      (&element_data->q_elem[0][0]), */
  /*                      (&element_data->q_elem[1][0]), */
  /*                      (&element_data->q_elem[2][0]), */
  /*                      vol_nodes */
  /* ); */


  /* DEBUG_PRINT_3ARR_DBL( */
                       /* (&Si_u[0][0]), */
                       /* (&Si_u[1][0]), */
                       /* (&Si_u[2][0]), */
                       /* vol_nodes */
  /* ); */
  
#ifdef D4EST_DEBUG
  if (debug_vecs != NULL && element_data->id == debug_vecs->elem_id){
    poisson_debug_vecs_set_lifteduflux(lifted_ustar_min_u, debug_vecs, d4est_ops);
    poisson_debug_vecs_set_q(element_data->q_elem, debug_vecs, d4est_ops);
    poisson_debug_vecs_set_Mq(Si_u, debug_vecs, d4est_ops);    
  }
  for (int i = 0; i < (P4EST_FACES); i++){
    for (int j = 0; j < (P4EST_DIM); j++){
      P4EST_FREE(lifted_ustar_min_u[i][j]);
    }
  }
#endif

  
  for (i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(Si_u[i]);
  }
  P4EST_FREE(vol_tmp);
  /* P4EST_FREE(ustar_minus_u); */
  P4EST_FREE(M_ustar_minus_u);

  /* linalg_fill_vec(&element_data->q_storage[0], 0., MAX_NODES); */
  /* linalg_copy_1st_to_2nd(element_data->q_elem, &(element_data->q_storage[0]), d4est_operators_get_nodes(dim, deg)); */
}


static
void poisson_compute_Au_elem(p4est_iter_volume_info_t* info, void* user_data)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t*) q->p.user_data;

  poisson_user_data_t* poisson_user_data = (poisson_user_data_t*) user_data;
  d4est_operators_t* d4est_ops = (d4est_operators_t*) poisson_user_data->d4est_ops;
  
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int faces = 2*dim;
  int face_nodes = d4est_operators_get_nodes(dim-1,deg);
  int vol_nodes = d4est_operators_get_nodes(dim,deg);

  double jacobian = element_data->jacobian;
  double surface_jacobian = element_data->surface_jacobian;
  double h = element_data->h;
  
  double* Au = element_data->Au_elem;
  linalg_fill_vec(Au, 0., vol_nodes);
  
  double* vol_tmp = P4EST_ALLOC(double, d4est_operators_get_nodes(dim,deg));
  double* vol_tmp2 = P4EST_ALLOC(double, d4est_operators_get_nodes(dim,deg));

  double n [(P4EST_DIM)];

  double* qstar_minus_qi = P4EST_ALLOC(double, faces*d4est_operators_get_nodes(dim-1,deg));
  double* M_qstar_minus_qi = P4EST_ALLOC(double, faces*d4est_operators_get_nodes(dim-1,deg));

#if D4EST_DEBUG
  poisson_debug_vecs_t* debug_vecs = (poisson_debug_vecs_t*)poisson_user_data->debug_vecs;
#endif
  
  /* compute vol_temp2[i] = S_i * q_i */
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    d4est_operators_apply_Dij(d4est_ops, element_data->q_elem[i], dim, deg, i, vol_tmp);
    linalg_vec_scale(2./h, vol_tmp, vol_nodes);

#ifdef D4EST_DEBUG
    if (debug_vecs != NULL && element_data->id == debug_vecs->elem_id){
      printf("2/h = %.20f\n", 2./h);
      DEBUG_PRINT_ARR_DBL(vol_tmp, vol_nodes);
    }
#endif
    
    d4est_operators_apply_Mij(d4est_ops, vol_tmp, dim, deg, vol_tmp2);
    linalg_vec_scale(jacobian, vol_tmp2, vol_nodes);
    linalg_vec_axpy(1.0, vol_tmp2, Au, vol_nodes);
  }

#ifdef D4EST_DEBUG
  double* lifted_qstar_min_q [P4EST_FACES];
  if (debug_vecs != NULL && element_data->id == debug_vecs->elem_id){
    poisson_debug_vecs_set_Mdivq(Au, debug_vecs, d4est_ops);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    lifted_qstar_min_q[f] = P4EST_ALLOC_ZERO(double, vol_nodes);
  }
#endif
  
/* #ifndef NDEBUG */
  /* util_print_matrix(Au, vol_nodes, 1, "Au = ", 0); */
/* #endif */
  
  /* compute surface_quadral[ n dot M*(q - q*) ] over boundary */
  int d,f;
  
  for (d = 0; d < (P4EST_DIM); d++){
    for (f = 0; f < faces; f++){
       /* d4est_operators_apply_slicer(element_data->q_elem[d], dim, f, deg, &qstar_minus_qi[f*face_nodes]); */
       /* linalg_vec_xpby( &(element_data->q_flux[d][f*face_nodes]), -1., &qstar_minus_qi[f*face_nodes], face_nodes); */
       d4est_operators_apply_Mij(d4est_ops, &element_data->qstar_min_q[d][f*face_nodes], dim - 1, deg, &M_qstar_minus_qi[f*face_nodes]);

       d4est_operators_get_normal(f, dim, &n[0]);
       linalg_vec_scale(n[d]*surface_jacobian, &M_qstar_minus_qi[f*face_nodes], face_nodes);

      /* lift face integral and dot with normal */
      d4est_operators_apply_LIFT(d4est_ops, &M_qstar_minus_qi[f*face_nodes], dim, deg, f, vol_tmp);

#ifdef D4EST_DEBUG
      linalg_vec_axpy(1., vol_tmp, lifted_qstar_min_q[f], vol_nodes);
#endif
      
      /* util_print_matrix(vol_tmp, vol_nodes, 1, "Mqstar_min_q lifted = ", 0); */
      linalg_vec_axpy(1., vol_tmp, Au, vol_nodes);
    }
  }

  /* Au *= -1 because matrix is negative definite! */
  linalg_vec_scale(-1., Au, vol_nodes);

  /* DEBUG_PRINT_ARR_DBL(Au, vol_nodes); */
  
#ifdef D4EST_DEBUG
  if (debug_vecs != NULL && element_data->id == debug_vecs->elem_id){
    poisson_debug_vecs_set_liftedqflux(lifted_qstar_min_q, debug_vecs, d4est_ops);
    poisson_debug_vecs_set_Au(Au, debug_vecs, d4est_ops);
  }
  for (int f = 0; f < (P4EST_FACES); f++){
    P4EST_FREE(lifted_qstar_min_q[f]);
  }
#endif
  /* util_print_matrix(element_data->u_elem, vol_nodes, 1, "u = ", 0); */
  /* for (d = 0; d < (P4EST_DIM); d++){ */
    /* util_print_matrix(element_data->q_elem[d], vol_nodes, 1, "q[d] = ", 0); */
  /* } */
  /* util_print_matrix(Au, vol_nodes, 1, "Au = ", 0); */
  
  P4EST_FREE(vol_tmp);
  P4EST_FREE(vol_tmp2);
  P4EST_FREE(qstar_minus_qi);
  P4EST_FREE(M_qstar_minus_qi);
}

static
void poisson_destroy_vecs(p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;
  
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(element_data->q_elem[i]);
    P4EST_FREE(element_data->qstar_min_q[i]);
    P4EST_FREE(element_data->du_elem[i]);
  }
  /* P4EST_FREE(element_data->u_storage); */
  P4EST_FREE(element_data->ustar_min_u);
}


/* static */
/* void poisson_build_rhs_with_strongbc_callback( */
/*                               p4est_iter_volume_info_t * info, */
/*                               void *user_data) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   element_data_t* elem_data = (element_data_t*) q->p.user_data; */

/*   poisson_user_data_t* poisson_user_data = (poisson_user_data_t*) user_data; */
/*   problem_data_t* prob_vecs = (problem_data_t*) poisson_user_data->problem_data; */
/*   d4est_operators_t* d4est_ops = (d4est_operators_t*) poisson_user_data->d4est_ops; */
  
/*   double* f_elem; */
/*   double* rhs_elem; */
/*   int stride = elem_data->stride; */
/*   int deg = elem_data->deg; */
/*   int dim = (P4EST_DIM); */

/*   f_elem = &(prob_vecs->f[stride]); */
/*   rhs_elem = &(prob_vecs->rhs[stride]); */

/*   d4est_operators_apply_Mij(d4est_ops, f_elem, dim, deg, rhs_elem); */

/*   /\* printf("stride = %d\n",stride); *\/ */
/*   /\* util_print_matrix(f_elem, d4est_operators_get_nodes(dim, deg), 1, "f_elem = ",0); *\/ */
/*   /\* util_print_matrix(rhs_elem, d4est_operators_get_nodes(dim, deg), 1, "f_elem = ",0); *\/ */
  
/*   /\* We will be multiplying Au by -1 since its negative definite so we must do the */
/*    * the same for the rhs vector *\/ */
/*   linalg_vec_scale(-1.*elem_data->jacobian, rhs_elem, d4est_operators_get_nodes(dim, deg)); */
/* } */

/* void */
/* poisson_build_rhs_with_strongbc */
/* ( */
/*  p4est_t* p4est, */
/*  p4est_ghost_t* ghost, */
/*  element_data_t* ghost_data, */
/*  problem_data_t* prob_vecs, */
/*  d4est_operators_t* d4est_ops */
/* ) */
/* { */
/*    */

/*   poisson_user_data_t poisson_user_data; */
/*   poisson_user_data.d4est_ops = d4est_ops; */
/*   poisson_user_data.problem_data = prob_vecs; */
  
/*   p4est_iterate(p4est, */
/* 		NULL, */
/* 		(void *)&poisson_user_data, */
/* 		poisson_build_rhs_with_strongbc_callback, */
/* 		NULL, */
/* #if (P4EST_DIM)==3 */
/*                  NULL, */
/* #endif */
/* 		NULL); */

/*   int local_nodes = prob_vecs->local_nodes; */
  
/*   double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes); */
/*   double* tmp = prob_vecs->u; */
/*   prob_vecs->u = u_eq_0;  */

/*   poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops); */
/*   linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes); */
  
/*   prob_vecs->u = tmp; */

/*   /\* zero out the boundary (strong enforcement) *\/ */
/*   /\* prob_vecs->u_at_bndry = zero_fcn; *\/ */
/*   /\* prob_vecs->scalar_flux_fcn_data.bndry_fcn = zero_fcn; *\/ */
/*   /\* prob_vecs->vector_flux_fcn_data.bndry_fcn = zero_fcn; *\/ */
  
/*   P4EST_FREE(u_eq_0); */
/*    */
/* } */

void poisson_apply_aij(
                         p4est_t* p4est,
                         p4est_ghost_t* ghost,
                         void* ghost_data,
                         problem_data_t* prob_vecs,
			 d4est_operators_t* d4est_ops,
                         d4est_geometry_t* d4est_geom
                         )
{
  /* p4est_ghost_t* ghost; */
  /* /\* element_data_t* ghost *\/_data; */
  poisson_user_data_t poisson_user_data;
  poisson_user_data.d4est_ops = d4est_ops;
  poisson_user_data.problem_data = prob_vecs;
 #if D4EST_DEBUG
  poisson_user_data.debug_vecs = NULL;
#endif
  
  compute_flux_user_data_t compute_flux_user_data;
  compute_flux_user_data.d4est_ops = d4est_ops;
  
  void* tmpptr = p4est->user_pointer;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &poisson_user_data,
		poisson_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data(p4est,ghost,ghost_data);
 
  compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->scalar_flux_fcn_data;
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
  		(void*)&poisson_user_data,
  		poisson_compute_q_elem,
  		NULL,
  #if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  

  
  /* P4EST_FREE(ghost_data); */
  /* ghost_data = P4EST_ALLOC(element_data_t, ghost->ghosts.elem_count); */
  p4est_ghost_exchange_data (p4est, ghost, ghost_data);
  
  
  compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->vector_flux_fcn_data;
  p4est->user_pointer = &compute_flux_user_data; 

  
  p4est_iterate (p4est,
  		 ghost,
  		 (void *) ghost_data,
  		 NULL,
                 compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  
  
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&poisson_user_data,
  		 poisson_compute_Au_elem,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  

  
  p4est_iterate(p4est,
		NULL,
  		(void*)&poisson_user_data,
		poisson_destroy_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  
  p4est->user_pointer = tmpptr;


  /* P4EST_FREE(ghost_data); */
  /* p4est_ghost_destroy(ghost); */
  /* printf("DONE POISSON CALL\n"); */
}


poisson_debug_vecs_t*
poisson_apply_aij_debug(
                         p4est_t* p4est,
                         p4est_ghost_t* ghost,
                         element_data_t* ghost_data,
                         problem_data_t* prob_vecs,
			 d4est_operators_t* d4est_ops,
                         int local_element_id
                         )
{
  /* p4est_ghost_t* ghost; */
  /* /\* element_data_t* ghost *\/_data; */
  poisson_user_data_t poisson_user_data;
  poisson_user_data.d4est_ops = d4est_ops;
  poisson_user_data.problem_data = prob_vecs;
#ifdef D4EST_DEBUG
  poisson_user_data.debug_vecs = poisson_debug_vecs_init
                                 (
                                  p4est,
                                  local_element_id
                                 );
  mpi_assert(poisson_user_data.debug_vecs != NULL);
#endif
  compute_flux_user_data_t compute_flux_user_data;
  compute_flux_user_data.d4est_ops = d4est_ops;

  
  void* tmpptr = p4est->user_pointer;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &poisson_user_data,
		poisson_init_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data(p4est,ghost,ghost_data);
 
  compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->scalar_flux_fcn_data;
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
  		(void*)&poisson_user_data,
  		poisson_compute_q_elem,
  		NULL,
  #if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  p4est_ghost_exchange_data (p4est, ghost, ghost_data);
  
  compute_flux_user_data.flux_fcn_ptrs = &prob_vecs->vector_flux_fcn_data;
  p4est->user_pointer = &compute_flux_user_data; 

  p4est_iterate (p4est,
  		 ghost,
  		 (void *) ghost_data,
  		 NULL,
                 compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  p4est_iterate (p4est,
  		 NULL,
                 (void*)&poisson_user_data,
  		 poisson_compute_Au_elem,
  		 NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  p4est_iterate(p4est,
		NULL,
  		(void*)&poisson_user_data,
		poisson_destroy_vecs,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);
  
  p4est->user_pointer = tmpptr;
#ifdef D4EST_DEBUG
  return poisson_user_data.debug_vecs;
#else
  return NULL;
#endif
  /* P4EST_FREE(ghost_data); */
  /* p4est_ghost_destroy(ghost); */
  /* printf("DONE POISSON CALL\n"); */
}

