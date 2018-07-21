#include <pXest.h>

static void
testd4est_mixed_element_data_init_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t *) q->p.user_data;
  p4est_connectivity_t* pconnect = info->p4est->connectivity;
  
  p4est_qcoord_t dq;

  int* testd4est_mixed_element_data_stride = ((int*)info->p4est->user_pointer);
  int set_deg = *((int*)user_data);
  
  if ( set_deg >= 0)
    elem_data->deg = set_deg;

  elem_data->id = info->quadid;
  elem_data->stride = *testd4est_mixed_element_data_stride;

  dq = P4EST_QUADRANT_LEN(q->level);
  p4est_qcoord_to_vertex(pconnect,
                         which_tree,
                         q->x,
                         q->y,
#if (P4EST_DIM)==3
                         q->z,
#endif
                         elem_data->xyz_corner);

  elem_data->h = dq/(double)P4EST_ROOT_LEN;
  elem_data->jacobian = util_dbl_pow_int(.5*elem_data->h, (P4EST_DIM) ); 
  elem_data->surface_jacobian = util_dbl_pow_int(.5*elem_data->h, (P4EST_DIM) - 1);
  *testd4est_mixed_element_data_stride += dgmath_get_nodes( (P4EST_DIM) , elem_data->deg);

}

typedef struct {

  grid_fcn_t init_fcn;
  double* vec;

} init_node_vec_user_data_t;

static void
testd4est_mixed_element_data_init_node_vec_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;

  init_node_vec_user_data_t* inv_user_data = (init_node_vec_user_data_t*) user_data;
  grid_fcn_t init_fcn = inv_user_data->init_fcn;
  double* vec = inv_user_data->vec;
    
  dgmath_jit_dbase_t* dgbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  
  int stride = elem_data->stride;
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );
  double h = elem_data->h;
  
  int i;
  double* x = dgmath_fetch_xyz_nd( dgbase, (P4EST_DIM),
                                   elem_data->deg,
                                   0);
  double xl = elem_data->xyz_corner[0];
  
  double* y = dgmath_fetch_xyz_nd( dgbase, (P4EST_DIM),
                                   elem_data->deg,
                                   1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd( dgbase, (P4EST_DIM),
                                   elem_data->deg,
                                   2);
  double zl = elem_data->xyz_corner[2];
#endif

  for (i = 0; i < volume_nodes; i++){
    vec[stride + i] = init_fcn
                      (
                       dgmath_rtox(x[i], xl, h),
                       dgmath_rtox(y[i], yl, h)
#if (P4EST_DIM)==3
                       ,
                       dgmath_rtox(z[i], zl, h)
#endif
                      );

#ifndef NDEBUG
    printf(" vec[stride + i], x, y, z = %f,%f,%f,%f",
          vec[stride+i],
          dgmath_rtox(x[i], xl, h),
          dgmath_rtox(y[i], yl, h),
#if (P4EST_DIM)==3
          dgmath_rtox(z[i], zl, h)
#else
          1.
#endif
         );
#endif
      }
}

static void
testd4est_mixed_element_data_store_local_estimator_in_corner_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  double* corner_vec = (double*)user_data;
  
  double est = elem_data->local_estimator;
  int quadid = elem_data->id;

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    corner_vec[quadid*(P4EST_CHILDREN) + i] = est;
  }
}

static void
testd4est_mixed_element_data_store_nodal_vec_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;

  double* nodal_vec = (double*)user_data;
  double* vertex_vec = (double*) info->p4est->user_pointer;
  int id = elem_data->id;
  int stride = elem_data->stride;
  
  int c;
  for (c = 0; c < (P4EST_CHILDREN); c++){
    int corner_to_node = dgmath_corner_to_node((P4EST_DIM), elem_data->deg, c);
    vertex_vec[id*(P4EST_CHILDREN) + c] = nodal_vec[stride + corner_to_node];
  }
}

void
testd4est_mixed_element_data_init_node_vec
(
 p4est_t* p4est,
 double* nodal_vec,
 grid_fcn_t init_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  init_node_vec_user_data_t inv_user_data;
  inv_user_data.vec = nodal_vec;
  inv_user_data.init_fcn = init_fcn;
  
  p4est_iterate(p4est,
		NULL,
		(void *) &inv_user_data,
		testd4est_mixed_element_data_init_node_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  /* reset pointer */
  p4est->user_pointer = tmp;
}

static void
testd4est_mixed_element_data_init_quadid_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t *) q->p.user_data;
  elem_data->id = info->quadid;
}

void
testd4est_mixed_element_data_init_quadid
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		testd4est_mixed_element_data_init_quadid_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
}

typedef struct {

  FILE* f;
  double* vec;


} print_node_vec_user_data_t;


static void
testd4est_mixed_element_data_print_node_vec_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)

{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  print_node_vec_user_data_t* pnv_user_data = (print_node_vec_user_data_t*)user_data;
  double* vec = pnv_user_data->vec;
  FILE* f = pnv_user_data->f; 
  int stride = elem_data->stride;
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );
  double h = elem_data->h;
  
  int i;
  double* x = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   elem_data->deg,
                                   0);
  double xl = elem_data->xyz_corner[0];
  
  double* y = dgmath_fetch_xyz_nd( dgmath_jit_dbase,(P4EST_DIM),
                                   elem_data->deg,
                                   1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd( dgmath_jit_dbase,(P4EST_DIM),
                                   elem_data->deg,
                                   2);
  double zl = elem_data->xyz_corner[2];
#endif

  
  if (!user_data){
    for (i = 0; i < volume_nodes; i++){
#if (P4EST_DIM)==3
      printf("%.15f, %.15f, %.15f | %.15f\n",
             dgmath_rtox(x[i], xl, h),
             dgmath_rtox(y[i], yl, h),
             dgmath_rtox(z[i], zl, h),
             vec[stride + i]
            );
#else
      printf("%.15f, %.15f | %.15f\n",
             dgmath_rtox(x[i], xl, h),
             dgmath_rtox(y[i], yl, h),
             vec[stride + i]
            );
#endif
    }
  }
  else {
    /* FILE * f = (FILE*)user_data; */
    for (i = 0; i < volume_nodes; i++){
#if (P4EST_DIM)==3
      fprintf(f, "%.15f, %.15f, %.15f, %.15f\n",
             dgmath_rtox(x[i], xl, h),
             dgmath_rtox(y[i], yl, h),
             dgmath_rtox(z[i], zl, h),
             vec[stride + i]
            );
#else
      fprintf(f, "%.15f, %.15f , %.15f\n",
             dgmath_rtox(x[i], xl, h),
             dgmath_rtox(y[i], yl, h),
             vec[stride + i]
            );
#endif
    }
  }
  
}

void
testd4est_mixed_element_data_print_node_vec
(
 p4est_t* p4est,
 double* nodal_vec,
 char* name_string,
 int mpi_rank,
 int save_to_file,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  print_node_vec_user_data_t pnv_user_data;
  pnv_user_data.vec = nodal_vec;
  pnv_user_data.f = NULL;
  
  if(save_to_file == 1){
    char save_as [200];
    sprintf(&save_as[0],"rank_%d_%s", mpi_rank, name_string);
    FILE* f = fopen(&save_as[0], "w");
    if (f == NULL){
      mpi_abort("Error opening file to write node vec!\n");
    }

    pnv_user_data.f = f;
    p4est_iterate(p4est,
                  NULL,
                  (void*)&pnv_user_data,
                  testd4est_mixed_element_data_print_node_vec_callback,
                  NULL,
#if (P4EST_DIM)==3
                  NULL,
#endif
                  NULL);

    fclose(f);
  }
  
  else {
    printf("%d: %s = \n", mpi_rank, name_string);
    p4est_iterate(p4est,
                  NULL,
                  (void*)&pnv_user_data,
                  testd4est_mixed_element_data_print_node_vec_callback,
                  NULL,
#if (P4EST_DIM)==3
                  NULL,
#endif
                  NULL);
  }
    
  /* reset pointer */
  p4est->user_pointer = tmp;
}

typedef struct {
  
  double* vec;
  double* norm_sqr;
  int* stride;

} compute_norm_user_data_t;

static
void
testd4est_mixed_element_data_compute_l2_norm_sqr_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  compute_norm_user_data_t* compute_norm_user_data
    = (compute_norm_user_data_t*) user_data;
  
  double* l2_norm_sqr = compute_norm_user_data->norm_sqr;
  double* vec = compute_norm_user_data->vec; 
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  
  int* stride = compute_norm_user_data->stride;
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );
    
  double* Mvec = P4EST_ALLOC(double, volume_nodes);
  dgmath_apply_Mij(dgmath_jit_dbase, &vec[*stride], (P4EST_DIM), elem_data->deg, Mvec);
  d4est_linalg_vec_scale(jacobian, Mvec, volume_nodes);

  /* if it's wanted */
  elem_data->local_estimator = d4est_linalg_vec_dot(&vec[*stride], Mvec, volume_nodes);

  /* printf("local_estimator = %f\n", elem_data->local_estimator); */
  *l2_norm_sqr += elem_data->local_estimator;

  *stride = *stride + volume_nodes;
  P4EST_FREE(Mvec);
}

static
void
testd4est_mixed_element_data_compute_l2_norm_sqr_no_local_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  compute_norm_user_data_t* compute_norm_user_data
    = (compute_norm_user_data_t*) user_data;
  
  double* l2_norm_sqr = compute_norm_user_data->norm_sqr;
  double* vec = compute_norm_user_data->vec; 
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  
  int* stride = compute_norm_user_data->stride;
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );
    
  double* Mvec = P4EST_ALLOC(double, volume_nodes);
  dgmath_apply_Mij(dgmath_jit_dbase, &vec[*stride], (P4EST_DIM), elem_data->deg, Mvec);
  d4est_linalg_vec_scale(jacobian, Mvec, volume_nodes);

  *l2_norm_sqr += d4est_linalg_vec_dot(&vec[*stride], Mvec, volume_nodes);;

  *stride = *stride + volume_nodes;
  P4EST_FREE(Mvec);
}

double
testd4est_mixed_element_data_compute_l2_norm_sqr_no_local
(
 p4est_t* p4est,
 double* nodal_vec,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  double l2_norm_sqr = 0.;
  compute_norm_user_data_t compute_norm_user_data;
  compute_norm_user_data.vec = nodal_vec;
  compute_norm_user_data.norm_sqr = &l2_norm_sqr; 
  compute_norm_user_data.stride = &stride; 
  
  p4est_iterate(p4est,
		NULL,
		(void *) &compute_norm_user_data,
		testd4est_mixed_element_data_compute_l2_norm_sqr_no_local_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
  
  /* reset pointer */
  p4est->user_pointer = tmp;
  return l2_norm_sqr;
}

/* double */
/* testd4est_mixed_element_data_compute_H1_norm_sqr */
/* ( */
/*  p4est_t* p4est, */
/*  double* nodal_vec */
/* ) */
/* { */
/*   void* tmp = p4est->user_pointer; */
/*   p4est->user_pointer = nodal_vec; */

/*   double h1_norm_sqr = 0.; */
  
/*   p4est_iterate(p4est, */
/* 		NULL, */
/* 		(void *) &(h1_norm_sqr), */
/* 		testd4est_mixed_element_data_compute_H1_norm_sqr_callback, */
/* 		NULL, */
/* #if (P4EST_DIM)==3 */
/*                  NULL,        */
/* #endif                 */
/* 		NULL); */


/*   /\* reset pointer *\/ */
/*   p4est->user_pointer = tmp; */
/*   return h1_norm_sqr; */
/* } */


static
void
testd4est_mixed_element_data_compute_DG_norm_sqr_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;

  compute_norm_user_data_t* compute_norm_user_data
    = (compute_norm_user_data_t*) user_data;
  
  double* dg_norm_sqr = compute_norm_user_data->norm_sqr;
  double* vec = compute_norm_user_data->vec; 
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  double surface_jacobian = elem_data->surface_jacobian;
  int stride = elem_data->stride;
  int dim = (P4EST_DIM);
  int deg = elem_data->deg;
  int faces = 2*dim;
  double h = elem_data->h;
  int volume_nodes = dgmath_get_nodes(dim, elem_data->deg );
  int face_nodes = dgmath_get_nodes(dim-1,deg);

  double* Mvec = P4EST_ALLOC(double, volume_nodes);
  double* dvec = P4EST_ALLOC(double, volume_nodes);
  double* Msigmavec_d = P4EST_ALLOC(double, face_nodes);
  double* sigmavec_d;
  
  int f,d;
  for (f = 0; f < faces; f++){
    for (d = 0; d < (P4EST_DIM); d++){
      sigmavec_d = &(elem_data->qstar_min_q[d][f*face_nodes]);
      dgmath_apply_Mij(dgmath_jit_dbase, sigmavec_d, dim - 1, deg, Msigmavec_d);
      d4est_linalg_vec_scale(surface_jacobian, Msigmavec_d, face_nodes);
      *dg_norm_sqr += d4est_linalg_vec_dot(sigmavec_d, Msigmavec_d, face_nodes);
    }      
  }
  
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_Dij(dgmath_jit_dbase, &vec[stride], dim, deg, i, dvec);
    d4est_linalg_vec_scale(2./h, dvec, dgmath_get_nodes(dim, deg));
    dgmath_apply_Mij(dgmath_jit_dbase, dvec, (P4EST_DIM), elem_data->deg, Mvec);
    d4est_linalg_vec_scale(jacobian, Mvec, volume_nodes);
    *dg_norm_sqr += d4est_linalg_vec_dot(dvec, Mvec, volume_nodes);
  }

  P4EST_FREE(Mvec);
  P4EST_FREE(Msigmavec_d);
  P4EST_FREE(dvec);
}


static void
testd4est_mixed_element_data_DG_norm_sqr_init_vecs_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  testd4est_mixed_element_data_t* testd4est_mixed_element_data = (testd4est_mixed_element_data_t *) q->p.user_data;
  double* nodal_vec = (double*) user_data;

  int dim = (P4EST_DIM);
  int deg = testd4est_mixed_element_data->deg;
  int i;
  
  for (i = 0; i < (P4EST_DIM); i++){
    testd4est_mixed_element_data->qstar_min_q[i] = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*dgmath_get_nodes(dim-1, deg));
  }

  testd4est_mixed_element_data->u_elem = &(nodal_vec[testd4est_mixed_element_data->stride]);

  d4est_util_copy_1st_to_2nd(
  testd4est_mixed_element_data->u_elem,
    &(testd4est_mixed_element_data->u_storage)[0],
    dgmath_get_nodes(dim, deg)
  );
}


double
testd4est_mixed_element_data_compute_l2_norm_sqr
(
 p4est_t* p4est,
 double* nodal_vec,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  double l2_norm_sqr = 0.;
  compute_norm_user_data_t compute_norm_user_data;
  compute_norm_user_data.vec = nodal_vec;
  compute_norm_user_data.norm_sqr = &l2_norm_sqr; 
  compute_norm_user_data.stride = &stride; 
    
  p4est_iterate(p4est,
		NULL,
		(void *) (&compute_norm_user_data),
		testd4est_mixed_element_data_compute_l2_norm_sqr_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
  
  /* reset pointer */
  p4est->user_pointer = tmp;
  return l2_norm_sqr;
}


void
testd4est_mixed_element_data_get_local_nodes_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  int* local_nodes = (int*) user_data;
  *local_nodes += dgmath_get_nodes( (P4EST_DIM),
                                   elem_data->deg);
}

int testd4est_mixed_element_data_get_local_nodes(p4est_t* p4est)
{
  int local_nodes = 0;
  
  p4est_iterate(p4est,
		NULL,
		(void *) (&local_nodes),
		testd4est_mixed_element_data_get_local_nodes_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
  
  return local_nodes;
}

int testd4est_mixed_element_data_get_local_matrix_nodes(p4est_t* p4est)
{
  int local_matrix_nodes = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        testd4est_mixed_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        local_matrix_nodes += volume_nodes*volume_nodes;
      }
    }
  return local_matrix_nodes;
}

/** 
 * Initialize jacobians,h,xyz,id
 * and elemental stride
 * 
 * @param p4est 
 * @param deg set to -1, if you do
 * not want to change the stored
 * degree
 */
void testd4est_mixed_element_data_init
(
 p4est_t* p4est,
 int deg
)
{
  int testd4est_mixed_element_data_stride = 0;
  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = (void*)&testd4est_mixed_element_data_stride;
  
  p4est_iterate(p4est,
		NULL,
		(void*) &deg,
		testd4est_mixed_element_data_init_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmpptr;
}

static void
testd4est_mixed_element_data_print_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  printf("**************\n");
  printf("Quad ID = %d\n", elem_data->id);
  printf("Stride = %d\n", elem_data->stride);
  printf("Degree = %d\n", elem_data->deg);
  printf("Estimator = %.25f\n", elem_data->local_estimator);
  printf("Predictor = %.25f\n", elem_data->local_predictor);
  printf("**************\n");
}

void
testd4est_mixed_element_data_print
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		testd4est_mixed_element_data_print_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}



/** 
 * TODO: iterate over each corner and avg values
 * from cells instead of just storing estimator in
 * each corner
 * 
 * @param p4est 
 * @param est_corner 
 * 
 * @return 
 */
void
testd4est_mixed_element_data_store_local_estimator_in_corner_array
(
 p4est_t* p4est,
 double* est_corner
)
{
  p4est_iterate(p4est,
		NULL,
		est_corner,
		testd4est_mixed_element_data_store_local_estimator_in_corner_array_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}

void
testd4est_mixed_element_data_store_nodal_vec_in_vertex_array
(
 p4est_t* p4est,
 double* nodal_vec,
 double* corner_vec
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = corner_vec;
  
  p4est_iterate(p4est,
		NULL,
		nodal_vec,
		testd4est_mixed_element_data_store_nodal_vec_in_vertex_array_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmp;
}


static void
testd4est_mixed_element_data_print_local_estimator_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  testd4est_mixed_element_data_t* elem_data = (testd4est_mixed_element_data_t*) q->p.user_data;
  printf("Quad %d: Local estimator %.20f, p = %d\n", elem_data->id, elem_data->local_estimator, elem_data->deg);
}


void
testd4est_mixed_element_data_print_local_estimator
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		testd4est_mixed_element_data_print_local_estimator_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}


static void
testd4est_mixed_element_data_copy_from_vec_to_storage_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  testd4est_mixed_element_data_t* testd4est_mixed_element_data = (testd4est_mixed_element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = testd4est_mixed_element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  d4est_util_copy_1st_to_2nd
    (
     &u[*stride],
     &(testd4est_mixed_element_data->u_storage)[0],
     volume_nodes
    );

  *stride += volume_nodes;
}

void
testd4est_mixed_element_data_copy_from_vec_to_storage(
 p4est_t* p4est,
 double* vec
 )
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = vec;
  int stride = 0;
  
  p4est_iterate(p4est,
		NULL,
		(void*)&stride,
		testd4est_mixed_element_data_copy_from_vec_to_storage_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmp;
}

static void
testd4est_mixed_element_data_copy_from_storage_to_vec_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  testd4est_mixed_element_data_t* testd4est_mixed_element_data = (testd4est_mixed_element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = testd4est_mixed_element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  d4est_util_copy_1st_to_2nd
    (
     &(testd4est_mixed_element_data->u_storage)[0],
     &u[*stride],
     volume_nodes
    );

  *stride += volume_nodes;
}

void
testd4est_mixed_element_data_copy_from_storage_to_vec
(
 p4est_t* p4est,
 double* vec
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = vec;
  int stride = 0;
  
  p4est_iterate(p4est,
		NULL,
		(void*)&stride,
		testd4est_mixed_element_data_copy_from_storage_to_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmp;
}


