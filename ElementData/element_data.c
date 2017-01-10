/* TURN DEBUG ON/OFF */
#define NDEBUG

#include <stdio.h>
#include <stdlib.h>
#include "../pXest/pXest.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/linalg.h"
#include "../dGMath/dgmath.h"
#include "../GridFunctions/grid_functions.h"
#include "../ElementData/element_data.h"
#include "../Flux/compute_flux.h"
#include "../Flux/dg_norm.h"

/* static int element_data_stride; */

typedef
struct {
  double *u;
  double *v;
  double* out;
  grid_fcn_ext_t f;
  int* stride;
  /* int nonlin_deg; */
} fofuv_data_t;

static void
element_data_init_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;
  element_data_t* elem_data = (element_data_t *) q->p.user_data;
  p4est_connectivity_t* pconnect = info->p4est->connectivity;
  
  p4est_qcoord_t dq;

  int* element_data_stride = ((int*)info->p4est->user_pointer);
  int set_deg = *((int*)user_data);
  /* if set_deg >= 0, use it as the degree of the element */
  /* if set_deg == -1, keep the degree as is */

  /* printf("deg, set_deg = %d,%d\n",elem_data->deg, set_deg); */
  
  if ( set_deg >= 0)
    elem_data->deg = set_deg;

  elem_data->id = info->quadid;
  elem_data->stride = *element_data_stride;

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
  *element_data_stride += dgmath_get_nodes( (P4EST_DIM) , elem_data->deg);

  /* if (elem_data->stride == 0) */
    /* dgmath_set_max_degree_used(elem_data->deg); */
}


/* static void */
/* element_data_apply_fof_uxyz_to_node_vec_callback */
/* ( */
/*  p4est_iter_volume_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   element_data_t* elem_data = (element_data_t*) q->p.user_data; */
/*   grid_uxy_fcn_t init_fcn = (grid_uxy_fcn_t)user_data; */

/*   double* vec = (double*)info->p4est->user_pointer; */
/*   int stride = elem_data->stride; */
/*   int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg ); */
/*   double h = elem_data->h; */
  
/*   int i; */
/*   double* x = dgmath_fetch_xyz_nd( (P4EST_DIM), */
/*                                    elem_data->deg, */
/*                                    0); */
/*   double xl = elem_data->xyz_corner[0]; */
  
/*   double* y = dgmath_fetch_xyz_nd( (P4EST_DIM), */
/*                                    elem_data->deg, */
/*                                    1); */
/*   double yl = elem_data->xyz_corner[1]; */
/* #if (P4EST_DIM)==3 */
/*   double* z = dgmath_fetch_xyz_nd( (P4EST_DIM), */
/*                                    elem_data->deg, */
/*                                    2); */
/*   double zl = elem_data->xyz_corner[2]; */
/* #endif */

/*   for (i = 0; i < volume_nodes; i++){ */
/*     vec[stride + i] = init_fcn */
/*                       ( */
/*                        dgmath_rtox(x[i], xl, h), */
/*                        dgmath_rtox(y[i], yl, h) */
/* #if (P4EST_DIM)==3 */
/*                        , */
/*                        dgmath_rtox(z[i], zl, h) */
/* #endif */
/*                       ); */

/* #ifndef NDEBUG */
/*     printf(" vec[stride + i], x, y, z = %f,%f,%f,%f", */
/*            vec[stride+i], */
/*            dgmath_rtox(x[i], xl, h), */
/*            dgmath_rtox(y[i], yl, h), */
/* #if (P4EST_DIM)==3 */
/*            dgmath_rtox(z[i], zl, h) */
/* #else */
/*            1. */
/* #endif */
/*           ); */
/* #endif */
/*   } */
  

/* } */

typedef struct {

  grid_fcn_t init_fcn;
  double* vec;

} init_node_vec_user_data_t;

static void
element_data_init_node_vec_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;

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
element_data_store_local_estimator_in_corner_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  double* corner_vec = (double*)user_data;
  
  double est = elem_data->local_estimator;
  int quadid = elem_data->id;

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    corner_vec[quadid*(P4EST_CHILDREN) + i] = est;
  }
}

static void
element_data_store_nodal_vec_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;

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
element_data_init_node_vec
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
		element_data_init_node_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,
#endif
		NULL);

  /* reset pointer */
  p4est->user_pointer = tmp;
}

static void
element_data_init_quadid_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t *) q->p.user_data;
  elem_data->id = info->quadid;
}

void
element_data_init_quadid
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		element_data_init_quadid_callback,
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
element_data_print_node_vec_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)

{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
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
element_data_print_node_vec
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
                  element_data_print_node_vec_callback,
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
                  element_data_print_node_vec_callback,
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
element_data_compute_l2_norm_sqr_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
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
  linalg_vec_scale(jacobian, Mvec, volume_nodes);

  /* if it's wanted */
  elem_data->local_estimator = linalg_vec_dot(&vec[*stride], Mvec, volume_nodes);

  /* printf("local_estimator = %f\n", elem_data->local_estimator); */
  *l2_norm_sqr += elem_data->local_estimator;

  *stride = *stride + volume_nodes;
  P4EST_FREE(Mvec);
}

/* static */
/* void */
/* element_data_compute_H1_norm_sqr_callback */
/* ( */
/*  p4est_iter_volume_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
  
/*   p4est_quadrant_t* q = info->quad; */
/*   element_data_t* elem_data = (element_data_t*) q->p.user_data; */
/*   double* h1_norm_sqr = (double*) user_data; */
/*   double* vec = (double*)info->p4est->user_pointer; */
/*   double jacobian = elem_data->jacobian; */

/*   int stride = elem_data->stride; */
/*   int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg ); */
/*   int dim = (P4EST_DIM); */
/*   int deg = elem_data->deg; */
/*   double h = elem_data->h; */
/*   double* Mvec = P4EST_ALLOC(double, volume_nodes); */
/*   double* dvec = P4EST_ALLOC(double, volume_nodes); */
/*   dgmath_apply_Mij(&vec[stride], (P4EST_DIM), elem_data->deg, Mvec); */
/*   linalg_vec_scale(jacobian, Mvec, volume_nodes); */

/*   /\* if it's wanted *\/ */

/*   elem_data->local_estimator = linalg_vec_dot(&vec[stride], Mvec, volume_nodes); */



/*   int i; */
/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     dgmath_apply_Dij(&vec[stride], dim, deg, i, dvec); */
/*     linalg_vec_scale(2./h, dvec, dgmath_get_nodes(dim, deg)); */
/*     dgmath_apply_Mij(dvec, (P4EST_DIM), elem_data->deg, Mvec); */
/*     linalg_vec_scale(jacobian, Mvec, volume_nodes); */
/*     elem_data->local_estimator += linalg_vec_dot(dvec, Mvec, volume_nodes); */
/*   } */
  
/*   *h1_norm_sqr += elem_data->local_estimator; */

/*   P4EST_FREE(Mvec); */
/*   P4EST_FREE(dvec); */
/* } */


static
void
element_data_compute_l2_norm_sqr_no_local_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
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
  linalg_vec_scale(jacobian, Mvec, volume_nodes);

  *l2_norm_sqr += linalg_vec_dot(&vec[*stride], Mvec, volume_nodes);;

  *stride = *stride + volume_nodes;
  P4EST_FREE(Mvec);
}

double
element_data_compute_l2_norm_sqr_no_local
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
		element_data_compute_l2_norm_sqr_no_local_callback,
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
/* element_data_compute_H1_norm_sqr */
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
/* 		element_data_compute_H1_norm_sqr_callback, */
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
element_data_compute_DG_norm_sqr_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;

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
      linalg_vec_scale(surface_jacobian, Msigmavec_d, face_nodes);
      *dg_norm_sqr += linalg_vec_dot(sigmavec_d, Msigmavec_d, face_nodes);
    }      
  }
  
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    dgmath_apply_Dij(dgmath_jit_dbase, &vec[stride], dim, deg, i, dvec);
    linalg_vec_scale(2./h, dvec, dgmath_get_nodes(dim, deg));
    dgmath_apply_Mij(dgmath_jit_dbase, dvec, (P4EST_DIM), elem_data->deg, Mvec);
    linalg_vec_scale(jacobian, Mvec, volume_nodes);
    *dg_norm_sqr += linalg_vec_dot(dvec, Mvec, volume_nodes);
  }

  P4EST_FREE(Mvec);
  P4EST_FREE(Msigmavec_d);
  P4EST_FREE(dvec);
}


static void
element_data_DG_norm_sqr_init_vecs_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;
  double* nodal_vec = (double*) user_data;

  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int i;
  
  for (i = 0; i < (P4EST_DIM); i++){
    element_data->qstar_min_q[i] = P4EST_ALLOC_ZERO(double, (P4EST_FACES)*dgmath_get_nodes(dim-1, deg));
  }

  element_data->u_elem = &(nodal_vec[element_data->stride]);

  linalg_copy_1st_to_2nd(
  element_data->u_elem,
    &(element_data->u_storage)[0],
    dgmath_get_nodes(dim, deg)
  );
}

static
void element_data_DG_norm_sqr_destroy_vecs_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;
  int i;
  for (i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(element_data->qstar_min_q[i]);
  }
}


/** 
 * The energy or dg norm, see
 *
 * @article
 * {
 * houston2008posteriori,
 * title={A posteriori error analysis of hp-version discontinuous Galerkin 
 * finite-element methods for second-order quasi-linear elliptic PDEs},
 * author={Houston, Paul and S{\"u}li, Endre and Wihler, Thomas P},
 * journal={IMA journal of numerical analysis},
 * volume={28},
 * number={2},
 * pages={245--273},
 * year={2008},
 * publisher={Oxford University Press}
 * }
 * 
 * @param p4est 
 * @param nodal_vec 
 * 
 * @return 
 */
double
element_data_compute_DG_norm_sqr
(
 p4est_t* p4est,
 double* nodal_vec,
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* ip_flux_params,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data
)
{ 
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = nodal_vec;

  double dg_norm_sqr = 0.;

  /* init q_flux in element_data to store [[u]], and copy node vec u to u_storage*/
  p4est_iterate(p4est,
		NULL,
		(void *) nodal_vec,
		element_data_DG_norm_sqr_init_vecs_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  /* exchange u_storage */
  p4est_ghost_exchange_data(p4est,ghost,ghost_data);  

  flux_fcn_ptrs_t dg_norm_sqr_flux_fcn_ptrs = dg_norm_ip_flux_dirichlet_fetch_fcns(
                                                                                   bndry_fcn,
                                                                                   ip_flux_params
                                                                                  );
  compute_flux_user_data_t compute_flux_user_data;
  compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  compute_flux_user_data.flux_fcn_ptrs = &dg_norm_sqr_flux_fcn_ptrs;
    
  /* p4est->user_pointer = &dg_norm_sqr_flux_fcn_ptrs; */
  p4est->user_pointer = &compute_flux_user_data;
  
  /* calculate sqrt(sigma)[[u]] on each face */
  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  compute_norm_user_data_t compute_norm_user_data;
  compute_norm_user_data.norm_sqr = &dg_norm_sqr;
  compute_norm_user_data.vec = nodal_vec;

  p4est->user_pointer = dgmath_jit_dbase;

  /* calcaulate dg norm sqr on each element */
  p4est_iterate(p4est,
  		NULL,
  		(void*)&compute_norm_user_data,
  		element_data_compute_DG_norm_sqr_callback,
  		NULL,
  #if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  /* free q_flux allocation */
  p4est_iterate(p4est,
		NULL,
		NULL,
		element_data_DG_norm_sqr_destroy_vecs_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
  
  /* reset pointer */
  p4est->user_pointer = tmp;
  return dg_norm_sqr;
}



double
element_data_compute_l2_norm_sqr
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
		element_data_compute_l2_norm_sqr_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
  
  /* reset pointer */
  p4est->user_pointer = tmp;
  return l2_norm_sqr;
}

static
void
element_data_integrate_au_andaddto_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  int stride = elem_data->stride;
  double deg = elem_data->deg;
  
  double** data = (double**) user_data;
  double a = *(data[0]);
  double* u = data[1];
  double* a_u = data[2];

  double* M_u = P4EST_ALLOC(
                            double,
                            dgmath_get_nodes
                            (
                             (P4EST_DIM),
                             deg
                            )
                           );

  dgmath_apply_Mij
    (
     dgmath_jit_dbase,
     &u[stride],
     (P4EST_DIM),
     deg,
     M_u
    );
  
  linalg_vec_axpy(
                  a*jacobian,
                  M_u,
                  &a_u[stride],
                  dgmath_get_nodes
                  (
                   (P4EST_DIM),
                   deg
                  )
                 );

  P4EST_FREE(M_u);
}

static
void
element_data_integrate_auv_andaddto_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  int stride = elem_data->stride;
  
  double** data = (double**) user_data;
  double a = *(data[0]);
  double* u = data[1];
  double* v = data[2];
  double* a_uv = data[3];

  int degH = elem_data->deg;
  int degh = 2*elem_data->deg;

  double* auv_restrict = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degH
                          )
                         );

  double* Mauv_restrict = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degH
                          )
                         );

  double* u_prolong = P4EST_ALLOC
                      (
                       double,
                       dgmath_get_nodes
                       (
                        (P4EST_DIM),
                        degh
                       )
                      );
  
  double* v_prolong = P4EST_ALLOC
                      (
                       double,
                       dgmath_get_nodes
                       (
                        (P4EST_DIM),
                        degh
                       )
                      );

  
  double* auv_prolong = P4EST_ALLOC
                      (
                       double,
                       dgmath_get_nodes
                       (
                        (P4EST_DIM),
                        degh
                       )
                      );

  dgmath_apply_p_prolong
    (
     dgmath_jit_dbase,
     &u[stride],
     degH,
     (P4EST_DIM),
     degh,
     u_prolong
    );

  dgmath_apply_p_prolong
    (
     dgmath_jit_dbase,
     &v[stride],
     degH,
     (P4EST_DIM),
     degh,
     v_prolong
    );

  int sqr_nodes = dgmath_get_nodes((P4EST_DIM), degh);
  int i;
  for (i = 0; i < sqr_nodes; i++)
    auv_prolong[i] = a*v_prolong[i]*u_prolong[i];
  
  dgmath_apply_p_restrict
    (
     dgmath_jit_dbase,
     auv_prolong,
     degh,
     (P4EST_DIM),
     degH,
     auv_restrict
    );

  dgmath_apply_Mij
    (
     dgmath_jit_dbase,
     auv_restrict,
     (P4EST_DIM),
     degH,
     Mauv_restrict
    );
  
  linalg_vec_axpy(jacobian, Mauv_restrict, &a_uv[stride], dgmath_get_nodes((P4EST_DIM),
                                                             degH));

  P4EST_FREE(auv_prolong);
  P4EST_FREE(Mauv_restrict);
  P4EST_FREE(u_prolong);
  P4EST_FREE(v_prolong);
  P4EST_FREE(auv_restrict);
}


/* static */
/* void */
/* element_data_integrate_fofuv_andaddto_callback */
/* ( */
/*  p4est_iter_volume_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   element_data_t* elem_data = (element_data_t*) q->p.user_data; */

/*   double jacobian = elem_data->jacobian; */
/*   double h = elem_data->h; */
/*   int stride = elem_data->stride; */

/*   fofuv_data_t* fofuv_data = (fofuv_data_t*) user_data; */
/*   double* u = fofuv_data->u; */
/*   double* v = fofuv_data->v; */
/*   double* out = fofuv_data->out; */
/*   grid_fcn_ext_t f = fofuv_data->f; */

/*   /\* printf("f(1.,1.,1.,1.) = %f \n", f(1.,1.,1.,1.)); *\/ */

/*   int degH = elem_data->deg; */
/*   int degh = fofuv_data->nonlin_deg; */

/*   /\* util_print_matrix(&u[stride], dgmath_get_nodes((P4EST_DIM), degH), 1, "u = ", 0); *\/ */
/*   /\* util_print_matrix(&v[stride], dgmath_get_nodes((P4EST_DIM), degH), 1, "v = ", 0); *\/ */
  
/*   if (degh == -1) */
/*     degh = degH; */

  
/*   double* fofu_restrict = P4EST_ALLOC */
/*                          ( */
/*                           double, */
/*                           dgmath_get_nodes */
/*                           ( */
/*                            (P4EST_DIM), */
/*                            degH */
/*                           ) */
/*                          ); */

/*   double* Mfofu_restrict = P4EST_ALLOC */
/*                          ( */
/*                           double, */
/*                           dgmath_get_nodes */
/*                           ( */
/*                            (P4EST_DIM), */
/*                            degH */
/*                           ) */
/*                          ); */


/*   double* u_prolong = P4EST_ALLOC */
/*                       ( */
/*                        double, */
/*                        dgmath_get_nodes */
/*                        ( */
/*                         (P4EST_DIM), */
/*                         degh */
/*                        ) */
/*                       ); */


/*   double* v_prolong = P4EST_ALLOC */
/*                       ( */
/*                        double, */
/*                        dgmath_get_nodes */
/*                        ( */
/*                         (P4EST_DIM), */
/*                         degh */
/*                        ) */
/*                       ); */

 
  
/*   double* fofu_prolong = P4EST_ALLOC */
/*                       ( */
/*                        double, */
/*                        dgmath_get_nodes */
/*                        ( */
/*                         (P4EST_DIM), */
/*                         degh */
/*                        ) */
/*                       ); */

/*   dgmath_apply_p_prolong */
/*     ( */
/*      &u[stride], */
/*      degH, */
/*      (P4EST_DIM), */
/*      degh, */
/*      u_prolong */
/*     ); */

/*   dgmath_apply_p_prolong */
/*     ( */
/*      &v[stride], */
/*      degH, */
/*      (P4EST_DIM), */
/*      degh, */
/*      v_prolong */
/*     ); */


/*   int nonlin_nodes = dgmath_get_nodes((P4EST_DIM), degh); */
/*   double* x = dgmath_fetch_xyz_nd((P4EST_DIM), degh, 0); */
/*   double xl = elem_data->xyz_corner[0]; */
/*   double* y = dgmath_fetch_xyz_nd((P4EST_DIM), degh, 1); */
/*   double yl = elem_data->xyz_corner[1]; */
/* #if (P4EST_DIM)==3 */
/*   double* z = dgmath_fetch_xyz_nd((P4EST_DIM), degh, 2); */
/*   double zl = elem_data->xyz_corner[2]; */
/* #endif */
  
/*   int i; */
/*   for (i = 0; i < nonlin_nodes; i++) */
/*     fofu_prolong[i] = f(dgmath_rtox(x[i], xl, h), */
/*                         dgmath_rtox(y[i], yl, h), */
/* #if (P4EST_DIM)==3 */
/*                         dgmath_rtox(z[i], zl, h), */
/* #endif */
/*                         u_prolong[i])*v_prolong[i]; */
  
/*   dgmath_apply_p_restrict */
/*     ( */
/*      fofu_prolong, */
/*      degh, */
/*      (P4EST_DIM), */
/*      degH, */
/*      fofu_restrict */
/*     ); */

/*   dgmath_apply_Mij */
/*     ( */
/*      fofu_restrict, */
/*      (P4EST_DIM), */
/*      degH, */
/*      Mfofu_restrict */
/*     ); */
  
/*   linalg_vec_axpy(jacobian, Mfofu_restrict, &out[stride], dgmath_get_nodes((P4EST_DIM), */
/*                                                              degH)); */

/*   free(Mfofu_restrict); */
/*   free(fofu_restrict); */
/*   free(u_prolong); */
/*   free(v_prolong); */
/*   free(fofu_prolong); */
/* } */

/* static */
/* void */
/* element_data_integrate_fofuv_andaddto_callback */
/* ( */
/*  p4est_iter_volume_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   element_data_t* elem_data = (element_data_t*) q->p.user_data; */
/*   dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer; */
  
/*   double jacobian = elem_data->jacobian; */
/*   double h = elem_data->h; */

/*   fofuv_data_t* fofuv_data = (fofuv_data_t*) user_data; */
/*   double* u = fofuv_data->u; */
/*   double* v = fofuv_data->v; */
/*   double* out = fofuv_data->out; */
/*   int* stride = fofuv_data->stride; */
/*   grid_fcn_ext_t f = fofuv_data->f; */

/*   int degH = elem_data->deg; */
  
/*   double* fofu = P4EST_ALLOC */
/*                          ( */
/*                           double, */
/*                           dgmath_get_nodes */
/*                           ( */
/*                            (P4EST_DIM), */
/*                            degH */
/*                           ) */
/*                          ); */

/*   double* Mfofu = P4EST_ALLOC */
/*                          ( */
/*                           double, */
/*                           dgmath_get_nodes */
/*                           ( */
/*                            (P4EST_DIM), */
/*                            degH */
/*                           ) */
/*                          ); */

/*   int volume_nodes = dgmath_get_nodes((P4EST_DIM), degH); */
/*   double* x = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 0); */
/*   double xl = elem_data->xyz_corner[0]; */
/*   double* y = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 1); */
/*   double yl = elem_data->xyz_corner[1]; */
/* #if (P4EST_DIM)==3 */
/*   double* z = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 2); */
/*   double zl = elem_data->xyz_corner[2]; */
/* #endif */
  
/*   int i; */
/*   for (i = 0; i < volume_nodes; i++) */
/*     fofu[i] = f(dgmath_rtox(x[i], xl, h), */
/*                 dgmath_rtox(y[i], yl, h), */
/* #if (P4EST_DIM)==3 */
/*                 dgmath_rtox(z[i], zl, h), */
/* #endif */
/*                 u[*stride + i])*v[*stride + i]; */


/*   dgmath_apply_Mij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      fofu, */
/*      (P4EST_DIM), */
/*      degH, */
/*      Mfofu */
/*     ); */
  
/*   linalg_vec_axpy( */
/*                   jacobian, */
/*                   Mfofu, */
/*                   &out[*stride], */
/*                   dgmath_get_nodes */
/*                   ( */
/*                    (P4EST_DIM), */
/*                    degH */
/*                   ) */
/*                  ); */

/*   *stride = *stride + volume_nodes; */

/*   P4EST_FREE(Mfofu); */
/*   P4EST_FREE(fofu); */
/* } */

static
void
element_data_integrate_fofuv_andaddto_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  double h = elem_data->h;

  fofuv_data_t* fofuv_data = (fofuv_data_t*) user_data;
  double* u = fofuv_data->u;
  double* v = fofuv_data->v;
  double* out = fofuv_data->out;
  int* stride = fofuv_data->stride;
  grid_fcn_ext_t f = fofuv_data->f;

  int degH = elem_data->deg;
  int degh = 2*degH;
  
  double* fofuv = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degH
                          )
                         );

  double* Mfofuv = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degH
                          )
                         );
  

  double* fofuv_prolong = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degh
                          )
                         );

   double* v_prolong = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degh
                          )
                         );



   double* u_prolong = P4EST_ALLOC
                         (
                          double,
                          dgmath_get_nodes
                          (
                           (P4EST_DIM),
                           degh
                          )
                         );

   dgmath_apply_p_prolong
     (
      dgmath_jit_dbase,
      &v[*stride],
      degH,
      (P4EST_DIM),
      degh,
      v_prolong
     );

   dgmath_apply_p_prolong
     (
      dgmath_jit_dbase,
      &u[*stride],
      degH,
      (P4EST_DIM),
      degh,
      u_prolong
     );

  int volume_nodes_H = dgmath_get_nodes((P4EST_DIM), degH);
  int volume_nodes_h = dgmath_get_nodes((P4EST_DIM), degh);
  
  double* x = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degh, 0);
  double xl = elem_data->xyz_corner[0];
  double* y = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degh, 1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degh, 2);
  double zl = elem_data->xyz_corner[2];
#endif
  
  int i;
  for (i = 0; i < volume_nodes_h; i++)
    fofuv_prolong[i] = f(dgmath_rtox(x[i], xl, h),
                dgmath_rtox(y[i], yl, h),
#if (P4EST_DIM)==3
                dgmath_rtox(z[i], zl, h),
#endif
                         u_prolong[i], NULL)*v_prolong[i];

  dgmath_apply_p_restrict
    (
     dgmath_jit_dbase,
     fofuv_prolong,
     degh,
     (P4EST_DIM),
     degH,
     fofuv
    );

  
  dgmath_apply_Mij
    (
     dgmath_jit_dbase,
     fofuv,
     (P4EST_DIM),
     degH,
     Mfofuv
    );
  
  linalg_vec_axpy(
                  jacobian,
                  Mfofuv,
                  &out[*stride],
                  dgmath_get_nodes
                  (
                   (P4EST_DIM),
                   degH
                  )
                 );

  *stride = *stride + volume_nodes_H;

  P4EST_FREE(Mfofuv);
  P4EST_FREE(fofuv);
  P4EST_FREE(fofuv_prolong);
  P4EST_FREE(u_prolong);
  P4EST_FREE(v_prolong);
}


static void
element_data_get_slice_data_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  slice_data_t* slice_data = (slice_data_t*) user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double h = elem_data->h;
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg);
  
  double* x = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 0);
  double xl = elem_data->xyz_corner[0];
  double* y = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 1);
  double yl = elem_data->xyz_corner[1];  
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 2);
  double zl = elem_data->xyz_corner[2];
#endif

  int i;
  for (i = 0; i < volume_nodes; i++){
    double xi = dgmath_rtox(x[i], xl, h);
    double yi = dgmath_rtox(y[i], yl, h);
#if (P4EST_DIM)==3
    double zi = dgmath_rtox(z[i], zl, h);
#endif
    /* printf("xi,yi,zi = %f,%f,%f\n",xi,yi,zi); */
    /* printf("scf(0.1,0.,0.0) = %d\n", slice_data->scf(.1,0.,0.)); */
    if (
        slice_data->scf( xi,
                         yi
#if (P4EST_DIM)==3
                         ,
                         zi
#endif
                       ) == 1
    )
      {
        /* just use function to count nodes in slice */
        if(slice_data->stride == -1){
          slice_data->slice_nodes++;
        }
        else{
          if(slice_data->xyz[0] != NULL)
            slice_data->xyz[0][slice_data->stride] = xi;
        
          if(slice_data->xyz[1] != NULL)
            slice_data->xyz[1][slice_data->stride] = yi;
        
#if (P4EST_DIM)==3
          if(slice_data->xyz[2] != NULL)
            slice_data->xyz[2][slice_data->stride] = zi;
#endif
          slice_data->slice[slice_data->stride] = slice_data->node_vec[elem_data->stride + i];
          slice_data->stride++;
        }
      }
  }
}

void
element_data_get_slice_data_count
(
 p4est_t* p4est,
 slice_data_t* slice_data
)
{
  slice_data->stride = -1;
  slice_data->slice_nodes = 0;

  p4est_iterate(p4est,
		NULL,
		(void *) (slice_data),
		element_data_get_slice_data_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);
}

void
element_data_get_slice_data
(
 p4est_t* p4est,
 slice_data_t* slice_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  slice_data->stride = 0;
  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;
  
  p4est_iterate(p4est,
		NULL,
		(void *) (slice_data),
		element_data_get_slice_data_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmpptr;
}

static
void
element_data_integrate_fofu_andaddto_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double jacobian = elem_data->jacobian;
  double h = elem_data->h;

  fofuv_data_t* fofuv_data = (fofuv_data_t*) user_data;
  double* u = fofuv_data->u;
  double* out = fofuv_data->out;
  int* stride = fofuv_data->stride;
  grid_fcn_ext_t f_fcn = fofuv_data->f;
 
  int degH = elem_data->deg;
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), degH);
  
  double* fofu = P4EST_ALLOC
                 (
                  double,
                  dgmath_get_nodes
                  (
                   (P4EST_DIM),
                   degH
                  )
                 );

  double* Mfofu = P4EST_ALLOC
                  (
                   double,
                   dgmath_get_nodes
                   (
                    (P4EST_DIM),
                    degH
                   )
                  );


  double* x = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 0);
  double xl = elem_data->xyz_corner[0];
  double* y = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), degH, 2);
  double zl = elem_data->xyz_corner[2];
#endif

  /* util_print_matrix(x, volume_nodes, 1, "x = ", 0); */

  /* printf("quadid = %d, h = %f\n", info->quadid, elem_data->h); */
  int i;
  for (i = 0; i < volume_nodes; i++){
    fofu[i] = f_fcn(dgmath_rtox(x[i], xl, h),
                        dgmath_rtox(y[i], yl, h),
#if (P4EST_DIM)==3
                        dgmath_rtox(z[i], zl, h),
#endif
                    u[*stride + i],
                   NULL);
    /* printf("f[%d] = %f, u[%d] =  %f, x,y,z = %f,%f,%f xl,yl,zl = %f,%f,%f\n", *stride + i, fofu[i], *stride + i, u[*stride + i], dgmath_rtox(x[i], xl, h), dgmath_rtox(y[i], yl, h), dgmath_rtox(z[i], zl, h), xl,yl,zl); */
  }

  dgmath_apply_Mij
    (
     dgmath_jit_dbase,
     fofu,
     (P4EST_DIM),
     degH,
     Mfofu
    );
  
  linalg_vec_axpy(jacobian, Mfofu, &out[*stride], dgmath_get_nodes((P4EST_DIM),
                                                             degH));


  *stride = *stride + volume_nodes;
  P4EST_FREE(Mfofu);
  P4EST_FREE(fofu);
}

void
element_data_integrate_fofuv_andaddto
(
 p4est_t* p4est,
 double* u,
 double* v,
 double* out,
 /* int nonlin_deg, */
 grid_fcn_ext_t f,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  int stride = 0;
  
  fofuv_data_t fofuv_data;
  fofuv_data.u = u;
  fofuv_data.v = v;
  fofuv_data.out = out;
  fofuv_data.stride = &stride;
  /* fofuv_data.nonlin_deg = nonlin_deg; */
  fofuv_data.f = f;

  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;
  
  if (v != NULL)
  p4est_iterate(p4est,
                NULL,
                (void *)(&fofuv_data),
                element_data_integrate_fofuv_andaddto_callback,
                NULL,
 #if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
  else
      p4est_iterate(p4est,
                NULL,
                (void *)(&fofuv_data),
                element_data_integrate_fofu_andaddto_callback,
                NULL,
 #if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);    

  p4est->user_pointer = tmpptr;
}

void
element_data_integrate_auv_andaddto
(
 p4est_t* p4est,
 double a,
 double* u,
 double* v,
 double* a_uv,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  double* data [4];
  data[0] = &a;
  data[1] = u;
  data[2] = v;
  data[3] = a_uv;

  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  p4est_iterate(p4est,
		NULL,
		(void *) (&data[0]),
		element_data_integrate_auv_andaddto_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);


  p4est->user_pointer = tmpptr;
}

void
element_data_integrate_au_andaddto
(
 p4est_t* p4est,
 double a,
 double* u,
 double* a_u,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  double* data [3];
  data[0] = &a;
  data[1] = u;
  data[2] = a_u;

  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;                 
  
  p4est_iterate(p4est,
		NULL,
		(void *) (&data[0]),
		element_data_integrate_au_andaddto_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmpptr;
}

void
element_data_get_local_nodes_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  int* local_nodes = (int*) user_data;
  *local_nodes += dgmath_get_nodes( (P4EST_DIM),
                                   elem_data->deg);
}

int element_data_get_local_nodes(p4est_t* p4est)
{
  int local_nodes = 0;
  
  p4est_iterate(p4est,
		NULL,
		(void *) (&local_nodes),
		element_data_get_local_nodes_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
  
  return local_nodes;
}

int element_data_get_local_matrix_nodes(p4est_t* p4est)
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
        element_data_t* ed = quad->p.user_data;
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
void element_data_init
(
 p4est_t* p4est,
 int deg
)
{
  int element_data_stride = 0;
  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = (void*)&element_data_stride;
  
  p4est_iterate(p4est,
		NULL,
		(void*) &deg,
		element_data_init_callback,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmpptr;
}

static void
element_data_print_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  printf("**************\n");
  printf("Quad ID = %d\n", elem_data->id);
  printf("Stride = %d\n", elem_data->stride);
  printf("Degree = %d\n", elem_data->deg);
  printf("Estimator = %.25f\n", elem_data->local_estimator);
  printf("Predictor = %.25f\n", elem_data->local_predictor);
  printf("**************\n");
}

static void
element_data_get_xyz_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  double h = elem_data->h;
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg);
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  
  double** xyz = (double**) user_data;
  
  double* x = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 0);
  double xl = elem_data->xyz_corner[0];
  double* y = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 1);
  double yl = elem_data->xyz_corner[1];  
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd(dgmath_jit_dbase, (P4EST_DIM), elem_data->deg, 2);
  double zl = elem_data->xyz_corner[2];
#endif

  int i;
  for (i = 0; i < volume_nodes; i++){
    xyz[0][i + elem_data->stride] = dgmath_rtox(x[i], xl, h);
    xyz[1][i + elem_data->stride] = dgmath_rtox(y[i], yl, h);
#if (P4EST_DIM)==3
    xyz[2][i + elem_data->stride] = dgmath_rtox(z[i], zl, h);
#endif
  }
}

void
element_data_get_xyz
(
 p4est_t* p4est,
 double* xyz [(P4EST_DIM)]
)
{
  p4est_iterate(p4est,
		NULL,
		&xyz[0],
		element_data_get_xyz_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}


void
element_data_print
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		element_data_print_callback,
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
element_data_store_local_estimator_in_corner_array
(
 p4est_t* p4est,
 double* est_corner
)
{
  p4est_iterate(p4est,
		NULL,
		est_corner,
		element_data_store_local_estimator_in_corner_array_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}

void
element_data_store_nodal_vec_in_vertex_array
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
		element_data_store_nodal_vec_in_vertex_array_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmp;
}


static void
element_data_print_local_estimator_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  printf("Quad %d: Local estimator %.20f, p = %d\n", elem_data->id, elem_data->local_estimator, elem_data->deg);
}


void
element_data_print_local_estimator
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		element_data_print_local_estimator_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}


static void
element_data_copy_from_vec_to_storage_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  linalg_copy_1st_to_2nd
    (
     &u[*stride],
     &(element_data->u_storage)[0],
     volume_nodes
    );

  *stride += volume_nodes;
}

void
element_data_copy_from_vec_to_storage(
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
		element_data_copy_from_vec_to_storage_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmp;
}

static void
element_data_copy_from_storage_to_vec_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  element_data_t* element_data = (element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  linalg_copy_1st_to_2nd
    (
     &(element_data->u_storage)[0],
     &u[*stride],
     volume_nodes
    );

  *stride += volume_nodes;
}

void
element_data_copy_from_storage_to_vec
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
		element_data_copy_from_storage_to_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmp;
}

typedef struct {

  grid_fcn_ext_t f_fcn;
  double* f;
  double* u;
  int* stride;

} compute_f_of_uxyz_user_data_t;

static void
element_data_compute_f_of_uxyz_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;
  
  compute_f_of_uxyz_user_data_t* fuxyz_user_data = (compute_f_of_uxyz_user_data_t*) user_data;
  grid_fcn_ext_t f_fcn = fuxyz_user_data->f_fcn;
  double* u = fuxyz_user_data->u;
  double* f = fuxyz_user_data->f;  
  int* stride = fuxyz_user_data->stride;
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );
  double h = elem_data->h;
  
  int i;
  double* x = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   elem_data->deg,
                                   0);
  double xl = elem_data->xyz_corner[0];
  
  double* y = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   elem_data->deg,
                                   1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   elem_data->deg,
                                   2);
  double zl = elem_data->xyz_corner[2];
#endif

  for (i = 0; i < volume_nodes; i++){
    f[*stride + i] = f_fcn
                      (
                       dgmath_rtox(x[i], xl, h),
                       dgmath_rtox(y[i], yl, h)
#if (P4EST_DIM)==3
                       ,
                       dgmath_rtox(z[i], zl, h)
#endif
                       ,
                       u[*stride + i],
                       NULL
                      );

    /* printf("f[%d] = %f, u[%d] = %f, x,y,z = %f,%f,%f\n", *stride + i, f[*stride + i], *stride + i, u[*stride + i], dgmath_rtox(x[i], xl, h), dgmath_rtox(y[i], yl, h), dgmath_rtox(z[i], zl, h)); */

      }
  
  *stride = *stride + volume_nodes;
}

void
element_data_compute_f_of_uxyz
(
 p4est_t* p4est,
 double* u,
 double* f,
 grid_fcn_ext_t f_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  
  compute_f_of_uxyz_user_data_t fuxyz_user_data;
  fuxyz_user_data.f_fcn = f_fcn;
  fuxyz_user_data.f = f;
  fuxyz_user_data.u = u;
  fuxyz_user_data.stride = &stride;
  
  p4est_iterate(p4est,
		NULL,
		(void*)&fuxyz_user_data,
		element_data_compute_f_of_uxyz_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);

  p4est->user_pointer = tmp;
}

typedef struct {
  double *u;
  double *Mu;
  int* stride;
}
apply_Mij_on_vec_user_data_t;


static
void element_data_apply_Mij_on_vec_callback(
                              p4est_iter_volume_info_t * info,
                              void *user_data)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  apply_Mij_on_vec_user_data_t* apply_Mij_on_vec_user_data = (apply_Mij_on_vec_user_data_t*)user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;

  int* stride = apply_Mij_on_vec_user_data->stride;
  double* u = apply_Mij_on_vec_user_data->u;
  double* Mu = apply_Mij_on_vec_user_data->Mu;                                       
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );  
  int deg = elem_data->deg;
  int dim = (P4EST_DIM);


  dgmath_apply_Mij(dgmath_jit_dbase, &u[*stride], dim, deg, &Mu[*stride]);
  linalg_vec_scale(elem_data->jacobian, &Mu[*stride], volume_nodes);
  *stride = *stride + volume_nodes;
}

void
element_data_apply_Mij_on_vec
(
 p4est_t* p4est,
 double* u,
 double* Mu,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  
  apply_Mij_on_vec_user_data_t apply_Mij_on_vec_user_data;
  apply_Mij_on_vec_user_data.u = u;
  apply_Mij_on_vec_user_data.Mu = Mu;
  apply_Mij_on_vec_user_data.stride = &stride;
  
  
  p4est_iterate(p4est,
		NULL,
		(void*)&apply_Mij_on_vec_user_data,
		element_data_apply_Mij_on_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmp;
}


/**
* proj_deltap tells us how high of a space we want to interpolate
* the vec to before we compute f of vec and Mij x f of vec
* p' = p + proj_deltap. If proj_deltap = 0, no interpolation/projection
* will occur.
*/


typedef struct {

  grid_fcn_ext_t f_fcn;
  double* u;
  double* Mu;
  int* stride;
  int proj_deltap;
  
} apply_Mij_on_f_of_vec_user_data_t;


static
void element_data_apply_Mij_on_f_of_vec_callback(
                              p4est_iter_volume_info_t * info,
                              void *user_data)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  apply_Mij_on_f_of_vec_user_data_t* apply_Mij_on_f_of_vec_user_data = (apply_Mij_on_f_of_vec_user_data_t*)user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;

  int* stride = apply_Mij_on_f_of_vec_user_data->stride;
  double* u = apply_Mij_on_f_of_vec_user_data->u;
  double* Mu = apply_Mij_on_f_of_vec_user_data->Mu;                                       
  grid_fcn_ext_t f_fcn = apply_Mij_on_f_of_vec_user_data->f_fcn;
  int proj_deltap = apply_Mij_on_f_of_vec_user_data->proj_deltap;

  int deg_old = elem_data->deg;
  int deg_new = elem_data->deg + proj_deltap;
  
  int volume_nodes_deg_old = dgmath_get_nodes( (P4EST_DIM), deg_old );  
  int volume_nodes_deg_new = dgmath_get_nodes( (P4EST_DIM), deg_new);

  double* u_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  double* f_of_u_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  double* M_f_of_u_interp = P4EST_ALLOC(double, volume_nodes_deg_new);

  dgmath_apply_p_prolong(dgmath_jit_dbase, &u[*stride], deg_old, (P4EST_DIM), deg_new, u_interp);

  double* x = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   0);
  double xl = elem_data->xyz_corner[0];
  
  double* y = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   2);
  double zl = elem_data->xyz_corner[2];
#endif

  for (int i = 0; i < volume_nodes_deg_new; i++){
    f_of_u_interp[i] = f_fcn
                      (
                       dgmath_rtox(x[i], xl, elem_data->h),
                       dgmath_rtox(y[i], yl, elem_data->h)
#if (P4EST_DIM)==3
                       ,
                       dgmath_rtox(z[i], zl, elem_data->h)
#endif
                       ,
                       u_interp[i],
                       NULL
                      );

  }
  
  dgmath_apply_Mij(dgmath_jit_dbase, f_of_u_interp, (P4EST_DIM), deg_new, M_f_of_u_interp);
  linalg_vec_scale(elem_data->jacobian, M_f_of_u_interp, volume_nodes_deg_new);
  dgmath_apply_p_prolong_transpose(dgmath_jit_dbase, M_f_of_u_interp, deg_new, (P4EST_DIM), deg_old, &Mu[*stride]);
  
  *stride = *stride + volume_nodes_deg_old;

  P4EST_FREE(M_f_of_u_interp);
  P4EST_FREE(f_of_u_interp);
  P4EST_FREE(u_interp);
}

void
element_data_apply_Mij_on_f_of_vec
(
 p4est_t* p4est,
 double* u,
 double* Mu,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 grid_fcn_ext_t f_fcn,
 int proj_deltap 
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  
  apply_Mij_on_f_of_vec_user_data_t apply_Mij_on_f_of_vec_user_data;
  apply_Mij_on_f_of_vec_user_data.u = u;
  apply_Mij_on_f_of_vec_user_data.Mu = Mu;
  apply_Mij_on_f_of_vec_user_data.stride = &stride;
  apply_Mij_on_f_of_vec_user_data.f_fcn = f_fcn;
  apply_Mij_on_f_of_vec_user_data.proj_deltap = proj_deltap;

  p4est_iterate(p4est,
		NULL,
		(void*)&apply_Mij_on_f_of_vec_user_data,
		element_data_apply_Mij_on_f_of_vec_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmp;
}



typedef struct {

  grid_fcn_ext_t f_fcn;
  double* vec1;
  double* vec2;
  double* M_f_of_vec1_x_vec2;
  int* stride;
  int proj_deltap;
  
} apply_Mij_on_f_of_vec1_x_vec2_user_data_t;


static
void element_data_apply_Mij_on_f_of_vec1_x_vec2_callback(
                              p4est_iter_volume_info_t * info,
                              void *user_data)
{
  p4est_quadrant_t* q = info->quad;
  element_data_t* elem_data = (element_data_t*) q->p.user_data;
  apply_Mij_on_f_of_vec1_x_vec2_user_data_t* apply_Mij_on_f_of_vec1_x_vec2_user_data = (apply_Mij_on_f_of_vec1_x_vec2_user_data_t*)user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer;

  int* stride = apply_Mij_on_f_of_vec1_x_vec2_user_data->stride;
  double* vec1 = apply_Mij_on_f_of_vec1_x_vec2_user_data->vec1;
  double* vec2 = apply_Mij_on_f_of_vec1_x_vec2_user_data->vec2;
  double* M_f_of_vec1_x_vec2= apply_Mij_on_f_of_vec1_x_vec2_user_data->M_f_of_vec1_x_vec2;                                       
  grid_fcn_ext_t f_fcn = apply_Mij_on_f_of_vec1_x_vec2_user_data->f_fcn;
  int proj_deltap = apply_Mij_on_f_of_vec1_x_vec2_user_data->proj_deltap;

  int deg_old = elem_data->deg;
  int deg_new = elem_data->deg + proj_deltap;
  
  int volume_nodes_deg_old = dgmath_get_nodes( (P4EST_DIM), deg_old );  
  int volume_nodes_deg_new = dgmath_get_nodes( (P4EST_DIM), deg_new);

  double* vec1_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  double* vec2_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  double* f_of_vec1_interp_x_vec2_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  double* M_f_of_vec1_interp_x_vec2_interp = P4EST_ALLOC(double, volume_nodes_deg_new);
  
  dgmath_apply_p_prolong(dgmath_jit_dbase, &vec1[*stride], deg_old, (P4EST_DIM), deg_new, vec1_interp);
  dgmath_apply_p_prolong(dgmath_jit_dbase, &vec2[*stride], deg_old, (P4EST_DIM), deg_new, vec2_interp);

  double* x = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   0);
  double xl = elem_data->xyz_corner[0];
  
  double* y = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   1);
  double yl = elem_data->xyz_corner[1];
#if (P4EST_DIM)==3
  double* z = dgmath_fetch_xyz_nd( dgmath_jit_dbase, (P4EST_DIM),
                                   deg_new,
                                   2);
  double zl = elem_data->xyz_corner[2];
#endif

  for (int i = 0; i < volume_nodes_deg_new; i++){
    f_of_vec1_interp_x_vec2_interp[i] = f_fcn
                      (
                       dgmath_rtox(x[i], xl, elem_data->h),
                       dgmath_rtox(y[i], yl, elem_data->h)
#if (P4EST_DIM)==3
                       ,
                       dgmath_rtox(z[i], zl, elem_data->h)
#endif
                       ,
                       vec1_interp[i],
                       NULL
                      )*vec2_interp[i];
  }
  
  dgmath_apply_Mij(dgmath_jit_dbase, f_of_vec1_interp_x_vec2_interp, (P4EST_DIM), deg_new, M_f_of_vec1_interp_x_vec2_interp);
  linalg_vec_scale(elem_data->jacobian, M_f_of_vec1_interp_x_vec2_interp, volume_nodes_deg_new);
  dgmath_apply_p_prolong_transpose(dgmath_jit_dbase, M_f_of_vec1_interp_x_vec2_interp, deg_new, (P4EST_DIM), deg_old, &M_f_of_vec1_x_vec2[*stride]);
  
  *stride = *stride + volume_nodes_deg_old;

  P4EST_FREE(M_f_of_vec1_interp_x_vec2_interp);
  P4EST_FREE(f_of_vec1_interp_x_vec2_interp);
  P4EST_FREE(vec2_interp);
  P4EST_FREE(vec1_interp);
}

void
element_data_apply_Mij_on_f_of_vec1_x_vec2
(
 p4est_t* p4est,
 double* vec1,
 double* vec2,
 double* M_f_of_vec1_x_vec2,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 grid_fcn_ext_t f_fcn,
 int proj_deltap 
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  int stride = 0;
  
  apply_Mij_on_f_of_vec1_x_vec2_user_data_t apply_Mij_on_f_of_vec1_x_vec2_user_data;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.vec1 = vec1;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.vec2 = vec2;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.M_f_of_vec1_x_vec2 = M_f_of_vec1_x_vec2;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.stride = &stride;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.f_fcn = f_fcn;
  apply_Mij_on_f_of_vec1_x_vec2_user_data.proj_deltap = proj_deltap;

  p4est_iterate(p4est,
		NULL,
		(void*)&apply_Mij_on_f_of_vec1_x_vec2_user_data,
		element_data_apply_Mij_on_f_of_vec1_x_vec2_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmp;
}

static
void element_data_compute_boundary_integral_callback
(
 p4est_iter_face_info_t * info,
 void *user_data
)
{
  boundary_integral_data_t* bi_data = (boundary_integral_data_t*) user_data;
  dgmath_jit_dbase_t * dgmath_jit_dbase = (dgmath_jit_dbase_t*) info->p4est->user_pointer;
  sc_array_t *sides = &(info->sides);

  if (sides->elem_count == 1)
    {
      p4est_iter_face_side_t *side;
      side = p4est_iter_fside_array_index_int(sides,0);
      *(bi_data->boundary_integral) +=
        bi_data->boundary_integral_fcn
        (
         (element_data_t *) side->is.full.quad->p.user_data,
         side->face,
         bi_data->ctx,
         dgmath_jit_dbase
        );

      /* printf("*(bi_data->boundary_integral) = %f\n", *(bi_data->boundary_integral)); */
    }
  
}

double
element_data_compute_boundary_integral
(
 p4est_t* p4est,
 boundary_integral_fcn_t boundary_integral_fcn,
 void* ctx,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = dgmath_jit_dbase;

  double boundary_integral = 0;
  boundary_integral_data_t bi_data;
  bi_data.boundary_integral_fcn = boundary_integral_fcn;
  bi_data.ctx = ctx;
  bi_data.boundary_integral = &boundary_integral;
  
 p4est_iterate(p4est,
               NULL,
               (void*) &bi_data,
               NULL,
               element_data_compute_boundary_integral_callback,
#if (P4EST_DIM)==3
               NULL,
#endif
               NULL);

 /* printf("boundary_integral = %f\n", boundary_integral); */
 
  p4est->user_pointer = tmp;
  return boundary_integral;
}

int
element_data_which_quadrant_of_root
(
 element_data_t* elem_data
)
{  
  double x0 = elem_data->xyz_corner[0];
  double y0 = elem_data->xyz_corner[1];
  double xh = x0 + .5*elem_data->h;
  double yh = y0 + .5*elem_data->h;
  double dx = xh - .5;
  double dy = yh - .5;
  
#if (P4EST_DIM)==3
  double z0 = elem_data->xyz_corner[2];
  double zh = z0 + .5*elem_data->h;
  double dz = zh - .5;
#else
  double dz = -1.;
#endif
  
  int octant = -1;
  if (dx < 0 && dy < 0 && dz < 0)
    octant = 0;
  else if (dx > 0 && dy < 0 && dz < 0)
    octant = 1;
  else if (dx < 0 && dy > 0 && dz < 0)
    octant = 2;
  else if (dx > 0 && dy > 0 && dz < 0)
    octant = 3;
  else if (dx < 0 && dy < 0 && dz > 0)
    octant = 4;
  else if (dx > 0 && dy < 0 && dz > 0)
    octant = 5;
  else if (dx < 0 && dy > 0 && dz > 0)
    octant = 6;
  else if (dx > 0 && dy > 0 && dz > 0)
    octant = 7;
  else
    mpi_abort("octant doesn't exist");

  return octant;
}


double
element_data_get_local_h_min
(
 p4est_t* p4est
)
{
  double h_min = 1.;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        double h = ((element_data_t*)(quad->p.user_data))->h; 
        if (q==0)
          h_min = h;
        else
          h_min = (h < h_min) ? h : h_min;
      }
    }

  return h_min;
}
/* int */
/* element_data_find_element_containing_global_id */
/* ( */
/*  p4est_t* p4est, */
/*  int global_id, */
/*  int local_nodes, */
/*  element_data_t* element_data */
/* ) */
/* { */
/*   int* local_nodes_array = P4EST_ALLOC_ZERO(int, p4est->mpisize); */
/*   int* local_strides_array = P4EST_ALLOC_ZERO(int, p4est->mpisize + 1); */

/*   sc_allgather */
/*     ( */
/*      &local_nodes, */
/*      1, */
/*      sc_MPI_INT, */
/*      local_nodes_array, */
/*      1, */
/*      sc_MPI_INT, */
/*      sc_MPI_COMM_WORLD */
/*     ); */

/*   int stride = 0; */
/*   int global_nodes = 0; */
/*   for(int i = 0; i < p4est->mpisize; i++){ */
/*     global_nodes += local_nodes_array[i]; */
/*     /\* printf("rank %d: local_nodes_array[%d] = %d\n", p4est->mpirank, i, local_nodes_array[i]); *\/ */
/*   } */
/*   /\* printf("global_nodes = %d\n", global_nodes); *\/ */
/*   for (int i = 0; i < p4est->mpisize + 1; i++){ */
/*     local_strides_array[i] = stride; */
/*     /\* printf("stride = %d\n", stride); *\/ */
/*     if (i < p4est->mpisize) */
/*       stride += local_nodes_array[i]; */
/*   } */

  
/*   /\* not on this process *\/ */
/*   if (global_id >= local_strides_array[p4est->mpirank] && global_id < local_strides_array[p4est->mpirank+1]){ */
/*     element_data = NULL; */
/*     return 0; */
/*   } */
/*   else { */
/*     p4est_iterate(p4est, */
/*                   NULL, */
/*                   (void*) &bi_data, */
/*                   NULL, */
/*                   element_data_compute_boundary_integral_callback, */
/* #if (P4EST_DIM)==3 */
/*                   NULL, */
/* #endif */
/*                   NULL); */

/*   } */

/*   P4EST_FREE(local_nodes_array); */
/*   P4EST_FREE(local_strides_array); */
/* } */

element_data_t*
element_data_get_element_data
(
 p4est_t* p4est,
 int local_element_id
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        element_data_t* ed = quad->p.user_data;
        if (ed->id == local_element_id)
          return ed;
      }
    }  
  return NULL;
}

void
element_data_get_level_range
(
 p4est_t* p4est,
 int* min_level,
 int* max_level
)
{

  int level;
  int min = P4EST_MAXLEVEL;
  int max = -1;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        level = quad->level;
        min = (level < min) ? level : min;
        max = (level > max) ? level : max;
      }
    }  

  *min_level = min;
  *max_level = max;
}

double
element_data_compute_l2_norm_error_no_local
(
 p4est_t* p4est,
 double* vec,
 int nodes,
 grid_fcn_t analytical_solution,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  double* vec_analytic = P4EST_ALLOC(double, nodes);
  element_data_init_node_vec(p4est,vec_analytic,analytical_solution,dgmath_jit_dbase);
  linalg_vec_axpy(-1., vec, vec_analytic, nodes);
  double err = element_data_compute_l2_norm_sqr_no_local(p4est,vec_analytic,dgmath_jit_dbase);
  return err;
}
