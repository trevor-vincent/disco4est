#include "../ElementData/curved_element_data.h"
#include "../GridFunctions/grid_functions.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/linalg.h"
#include <curved_compute_flux.h>
#include <curved_dg_norm.h>
#include <ip_flux_aux.h>
#include <d4est_geometry.h>
#include <sc_reduce.h>

typedef struct {
  grid_fcn_t init_fcn;
  double* vec;
  int* stride;
} init_node_vec_user_data_t;

typedef struct {
  double* vec;
  double* norm_sqr;
  int* stride;
} compute_norm_user_data_t;

/* static void */
/* curved_element_data_init_node_vec_callback */
/* ( */
/*  p4est_iter_volume_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data; */
/*   init_node_vec_user_data_t* inv_user_data = (init_node_vec_user_data_t*) user_data; */
/*   grid_fcn_t init_fcn = inv_user_data->init_fcn; */
/*   double* vec = inv_user_data->vec; */
/*   /\* dgmath_jit_dbase_t* dgbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer; *\/ */
  
/*   int* stride = inv_user_data->stride; */
/*   int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg ); */

/*   int i; */
/*   for (i = 0; i < volume_nodes; i++){ */
/*     vec[*stride + i] = init_fcn */
/*                        ( */
/*                         elem_data->xyz[0][i], */
/*                         elem_data->xyz[1][i] */
/* #if (P4EST_DIM)==3 */
/*                         , */
/*                         elem_data->xyz[2][i] */
/* #endif */
/*                        ); */
/*   } */
/*   *stride += volume_nodes; */
/* } */

/* must be run after the degrees have been set by curved_element_data_init for example */
void
curved_element_data_set_degrees
(
 p4est_t* p4est,
 int (*set_deg_fcn)(int, int, int)
)
{
  int k = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        k++;
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        ed->deg = set_deg_fcn(tt,k,ed->deg);
      }
    }    
}


void
curved_element_data_init_node_vec
(
 p4est_t* p4est,
 double* node_vec,
 grid_fcn_t init_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom
)
{

  double* xyz_temp [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    xyz_temp[d] = P4EST_ALLOC(double, dgmath_get_nodes((P4EST_DIM), (MAX_DEGREE)));
  }
  

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;        
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);

        d4est_geometry_compute_xyz
          (
           dgmath_jit_dbase,
           d4est_geom,
           tt,
           ed->deg,
           LOBATTO,
           ed->q,
           ed->dq,
           xyz_temp
          );

        
        for (int i = 0; i < volume_nodes; i++){
          node_vec[ed->nodal_stride + i] = init_fcn(xyz_temp[0][i],
                                                    xyz_temp[1][i]
#if (P4EST_DIM)==3
 ,  
                                                    xyz_temp[2][i]
#endif
                                                   );

          /* printf(" %d %.15f %.15f %.15f %.15f\n",ed->id, xyz_temp[0][i], xyz_temp[1][i], xyz_temp[2][i], node_vec[ed->nodal_stride + i]); */
        }

        

        
      }
    }


  for (int d = 0; d < (P4EST_DIM); d++) {
    P4EST_FREE(xyz_temp[d]);
  }
  
}


typedef struct {

  dgmath_jit_dbase_t* dgmath_jit_dbase;
  /* p4est_geometry_t* p4est_geometry; */
  d4est_geometry_t* d4est_geometry;
  geometric_factors_t* geometric_factors;
  int sqr_nodal_stride;
  int sqr_trace_stride;
  int nodal_stride;
  int integ_stride;

  int id_stride;
  int deg;
  int deg_integ;
  /* integ_type_t integ_type; */ /* deprecated */
  int local_nodes;
  int local_sqr_nodes;
  int local_sqr_trace_nodes;
  int local_nodes_integ;
  
} curved_element_data_init_ctx_t;


/* static void */
/* curved_element_data_destroy_callback */
/* ( */
/*  p4est_iter_volume_info_t * info, */
/*  void *user_data */
/* ) */
/* {  */
/*   p4est_quadrant_t* q = info->quad; */
/*   curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data; */
/*   /\* p4est_connectivity_t* pconnect = info->p4est->connectivity; *\/ */

/*   int i,j; */
/*   P4EST_FREE(elem_data->J); */
/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     P4EST_FREE(elem_data->xyz[i]); */
/*     for (j = 0; j < (P4EST_DIM); j++){ */
/*       P4EST_FREE(elem_data->xyz_rst[i][j]); */
/*     } */
/*   } */
/* } */

geometric_factors_t*
geometric_factors_init
(
 p4est_t* p4est
)
{
  geometric_factors_t* geometric_factors = P4EST_ALLOC(geometric_factors_t, 1);
  /* geometric_factors->J = NULL; */
  geometric_factors->J_integ = NULL; 
  geometric_factors->xyz = NULL;
  geometric_factors->xyz_integ = NULL;
  /* geometric_factors->xyz_rst = NULL; */
  geometric_factors->xyz_rst_integ = NULL;
  geometric_factors->xyz_rst_Lobatto_integ = NULL;
  /* geometric_factors->rst_xyz = NULL; */
  geometric_factors->rst_xyz_integ = NULL;
  /* geometric_factors->invM = NULL; */
  /* geometric_factors->invMface = NULL; */

  return geometric_factors;
}

void
geometric_factors_reinit
(
 p4est_t* p4est,
 geometric_factors_t* geometric_factors,
 curved_element_data_local_sizes_t local_sizes
 /* int local_nodes, */
 /* int local_nodes_integ, */
 /* int local_sqr_nodes, */
 /* int local_sqr_trace_nodes */
)
{

  int local_nodes = local_sizes.local_nodes;
  int local_nodes_integ = local_sizes.local_nodes_integ;
  int local_sqr_nodes = local_sizes.local_sqr_nodes;
  int local_sqr_trace_nodes = local_sizes.local_sqr_trace_nodes;
  /* int local_sqr_nodes_invM = local_sizes.local_sqr_nodes_invM; */
    
  int vector_nodes = local_nodes*(P4EST_DIM); 
  /* int matrix_nodes = local_nodes*(P4EST_DIM)*(P4EST_DIM); */
  
  /* geometric_factors->J = P4EST_REALLOC(geometric_factors->J,double,local_nodes); */
  geometric_factors->xyz = P4EST_REALLOC(geometric_factors->xyz,double,vector_nodes);
  geometric_factors->xyz_integ = P4EST_REALLOC(geometric_factors->xyz_integ, double, (P4EST_DIM)*local_nodes_integ);
  /* geometric_factors->xyz_rst = P4EST_REALLOC(geometric_factors->xyz_rst,double,matrix_nodes);   */
  /* geometric_factors->rst_xyz = P4EST_REALLOC(geometric_factors->rst_xyz,double,matrix_nodes); */

  int matrix_nodes_integ = local_nodes_integ*(P4EST_DIM)*(P4EST_DIM);
  geometric_factors->J_integ = P4EST_REALLOC(geometric_factors->J_integ,double,local_nodes_integ);
  geometric_factors->xyz_rst_integ = P4EST_REALLOC(geometric_factors->xyz_rst_integ,double,matrix_nodes_integ);
  geometric_factors->xyz_rst_Lobatto_integ = P4EST_REALLOC(geometric_factors->xyz_rst_Lobatto_integ,double,matrix_nodes_integ);
  geometric_factors->rst_xyz_integ = P4EST_REALLOC(geometric_factors->rst_xyz_integ,double,matrix_nodes_integ);
  /* geometric_factors->invM = P4EST_REALLOC(geometric_factors->invM, double, local_sqr_nodes_invM); */
  /* geometric_factors->invMface = P4EST_REALLOC(geometric_factors->invMface, double, local_sqr_trace_nodes); */

}

void
geometric_factors_destroy
(
 geometric_factors_t* geometric_factors
)
{
  /* P4EST_FREE(geometric_factors->J); */
  P4EST_FREE(geometric_factors->J_integ);
  P4EST_FREE(geometric_factors->xyz);
  P4EST_FREE(geometric_factors->xyz_integ);
  /* P4EST_FREE(geometric_factors->xyz_rst); */
  P4EST_FREE(geometric_factors->xyz_rst_integ);
  P4EST_FREE(geometric_factors->xyz_rst_Lobatto_integ);
  /* P4EST_FREE(geometric_factors->rst_xyz); */
  P4EST_FREE(geometric_factors->rst_xyz_integ);
  /* P4EST_FREE(geometric_factors->invM); */
  /* P4EST_FREE(geometric_factors->invMface); */
  P4EST_FREE(geometric_factors);
}

static void
curved_element_data_print_node_vec_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;
  double* vec = (double*)user_data;
  int volume_nodes = dgmath_get_nodes(
                                      (P4EST_DIM),
                                      elem_data->deg
  );
  int i;
  for (i = 0; i < volume_nodes; i++){
#if (P4EST_DIM)==3
    printf("x, y, z, vec = %f, %f, %f, %f\n",
           elem_data->xyz[0][i],
           elem_data->xyz[1][i],
           elem_data->xyz[2][i],
           vec[elem_data->nodal_stride + i]
          );
#else
    printf("x, y, vec = %f, %f, %f\n",
           elem_data->xyz[0][i],
           elem_data->xyz[1][i],
           vec[elem_data->nodal_stride + i]
          );
#endif  
  }
}

typedef struct {

  double** vecs;
  int num_vecs;
  int* elems;
  int num_elems;
  
} debug_print_node_vecs_t;

static void
curved_element_data_debug_print_node_vecs_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;
  debug_print_node_vecs_t* dpnv = (debug_print_node_vecs_t*)user_data;
  
  int volume_nodes = dgmath_get_nodes(
                                      (P4EST_DIM),
                                      elem_data->deg
  );
  int stride = elem_data->nodal_stride;

  /* check if in element list */
  int i,j;
  int found = 0;
  for (i = 0; i < dpnv->num_elems; i++){
    if (dpnv->elems[i] == elem_data->id){
      found = 1;
      break;
    }
  }
  
  if (found == 1){
    printf("[Element Info]: \n");
    printf("[id]: %d\n", elem_data->id);
    printf("[stride]: %d\n", elem_data->nodal_stride);
    printf("[diam]: %f\n", elem_data->diam);       
    printf("[deg]: %d\n", elem_data->deg);

    printf("i x y z ");
    for (i = 0; i < dpnv->num_vecs; i++){
      printf("v%d ", i);
    }
    /* printf("jac \n"); */

    for (i = 0; i < volume_nodes; i++){
#if (P4EST_DIM)==3
      printf("%d %.10f %.10f %.10f ", i, elem_data->xyz[0][i], elem_data->xyz[1][i], elem_data->xyz[2][i]);
#else
      printf("%d %.10f %.10f", i, elem_data->xyz[0][i], elem_data->xyz[1][i]);
#endif
      for (j = 0; j < dpnv->num_vecs; j++){
        printf("%.15f ", dpnv->vecs[j][i + stride]);
      }
      /* printf("%.15f \n", elem_data->J[i]); */
    }
  }
}

void
curved_element_data_print_node_vec
(
 p4est_t* p4est,
 double* vec
)
{
  p4est_iterate(
                p4est,
                NULL,
                (void*) vec,
                curved_element_data_print_node_vec_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL
  );
}

void
curved_element_data_debug_print_node_vecs
(
 p4est_t* p4est,
 double** vecs,
 int num_vecs,
 int* elems,
 int num_elems
)
{
  debug_print_node_vecs_t dpnv;
  dpnv.vecs = vecs;
  dpnv.num_vecs = num_vecs;
  dpnv.elems = elems;
  dpnv.num_elems = num_elems;
  
  p4est_iterate(
                p4est,
                NULL,
                (void*) &dpnv,
                curved_element_data_debug_print_node_vecs_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL
  );
}


double
curved_element_data_compute_diam
(
 double* xyz [(P4EST_DIM)],
 int deg,
 diam_compute_option_t option
)
{
  double diam = 0.;
    
  /* Use an approximate method to calculate diam: iterate through corners of element*/
  if (option == DIAM_APPROX || option == DIAM_APPROX_CUBE){
    for (int i = 0; i < (P4EST_CHILDREN); i++){
      for (int j = 0; j < (P4EST_CHILDREN); j++){
        int corner_node_i = dgmath_corner_to_node((P4EST_DIM), deg, i);
        int corner_node_j = dgmath_corner_to_node((P4EST_DIM), deg, j);
        double diam_temp = 0;
        for (int d = 0; d < (P4EST_DIM); d++){
          double diam_dx = xyz[d][corner_node_i] - xyz[d][corner_node_j];
          diam_temp += diam_dx*diam_dx;
        }
        diam_temp = sqrt(diam_temp);      
        diam = (diam_temp > diam) ? diam_temp : diam;
      }
    }

    if (option == DIAM_APPROX_CUBE){
      diam *= 1./sqrt(3.);
    }
    
  }
  else {
    int volume_nodes = dgmath_get_nodes((P4EST_DIM), deg);
    for (int i = 0; i < volume_nodes; i++){
      for (int j = 0; j < volume_nodes; j++){
        double diam_temp = 0;
        for (int d = 0; d < (P4EST_DIM); d++){
          double diam_dx = xyz[d][i] - xyz[d][j];
          diam_temp += diam_dx*diam_dx;
        }
        diam_temp = sqrt(diam_temp);      
        diam = (diam_temp > diam) ? diam_temp : diam;
      }
    }
  }

  return diam;
}



int
curved_element_data_is_it_on_boundary
(
 p4est_t* p4est,
 p4est_quadrant_t* q,
 int which_tree
)
{
  p4est_qcoord_t      dh, xyz_temp;
  p4est_connectivity_t *conn = p4est->connectivity;
  int fbsum = 0;
  for (int face = 0; face < (P4EST_FACES); face++){
    if (conn->tree_to_tree[P4EST_FACES * which_tree + face] != which_tree ||
        (int) conn->tree_to_face[P4EST_FACES * which_tree + face] != face) {
    }
    else {
      dh = P4EST_LAST_OFFSET (q->level);
      switch (face / 2) {
      case 0:
        xyz_temp = q->x;
        break;
      case 1:
        xyz_temp = q->y;
        break;
#ifdef P4_TO_P8
      case 2:
        xyz_temp = q->z;
        break;
#endif
      default:
        SC_ABORT_NOT_REACHED ();
        break;
      }
      fbsum += (xyz_temp == ((face & 0x01) ? dh : 0));
    }
    if (fbsum > 0) break;
  }
  
  return (fbsum > 0);
}

double
curved_element_data_compute_element_volume
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int deg_GL,
 double* jac_GL
)
{
  int volume_nodes_GL = dgmath_get_nodes((P4EST_DIM), deg_GL);
  double* MJac = P4EST_ALLOC(double, volume_nodes_GL);
  double* integ_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, deg_GL);
  double volume = 0;
  
#if (P4EST_DIM)==3
  linalg_kron_vec_o_vec_o_vec_dot_x(integ_weights, jac_GL, deg_GL + 1, MJac);
#elif (P4EST_DIM)==2
  linalg_kron_vec_o_vec_dot_x(integ_weights, jac_GL, deg_GL + 1, MJac);
#else
  mpi_abort("only DIM=2 or DIM=3");
#endif
  
  for (int i = 0; i < volume_nodes_GL; i++){
    volume += MJac[i];
  }
  
  P4EST_FREE(MJac);
  return volume;
}

double
curved_element_data_compute_element_face_area
(
 curved_element_data_t* elem_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* geom,
 int face,
 int deg
)
{
  int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, deg);
  double* n_on_face [(P4EST_DIM)];
  double* xyz_on_face [(P4EST_DIM)];
  for (int i = 0; i < (P4EST_DIM); i++){
    n_on_face[i] = NULL;
    xyz_on_face[i] = P4EST_ALLOC(double,face_nodes);
  }

  double* sj_on_face = P4EST_ALLOC(double, face_nodes);

  d4est_geometry_compute_geometric_data_on_mortar(
                                                  elem_data->tree,
                                                  elem_data->q,
                                                  elem_data->dq,
                                                  1,
                                                  1,
                                                  &deg,
                                                  face,
                                                  NULL,
                                                  sj_on_face,
                                                  NULL,
                                                  NULL,
                                                  GAUSS,
                                                  geom,
                                                  dgmath_jit_dbase
                                                 );
                                                  
  double* Msj = P4EST_ALLOC(double, face_nodes);
  double* integ_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, deg);
  double area = 0;
  
#if (P4EST_DIM)==3
  linalg_kron_vec_o_vec_dot_x(integ_weights, sj_on_face, deg + 1, Msj);
#elif (P4EST_DIM)==2
  linalg_kron_vec_dot_x(integ_weights, sj_on_face, deg + 1, Msj);
#else
  mpi_abort("only DIM=2 or DIM=3");
#endif
  
  for (int i = 0; i < face_nodes; i++){
    area += Msj[i];
  }
  P4EST_FREE(Msj);
  
  P4EST_FREE(sj_on_face);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(xyz_on_face[i]);
  }

  return area;
}

curved_element_data_local_sizes_t
curved_element_data_compute_strides_and_sizes
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 curved_element_data_user_fcn_t user_fcn,
 void* user_ctx
)
{
  /* sizes */
  int local_nodes = 0;
  int local_sqr_nodes = 0;
  int local_sqr_trace_nodes = 0;
  int local_nodes_integ = 0;
  /* int local_sqr_nodes_invM = 0; */

  /* strides */
  int sqr_nodal_stride = 0;
  int sqr_trace_stride = 0;
  int nodal_stride = 0;
  int integ_stride = 0;
  int id_stride = 0;  
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* elem_data = (curved_element_data_t*)(quad->p.user_data);
        elem_data->tree = tt;
        elem_data->tree_quadid = q;
        elem_data->dq = P4EST_QUADRANT_LEN(quad->level);
        elem_data->q[0] = quad->x;
        elem_data->q[1] = quad->y;
#if (P4EST_DIM)==3  
        elem_data->q[2] = quad->z;
#endif        
        /* user_fcn should set degree, 
           or the degree will be assumed to be set */
        
        user_fcn(elem_data, user_ctx);
        
        mpi_assert(elem_data->deg > 0 && elem_data->deg_integ > 0);
        
        elem_data->id = id_stride;
        elem_data->sqr_nodal_stride = sqr_nodal_stride;
        elem_data->sqr_trace_stride = sqr_trace_stride;
        elem_data->nodal_stride = nodal_stride;
        elem_data->integ_stride = integ_stride;
        
        int nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg);
        int nodes_integ = dgmath_get_nodes((P4EST_DIM), elem_data->deg_integ);
        int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, elem_data->deg);
        local_nodes += nodes;
        local_sqr_nodes += nodes*nodes;
        local_sqr_trace_nodes += (P4EST_FACES)*face_nodes*face_nodes;
        local_nodes_integ += dgmath_get_nodes((P4EST_DIM), elem_data->deg_integ);

        /* if (elem_data->deg < elem_data->deg_integ){ */
          /* local_sqr_nodes_invM += nodes*nodes; */
        /* } */
        
        sqr_nodal_stride += nodes*nodes;
        sqr_trace_stride += face_nodes*face_nodes*(P4EST_FACES);
        nodal_stride += nodes;
        integ_stride += nodes_integ;
        id_stride += 1;
      }
    }

  curved_element_data_local_sizes_t local_sizes;
  local_sizes.local_nodes = local_nodes;
  local_sizes.local_sqr_nodes = local_sqr_nodes;
  local_sizes.local_sqr_trace_nodes = local_sqr_trace_nodes;
  local_sizes.local_nodes_integ = local_nodes_integ;
  /* local_sizes.local_sqr_nodes_invM = local_sqr_nodes_invM; */
  return local_sizes;
  
}



void
curved_element_data_init_new
(
 p4est_t* p4est,
 geometric_factors_t* geometric_factors,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 curved_element_data_user_fcn_t user_fcn,
 void* user_ctx,
 int compute_geometric_data,
 int set_geometric_aliases
)
{
  curved_element_data_local_sizes_t local_sizes
    = curved_element_data_compute_strides_and_sizes(p4est, dgmath_jit_dbase, d4est_geometry, user_fcn, user_ctx);

  if (compute_geometric_data){
    geometric_factors_reinit
      (
       p4est,
       geometric_factors,
       local_sizes
       /* local_sizes.local_nodes, */
       /* local_sizes.local_nodes_integ, */
       /* local_sizes.local_sqr_nodes, */
       /* local_sizes.local_sqr_trace_nodes */
      );
  }
  
  /* int invM_stride = 0; */
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        curved_element_data_t* elem_data = (curved_element_data_t*)(quad->p.user_data);

        if (set_geometric_aliases){
          elem_data->J_integ = &geometric_factors->J_integ[elem_data->integ_stride];  
          for (int i = 0; i < (P4EST_DIM); i++){
            elem_data->xyz[i] = &geometric_factors->xyz[i*local_sizes.local_nodes + elem_data->nodal_stride];
            elem_data->xyz_integ[i] = &geometric_factors->xyz_integ[i*local_sizes.local_nodes_integ + elem_data->integ_stride];
            for (int j = 0; j < (P4EST_DIM); j++){
              elem_data->xyz_rst_integ[i][j] = &geometric_factors->xyz_rst_integ[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_integ + elem_data->integ_stride];
              elem_data->rst_xyz_integ[i][j] = &geometric_factors->rst_xyz_integ[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_integ + elem_data->integ_stride];
            }
          }
        }
        
        int volume_nodes_integ = dgmath_get_nodes((P4EST_DIM), elem_data->deg_integ);
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg);

        if (compute_geometric_data && set_geometric_aliases){
        d4est_geometry_compute_xyz
          (
           dgmath_jit_dbase,
           d4est_geometry,
           tt,
           elem_data->deg,
           LOBATTO,
           elem_data->q,
           elem_data->dq,
           elem_data->xyz
          );

        d4est_geometry_compute_xyz
          (
           dgmath_jit_dbase,
           d4est_geometry,
           tt,
           elem_data->deg_integ,
           GAUSS,
           elem_data->q,
           elem_data->dq,
           elem_data->xyz_integ
          );

        d4est_geometry_compute_dxyz_drst
          (
           elem_data->tree,
           elem_data->q,
           elem_data->dq,
           elem_data->deg_integ,           
           GAUSS,
           d4est_geometry,
           dgmath_jit_dbase,
           elem_data->xyz_rst_integ
          );

        d4est_geometry_compute_jacobian_and_drst_dxyz
          (
           elem_data->xyz_rst_integ,
           elem_data->J_integ,
           elem_data->rst_xyz_integ,
           volume_nodes_integ
          );
        }
        /* if(elem_data->deg == elem_data->deg_integ) */
        /*   elem_data->invM == NULL; */
        /* else{ */
        /*   elem_data->invM = &geometric_factors->invM[invM_stride]; */
        /*   dgmath_compute_curvedInverseGaussMass */
        /*     ( */
        /*      dgmath_jit_dbase, */
        /*      elem_data->deg, */
        /*      elem_data->deg_integ, */
        /*      (P4EST_DIM), */
        /*      elem_data->J_integ, */
        /*      elem_data->invM */
        /*     ); */
        /*   invM_stride += volume_nodes*volume_nodes; */
        /* } */

        if (set_geometric_aliases){
        elem_data->volume
          = curved_element_data_compute_element_volume
          (
           dgmath_jit_dbase,
           elem_data->deg_integ,
           elem_data->J_integ
          );
        
        for (int face = 0; face < (P4EST_FACES); face++){
          elem_data->surface_area[face] = curved_element_data_compute_element_face_area(elem_data,dgmath_jit_dbase, d4est_geometry,face, elem_data->deg_integ);
        }
        elem_data->diam = curved_element_data_compute_diam(elem_data->xyz, elem_data->deg, DIAM_APPROX_CUBE);
        }
        
      }
    }
  
}

/* void */
/* curved_element_data_init */
/* ( */
/*  p4est_t* p4est, */
/*  geometric_factors_t* geometric_factors, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  /\* p4est_geometry_t* p4est_geometry, *\/ */
/*  d4est_geometry_t* d4est_geometry, */
/*  int deg, */
/*  /\* int* deg_integ_diff, *\/ */
/*  int deg_integ */
/* ) */
/* {   */
/*   curved_element_data_init_ctx_t ctx; */
/*   /\* ctx.p4est_geometry = p4est_geometry; *\/ */
/*   ctx.d4est_geometry = d4est_geometry; */
/*   ctx.dgmath_jit_dbase = dgmath_jit_dbase; */
/*   ctx.geometric_factors = geometric_factors; */
/*   ctx.deg = deg; */
/*   ctx.deg_integ = deg_integ; */
/*   /\* ctx.deg_integ_diff = deg_integ_diff; *\/ */
/*   ctx.id_stride = 0; */
/*   /\* ctx.integ_type = integ_type; *\/ */

/*   int local_nodes = 0.; */
/*   int local_nodes_integ = 0.; */
/*   int local_sqr_nodes = 0.; */
/*   int local_sqr_trace_nodes = 0.; */
/*   if (deg == -1){ */
/*     for (p4est_topidx_t tt = p4est->first_local_tree; */
/*          tt <= p4est->last_local_tree; */
/*          ++tt) */
/*       { */
/*         p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*         sc_array_t* tquadrants = &tree->quadrants; */
/*         int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*         for (int q = 0; q < Q; ++q) { */
/*           p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*           curved_element_data_t* elem_data = (curved_element_data_t*)(quad->p.user_data);  */
/*           /\* elem_data->deg_integ = elem_data->deg + deg_integ_diff[tt]; *\/ */
/*           elem_data->deg_integ = deg_integ; */
/*           local_nodes_integ += dgmath_get_nodes((P4EST_DIM), elem_data->deg_integ); */
/*           int nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg); */
/*           int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, elem_data->deg); */
/*           local_nodes += nodes; */
/*           local_sqr_nodes += nodes*nodes; */
/*           local_sqr_trace_nodes += (P4EST_FACES)*face_nodes*face_nodes; */
/*         } */
/*       }     */
/*   } */
/*   else if (deg > 0){ */
/*     int nodes = dgmath_get_nodes((P4EST_DIM), deg); */
/*     int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, deg); */
/*     local_nodes = p4est->local_num_quadrants*nodes; */
/*     local_sqr_nodes = p4est->local_num_quadrants*nodes*nodes;  */
/*     local_sqr_trace_nodes = p4est->local_num_quadrants*face_nodes*face_nodes*(P4EST_FACES);  */
/*     for (p4est_topidx_t tt = p4est->first_local_tree; */
/*          tt <= p4est->last_local_tree; */
/*          ++tt) */
/*       { */
/*         p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*         sc_array_t* tquadrants = &tree->quadrants; */
/*         int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*         for (int q = 0; q < Q; ++q) { */
/*           /\* local_nodes_integ += dgmath_get_nodes((P4EST_DIM), deg_integ_diff[tt] + deg); *\/ */
/*           local_nodes_integ += dgmath_get_nodes((P4EST_DIM), deg_integ); */
/*         } */
/*       } */
/*   } */
/*   else */
/*     mpi_abort("[D4EST_ERROR]: curved_element_data_init, deg bad"); */

/*   ctx.local_nodes = local_nodes; */
/*   ctx.local_nodes_integ = local_nodes_integ; */
/*   ctx.local_sqr_nodes = local_sqr_nodes; */
/*   ctx.local_sqr_trace_nodes = local_sqr_trace_nodes; */
/*   ctx.nodal_stride = 0; */
/*   ctx.sqr_nodal_stride = 0; */
/*   ctx.sqr_trace_stride = 0; */
/*   ctx.integ_stride = 0; */
/*   /\* ctx.nodal_vector_stride = 0; *\/ */
/*   /\* ctx.nodal_matrix_stride = 0; *\/ */

/*   /\* reinitialize geometric factors *\/ */

/*   /\* printf("local_nodes_integ = %d\n", local_nodes_integ); *\/ */
/*   geometric_factors_reinit */
/*     ( */
/*      p4est, */
/*      geometric_factors, */
/*      /\* deg, *\/ */
/*      /\* local_nodes, *\/ */
/*      /\* local_nodes_integ, *\/ */
/*      /\* local_sqr_nodes, *\/ */
/*      /\* local_sqr_trace_nodes *\/ */
/*     ); */
  
/*   p4est_iterate(p4est, */
/*                 NULL, */
/*                 (void*) &ctx, */
/*                 curved_element_data_init_callback, */
/*                 NULL, */
/* #if (P4EST_DIM)==3 */
/*                 NULL,        */
/* #endif                 */
/*                 NULL); */

/* } */


/* void */
/* curved_element_data_destroy */
/* ( */
/*  p4est_t* p4est */
/* ) */
/* {  */
/*   p4est_iterate(p4est, */
/* 		NULL, */
/* 		NULL, */
/* 		curved_element_data_destroy_callback, */
/* 		NULL, */
/* #if (P4EST_DIM)==3 */
/*                  NULL,        */
/* #endif                 */
/* 		NULL); */
/* } */


typedef struct {
  double *u;
  double *Mu;
  int* stride;
}
  apply_Mij_on_vec_user_data_t;

/* static */
/* void curved_element_data_apply_Mij_on_vec_callback( */
/*                               p4est_iter_volume_info_t * info, */
/*                               void *user_data) */
/* { */
/*   p4est_quadrant_t* q = info->quad; */
/*   curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data; */
/*   apply_Mij_on_vec_user_data_t* apply_Mij_on_vec_user_data = (apply_Mij_on_vec_user_data_t*)user_data; */
/*   /\* dgmath_jit_dbase_t* dgmath_jit_dbase = (dgmath_jit_dbase_t*)info->p4est->user_pointer; *\/ */

/*   int* stride = apply_Mij_on_vec_user_data->stride; */
/*   /\* double* u = apply_Mij_on_vec_user_data->u; *\/ */
/*   /\* double* Mu = apply_Mij_on_vec_user_data->Mu;                                        *\/ */
/*   int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg );   */
/*   /\* int deg = elem_data->deg; *\/ */
/*   /\* int dim = (P4EST_DIM); *\/ */

/* #ifdef GEOM_DEALIASING */
/*   dgmath_compute_M_vec_wo_aliasing */
/*     ( */
/*      dgmath_jit_dbase, */
/*      elem_data->xyz_rst, */
/*      &u[*stride], */
/*      elem_data->deg, */
/*      &Mu[*stride], */
/*      2*elem_data->deg */
/*     ); */
/* #else */
/*   mpi_abort("NOT SUPPORTED ANYMORE"); */
/*   /\* double* Ju = P4EST_ALLOC(double, volume_nodes);  *\/ */
/*   /\* linalg_component_mult(&u[*stride], &elem_data->J[0], Ju, volume_nodes); *\/ */
/*   /\* dgmath_apply_Mij(dgmath_jit_dbase, Ju, (P4EST_DIM), elem_data->deg, &Mu[*stride]); *\/ */
/*   /\* P4EST_FREE(Ju); *\/ */
/* #endif */
/*   *stride = *stride + volume_nodes;   */
/* } */

/* void */
/* curved_element_data_apply_Mij_on_vec */
/* ( */
/*  p4est_t* p4est, */
/*  double* u, */
/*  double* Mu, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
/*   void* tmp = p4est->user_pointer; */
/*   p4est->user_pointer = dgmath_jit_dbase; */

/*   int stride = 0; */
  
/*   apply_Mij_on_vec_user_data_t apply_Mij_on_vec_user_data; */
/*   apply_Mij_on_vec_user_data.u = u; */
/*   apply_Mij_on_vec_user_data.Mu = Mu; */
/*   apply_Mij_on_vec_user_data.stride = &stride; */
  
/*   p4est_iterate( */
/*                 p4est, */
/*  		NULL, */
/* 		(void*)&apply_Mij_on_vec_user_data, */
/* 		curved_element_data_apply_Mij_on_vec_callback, */
/* 		NULL, */
/* #if (P4EST_DIM)==3 */
/*                 NULL, */
/* #endif */
/* 		NULL */
/*   ); */

/*   p4est->user_pointer = tmp; */
/* } */

/* double */
/* curved_element_data_compute_l2_norm_sqr_no_local */
/* ( */
/*  p4est_t* p4est, */
/*  double* vec, */
/*  int local_nodes, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
/*   double* Mvec = P4EST_ALLOC(double, local_nodes); */
/*   curved_element_data_apply_Mij_on_vec */
/*     ( */
/*      p4est, */
/*      vec, */
/*      Mvec, */
/*      dgmath_jit_dbase */
/*     ); */
/*   double vTMv = linalg_vec_dot(vec, Mvec, local_nodes); */
/*   P4EST_FREE(Mvec); */
/*   return vTMv; */
/* } */

double
curved_element_data_compute_l2_norm_sqr
(
 p4est_t* p4est,
 double* nodal_vec,
 int local_nodes,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 norm_storage_option_t store_local
)
{
  double* Mvec = P4EST_ALLOC(double, local_nodes);
  double l2_norm_sqr = 0.;  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        dgmath_apply_curvedGaussMass(
                                     dgmath_jit_dbase,
                                     &nodal_vec[ed->nodal_stride],
                                     ed->deg,
                                     ed->J_integ,
                                     ed->deg_integ,
                                     (P4EST_DIM),
                                     &Mvec[ed->nodal_stride]
                                    );
        
        double norm2 = linalg_vec_dot(&nodal_vec[ed->nodal_stride], &Mvec[ed->nodal_stride], volume_nodes);
        
        if (store_local == STORE_LOCALLY){
          ed->local_estimator = norm2;
        }
        
        l2_norm_sqr += norm2;
      }
    }
  P4EST_FREE(Mvec);
  return l2_norm_sqr;
}


double
curved_element_data_compute_dg_norm_sqr
(
 p4est_t* p4est,
 double* nodal_vec,
 int local_nodes,
 ip_flux_params_t* ip_flux_params,
 d4est_geometry_t* d4est_geom,
 p4est_ghost_t* ghost,
 void* ghost_data,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  curved_element_data_copy_from_vec_to_storage
    (
     p4est,
     nodal_vec
    );
  
  double dg_norm_sqr = 0.;  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), ed->deg_integ);
        int volume_nodes_Lobatto = dgmath_get_nodes((P4EST_DIM), ed->deg);
        double* dnodal_vec [(P4EST_DIM)]; D4EST_ALLOC_DIM_VEC(dnodal_vec, volume_nodes_Gauss);

        curved_element_compute_derivative_on_Gauss
          (
           &nodal_vec[ed->nodal_stride],
           ed->rst_xyz_integ,
           dnodal_vec,
           ed->deg,
           ed->deg_integ,
           dgmath_jit_dbase
          );
  
        
        for (int d = 0; d < (P4EST_DIM); d++)
          {
            dg_norm_sqr += dgmath_Gauss_quadrature
                           (
                            dgmath_jit_dbase,
                            dnodal_vec[d],
                            dnodal_vec[d],
                            ed->J_integ,
                            ed->deg_integ,
                            (P4EST_DIM)
                           );
          }

        D4EST_FREE_DIM_VEC(dnodal_vec);
      }
    }

  curved_dg_norm_params_t curved_dg_params;
  curved_dg_params.ip_flux_params = ip_flux_params;
  curved_dg_params.dg_norm_face_term = 0.;
  
  curved_flux_fcn_ptrs_t flux_fcn_ptrs = curved_dg_norm_fetch_fcns(&curved_dg_params);

  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.dgmath_jit_dbase = dgmath_jit_dbase;
  curved_compute_flux_user_data.geom = d4est_geom;
  curved_compute_flux_user_data.flux_fcn_ptrs = &flux_fcn_ptrs;
  
  void* tmpptr = p4est->user_pointer;
  p4est->user_pointer = &curved_compute_flux_user_data;
  
  p4est_iterate(p4est,
		ghost,
		ghost_data,
		NULL,
		curved_compute_flux_on_local_elements,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmpptr;

  printf("dg_norm_sqr before face term = %.25f\n", dg_norm_sqr);
  
  dg_norm_sqr += curved_dg_params.dg_norm_face_term;

  printf("dg_norm_sqr after face term = %.25f\n", dg_norm_sqr);
  return dg_norm_sqr;
}

void
curved_element_data_get_local_nodes_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;
  int* local_nodes = (int*) user_data;
  *local_nodes += dgmath_get_nodes( (P4EST_DIM),
                                    elem_data->deg);
}

int curved_element_data_get_local_nodes(p4est_t* p4est)
{
  int local_nodes = 0;
  
  p4est_iterate(p4est,
                NULL,
                (void *) (&local_nodes),
                curved_element_data_get_local_nodes_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);
  
  return local_nodes;
}

static void
curved_element_data_copy_from_vec_to_storage_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* curved_element_data = (curved_element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = curved_element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  linalg_copy_1st_to_2nd
    (
     &u[*stride],
     &(curved_element_data->u_storage)[0],
     volume_nodes
    );

  *stride += volume_nodes;
}

void
curved_element_data_copy_from_vec_to_storage(
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
                curved_element_data_copy_from_vec_to_storage_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);

  p4est->user_pointer = tmp;
}

static void
curved_element_data_copy_from_storage_to_vec_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t *q = info->quad;
  curved_element_data_t* curved_element_data = (curved_element_data_t *) q->p.user_data;
  double* u = (double*) info->p4est->user_pointer;
  int* stride = (int*) user_data;
  
  int dim = (P4EST_DIM);
  int deg = curved_element_data->deg;
  int volume_nodes = dgmath_get_nodes(dim,deg);
  
  linalg_copy_1st_to_2nd
    (
     &(curved_element_data->u_storage)[0],
     &u[*stride],
     volume_nodes
    );

  *stride += volume_nodes;
}


void
curved_element_data_copy_from_storage_to_vec
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
                curved_element_data_copy_from_storage_to_vec_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);

  p4est->user_pointer = tmp;
}

static void
curved_element_data_store_element_scalar_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;
  double* vertex_array_vec = (double*)user_data;
  double (*get_local_scalar_fcn)(curved_element_data_t*) = (double (*)(curved_element_data_t*))info->p4est->user_pointer;  

  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid; 
  p4est_tree_t       *tree;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (info->p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = (P4EST_CHILDREN) * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in */

  double est = get_local_scalar_fcn(elem_data);
  /* int quadid = elem_data->id; */

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    vertex_array_vec[arrayoffset + i] = est;
  }
}

void
curved_element_data_store_element_scalar_in_vertex_array
(
 p4est_t* p4est,
 double* vertex_array,
 double (*get_local_scalar_fcn)(curved_element_data_t*)
)
{
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = (void*)get_local_scalar_fcn;
  
  p4est_iterate(p4est,
                NULL,
                (void*)vertex_array,
                curved_element_data_store_element_scalar_in_vertex_array_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);

  p4est->user_pointer = tmp;
}


static void
curved_element_data_store_nodal_vec_in_vertex_array_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;

  double* nodal_vec = (double*)user_data;
  double* vertex_vec = (double*) info->p4est->user_pointer;
  /* int id = elem_data->id; */
  int stride = elem_data->nodal_stride;

  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid; 
  p4est_tree_t       *tree;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (info->p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = (P4EST_CHILDREN) * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in */
  
  int c;
  for (c = 0; c < (P4EST_CHILDREN); c++){
    int corner_to_node = dgmath_corner_to_node((P4EST_DIM), elem_data->deg, c);
    vertex_vec[arrayoffset + c] = nodal_vec[stride + corner_to_node];
  }
}

typedef struct{

  int id;
  double r;
  double vec_norm;
  double xh;
  double yh;
  double zh;
  int treeid;
  
} debug_spheresym_t;

typedef struct {

  dgmath_jit_dbase_t* dgbase;
  double* vec;

} dgbase_and_vec_t;


/* Return value meaning */
/* <0 The element pointed by p1 goes before the element pointed by p2 */
/* 0  The element pointed by p1 is equivalent to the element pointed by p2 */
/* >0 The element pointed by p1 goes after the element pointed by p2 */
static
int curved_element_data_debug_sphere_compare
(
 const void *p1,
 const void *p2
)
{
  debug_spheresym_t *e1 = (debug_spheresym_t*)p1;
  debug_spheresym_t *e2 = (debug_spheresym_t*)p2;
  double r1 = e1->r;
  double r2 = e2->r;
  return (r1 > r2) - (r1 < r2);
}


static void
curved_element_data_debug_spheresym_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  mpi_abort("deprecated");

  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;

  dgbase_and_vec_t* dgbase_and_vec = (dgbase_and_vec_t*)user_data;
  debug_spheresym_t* debug = (debug_spheresym_t*) info->p4est->user_pointer;

  double* vec = dgbase_and_vec->vec;
  dgmath_jit_dbase_t* dgbase = dgbase_and_vec->dgbase;
  
  int id = elem_data->id;
  int stride = elem_data->nodal_stride;
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), elem_data->deg);
  
  double* Jvec = P4EST_ALLOC(double, volume_nodes);

  int i;
  for (i = 0; i < volume_nodes; i++){
    /* Jvec[i] = elem_data->J[i]*vec[stride + i]; */
  }

  double* Mvec = P4EST_ALLOC(double, volume_nodes);
  dgmath_apply_Mij(dgbase, Jvec, (P4EST_DIM), elem_data->deg, Mvec);
  debug[id].vec_norm = sqrt(linalg_vec_dot(&vec[stride], Mvec, volume_nodes));
  debug[id].id = id;

  double x_0 = elem_data->xyz[0][0];
  double y_0 = elem_data->xyz[1][0];
  double z_0 = elem_data->xyz[2][0];

  double x_1;
  double y_1;
  double z_1;

  for (i = 0; i < (P4EST_CHILDREN); i++){
    int corner_node_i = dgmath_corner_to_node((P4EST_DIM), elem_data->deg, i);
    x_1 = elem_data->xyz[0][corner_node_i];
    y_1 = elem_data->xyz[1][corner_node_i];
    z_1 = elem_data->xyz[2][corner_node_i];
    if (fabs(x_1 - x_0) > 0 && fabs(y_1 - y_0) > 0 && fabs(z_1 - z_0) > 0){
      break;
    }
  }
  double x_half = .5*(x_1 + x_0);
  double y_half = .5*(y_1 + y_0);
  double z_half = .5*(z_1 + z_0);
  debug[id].r = sqrt(x_half*x_half + y_half*y_half + z_half*z_half);
  debug[id].xh = x_half;
  debug[id].yh = y_half;
  debug[id].zh = z_half;
  debug[id].treeid = info->treeid;
  
  P4EST_FREE(Mvec);
  P4EST_FREE(Jvec);
}


void
curved_element_data_debug_spheresym
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgbase,
 double* vec
)
{
  debug_spheresym_t* debug_sphere = P4EST_ALLOC(debug_spheresym_t, p4est->local_num_quadrants);
  dgbase_and_vec_t dgbase_and_vec;
  dgbase_and_vec.dgbase = dgbase;
  dgbase_and_vec.vec = vec;
  
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = debug_sphere;
  
  p4est_iterate(p4est,
                NULL,
                (void*)&dgbase_and_vec,
                curved_element_data_debug_spheresym_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
                NULL);

  p4est->user_pointer = tmp;
  qsort(debug_sphere, p4est->local_num_quadrants, sizeof(debug_spheresym_t), curved_element_data_debug_sphere_compare);

  int i;
  double r_temp = -1;
  double norm_temp = -1;
  for(i = 0; i < p4est->local_num_quadrants; i++){
    if (fabs(debug_sphere[i].r - r_temp) > 0) {
      r_temp = debug_sphere[i].r;
      norm_temp = debug_sphere[i].vec_norm;
    }
    else {
      if (debug_sphere[i].vec_norm > 10*norm_temp || debug_sphere[i].vec_norm < norm_temp/10.){
        printf("Problem Here\n");
      }
    }
    printf("r, id, treeid, x, y, z, norm  = %.10f, %d,%d,%.10f, %.10f, %.10f, %.10f\n", debug_sphere[i].r, debug_sphere[i].id, debug_sphere[i].treeid, debug_sphere[i].xh, debug_sphere[i].yh, debug_sphere[i].zh, debug_sphere[i].vec_norm);
  }
  
  P4EST_FREE(debug_sphere);
}

void
curved_element_data_store_nodal_vec_in_vertex_array
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
                curved_element_data_store_nodal_vec_in_vertex_array_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,
#endif
                NULL);

  p4est->user_pointer = tmp;
}

 
void
curved_element_data_reorient_f_p_elements_to_f_m_order
(
 curved_element_data_t** e_p,
 int face_dim,
 int f_m,
 int f_p,
 int o,
 int faces_p,
 curved_element_data_t* e_p_oriented [(P4EST_HALF)]
)
{

  if (faces_p == 1){
    e_p_oriented[0] = e_p[0];
    return;
  }
  for (int i = 0; i < faces_p; i++){
    int inew = dgmath_reorient_face_order(face_dim, f_m, f_p, o, i);
    e_p_oriented[i] = e_p[inew];
  }
}
 
/* only for serial use */
int
curved_element_data_debug_find_node
(
 p4est_t* p4est,
 int node
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
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        if(node >= ed->nodal_stride && node < (ed->nodal_stride + volume_nodes))
          return ed->id;
      }
    }
  return -1;
}

int
curved_element_data_global_node_to_local_node
(
 p4est_t* p4est,
 int global_node
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
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        if(global_node >= ed->nodal_stride && global_node < (ed->nodal_stride + volume_nodes))
          return global_node - ed->nodal_stride;
      }
    }
  return -1;
}


curved_element_data_t*
curved_element_data_get_element_data
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
        curved_element_data_t* ed = quad->p.user_data;
        if (ed->id == local_element_id)
          return ed;
      }
    }  
  return NULL;
}

int
curved_element_data_count_boundary_elements
(
 p4est_t* p4est
)
{
  int belems = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        belems += ed->on_bdry;
      }
    }  
  return belems;
}

void
curved_element_data_print_element_data_debug
(
 p4est_t* p4est
)
{
/* #ifndef D4EST_DEBUG */
/*   mpi_abort("compile with the debug flag if you want to print curved element data"); */
/* #endif */

  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;

        printf("** Element %d **\n", ed->id);
        printf("deg, deg_integ = %d, %d\n", ed->deg, ed->deg_integ);
        printf("tree, tree_quadid = %d, %d\n", ed->tree, ed->tree_quadid);

        
        int volume_nodes = dgmath_get_nodes( (P4EST_DIM), ed->deg );
       
        
#if (P4EST_DIM)==2 
        printf("q = %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->dq);
        DEBUG_PRINT_2ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes);
        /* DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes); */
#elif (P4EST_DIM)==3
        printf("q = %d, %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->q[2], ed->dq);
        DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], ed->xyz[2], volume_nodes);
#else
        mpi_abort("DIM = 2 or 3");
#endif
        
/* #else */
        /* mpi_abort("DEBUG flag must be set"); */
/* #endif */
      }
    }  
}

static void
curved_element_data_print_local_estimator_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  curved_element_data_t* elem_data = (curved_element_data_t*) q->p.user_data;
  printf("Quad %d: Local estimator %.20f, p = %d\n", elem_data->id, elem_data->local_estimator, elem_data->deg);
}


void
curved_element_data_print_local_estimator
(
 p4est_t* p4est
)
{
  p4est_iterate(p4est,
		NULL,
		NULL,
		curved_element_data_print_local_estimator_callback,
		NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
		NULL);
}


void
curved_element_data_print_number_of_elements_per_tree
(
 p4est_t* p4est
)
{
  int num_trees = p4est->connectivity->num_trees;
  int* elements_per_tree_local = P4EST_ALLOC_ZERO(int, num_trees);
  int* elements_per_tree_global = P4EST_ALLOC_ZERO(int, num_trees);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      elements_per_tree_local[tt] = Q;
    }

    sc_reduce
      (
       &elements_per_tree_local[0],
       &elements_per_tree_global[0],
       num_trees,
       sc_MPI_INT,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    if (p4est->mpirank == 0){
      for (int i = 0; i < num_trees; i++){
        printf(" Tree %d: Number of Elements = %d\n", i, elements_per_tree_global[i]);
      }
    }    
  
  P4EST_FREE(elements_per_tree_local);
  P4EST_FREE(elements_per_tree_global);
}

void curved_element_data_apply_fofufofvlilj_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 double* vec,
 double* u,
 double* v,
 curved_element_data_t* elem_data,
 int deg_Gauss,
 int dim,
 double* Mvec,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{
  int deg_Lobatto = elem_data->deg;
  
  if (deg_Gauss == elem_data->deg_integ)
    {
      dgmath_apply_fofufofvlilj_Gaussnodes
        (
         dgmath_jit_dbase,
         vec,
         u,
         v,
         deg_Lobatto,
         elem_data->J_integ,
         elem_data->xyz_integ,
         deg_Gauss,
         dim,
         Mvec,
         fofu_fcn,
         fofu_ctx,
         fofv_fcn,
         fofv_ctx
        );
    }
  
  else if (deg_Gauss > elem_data->deg_integ){


    int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss);
    double* J_integ = P4EST_ALLOC(double, volume_nodes_Gauss);
    double* xyz_integ [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_integ, volume_nodes_Gauss);
    double* xyz_rst_integ [(P4EST_DIM)][(P4EST_DIM)];
    double* rst_xyz_integ [(P4EST_DIM)][(P4EST_DIM)];
    D4EST_ALLOC_DBYD_MAT(xyz_rst_integ, volume_nodes_Gauss);
    D4EST_ALLOC_DBYD_MAT(rst_xyz_integ, volume_nodes_Gauss);
    
    d4est_geometry_compute_xyz
      (
       dgmath_jit_dbase,
       d4est_geometry,
       elem_data->tree,
       deg_Gauss,
       GAUSS,
       elem_data->q,
       elem_data->dq,
       xyz_integ
      );


    d4est_geometry_compute_dxyz_drst
      (
       elem_data->tree,
       elem_data->q,
       elem_data->dq,
       elem_data->deg_integ,       
       GAUSS,
       d4est_geometry,
       dgmath_jit_dbase,
       elem_data->xyz_rst_integ
      );
    
    d4est_geometry_compute_jacobian_and_drst_dxyz
      (
       elem_data->xyz_rst_integ,
       elem_data->J_integ,
       elem_data->rst_xyz_integ,
       volume_nodes_Gauss
      );
    
    dgmath_apply_fofufofvlilj_Gaussnodes
      (
       dgmath_jit_dbase,
       vec,
       u,
       v,
       deg_Lobatto,
       J_integ,
       xyz_integ,
       deg_Gauss,
       dim,
       Mvec,
       fofu_fcn,
       fofu_ctx,
       fofv_fcn,
       fofv_ctx
      );

    P4EST_FREE(J_integ);
    D4EST_FREE_DIM_VEC(xyz_integ);
    D4EST_FREE_DBYD_MAT(xyz_rst_integ);
    D4EST_FREE_DBYD_MAT(rst_xyz_integ);
  }
  else {
    mpi_abort("deg_Lobatto == elem_data->deg && deg_Gauss >= elem_data->deg_integ");
  }
}


void curved_element_data_apply_fofufofvlj_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 double* u,
 double* v,
 curved_element_data_t* elem_data,
 int deg_Gauss,
 int dim,
 double* out,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{

  int deg_Lobatto = elem_data->deg;
  
  if (deg_Gauss == elem_data->deg_integ)
    {
      /* int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss); */
      /* DEBUG_PRINT_ARR_DBL(elem_data->J_integ, volume_nodes_Gauss); */
      dgmath_apply_fofufofvlj_Gaussnodes
        (
         dgmath_jit_dbase,
         u,
         v,
         deg_Lobatto,
         elem_data->J_integ,
         elem_data->xyz_integ,
         deg_Gauss,
         dim,
         out,
         fofu_fcn,
         fofu_ctx,
         fofv_fcn,
         fofv_ctx
        );

      
    }
  else if (deg_Gauss > elem_data->deg_integ){

    int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss);
    double* J_integ = P4EST_ALLOC(double, volume_nodes_Gauss);
    double* xyz_integ [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_integ, volume_nodes_Gauss);
    double* xyz_rst_integ [(P4EST_DIM)][(P4EST_DIM)];
    double* rst_xyz_integ [(P4EST_DIM)][(P4EST_DIM)];
    D4EST_ALLOC_DBYD_MAT(xyz_rst_integ, volume_nodes_Gauss);
    D4EST_ALLOC_DBYD_MAT(rst_xyz_integ, volume_nodes_Gauss);
    
    
    d4est_geometry_compute_xyz
      (
       dgmath_jit_dbase,
       d4est_geometry,
       elem_data->tree,
       deg_Gauss,
       GAUSS,
       elem_data->q,
       elem_data->dq,
       xyz_integ
      );

    d4est_geometry_compute_dxyz_drst
      (
       elem_data->tree,
       elem_data->q,
       elem_data->dq,
       elem_data->deg_integ,       
       GAUSS,
       d4est_geometry,
       dgmath_jit_dbase,
       elem_data->xyz_rst_integ
      );

    
    d4est_geometry_compute_jacobian_and_drst_dxyz
      (
       elem_data->xyz_rst_integ,
       elem_data->J_integ,
       elem_data->rst_xyz_integ,
       volume_nodes_Gauss
      );
    

    dgmath_apply_fofufofvlj_Gaussnodes
      (
       dgmath_jit_dbase,
       u,
       v,
       deg_Lobatto,
       J_integ,
       xyz_integ,
       deg_Gauss,
       dim,
       out,
       fofu_fcn,
       fofu_ctx,
       fofv_fcn,
       fofv_ctx
      );

    P4EST_FREE(J_integ);
    D4EST_FREE_DIM_VEC(xyz_integ);
    D4EST_FREE_DBYD_MAT(xyz_rst_integ);
    D4EST_FREE_DBYD_MAT(rst_xyz_integ);
  }
  else {
    mpi_abort("deg_Lobatto == elem_data->deg && deg_Gauss >= elem_data->deg_integ");
  }
}


void curved_element_data_form_fofufofvlilj_matrix_Gaussnodes
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 double* u,
 double* v,
 curved_element_data_t* elem_data,
 int deg_Gauss,
 int dim,
 double* mat,
 grid_fcn_ext_t fofu_fcn,
 void* fofu_ctx,
 grid_fcn_ext_t fofv_fcn,
 void* fofv_ctx
)
{
  int deg_Lobatto = elem_data->deg;
  
  if (deg_Gauss == elem_data->deg_integ)
    {
      dgmath_form_fofufofvlilj_matrix_Gaussnodes
        (
         dgmath_jit_dbase,
         u,
         v,
         deg_Lobatto,
         elem_data->xyz_integ,
         elem_data->J_integ,
         deg_Gauss,
         dim,
         mat,
         fofu_fcn,
         fofu_ctx,
         fofv_fcn,
         fofv_ctx
        );
    }
  
  else if (deg_Gauss > elem_data->deg_integ){


    int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss);
    double* J_integ = P4EST_ALLOC(double, volume_nodes_Gauss);
    double* xyz_integ [(P4EST_DIM)];
    D4EST_ALLOC_DIM_VEC(xyz_integ, volume_nodes_Gauss);
    double* xyz_rst_integ [(P4EST_DIM)][(P4EST_DIM)];
    double* rst_xyz_integ [(P4EST_DIM)][(P4EST_DIM)];
    D4EST_ALLOC_DBYD_MAT(xyz_rst_integ, volume_nodes_Gauss);
    D4EST_ALLOC_DBYD_MAT(rst_xyz_integ, volume_nodes_Gauss);
    
    d4est_geometry_compute_xyz
      (
       dgmath_jit_dbase,
       d4est_geometry,
       elem_data->tree,
       deg_Gauss,
       GAUSS,
       elem_data->q,
       elem_data->dq,
       xyz_integ
      );


d4est_geometry_compute_dxyz_drst
      (
       elem_data->tree,
       elem_data->q,
       elem_data->dq,
       elem_data->deg_integ,       
       GAUSS,
       d4est_geometry,
       dgmath_jit_dbase,
       elem_data->xyz_rst_integ
      );

    
    d4est_geometry_compute_jacobian_and_drst_dxyz
      (
       elem_data->xyz_rst_integ,
       elem_data->J_integ,
       elem_data->rst_xyz_integ,
       volume_nodes_Gauss
      );
    
      dgmath_form_fofufofvlilj_matrix_Gaussnodes
      (
       dgmath_jit_dbase,
       u,
       v,
       deg_Lobatto,
       xyz_integ,
       J_integ,
       deg_Gauss,
       dim,
       mat,
       fofu_fcn,
       fofu_ctx,
       fofv_fcn,
       fofv_ctx
      );

    P4EST_FREE(J_integ);
    D4EST_FREE_DIM_VEC(xyz_integ);
    D4EST_FREE_DBYD_MAT(xyz_rst_integ);
    D4EST_FREE_DBYD_MAT(rst_xyz_integ);
  }
  else {
    mpi_abort("deg_Lobatto == elem_data->deg && deg_Gauss >= elem_data->deg_integ");
  }
}

int curved_element_data_get_local_matrix_nodes(p4est_t* p4est){

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
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);
        local_matrix_nodes += volume_nodes*volume_nodes;
      }
    }
  return local_matrix_nodes;
}

void curved_element_data_apply_curvedGaussStiff
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 curved_element_data_t* elem_data,
 double* vec,
 double* stiff_vec
)
{
  mpi_assert(elem_data->deg_stiffness > 0);
  
  if (elem_data->deg_stiffness == elem_data->deg_integ)
    {
      dgmath_apply_curvedGaussStiff
        (
         dgmath_jit_dbase,
         vec,
         elem_data->deg,
         elem_data->J_integ,
         elem_data->rst_xyz_integ,
         elem_data->deg_integ,
         (P4EST_DIM),
         stiff_vec
        );     
    }
  
  else if (elem_data->deg_stiffness > elem_data->deg_integ){

    int volume_nodes_stiffness = dgmath_get_nodes((P4EST_DIM), elem_data->deg_stiffness);
    double* J_stiffness = P4EST_ALLOC(double, volume_nodes_stiffness);
    double* xyz_rst_stiffness [(P4EST_DIM)][(P4EST_DIM)];
    double* rst_xyz_stiffness [(P4EST_DIM)][(P4EST_DIM)];
    D4EST_ALLOC_DBYD_MAT(xyz_rst_stiffness, volume_nodes_stiffness);
    D4EST_ALLOC_DBYD_MAT(rst_xyz_stiffness, volume_nodes_stiffness);
    
 
/*     curved_element_data_compute_dxyz_drst */
/*       ( */
/*        dgmath_jit_dbase, */
/*        elem_data->q, */
/*        elem_data->dq, */
/*        elem_data->tree, */
/*        d4est_geometry, */
/*        elem_data->deg_stiffness, */
/*        1, */
/*        xyz_rst_stiffness, */
/*        (double* [(P4EST_DIM)]){NULL, NULL */
/* #if (P4EST_DIM)==3 */
/*            , NULL */
/* #endif */
/*            } */
/*       ); */

/*     curved_element_data_compute_J_and_rst_xyz */
/*       ( */
/*        xyz_rst_stiffness, */
/*        J_stiffness, */
/*        rst_xyz_stiffness, */
/*        volume_nodes_stiffness */
/*       ); */

d4est_geometry_compute_dxyz_drst
      (
       elem_data->tree,
       elem_data->q,
       elem_data->dq,
       elem_data->deg_stiffness,       
       GAUSS,
       d4est_geometry,
       dgmath_jit_dbase,
       xyz_rst_stiffness
      );

   
    
    d4est_geometry_compute_jacobian_and_drst_dxyz
      (
       xyz_rst_stiffness,
       J_stiffness,
       rst_xyz_stiffness,
       volume_nodes_stiffness
      );


    
      dgmath_apply_curvedGaussStiff
        (
         dgmath_jit_dbase,
         vec,
         elem_data->deg,
         J_stiffness,
         rst_xyz_stiffness,
         elem_data->deg_stiffness,
         (P4EST_DIM),
         stiff_vec
        );



    P4EST_FREE(J_stiffness);
    D4EST_FREE_DBYD_MAT(xyz_rst_stiffness);
    D4EST_FREE_DBYD_MAT(rst_xyz_stiffness);
  }
  else {
    mpi_abort("deg_stiffness >= elem_data->deg_integ");
  }
}

void
curved_element_compute_derivative_on_Gauss
(
 double* vec,
 double* rst_xyz_Gauss [(P4EST_DIM)][(P4EST_DIM)],
 double* dvec [(P4EST_DIM)],
 int deg_Lobatto,
 int deg_Gauss,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{

  int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss);
  int volume_nodes_Lobatto = dgmath_get_nodes((P4EST_DIM), deg_Lobatto);

  double* dvec_di_prolonged = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* dvec_di_Gauss = P4EST_ALLOC(double, volume_nodes_Gauss);
  double* dvec_di_Lobatto = P4EST_ALLOC(double, volume_nodes_Lobatto);

  for (int j = 0; j < (P4EST_DIM); j++){
    for (int k = 0; k < volume_nodes_Gauss; k++){
      dvec[j][k] = 0.;
    }
  }
  
  for (int i = 0; i < (P4EST_DIM); i++){

    dgmath_apply_Dij(dgmath_jit_dbase, vec, (P4EST_DIM), deg_Lobatto, i, dvec_di_Lobatto);
    
    dgmath_apply_p_prolong(dgmath_jit_dbase,
                           dvec_di_Lobatto,
                           deg_Lobatto,
                           (P4EST_DIM),
                           deg_Gauss,
                           dvec_di_prolonged);

    dgmath_interp_GLL_to_GL(dgmath_jit_dbase,
                            dvec_di_prolonged,
                            deg_Gauss,
                            deg_Gauss,
                            dvec_di_Gauss,
                            (P4EST_DIM));


    for (int j = 0; j < (P4EST_DIM); j++){
      for (int k = 0; k < volume_nodes_Gauss; k++){
       dvec[j][k] += rst_xyz_Gauss[i][j][k]*dvec_di_Gauss[k];
      }
    }
  }

  P4EST_FREE(dvec_di_prolonged);
  P4EST_FREE(dvec_di_Gauss);
  P4EST_FREE(dvec_di_Lobatto);
}
