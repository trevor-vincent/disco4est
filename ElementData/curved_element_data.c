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
 p4est_geometry_t* p4est_geom
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
        curved_element_data_compute_xyz(dgmath_jit_dbase,p4est_geom,tt,ed->deg,LOBATTO,ed->q,ed->dq,xyz_temp);

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
curved_element_data_init_callback
(
 p4est_iter_volume_info_t * info,
 void *user_data
)
{
  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;
  curved_element_data_t* elem_data = (curved_element_data_t *) q->p.user_data;
  curved_element_data_init_ctx_t* ctx = (curved_element_data_init_ctx_t*)user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = ctx->dgmath_jit_dbase;
  geometric_factors_t* geometric_factors = ctx->geometric_factors;
  d4est_geometry_t* d4est_geometry = ctx->d4est_geometry;
  p4est_geometry_t* p4est_geometry = ctx->d4est_geometry->p4est_geom;

  elem_data->tree = which_tree;
  
  if ( ctx->deg >= 0){
    elem_data->deg = ctx->deg;
    /* elem_data->deg_integ = ctx->deg + ctx->deg_integ_diff[which_tree]; */
    elem_data->deg_integ = ctx->deg_integ;//ctx->deg + ctx->deg_integ_diff[which_tree];
  }
  else {
    /* elem_data->deg_integ = elem_data->deg + ctx->deg_integ_diff[which_tree]; */
    elem_data->deg_integ = ctx->deg_integ;//elem_data->deg + ctx->deg_integ_diff[which_tree];
  }
  
  p4est_qcoord_t dq;
  elem_data->id = ctx->id_stride;
  elem_data->sqr_nodal_stride = ctx->sqr_nodal_stride;
  elem_data->sqr_trace_stride = ctx->sqr_trace_stride;
  elem_data->nodal_stride = ctx->nodal_stride;
  elem_data->integ_stride = ctx->integ_stride;
  
  int volume_nodes = dgmath_get_nodes( (P4EST_DIM), elem_data->deg);
  int volume_nodes_integ = dgmath_get_nodes( (P4EST_DIM), elem_data->deg_integ);

  elem_data->J_integ = &geometric_factors->J_integ[ctx->integ_stride];
  
  int i,j;
  for (i = 0; i < (P4EST_DIM); i++){
    elem_data->xyz[i] = &geometric_factors->xyz[i*ctx->local_nodes + ctx->nodal_stride];
    elem_data->xyz_integ[i] = &geometric_factors->xyz_integ[i*ctx->local_nodes_integ + ctx->integ_stride];
    for (j = 0; j < (P4EST_DIM); j++){
      elem_data->xyz_rst_integ[i][j] = &geometric_factors->xyz_rst_integ[(i*(P4EST_DIM) + j)*ctx->local_nodes_integ + ctx->integ_stride];
      elem_data->xyz_rst_Lobatto_integ[i][j] = &geometric_factors->xyz_rst_Lobatto_integ[(i*(P4EST_DIM) + j)*ctx->local_nodes_integ + ctx->integ_stride];
      elem_data->rst_xyz_integ[i][j] = &geometric_factors->rst_xyz_integ[(i*(P4EST_DIM) + j)*ctx->local_nodes_integ + ctx->integ_stride];
    }
  }
  
  dq = P4EST_QUADRANT_LEN(q->level);
  
  double* r = dgmath_fetch_xyz_nd( dgmath_jit_dbase,
                                   (P4EST_DIM),
                                   elem_data->deg,
                                   0);
  
  double* s = dgmath_fetch_xyz_nd(dgmath_jit_dbase,
                                  (P4EST_DIM),
                                  elem_data->deg,
                                  1);
#if (P4EST_DIM) == 3
  double* t = dgmath_fetch_xyz_nd( dgmath_jit_dbase,
                                   (P4EST_DIM),
                                   elem_data->deg,
                                   2);
#endif

  elem_data->dq = dq;

  elem_data->q[0] = q->x;
  elem_data->q[1] = q->y;
#if (P4EST_DIM)==3  
  elem_data->q[2] = q->z;
#endif
  
  double abc [3]; /* [0,1]**DIM */
  double xyz [3]; /* curvilinear coordinates */
  /* int d,d1; */
  elem_data->debug_flag = 0;
  for (i = 0; i < volume_nodes; i++){
    abc[0] = dgmath_rtox(r[i], (double)q->x, (double)dq)/(double)P4EST_ROOT_LEN;
    abc[1] = dgmath_rtox(s[i], (double)q->y, (double)dq)/(double)P4EST_ROOT_LEN;
#if (P4EST_DIM) == 3 
    abc[2] = dgmath_rtox(t[i], (double)q->z, (double)dq)/(double)P4EST_ROOT_LEN;
#endif        
    p4est_geometry->X(p4est_geometry, which_tree, abc, xyz);

    for (int d = 0; d < (P4EST_DIM); d++){
      elem_data->xyz[d][i] = xyz[d];
    }
  }

  /* DEBUG_PRINT_ARR_DBL(elem_data->xyz[, volume_nodes); */
  
  double diam;
  double diam_dx;
  int corner_node_i, corner_node_j;
  elem_data->diam = -1.;
  /* loop through all possible corner combinations to get biggest diameter */
  /* This is done because corner numbers change between tree coordinate systems */
  /* So in a cubed sphere example, the diameter for elements with the same radial
   * coordinate, but different trees will not have the same diameter unless we
   * we pass through all possible combinations and find the biggest */
  for (i = 0; i < (P4EST_CHILDREN); i++){
    for (j = 0; j < (P4EST_CHILDREN); j++){
      corner_node_i = dgmath_corner_to_node((P4EST_DIM), elem_data->deg, i);
      corner_node_j = dgmath_corner_to_node((P4EST_DIM), elem_data->deg, j);
      diam = 0;
      for (int d = 0; d < (P4EST_DIM); d++){
        diam_dx = elem_data->xyz[d][corner_node_i] - elem_data->xyz[d][corner_node_j];
        diam += diam_dx*diam_dx;
      }
      diam = sqrt(diam);      
      elem_data->diam = (diam > elem_data->diam) ? diam : elem_data->diam;
    }
  }


  p4est_qcoord_t      dh, xyz_temp;
  p4est_connectivity_t *conn = info->p4est->connectivity;
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
  elem_data->on_bdry = (fbsum > 0);

  curved_element_data_compute_dxyz_drst
    (
     dgmath_jit_dbase,
     elem_data->q,
     elem_data->dq,
     elem_data->tree,
     d4est_geometry->p4est_geom,
     elem_data->deg_integ,
     1,
     elem_data->xyz_rst_integ,
     elem_data->xyz_integ
    );

  for (int i = 0; i < volume_nodes_integ; i++){
    
    double xr = elem_data->xyz_rst_integ[0][0][i];
    double xs = elem_data->xyz_rst_integ[0][1][i];
#if (P4EST_DIM)==3
    double xt = elem_data->xyz_rst_integ[0][2][i];
#endif
    
    double yr = elem_data->xyz_rst_integ[1][0][i];
    double ys = elem_data->xyz_rst_integ[1][1][i];
#if (P4EST_DIM)==3
    double yt = elem_data->xyz_rst_integ[1][2][i];
    
    double zr = elem_data->xyz_rst_integ[2][0][i];
    double zs = elem_data->xyz_rst_integ[2][1][i];
    double zt = elem_data->xyz_rst_integ[2][2][i];
#endif
    
    double* rx = &elem_data->rst_xyz_integ[0][0][i];
    double* ry = &elem_data->rst_xyz_integ[0][1][i];
#if (P4EST_DIM)==3
    double* rz = &elem_data->rst_xyz_integ[0][2][i];
#endif
    double* sx = &elem_data->rst_xyz_integ[1][0][i];
    double* sy = &elem_data->rst_xyz_integ[1][1][i];
#if (P4EST_DIM)==3
    double* sz = &elem_data->rst_xyz_integ[1][2][i];
    
    double* tx = &elem_data->rst_xyz_integ[2][0][i];
    double* ty = &elem_data->rst_xyz_integ[2][1][i];
    double* tz = &elem_data->rst_xyz_integ[2][2][i];
#endif
    double* J = &(elem_data->J_integ[i]);

#if (P4EST_DIM) == 3
    *J = xr*(ys*zt-zs*yt)
                - yr*(xs*zt-zs*xt)
                + zr*(xs*yt-ys*xt);
    *rx =  (ys*zt - zs*yt)/(*J);
    *ry = -(xs*zt - zs*xt)/(*J);
    *rz =  (xs*yt - ys*xt)/(*J);
    *sx = -(yr*zt - zr*yt)/(*J);
    *sy =  (xr*zt - zr*xt)/(*J);
    *sz = -(xr*yt - yr*xt)/(*J);
    *tx =  (yr*zs - zr*ys)/(*J);
    *ty = -(xr*zs - zr*xs)/(*J);
    *tz =  (xr*ys - yr*xs)/(*J);
#elif (P4EST_DIM) == 2
    *J = -xs*yr + xr*ys;
    *rx = ys/(*J);
    *sx =-yr/(*J);
    *ry =-xs/(*J);
    *sy = xr/(*J);
#else
    mpi_abort("DIM must be 2 or 3");
#endif  
  }

  /* if(elem_data->deg == elem_data->deg_integ) */
    /* elem_data->invM == NULL; */
  /* else{ */
  /* elem_data->invM = &geometric_factors->invM[ctx->sqr_nodal_stride]; */
  /* dgmath_compute_curvedInverseGaussMass */
  /*   ( */
  /*    dgmath_jit_dbase, */
  /*    elem_data->deg, */
  /*    elem_data->deg_integ, */
  /*    (P4EST_DIM), */
  /*    elem_data->J_integ, */
  /*    elem_data->invM */
  /*   ); */
  /* } */
  int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, elem_data->deg);  
  double* MJac = P4EST_ALLOC(double, volume_nodes_integ);
  elem_data->volume = 0;
  double* integ_weights = NULL;
  integ_weights = dgmath_fetch_GL_weights_1d(dgmath_jit_dbase, elem_data->deg_integ);
  /* DEBUG_PRINT_ARR_DBL(elem_data->J_integ, volume_nodes_integ); */
  
#if (P4EST_DIM)==3
  linalg_kron_vec_o_vec_o_vec_dot_x(integ_weights, elem_data->J_integ, elem_data->deg_integ + 1, MJac);
#elif (P4EST_DIM)==2
  linalg_kron_vec_o_vec_dot_x(integ_weights, elem_data->J_integ, elem_data->deg_integ + 1, MJac);
#else
  mpi_abort("only DIM=2 or DIM=3");
#endif  
  for (int i = 0; i < volume_nodes_integ; i++){
    elem_data->volume += MJac[i];
  }
  
  P4EST_FREE(MJac);

  ctx->sqr_nodal_stride += volume_nodes*volume_nodes;
  ctx->sqr_trace_stride += face_nodes*face_nodes*(P4EST_FACES);
  ctx->nodal_stride += volume_nodes;
  ctx->integ_stride += volume_nodes_integ;
  /* ctx->nodal_vector_stride += (P4EST_DIM)*volume_nodes; */
  /* ctx->nodal_matrix_stride += (P4EST_DIM)*(P4EST_DIM)*volume_nodes; */
  ctx->id_stride += 1;
}

/* typedef enum {INTERP_X_ON_LOBATTO, COMPUTE_DX_ON_LOBATTO, COMPUTE_DX_ON_GAUSS} dxdr_method_t; */

void
curved_element_data_compute_dxyz_drst
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 p4est_geometry_t* p4est_geom,
 int deg,
 int interp_to_Gauss, /* interp Lobatto values to Gauss values */
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* xyz_store [(P4EST_DIM)]
)
{
  double* xyz [(P4EST_DIM)];
  int volume_nodes = dgmath_get_nodes((P4EST_DIM),deg);
  if (xyz_store[0] == NULL){
    for (int i = 0; i < (P4EST_DIM); i++){
      xyz[i] = P4EST_ALLOC(double, volume_nodes);
    }
  }
  else {
    for (int i = 0; i < (P4EST_DIM); i++){
      xyz[i] = xyz_store[i];
    }    
  }

  double rst [(P4EST_DIM)];
  double xyz_i [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];
  double abc [(P4EST_DIM)];

  dgmath_rst_t rst_points
    = dgmath_get_rst_points(dgmath_jit_dbase, deg, (P4EST_DIM), LOBATTO);
  
  for (int i = 0; i < volume_nodes; i++){
    rst[0] = rst_points.r[i];
    rst[1] = rst_points.s[i];
#if (P4EST_DIM)==3
    rst[2] = rst_points.t[i];
#endif

    /* if (dxdr_method == INTERP_X_ON_LOBATTO){ */
      abc[0] = dgmath_rtox(rst_points.r[i], (double)q0[0], (double)dq)/(double)P4EST_ROOT_LEN;
      abc[1] = dgmath_rtox(rst_points.s[i], (double)q0[1], (double)dq)/(double)P4EST_ROOT_LEN;
#if (P4EST_DIM) == 3 
      abc[2] = dgmath_rtox(rst_points.t[i], (double)q0[2], (double)dq)/(double)P4EST_ROOT_LEN;
#endif        
      p4est_geom->X(p4est_geom, which_tree, abc, xyz_i);

      for (int d = 0; d < (P4EST_DIM); d++){
        xyz[d][i] = xyz_i[d];
      }
  }

  double* tmp = P4EST_ALLOC(double, volume_nodes);
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      dgmath_apply_Dij(dgmath_jit_dbase, &xyz[d][0], (P4EST_DIM), deg, d1, &dxyz_drst[d][d1][0]);
    }
  }

  if (interp_to_Gauss){
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int d1 = 0; d1 < (P4EST_DIM); d1++){
        dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dxyz_drst[d][d1][0], deg, deg, tmp, (P4EST_DIM));
        linalg_copy_1st_to_2nd(tmp, &dxyz_drst[d][d1][0], volume_nodes);
      }
    }
  }
  P4EST_FREE(tmp);
  if (xyz_store[0] == NULL){
    for (int i = 0; i < (P4EST_DIM); i++){
      P4EST_FREE(xyz[i]);
    }
  } 
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
  
  curved_element_data_compute_mortar_normal_and_sj_using_face_data
    (
     &elem_data,
     1,
     1,
     &deg,
     face,
     INTERP_X_ON_LOBATTO,
     1,
     n_on_face,
     sj_on_face,
     geom,
     dgmath_jit_dbase,
     xyz_on_face
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
curved_element_data_compute_xyz
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_geometry_t* p4est_geometry,
 int which_tree,
 int deg,
 quadrature_type_t type,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 double* xyz [(P4EST_DIM)]
)
{  
  dgmath_rst_t rst_points = dgmath_get_rst_points(dgmath_jit_dbase,
                                                  deg,
                                                  (P4EST_DIM),
                                                  type);
  
  double* rst [3] = {rst_points.r, rst_points.s, NULL};
#if (P4EST_DIM)==3
  rst[2] = rst_points.t;
#endif

  int volume_nodes = dgmath_get_nodes((P4EST_DIM), deg);
  
  double abc_i [3]; 
  double xyz_i [3];
  for (int i = 0; i < volume_nodes; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      abc_i[d] = dgmath_rtox(rst[d][i], (double)q[d], (double)dq)/(double)(P4EST_ROOT_LEN);
    }
    p4est_geometry->X(p4est_geometry, which_tree, abc_i, xyz_i);
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d][i] = xyz_i[d];
    }
  }
}

void
curved_element_data_init_new
(
 p4est_t* p4est,
 geometric_factors_t* geometric_factors,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geometry,
 curved_element_data_user_fcn_t user_fcn,
 void* user_ctx
)
{
  curved_element_data_local_sizes_t local_sizes
    = curved_element_data_compute_strides_and_sizes(p4est, dgmath_jit_dbase, d4est_geometry, user_fcn, user_ctx);

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

        elem_data->J_integ = &geometric_factors->J_integ[elem_data->integ_stride];  
        for (int i = 0; i < (P4EST_DIM); i++){
          elem_data->xyz[i] = &geometric_factors->xyz[i*local_sizes.local_nodes + elem_data->nodal_stride];
          elem_data->xyz_integ[i] = &geometric_factors->xyz_integ[i*local_sizes.local_nodes_integ + elem_data->integ_stride];
          for (int j = 0; j < (P4EST_DIM); j++){
            elem_data->xyz_rst_integ[i][j] = &geometric_factors->xyz_rst_integ[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_integ + elem_data->integ_stride];
            elem_data->rst_xyz_integ[i][j] = &geometric_factors->rst_xyz_integ[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_integ + elem_data->integ_stride];
          }
        }

        int volume_nodes_integ = dgmath_get_nodes((P4EST_DIM), elem_data->deg_integ);
        int volume_nodes= dgmath_get_nodes((P4EST_DIM), elem_data->deg);
        
        curved_element_data_compute_xyz
          (
           dgmath_jit_dbase,
           d4est_geometry->p4est_geom,
           tt,
           elem_data->deg,
           LOBATTO,
           elem_data->q,
           elem_data->dq,
           elem_data->xyz
          );

        curved_element_data_compute_xyz
          (
           dgmath_jit_dbase,
           d4est_geometry->p4est_geom,
           tt,
           elem_data->deg_integ,
           GAUSS,
           elem_data->q,
           elem_data->dq,
           elem_data->xyz_integ
          );

        /* printf("\nElement %d\n", elem_data->id); */
        /* DEBUG_PRINT_3ARR_DBL(elem_data->xyz_integ[0], elem_data->xyz_integ[1], elem_data->xyz_integ[2], volume_nodes_integ); */
        
        curved_element_data_compute_dxyz_drst
          (
           dgmath_jit_dbase,
           elem_data->q,
           elem_data->dq,
           elem_data->tree,
           d4est_geometry->p4est_geom,
           elem_data->deg_integ,
           1,
           elem_data->xyz_rst_integ,
           (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
               , NULL
#endif
               }

          );

        curved_element_data_compute_J_and_rst_xyz
          (
           elem_data->xyz_rst_integ,
           elem_data->J_integ,
           elem_data->rst_xyz_integ,
           volume_nodes_integ
          );
        
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
 
void
curved_element_data_compute_mortar_normal_and_sj_using_face_data
(
 curved_element_data_t** e,
 int num_faces_side,
 int num_faces_mortar,
 int* deg_mortar,
 int face_side,
 dxdr_method_t dxdr_method,
 int interp_to_Gauss,
 double* n [(P4EST_DIM)],
 double* sj,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* xyz_storage [(P4EST_DIM)]
)
{

  double* xyz [(P4EST_DIM)];
  if (xyz_storage[0] != NULL){
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d] = xyz_storage[d];
    }
  }
  else {
    int total_face_mortar_nodes = 0;
    for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
      total_face_mortar_nodes += dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]);
    }
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d] = P4EST_ALLOC(double, total_face_mortar_nodes);
    }
  }

  /* Calculate the four "0" corners of 
   * the mortar faces. In the case that
   * there is only one mortar face, these
   * will be the four corners of that face
   */
  
  p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)];
  
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face_side][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = dgmath_is_child_left_or_right(c, d);
      q0[j][d] = e[0]->q[d] + cd*e[0]->dq;
    }
  }

  /* Calculate the vectors that span the face 
   * there will be one in 2-D and two in 3-d */
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]);
      if (num_faces_side != num_faces_mortar)
        dqa[dir][d] /= 2;
    }
  }

  if (num_faces_side != num_faces_mortar){
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int c = 0; c < (P4EST_HALF); c++){
        q0[c][d] = q0[0][d];
        for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
          int cd = dgmath_is_child_left_or_right(c, dir);
          q0[c][d] += cd*dqa[dir][d];
        }
      }
    }
  }
  
  double* a [((P4EST_DIM)-1)];
  /* double* xyz [(P4EST_DIM)]; */
  double* dxda [(P4EST_DIM)][((P4EST_DIM)-1)];
  double dxyz_drs_i [(P4EST_DIM)][((P4EST_DIM)-1)];
  double abc [] = {0.,0.,0.};
  double xyz_i [] = {0.,0.,0.};
  int face_mortar_nodal_stride = 0;
  dgmath_rst_t rst_points;
  
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){


    int face_mortar_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]);
    /* compute the LGL nodes in the directions of the face_mortar vectors */

    double* tmp = P4EST_ALLOC(double,face_mortar_nodes);
    for (int d = 0; d < (P4EST_DIM); d++){
      /* xyz[d] = P4EST_ALLOC(double, face_mortar_nodes); */
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
        dxda[d][dir] = P4EST_ALLOC(double, face_mortar_nodes);
    }
     
    if (dxdr_method == COMPUTE_DX_ON_GAUSS)
      rst_points = dgmath_get_rst_points(dgmath_jit_dbase, deg_mortar[face_mortar], (P4EST_DIM)-1, GAUSS);
    else
      rst_points = dgmath_get_rst_points(dgmath_jit_dbase, deg_mortar[face_mortar], (P4EST_DIM)-1, LOBATTO);
      
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      if (dir == 0)
        a[dir] = rst_points.r;
      else
        a[dir] = rst_points.s;
    }

    
    for (int i = 0; i < face_mortar_nodes; i++){
      if (xyz[0] != NULL){
        for (int d = 0; d < (P4EST_DIM); d++){
          /* get "0" corner of this face_mortar */
          abc[d] = (double)q0[face_mortar][d];
       
          for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
            /* add a fraction of the face_mortar vector in direction dir
           * corresponding to the placement of the LGL node */
            double da = (a[dir][i] + 1.)/2.;
            abc[d] += da*((double)dqa[dir][d]);
            /* rs[dir] = a[dir][i]; */
            /* printf("abc[%d] = %f\n",d, abc[d]); */
          }
        
          abc[d] /= (double)(P4EST_ROOT_LEN);
        }
        /* convert vertex coords to physical coords */
        d4est_geom->p4est_geom->X(d4est_geom->p4est_geom, e[0]->tree, abc, xyz_i);     
        for (int d = 0; d < (P4EST_DIM); d++){
          xyz[d][i] = xyz_i[d];
        }
      }
      if (dxdr_method != INTERP_X_ON_LOBATTO){
        double rs [(P4EST_DIM)-1];
        for (int dir = 0; dir < (P4EST_DIM)-1; dir++){
          rs[dir] = a[dir][i];
        }
        d4est_geom->dxdr_face(d4est_geom->p4est_geom, e[0]->tree, rs, dxyz_drs_i, q0[face_mortar], dqa);      
        for (int d = 0; d < (P4EST_DIM); d++){
          for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
            dxda[d][dir][i] = dxyz_drs_i[d][dir];
          }
        }
      }
    }

    /* compute the tangent vectors in direction(s) "dir" */
    if (dxdr_method == INTERP_X_ON_LOBATTO){
      for (int d = 0; d < (P4EST_DIM); d++){
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
          dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], dir, dxda[d][dir]);
        }
      }
    }

    if (interp_to_Gauss && dxdr_method != COMPUTE_DX_ON_GAUSS){
      for (int d = 0; d < (P4EST_DIM); d++){
        if (xyz[0] != NULL){
          dgmath_interp_GLL_to_GL(dgmath_jit_dbase, xyz[d], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1);
          linalg_copy_1st_to_2nd(tmp, xyz[d], face_mortar_nodes);
        }
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
          dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxda[d][dir], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1);
          linalg_copy_1st_to_2nd(tmp, dxda[d][dir], face_mortar_nodes);
        }
      }
    }
    /* get the normal by taking the cross product of the tangent vectors
     * in 2-d, we take the cross product of the tangent vector and zhat*/
    for (int i = 0; i < face_mortar_nodes; i++){
      double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}};
      double n_i [] = {0.,0.,0.};
      for (int d = 0; d < (P4EST_DIM); d++)
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
          vecs[dir][d] = dxda[d][dir][i];

      linalg_cross_prod
        (
         vecs[0][0],
         vecs[0][1],
         vecs[0][2],
         vecs[1][0],
         vecs[1][1],
         vecs[1][2],
         &(n_i[0]),
         &(n_i[1]),
         &(n_i[2])
        );

      sj[i + face_mortar_nodal_stride] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) */
        if (face_side == 0 || face_side == 3 || face_side == 4){
          n_i[d] *= -1.;
        }
        sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d];
      }
      sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]);
      if (n[0] != NULL){
        for (int d = 0; d < (P4EST_DIM); d++){
          n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride];
        }
      }
    }

    face_mortar_nodal_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar[face_mortar]);
    for (int d = 0; d < (P4EST_DIM); d++){
      /* P4EST_FREE(xyz[d]); */
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
        P4EST_FREE(dxda[d][dir]);
    }
    P4EST_FREE(tmp);
  }

  if (xyz_storage[0] == NULL){
    for (int d = 0; d < (P4EST_DIM); d++){
      P4EST_FREE(xyz[d]);
    }
  }
  
}

/* void */
/* curved_element_data_compute_derivative_at_Gauss_nodes_on_face */
/* ( */
/*  curved_element_data_t* elem_data, */
/*  int face_side, */
/*  int deg_Gauss, */
/*  double* u, */
/*  double* du_face [(P4EST_DIM)], */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  p4est_geometry_t* p4est_geometry */
/* ) */
/* { */
/*   int volume_nodes_Lobatto = dgmath_get_nodes((P4EST_DIM), elem_data->deg); */
/*   int face_nodes_Lobatto = dgmath_get_nodes((P4EST_DIM)-1, elem_data->deg); */
/*   int volume_nodes_Gauss = dgmath_get_nodes((P4EST_DIM), deg_Gauss); */
/*   int face_nodes_Gauss = dgmath_get_nodes((P4EST_DIM)-1, deg_Gauss); */

/*   double* xyz[(P4EST_DIM)]; */
/*   double* du_di = P4EST_ALLOC(double, volume_nodes_Lobatto); */
/*   double* du_di_face = P4EST_ALLOC(double, face_nodes_Lobatto); */
/*   double* du_di_face_Gauss = P4EST_ALLOC(double, face_nodes_Gauss); */
/*   double* xyz_rst_d_d1 = P4EST_ALLOC(double, volume_nodes_Lobatto); */
/*   double* xyz_rst_face_d_d1 = P4EST_ALLOC(double, face_nodes_Lobatto); */
/*   double* xyz_rst_face_Gauss [(P4EST_DIM)][(P4EST_DIM)]; */
/*   for (int i = 0; i < (P4EST_DIM); i++) { */
/*     for (int j = 0; j < (P4EST_DIM); j++) { */
/*       xyz_rst_face_Gauss[i][j] = P4EST_ALLOC(double, face_nodes_Gauss); */
/*     } */
/*   } */
  
/*   double* r = dgmath_fetch_xyz_nd( dgmath_jit_dbase, */
/*                                   (P4EST_DIM), */
/*                                    deg_Gauss, */
/*                                    0); */
  
/*   double* s = dgmath_fetch_xyz_nd(dgmath_jit_dbase, */
/*                                   (P4EST_DIM), */
/*                                    deg_Gauss, */
/*                                    1); */
/* #if (P4EST_DIM) == 3 */
/*   double* t = dgmath_fetch_xyz_nd( dgmath_jit_dbase, */
/*                                    (P4EST_DIM), */
/*                                    deg_Gauss, */
/*                                    2); */
/* #endif */

/*   p4est_qcoord_t dq = elem_data->dq; */
/*   p4est_qcoord_t qx = elem_data->q[0]; */
/*   p4est_qcoord_t qy = elem_data->q[1]; */
/*   p4est_qcoord_t qz; */
/* #if (P4EST_DIM)==3 */
/*   qz = elem_data->q[2]; */
/* #endif */
  
/*   double abc [3]; /\* [0,1]**DIM *\/ */
/*   double xyz [3]; /\* curvilinear coordinates *\/ */
/*   int d,d1; */
/*   for (int i = 0; i < volume_nodes_Gauss; i++){ */
/*     abc[0] = dgmath_rtox(r[i], (double)qx, (double)dq)/(double)P4EST_ROOT_LEN; */
/*     abc[1] = dgmath_rtox(s[i], (double)qy, (double)dq)/(double)P4EST_ROOT_LEN; */
/* #if (P4EST_DIM) == 3  */
/*     abc[2] = dgmath_rtox(t[i], (double)qz, (double)dq)/(double)P4EST_ROOT_LEN; */
/* #endif         */
/*     p4est_geometry->X(p4est_geometry, elem_data->tree, abc, xyz); */

/*     for (d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d][i] = xyz[d]; */
/*     } */
/*   } */
  
/*   for (d = 0; d < (P4EST_DIM); d++){ */
/*     for (d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       dgmath_apply_Dij(dgmath_jit_dbase, &elem_data->xyz[d][0], (P4EST_DIM), elem_data->deg, d1, xyz_rst_d_d1); */
/*       dgmath_apply_slicer(dgmath_jit_dbase, xyz_rst_d_d1, (P4EST_DIM), face_side, elem_data->deg, xyz_rst_face_d_d1);       */
/*       dgmath_interp_GLL_to_GL(dgmath_jit_dbase, xyz_rst_face_d_d1, elem_data->deg, elem_data->deg + elem_data->deg_Gauss_diff, &xyz_rst_face_Gauss[d][d1][0], (P4EST_DIM)-1); */
/*     }     */
/*   } */

/*   for (int i = 0; i < face_nodes; i++){ */
/*     double xr = xyz_rst_face_Gauss[0][0][i]; */
/*     double xs = xyz_rst_face_Gauss[0][1][i]; */
/* #if (P4EST_DIM)==3 */
/*     double xt = xyz_rst_face_Gauss[0][2][i]; */
/* #endif */
    
/*     double yr = xyz_rst_face_Gauss[1][0][i]; */
/*     double ys = xyz_rst_face_Gauss[1][1][i]; */
/* #if (P4EST_DIM)==3 */
/*     double yt = xyz_rst_face_Gauss[1][2][i]; */
    
/*     double zr = xyz_rst_face_Gauss[2][0][i]; */
/*     double zs = xyz_rst_face_Gauss[2][1][i]; */
/*     double zt = xyz_rst_face_Gauss[2][2][i]; */
/* #endif */
    
/*     double* rx = &rst_xyz_face_Gauss[0][0][i]; */
/*     double* ry = &rst_xyz_face_Gauss[0][1][i]; */
/* #if (P4EST_DIM)==3 */
/*     double* rz = &rst_xyz_face_Gauss[0][2][i]; */
/* #endif */
/*     double* sx = &rst_xyz_face_Gauss[1][0][i]; */
/*     double* sy = &rst_xyz_face_Gauss[1][1][i]; */
/* #if (P4EST_DIM)==3 */
/*     double* sz = &rst_xyz_face_Gauss[1][2][i]; */
    
/*     double* tx = &rst_xyz_face_Gauss[2][0][i]; */
/*     double* ty = &rst_xyz_face_Gauss[2][1][i]; */
/*     double* tz = &rst_xyz_face_Gauss[2][2][i]; */
/* #endif */

/* #if (P4EST_DIM) == 3 */
/*     double J = xr*(ys*zt-zs*yt) */
/*          - yr*(xs*zt-zs*xt) */
/*          + zr*(xs*yt-ys*xt); */
/*     *rx =  (ys*zt - zs*yt)/(J); */
/*     *ry = -(xs*zt - zs*xt)/(J); */
/*     *rz =  (xs*yt - ys*xt)/(J); */
/*     *sx = -(yr*zt - zr*yt)/(J); */
/*     *sy =  (xr*zt - zr*xt)/(J); */
/*     *sz = -(xr*yt - yr*xt)/(J); */
/*     *tx =  (yr*zs - zr*ys)/(J); */
/*     *ty = -(xr*zs - zr*xs)/(J); */
/*     *tz =  (xr*ys - yr*xs)/(J); */
/* #elif (P4EST_DIM) == 2 */
/*     double J = -xs*yr + xr*ys; */
/*     *rx = ys/(J); */
/*     *sx =-yr/(J); */
/*     *ry =-xs/(J); */
/*     *sy = xr/(J); */
/* #else */
/*     mpi_abort("DIM must be 2 or 3"); */
/* #endif   */
/*   } */

/*   for (i = 0; i < (P4EST_DIM); i++){ */
/*     dgmath_apply_Dij(dgmath_jit_dbase, u, dim, deg, i, du_di); */
/*     dgmath_apply_slicer(dgmath_jit_dbase, du_di, (P4EST_DIM), face_side, elem_data->deg, du_di_face); */
/*     dgmath_interp_GLL_to_GL(dgmath_jit_dbase, du_di_face, elem_data->deg, deg_face, du_di_face_Gauss, (P4EST_DIM)-1); */
/*     for (j = 0; j < (P4EST_DIM); j++){ */
/*       for (k = 0; k < face_nodes_Gauss; k++){ */
/*         du_face[j][k] += rst_xyz_face_Gauss[i][j][k]*du_di_face_Gauss[k]; */
/*       } */
/*     }     */
/*   } */

/*   for (int i = 0; i < (P4EST_DIM); i++) { */
/*     for (int j = 0; j < (P4EST_DIM); j++) {   */
/*       P4EST_FREE(xyz_rst_face_Gauss[i][j]); */
/*     } */
/*   }   */
/*   P4EST_FREE(xyz_rst_face_d_d1); */
/*   P4EST_FREE(xyz_rst_d_d1); */
/*   P4EST_FREE(du_di_face_Gauss); */
/*   P4EST_FREE(du_di_face); */
/*   P4EST_FREE(du_di); */
/* } */

/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Lobatto, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  p4est_geometry_t* p4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */
  
/*   double* a [((P4EST_DIM)-1)]; */
/*   double* xyz [(P4EST_DIM)]; */
/*   double* dxda_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   double abc [] = {0.,0.,0.}; */
/*   double xyz_i [] = {0.,0.,0.}; */
/*   int face_mortar_nodal_stride = 0; */

/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */
/*     /\* compute the LGL nodes in the directions of the face_mortar vectors *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dxda_Lobatto[d][dir] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*         dxda_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*       } */
/*     } */

    
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_xyz_nd(dgmath_jit_dbase, */
/*                                     (P4EST_DIM) - 1, */
/*                                     deg_mortar_Lobatto[face_mortar], */
/*                                     dir);  */
/*     } */

    
/*     for (int i = 0; i < face_mortar_Lobatto_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*           /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
/*       /\* convert vertex coords to physical coords *\/ */
/*       p4est_geom->X(p4est_geom, e[0]->tree, abc, xyz_i); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         xyz[d][i] = xyz_i[d]; */
/*         /\* printf("xyz_i[d] = %f\n", xyz_i[d]); *\/ */
/*       } */
/*     } */

/*     /\* compute the tangent vectors in direction(s) "dir" *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++) */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar_Lobatto[face_mortar], dir, dxda_Lobatto[d][dir]); */
/*         dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxda_Lobatto[d][dir], deg_mortar_Lobatto[face_mortar], deg_mortar_Gauss[face_mortar], dxda_Gauss[d][dir], (P4EST_DIM)-1); */

/*         /\* dgmath_apply_GaussDij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], deg_mortar[face_mortar], dir, dxda_Gauss[d][dir]); *\/ */
        
/*       } */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_nodal_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride]; */
/*       } */
/*     } */

/*     face_mortar_nodal_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face_mortar]); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz[d]); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_Gauss[d][dir]); */
/*         P4EST_FREE(dxda_Lobatto[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */



/* void */
/* curved_element_data_compute_mortar_normal_and_sj */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  dxdr_method_t dxdr_method, */
/*  int deg, */
/*  int interp_to_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  p4est_geometry_t* p4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */
  
/*   double* a [((P4EST_DIM)-1)]; */
/*   double* xyz [(P4EST_DIM)]; */
/*   double* dxda_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   double abc [] = {0.,0.,0.}; */
/*   double xyz_i [] = {0.,0.,0.}; */
/*   int face_mortar_nodal_stride = 0; */

/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */
/*     /\* compute the LGL nodes in the directions of the face_mortar vectors *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dxda_Lobatto[d][dir] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*         dxda_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*       } */
/*     } */

    
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_xyz_nd(dgmath_jit_dbase, */
/*                                     (P4EST_DIM) - 1, */
/*                                     deg_mortar_Lobatto[face_mortar], */
/*                                     dir);  */
/*     } */

    
/*     for (int i = 0; i < face_mortar_Lobatto_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*         /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
/*       /\* convert vertex coords to physical coords *\/ */
/*       p4est_geom->X(p4est_geom, e[0]->tree, abc, xyz_i); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         xyz[d][i] = xyz_i[d]; */
/*         /\* printf("xyz_i[d] = %f\n", xyz_i[d]); *\/ */
/*       } */
/*     } */

/*     /\* compute the tangent vectors in direction(s) "dir" *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++) */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar_Lobatto[face_mortar], dir, dxda_Lobatto[d][dir]); */
/*         dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxda_Lobatto[d][dir], deg_mortar_Lobatto[face_mortar], deg_mortar_Gauss[face_mortar], dxda_Gauss[d][dir], (P4EST_DIM)-1); */

/*         /\* dgmath_apply_GaussDij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], deg_mortar[face_mortar], dir, dxda_Gauss[d][dir]); *\/ */
        
/*       } */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_nodal_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride]; */
/*       } */
/*     } */

/*     face_mortar_nodal_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face_mortar]); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz[d]); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_Gauss[d][dir]); */
/*         P4EST_FREE(dxda_Lobatto[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */
    


/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_use_deriv */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Lobatto, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  d4est_geometry_t* d4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */

/*   double dqa_norm = 0.; */
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     double dqad = (double)dqa[0][d]/(double)(P4EST_ROOT_LEN); */
/*     dqa_norm += dqad*dqad; */
/*   } */
/*   dqa_norm = sqrt(dqa_norm); */

  
/*   double* a [((P4EST_DIM)-1)]; */
/*   double* xyz [(P4EST_DIM)]; */
/*   double* dxda_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   double abc [] = {0.,0.,0.}; */
/*   double xyz_i [] = {0.,0.,0.}; */
/*   int face_mortar_nodal_stride = 0; */

/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */
/*     /\* compute the LGL nodes in the directions of the face_mortar vectors *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dxda_Lobatto[d][dir] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*         dxda_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*       } */
/*     } */

    
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_xyz_nd(dgmath_jit_dbase, */
/*                                     (P4EST_DIM) - 1, */
/*                                     deg_mortar_Lobatto[face_mortar], */
/*                                     dir);  */
/*     } */

    
/*     for (int i = 0; i < face_mortar_Lobatto_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*         /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
/*       /\* convert vertex coords to physical coords *\/ */
/*       p4est_geom->X(p4est_geom, e[0]->tree, abc, xyz_i); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         xyz[d][i] = xyz_i[d]; */
/*         /\* printf("xyz_i[d] = %f\n", xyz_i[d]); *\/ */
/*       } */

/*       d4est_geom->dxda */
      
/*     } */

/*     /\* compute the tangent vectors in direction(s) "dir" *\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++) */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar_Lobatto[face_mortar], dir, dxda_Lobatto[d][dir]); */
/*         dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxda_Lobatto[d][dir], deg_mortar_Lobatto[face_mortar], deg_mortar_Gauss[face_mortar], dxda_Gauss[d][dir], (P4EST_DIM)-1); */

/*         /\* dgmath_apply_GaussDij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], deg_mortar[face_mortar], dir, dxda_Gauss[d][dir]); *\/ */
        
/*       } */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_nodal_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride]; */
/*       } */
/*     } */

/*     face_mortar_nodal_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face_mortar]); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz[d]); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_Gauss[d][dir]); */
/*         P4EST_FREE(dxda_Lobatto[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */

/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_no_interp */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  d4est_geometry_t* d4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase, */
/*  int store_xyz, */
/*  double* xyz_on_mortar_Gauss [(P4EST_DIM)] */
/* ) */
/* { */
/*   /\* printf("\n\nNEW NORMAL CALCULATION\n\n"); *\/ */

  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
/*   double dqa_norm = 0.; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */

/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     double dqad = (double)dqa[0][d]/(double)(P4EST_ROOT_LEN); */
/*     dqa_norm += dqad*dqad; */
/*   } */
/*   dqa_norm = sqrt(dqa_norm); */
  
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][0]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][1]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL(dqa_norm);  *\/ */

  
/*   double* a [((P4EST_DIM)-1)]; */
/*   double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   double abc [] = {0.,0.,0.}; */
/*   int face_mortar_Gauss_stride = 0; */

/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     /\* int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); *\/ */
/*     /\* compute the LGL nodes in the directions of the face_mortar vectors x*\/ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dxda_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*       } */
/*     } */

    
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, */
/*                                           (P4EST_DIM) - 1, */
/*                                           deg_mortar_Gauss[face_mortar], */
/*                                           dir);  */
/*     } */

    
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*         /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
      
/*       double dxyz_dabc [(P4EST_DIM)][(P4EST_DIM)]; */

/*       d4est_geom->dxda(d4est_geom->p4est_geom, e[0]->tree, abc, dxyz_dabc);     */

/*       /\* printf("dxyz_dabc[0][0] = %f\n", dxyz_dabc[0][0]); *\/ */
      
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           dxda_Gauss[d][dir][i] = 0.; */
/*           for (int k = 0; k < (P4EST_DIM); k++){ */
/*             //directional derivative on face Gauss points in direction of dqa[dir] */
/*             dxda_Gauss[d][dir][i] += dxyz_dabc[d][k]*((double)dqa[dir][k]/(double)(P4EST_ROOT_LEN))/dqa_norm; */
/*             //jacobian between logical points and GL points */
/*             /\* printf("dxda_Gauss[d][dir][i] = %f\n",dxda_Gauss[d][dir][i]); *\/ */
/*           } */
/*           dxda_Gauss[d][dir][i] *= dqa_norm; */
/*           dxda_Gauss[d][dir][i] *= .5; */
/*         } */
/*       } */

/*       if (store_xyz){         */
/*         double xyz_i [(P4EST_DIM)]; */
/*         d4est_geom->p4est_geom->X(d4est_geom->p4est_geom, e[0]->tree, abc, xyz_i); */
/*         for (int d = 0; d < (P4EST_DIM); d++) { */
/*           xyz_on_mortar_Gauss[d][i+face_mortar_Gauss_stride] = xyz_i[d]; */
/*         } */
/*       }            */
/*     } */

/*     /\* DEBUG_PRINT_2ARR_DBL(dxda_Gauss[0][0], dxda_Gauss[1][0], face_mortar_Gauss_nodes); *\/ */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_Gauss_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_Gauss_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_Gauss_stride] = sqrt(sj[i + face_mortar_Gauss_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_Gauss_stride] = n_i[d]/sj[i + face_mortar_Gauss_stride]; */
/*       } */
/*     } */

/*     face_mortar_Gauss_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_Gauss[face_mortar]); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_Gauss[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */

/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_using_interp */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Lobatto, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  p4est_geometry_t* p4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*     } */
/*   } */
  
/*   double* a [((P4EST_DIM)-1)]; */
/*   double* xyz [(P4EST_DIM)]; */
/*   double* dxda_side_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_mortar_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_mortar_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   double abc [] = {0.,0.,0.}; */
/*   double xyz_i [] = {0.,0.,0.}; */



/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     xyz[d] = P4EST_ALLOC(double, face_side_Lobatto_nodes); */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dxda_side_Lobatto[d][dir] = P4EST_ALLOC(double, face_side_Lobatto_nodes); */
/*       dxda_mortar_Lobatto[d][dir] = P4EST_ALLOC(double, face_side_Lobatto_nodes); */
/*       dxda_mortar_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*     } */
/*   } */

/*   int face_side_Lobatto_stride = 0; */
/*   for (int face_side = 0; face_side < faces_side; face_side++){ */
      
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_xyz_nd(dgmath_jit_dbase, */
/*                                     (P4EST_DIM) - 1, */
/*                                     deg_side_Lobatto[face_side], */
/*                                     dir);  */
/*     } */


/*     for (int i = 0; i < face_side_Lobatto_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_side][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*           /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
/*       /\* convert vertex coords to physical coords *\/ */
/*       p4est_geom->X(p4est_geom, e[0]->tree, abc, xyz_i); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         xyz[d][i] = xyz_i[d]; */
/*         /\* printf("xyz_i[d] = %f\n", xyz_i[d]); *\/ */
/*       } */
/*     } */

/*     /\* compute the tangent vectors in direction(s) "dir" *\/ */
/*     dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar_Lobatto[face_side], dir, &dxda_side_Lobatto[d][dir][face_side_Lobatto_stride]); */
/*     face_side_Lobatto_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_side_Lobatto[face_side]); */
/*   } */

/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dgmath_project_side_onto_mortar_space(dgmath_jit_dbase, dxda_side_Lobatto[d][dir], faces_side, deg_side_Lobatto, dxda_mortar_Lobatto[d][dir], faces_mortar, deg_mortar_Lobatto);         */
/*     } */
/*   } */

/*   int face_mortar_Gauss_stride = 0; */
/*   int face_mortar_Lobatto_stride = 0; */
/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */

/*     dgmath_interp_GLL_to_GL */
/*       ( */
/*        dgmath_jit_dbase, */
/*        &dxda_mortar_Lobatto[d][dir][face_mortar_Lobatto_stride], */
/*        deg_mortar_Lobatto[face_mortar], */
/*        deg_mortar_Gauss[face_mortar], */
/*        dxda_Gauss[d][dir], */
/*        (P4EST_DIM)-1 */
/*       ); */
    
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_Gauss_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_Gauss_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_Gauss_stride] = sqrt(sj[i + face_mortar_Gauss_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_Gauss_stride] = n_i[d]/sj[i + face_mortar_Gauss_stride]; */
/*       } */
/*     } */
      
/*     face_mortar_Gauss_stride += face_mortar_Gauss_nodes; */
/*     face_mortar_Lobatto_stride += face_mortar_Lobatto_nodes; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz[d]); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_mortar_Gauss[d][dir]); */
/*         P4EST_FREE(dxda_side_Lobatto[d][dir]); */
/*         P4EST_FREE(dxda_mortar_Lobatto[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */


/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_using_interp_cheap */
/* ( */
/*  curved_element_data_t** e_m, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Lobatto, */
/*  int* deg_mortar_Gauss, */
/*  int f_m, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  p4est_geometry_t* p4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
  
/*   /\* double* a [((P4EST_DIM)-1)]; *\/ */
/*   double* xyz_face [(P4EST_DIM)]; */
/*   double* dxda_side_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_mortar_Lobatto [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double* dxda_mortar_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; */

/*   int deg_side_Lobatto [P4EST_HALF]; */
/*   int face_mortar_Gauss_nodes = 0; */
/*   int face_mortar_Lobatto_nodes = 0; */
/*   int face_side_Lobatto_nodes = 0; */
  
/*   for (int face_side = 0; face_side < num_faces_side; face_side++){ */
/*     deg_side_Lobatto[face_side] = e_m[face_side]->deg; */
/*     face_side_Lobatto_nodes += dgmath_get_nodes((P4EST_DIM) - 1, deg_side_Lobatto[face_side]); */
/*   } */
  
/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
/*     face_mortar_Gauss_nodes += dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     face_mortar_Lobatto_nodes += dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */
/*   } */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     xyz_face[d] = P4EST_ALLOC(double, face_side_Lobatto_nodes); */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dxda_side_Lobatto[d][dir] = P4EST_ALLOC(double, face_side_Lobatto_nodes); */
/*       dxda_mortar_Lobatto[d][dir] = P4EST_ALLOC(double, face_mortar_Lobatto_nodes); */
/*       dxda_mortar_Gauss[d][dir] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*     } */
/*   } */


/*   /\* printf("NEW SIDE\n"); *\/ */
/*   int face_side_Lobatto_stride = 0; */
/*   for (int face_side = 0; face_side < num_faces_side; face_side++){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       dgmath_apply_slicer(dgmath_jit_dbase, e_m[face_side]->xyz[d], (P4EST_DIM), f_m, e_m[face_side]->deg, xyz_face[d]); */
/*       /\* DEBUG_PRINT_ARR_DBL(xyz_face[d], dgmath_get_nodes((P4EST_DIM)-1, e_m[face_side]->deg)); *\/ */
/*       for (int dir = 0; dir < (P4EST_DIM)-1; dir++){ */
/*         dgmath_apply_Dij(dgmath_jit_dbase, xyz_face[d], ((P4EST_DIM))-1, deg_side_Lobatto[face_side], dir, &dxda_side_Lobatto[d][dir][face_side_Lobatto_stride]); */
/*         /\* printf("face_side = %d, d, dir = %d,%d\n", face_side,d, dir); *\/ */
/*         /\* double* aliastmp = &dxda_side_Lobatto[d][dir][face_side_Lobatto_stride]; *\/ */
/*         /\* DEBUG_PRINT_ARR_DBL(aliastmp, face_side_Lobatto_nodes); *\/ */
/*       } */
/*     } */
/*     face_side_Lobatto_stride += dgmath_get_nodes((P4EST_DIM)-1, e_m[face_side]->deg); */
/*   } */



/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dgmath_project_side_onto_mortar_space(dgmath_jit_dbase, */
/*                                             dxda_side_Lobatto[d][dir], */
/*                                             num_faces_side, */
/*                                             &deg_side_Lobatto[0], */
/*                                             dxda_mortar_Lobatto[d][dir], */
/*                                             num_faces_mortar, */
/*                                             deg_mortar_Lobatto); */
/*       /\* printf("d,dir = %d,%d\n", d, dir); *\/ */
/*       /\* DEBUG_PRINT_ARR_DBL(dxda_mortar_Lobatto[d][dir], face_mortar_Lobatto_nodes); *\/ */
/*     } */
/*   } */

/*   int face_mortar_Gauss_stride = 0; */
/*   int face_mortar_Lobatto_stride = 0; */
/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */

/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         dgmath_interp_GLL_to_GL */
/*           ( */
/*            dgmath_jit_dbase, */
/*            &dxda_mortar_Lobatto[d][dir][face_mortar_Lobatto_stride], */
/*            deg_mortar_Lobatto[face_mortar], */
/*            deg_mortar_Gauss[face_mortar], */
/*            dxda_mortar_Gauss[d][dir], */
/*            (P4EST_DIM)-1 */
/*           ); */
/*       } */
/*     } */
    
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */
/*     int face_mortar_Lobatto_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Lobatto[face_mortar]); */

/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda_mortar_Gauss[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_Gauss_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (f_m == 0 || f_m == 3 || f_m == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_Gauss_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_Gauss_stride] = sqrt(sj[i + face_mortar_Gauss_stride]); */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][i + face_mortar_Gauss_stride] = n_i[d]/sj[i + face_mortar_Gauss_stride]; */
/*       } */
/*     } */
      
/*     face_mortar_Gauss_stride += face_mortar_Gauss_nodes; */
/*     face_mortar_Lobatto_stride += face_mortar_Lobatto_nodes; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz_face[d]); */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*         P4EST_FREE(dxda_mortar_Gauss[d][dir]); */
/*         P4EST_FREE(dxda_side_Lobatto[d][dir]); */
/*         P4EST_FREE(dxda_mortar_Lobatto[d][dir]); */
/*       } */
/*     } */
/*   } */
/* } */

void
curved_element_data_compute_drst_dxyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)],
 int nodes
)
{
 for (int i = 0; i < nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
#endif
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];
#if (P4EST_DIM)==3
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif
    
    double* rx = &drst_dxyz[0][0][i];
    double* ry = &drst_dxyz[0][1][i];
#if (P4EST_DIM)==3
    double* rz = &drst_dxyz[0][2][i];
#endif
    double* sx = &drst_dxyz[1][0][i];
    double* sy = &drst_dxyz[1][1][i];
#if (P4EST_DIM)==3
    double* sz = &drst_dxyz[1][2][i];
    
    double* tx = &drst_dxyz[2][0][i];
    double* ty = &drst_dxyz[2][1][i];
    double* tz = &drst_dxyz[2][2][i];
#endif

#if (P4EST_DIM) == 3
    double J = xr*(ys*zt-zs*yt)
                 - yr*(xs*zt-zs*xt)
                 + zr*(xs*yt-ys*xt);
    *rx =  (ys*zt - zs*yt)/(J);
    *ry = -(xs*zt - zs*xt)/(J);
    *rz =  (xs*yt - ys*xt)/(J);
    *sx = -(yr*zt - zr*yt)/(J);
    *sy =  (xr*zt - zr*xt)/(J);
    *sz = -(xr*yt - yr*xt)/(J);
    *tx =  (yr*zs - zr*ys)/(J);
    *ty = -(xr*zs - zr*xs)/(J);
    *tz =  (xr*ys - yr*xs)/(J);
#elif (P4EST_DIM) == 2
    double J = -xs*yr + xr*ys;
    *rx = ys/(J);
    *sx =-yr/(J);
    *ry =-xs/(J);
    *sy = xr/(J);
#else
    mpi_abort("DIM must be 2 or 3");
#endif 
  }  
}



void
curved_element_data_compute_J_and_rst_xyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* jac,
 double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)],
 int volume_nodes
)
{

  for (int i = 0; i < volume_nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
#endif
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];
#if (P4EST_DIM)==3
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif
    
    double* rx = &drst_dxyz[0][0][i];
    double* ry = &drst_dxyz[0][1][i];
#if (P4EST_DIM)==3
    double* rz = &drst_dxyz[0][2][i];
#endif
    double* sx = &drst_dxyz[1][0][i];
    double* sy = &drst_dxyz[1][1][i];
#if (P4EST_DIM)==3
    double* sz = &drst_dxyz[1][2][i];
    
    double* tx = &drst_dxyz[2][0][i];
    double* ty = &drst_dxyz[2][1][i];
    double* tz = &drst_dxyz[2][2][i];
#endif

    double* J = &jac[i];
    
#if (P4EST_DIM) == 3
    *J = xr*(ys*zt-zs*yt)
                - yr*(xs*zt-zs*xt)
                + zr*(xs*yt-ys*xt);
    *rx =  (ys*zt - zs*yt)/(*J);
    *ry = -(xs*zt - zs*xt)/(*J);
    *rz =  (xs*yt - ys*xt)/(*J);
    *sx = -(yr*zt - zr*yt)/(*J);
    *sy =  (xr*zt - zr*xt)/(*J);
    *sz = -(xr*yt - yr*xt)/(*J);
    *tx =  (yr*zs - zr*ys)/(*J);
    *ty = -(xr*zs - zr*xs)/(*J);
    *tz =  (xr*ys - yr*xs)/(*J);
#elif (P4EST_DIM) == 2
    *J = -xs*yr + xr*ys;
    *rx = ys/(*J);
    *sx =-yr/(*J);
    *ry =-xs/(*J);
    *sy = xr/(*J);
#else
    mpi_abort("DIM must be 2 or 3");
#endif 
  }  
}




 
/* void */
/* curved_element_data_compute_J_and_rst_xyz */
/* ( */
/*  double dxyz_drst [(P4EST_DIM)][(P4EST_DIM)], */
/*  double* J, */
/*  double drst_dxyz [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*     double xr = dxyz_drst[0][0]; */
/*     double xs = dxyz_drst[0][1]; */
/* #if (P4EST_DIM)==3 */
/*     double xt = dxyz_drst[0][2]; */
/* #endif */
    
/*     double yr = dxyz_drst[1][0]; */
/*     double ys = dxyz_drst[1][1]; */
/* #if (P4EST_DIM)==3 */
/*     double yt = dxyz_drst[1][2]; */
    
/*     double zr = dxyz_drst[2][0]; */
/*     double zs = dxyz_drst[2][1]; */
/*     double zt = dxyz_drst[2][2]; */
/* #endif */
    
/*     double* rx = &drst_dxyz[0][0]; */
/*     double* ry = &drst_dxyz[0][1]; */
/* #if (P4EST_DIM)==3 */
/*     double* rz = &drst_dxyz[0][2]; */
/* #endif */
/*     double* sx = &drst_dxyz[1][0]; */
/*     double* sy = &drst_dxyz[1][1]; */
/* #if (P4EST_DIM)==3 */
/*     double* sz = &drst_dxyz[1][2]; */
    
/*     double* tx = &drst_dxyz[2][0]; */
/*     double* ty = &drst_dxyz[2][1]; */
/*     double* tz = &drst_dxyz[2][2]; */
/* #endif */

/* #if (P4EST_DIM) == 3 */
/*     *J = xr*(ys*zt-zs*yt) */
/*          - yr*(xs*zt-zs*xt) */
/*          + zr*(xs*yt-ys*xt); */
/*     *rx =  (ys*zt - zs*yt)/(*J); */
/*     *ry = -(xs*zt - zs*xt)/(*J); */
/*     *rz =  (xs*yt - ys*xt)/(*J); */
/*     *sx = -(yr*zt - zr*yt)/(*J); */
/*     *sy =  (xr*zt - zr*xt)/(*J); */
/*     *sz = -(xr*yt - yr*xt)/(*J); */
/*     *tx =  (yr*zs - zr*ys)/(*J); */
/*     *ty = -(xr*zs - zr*xs)/(*J); */
/*     *tz =  (xr*ys - yr*xs)/(*J); */
/* #elif (P4EST_DIM) == 2 */
/*     *J = -xs*yr + xr*ys; */
/*     *rx = ys/(*J); */
/*     *sx =-yr/(*J); */
/*     *ry =-xs/(*J); */
/*     *sy = xr/(*J); */
/* #else */
/*     mpi_abort("DIM must be 2 or 3"); */
/* #endif  */
/* } */



/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_no_interp_usingJ */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  d4est_geometry_t* d4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
/*   /\* printf("\n\nNEW NORMAL CALCULATION\n\n"); *\/ */

  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
/*   double dqa_norm = 0.; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */

/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     double dqad = (double)dqa[0][d]/(double)(P4EST_ROOT_LEN); */
/*     dqa_norm += dqad*dqad; */
/*   } */
/*   dqa_norm = sqrt(dqa_norm); */
  
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][0]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][1]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL(dqa_norm);  *\/ */

  
/*   double* a [((P4EST_DIM)-1)]; */
/*   /\* double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; *\/ */

/*   double abc [] = {0.,0.,0.}; */
/*   int face_mortar_Gauss_stride = 0; */

/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */

/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, */
/*                                           (P4EST_DIM) - 1, */
/*                                           deg_mortar_Gauss[face_mortar], */
/*                                           dir);  */
/*     } */

    
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*           /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
      
/*       double dxyz_dabc [(P4EST_DIM)][(P4EST_DIM)]; */
/*       double drst_dxyz [(P4EST_DIM)][(P4EST_DIM)]; */
/*       double J; */


/*       d4est_geom->dxda(d4est_geom->p4est_geom, e[0]->tree, abc, dxyz_dabc);     */
/*       for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*         for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*           dxyz_dabc[d1][d2] *= dqa_norm; */
/*           dxyz_dabc[d1][d2] *= .5; */
/*         } */
/*       } */
        
/*       curved_element_data_compute_J_and_rst_xyz(dxyz_dabc, &J, drst_dxyz);     */
      
/*       int i0 = -1; */
/*       if (face_side == 0 || face_side == 1){ */
/*         i0 = 0; */
/*       } */
/*       else if (face_side == 2 || face_side == 3){ */
/*         i0 = 1; */
/*       } */
/*       else if (face_side == 4 || face_side == 5){ */
/*         i0 = 2; */
/*       } */
/*       else { */
/*         mpi_abort("face_side must be < 6\n"); */
/*       } */
  
/*       double sgn = (face_side == 0 || face_side == 2 || face_side == 4) ? -1. : 1.; */
  
/*       int ks = i + face_mortar_Gauss_stride; */
/*       sj[ks] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][ks] = sgn*drst_dxyz[i0][d]*J; */
/*         sj[ks] += n[d][ks]*n[d][ks]; */
/*       } */
/*       sj[ks] = sqrt(sj[ks]); */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         n[d][ks] /= sj[ks]; */
/*     } */
/*     face_mortar_Gauss_stride += face_mortar_Gauss_nodes;     */
/*   } */
/* } */

/* void */
/* curved_element_data_compute_mortar_normal_and_sj_using_face_data_at_Gauss_nodes_usingJ_Lobatto */
/* ( */
/*  curved_element_data_t** e, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar_Gauss, */
/*  int face_side, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  d4est_geometry_t* d4est_geom, */
/*  dgmath_jit_dbase_t* dgmath_jit_dbase */
/* ) */
/* { */
/*   /\* printf("\n\nNEW NORMAL CALCULATION\n\n"); *\/ */

  
/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = dgmath_is_child_left_or_right(c, d); */
/*       q0[j][d] = e[0]->q[d] + cd*e[0]->dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
/*   double dqa_norm = 0.; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = dgmath_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */

/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     double dqad = (double)dqa[0][d]/(double)(P4EST_ROOT_LEN); */
/*     dqa_norm += dqad*dqad; */
/*   } */
/*   dqa_norm = sqrt(dqa_norm); */
  
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][0]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL((double)dqa[0][1]/(double)(P4EST_ROOT_LEN));  *\/ */
/*   /\* DEBUG_PRINT_DBL(dqa_norm);  *\/ */

  
/*   double* a [((P4EST_DIM)-1)]; */
/*   /\* double* dxda_Gauss [(P4EST_DIM)][((P4EST_DIM)-1)]; *\/ */

/*   double abc [] = {0.,0.,0.}; */
/*   int face_mortar_Gauss_stride = 0; */

/*   double* dxyz_drst_Lobatto [(P4EST_DIM)][(P4EST_DIM)]; */
/*   double* dxyz_drst_Gauss [(P4EST_DIM)][(P4EST_DIM)]; */
  
/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
  
/*     int face_mortar_Gauss_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_Gauss[face_mortar]); */

/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       a[dir] =  dgmath_fetch_Gauss_xyz_nd(dgmath_jit_dbase, */
/*                                           (P4EST_DIM) - 1, */
/*                                           deg_mortar_Gauss[face_mortar], */
/*                                           dir);  */
/*     } */

    
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         dxyz_drst_Lobatto[d1][d2] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*         dxyz_drst_Gauss[d1][d2] = P4EST_ALLOC(double, face_mortar_Gauss_nodes); */
/*       } */
/*     } */
    
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* get "0" corner of this face_mortar *\/ */
/*         abc[d] = (double)q0[face_mortar][d]; */
       
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*           double da = (a[dir][i] + 1.)/2.; */
/*           abc[d] += da*((double)dqa[dir][d]); */
/*           /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*         } */
/*         abc[d] /= (double)(P4EST_ROOT_LEN); */
/*       } */
      
/*       double dxyz_dabc [(P4EST_DIM)][(P4EST_DIM)]; */



/*       d4est_geom->dxda(d4est_geom->p4est_geom, e[0]->tree, abc, dxyz_dabc);     */
/*       for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*         for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*           dxyz_drst_Lobatto[d1][d2][i] = dxyz_dabc[d1][d2]; */
/*           dxyz_drst_Lobatto[d1][d2][i] *= dqa_norm; */
/*           dxyz_drst_Lobatto[d1][d2][i] *= .5; */
/*         } */
/*       } */
/*     } */

/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxyz_drst_Lobatto[d1][d2], deg_mortar_Gauss[face_mortar], deg_mortar_Gauss[face_mortar], dxyz_drst_Gauss[d1][d2], (P4EST_DIM)-1); */
/*       } */
/*     } */
    
/*     for (int i = 0; i < face_mortar_Gauss_nodes; i++){ */

/*       double drst_dxyz [(P4EST_DIM)][(P4EST_DIM)]; */
/*       double dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]; */
/*       double J; */
      
/*       for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*         for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*           dxyz_drst[d1][d2] = dxyz_drst_Gauss[d1][d2][i]; */
/*         } */
/*       } */
/*       curved_element_data_compute_J_and_rst_xyz(dxyz_drst, &J, drst_dxyz);     */
      
/*       int i0 = -1; */
/*       if (face_side == 0 || face_side == 1){ */
/*         i0 = 0; */
/*       } */
/*       else if (face_side == 2 || face_side == 3){ */
/*         i0 = 1; */
/*       } */
/*       else if (face_side == 4 || face_side == 5){ */
/*         i0 = 2; */
/*       } */
/*       else { */
/*         mpi_abort("face_side must be < 6\n"); */
/*       } */
  
/*       double sgn = (face_side == 0 || face_side == 2 || face_side == 4) ? -1. : 1.; */
  
/*       int ks = i + face_mortar_Gauss_stride; */
/*       sj[ks] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         n[d][ks] = sgn*drst_dxyz[i0][d]*J; */
/*         sj[ks] += n[d][ks]*n[d][ks]; */
/*       } */
/*       sj[ks] = sqrt(sj[ks]); */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         n[d][ks] /= sj[ks]; */
/*     } */

/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         P4EST_FREE(dxyz_drst_Lobatto[d1][d2]); */
/*         P4EST_FREE(dxyz_drst_Gauss[d1][d2]); */
/*       } */
/*     } */
    
/*     face_mortar_Gauss_stride += face_mortar_Gauss_nodes;     */
/*   } */
/* } */

void
curved_element_data_compute_surface_jacobian_and_normal_from_rst_xyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* n [(P4EST_DIM)],
 double* sj,
 int face,
 int deg,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), deg);
  double* J = P4EST_ALLOC(double, volume_nodes);
  double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_MAT(drst_dxyz, P4EST_DIM, P4EST_DIM, volume_nodes);
  
  curved_element_data_compute_J_and_rst_xyz
    (
     dxyz_drst,
     J,
     drst_dxyz,
     volume_nodes
    );

  curved_element_data_compute_surface_jacobian_and_normal
    (
     drst_dxyz,
     J,
     n,
     sj,
     (P4EST_DIM),
     face,
     deg,
     dgmath_jit_dbase
    );

  D4EST_FREE_MAT(drst_dxyz, P4EST_DIM, P4EST_DIM);
  P4EST_FREE(J);  
}


void
curved_element_data_compute_surface_jacobian_and_normal_from_rst_xyz_interp_to_Gauss
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* n [(P4EST_DIM)],
 double* sj,
 int face,
 int deg,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  int face_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg);
  double* J_on_face_Gauss;
  double* drst_dxyz_on_face_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  double* dxyz_drst_on_face [(P4EST_DIM)][(P4EST_DIM)];
  double* dxyz_drst_on_face_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  
  J_on_face_Gauss = P4EST_ALLOC(double, face_nodes);
  D4EST_ALLOC_MAT(drst_dxyz_on_face_Gauss, P4EST_DIM, P4EST_DIM, face_nodes);
  D4EST_ALLOC_MAT(dxyz_drst_on_face_Gauss, P4EST_DIM, P4EST_DIM, face_nodes);
  D4EST_ALLOC_MAT(dxyz_drst_on_face, P4EST_DIM, P4EST_DIM, face_nodes);
  
  int i0 = -1;
  
  if (face == 0 || face == 1){
    i0 = 0;
  }
  else if (face == 2 || face == 3){
    i0 = 1;
  }
  else if (face == 4 || face == 5){
    i0 = 2;
  }
  else {
    mpi_abort("face must be < 6\n");
  }
  

  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         dxyz_drst[i][j],
         (P4EST_DIM),
         face,
         deg,
         dxyz_drst_on_face[i][j]
        );

      dgmath_interp_GLL_to_GL(dgmath_jit_dbase,
                              dxyz_drst_on_face[i][j],
                              deg,
                              deg,
                              dxyz_drst_on_face_Gauss[i][j],
                              (P4EST_DIM)-1);
    }

  curved_element_data_compute_J_and_rst_xyz
    (
     dxyz_drst_on_face_Gauss,
     J_on_face_Gauss,
     drst_dxyz_on_face_Gauss,
     face_nodes
    );

    
  double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
  
  for (int i = 0; i < face_nodes; i++){
    sj[i] = 0.;
    for (int d = 0; d < (P4EST_DIM); d++){
      n[d][i] = sgn*drst_dxyz_on_face_Gauss[i0][d][i]*J_on_face_Gauss[i];
      /* if (n[d][i] == 0) */
      /* continue; */
      sj[i] += n[d][i]*n[d][i];
    }
    sj[i] = sqrt(sj[i]);
    for (int d = 0; d < (P4EST_DIM); d++)
      n[d][i] /= sj[i];
  }

  D4EST_FREE_MAT(drst_dxyz_on_face_Gauss, P4EST_DIM, P4EST_DIM);
  D4EST_FREE_MAT(dxyz_drst_on_face_Gauss, P4EST_DIM, P4EST_DIM);
  D4EST_FREE_MAT(dxyz_drst_on_face, P4EST_DIM, P4EST_DIM);
  P4EST_FREE(J_on_face_Gauss);
}


void
curved_element_data_compute_surface_jacobian_and_normal
(
 double* rst_xyz [(P4EST_DIM)][(P4EST_DIM)],
 double* J,
 double* n [(P4EST_DIM)],
 double* sj,
 int dim,
 int face,
 int deg,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{
  int i0 = -1;
  
  if (face == 0 || face == 1){
    i0 = 0;
  }
  else if (face == 2 || face == 3){
    i0 = 1;
  }
  else if (face == 4 || face == 5){
    i0 = 2;
  }
  else {
    mpi_abort("face must be < 6\n");
  }
  
  int face_nodes = dgmath_get_nodes(dim - 1, deg);
  double* rst_xyz_on_face [3][3];
  double* J_on_face = P4EST_ALLOC(double, face_nodes);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++){
      rst_xyz_on_face[i][j] = P4EST_ALLOC(double, face_nodes);
      dgmath_apply_slicer
        (
         dgmath_jit_dbase,
         rst_xyz[i][j],
         dim,
         face,
         deg,
         rst_xyz_on_face[i][j]
        );
    }

  dgmath_apply_slicer
    (
     dgmath_jit_dbase,
     J,
     dim,
     face,
     deg,
     J_on_face
    );
    
  double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
  
  for (int i = 0; i < face_nodes; i++){
    sj[i] = 0.;
    for (int d = 0; d < dim; d++){
      n[d][i] = sgn*rst_xyz_on_face[i0][d][i]*J_on_face[i];
      /* if (n[d][i] == 0) */
      /* continue; */
      sj[i] += n[d][i]*n[d][i];
    }
    sj[i] = sqrt(sj[i]);
    for (int d = 0; d < dim; d++)
      n[d][i] /= sj[i];
  }

  P4EST_FREE(J_on_face);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      P4EST_FREE(rst_xyz_on_face[i][j]);
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

void
curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
(
 curved_element_data_t** e,
 int num_faces_side,
 int num_faces_mortar,
 int* deg_mortar_integ,
 int face,
 double* drst_dxyz_on_mortar_Gauss [(P4EST_DIM)][(P4EST_DIM)],
 double* sj_on_mortar_Gauss,
 double* n_on_mortar_Gauss [(P4EST_DIM)],
 p4est_geometry_t* p4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* j_div_sj_mortar_Gauss
)
{
  double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)];
  double* dxyz_drst_on_face_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  int max_deg = 0;
  for (int i = 0; i < num_faces_mortar; i++){
    max_deg = (deg_mortar_integ[i] > max_deg) ? deg_mortar_integ[i] : max_deg;
  }
  int volume_nodes_max = dgmath_get_nodes((P4EST_DIM), max_deg);
  int face_nodes_max = dgmath_get_nodes((P4EST_DIM)-1, max_deg);
  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      dxyz_drst[i][j] = P4EST_ALLOC(double, volume_nodes_max);
      dxyz_drst_on_face_Gauss[i][j] = P4EST_ALLOC(double, face_nodes_max);
    }

  double* temp = P4EST_ALLOC(double, volume_nodes_max);
  double* J_on_face_Gauss = P4EST_ALLOC(double, face_nodes_max);


  
  p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)];
  
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = dgmath_is_child_left_or_right(c, d);
      q0[j][d] = e[0]->q[d] + cd*e[0]->dq;
    }
  }
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
        dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]);
        if (num_faces_side != num_faces_mortar)
          dqa[dir][d] /= 2;
      }
    }    

    p4est_qcoord_t dq0mf0 [(P4EST_DIM)];
    for (int d = 0; d < (P4EST_DIM); d++){
      dq0mf0[d] = (q0[0][d] - e[0]->q[d])/2;
    }
    
    for (int d = 0; d < (P4EST_DIM); d++){
        for (int c = 0; c < (P4EST_HALF); c++){
          q0[c][d] = e[0]->q[d];
          if (num_faces_side != num_faces_mortar)
            q0[c][d] += dq0mf0[d];
          for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
            int cd = dgmath_is_child_left_or_right(c, dir);
            q0[c][d] += cd*dqa[dir][d];
          }
        }
    }

    p4est_qcoord_t q [(P4EST_DIM)];
    p4est_qcoord_t mortar_dq = (num_faces_side == num_faces_mortar) ? e[0]->dq : e[0]->dq/2;
 
  int face_mortar_integ_stride = 0;
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
  
    int face_mortar_integ_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_integ[face_mortar]);
    /* double* xyz [(P4EST_DIM)]; */
    for (int d = 0; d < (P4EST_DIM); d++){
      q[d] = q0[face_mortar][d];
    }

    curved_element_data_compute_dxyz_drst
      (
       dgmath_jit_dbase,
       q,
       mortar_dq,
       e[0]->tree,
       p4est_geom,
       deg_mortar_integ[face_mortar],
       0,
       dxyz_drst,
       (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
           , NULL
#endif
           }
      );

    double* drst_dxyz_on_face_Gauss[(P4EST_DIM)][(P4EST_DIM)];
    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
        dgmath_apply_slicer(dgmath_jit_dbase,dxyz_drst[i][j],
                            (P4EST_DIM),
                            face,
                            deg_mortar_integ[face_mortar],
                            temp);

        dgmath_interp_GLL_to_GL
          (
           dgmath_jit_dbase,
           temp,
           deg_mortar_integ[face_mortar],
           deg_mortar_integ[face_mortar],
           dxyz_drst_on_face_Gauss[i][j],
           (P4EST_DIM)-1
          );

        drst_dxyz_on_face_Gauss[i][j] = &drst_dxyz_on_mortar_Gauss[i][j][face_mortar_integ_stride];
      }
    }
    
    if (sj_on_mortar_Gauss != NULL){
      curved_element_data_compute_J_and_rst_xyz(dxyz_drst_on_face_Gauss, J_on_face_Gauss, drst_dxyz_on_face_Gauss, face_mortar_integ_nodes);

      int i0 = -1; 
      if (face == 0 || face == 1){
        i0 = 0;
      }
      else if (face == 2 || face == 3){
        i0 = 1;
      }
      else if (face == 4 || face == 5){
        i0 = 2;
      }
      else {
        mpi_abort("face must be < 6\n");
      }
      double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
      for (int i = 0; i < face_mortar_integ_nodes; i++){
        sj_on_mortar_Gauss[face_mortar_integ_stride + i] = 0.;
        int is = face_mortar_integ_stride + i;
        for (int d = 0; d < (P4EST_DIM); d++){
          n_on_mortar_Gauss[d][is] = sgn*drst_dxyz_on_face_Gauss[i0][d][i]*J_on_face_Gauss[i];
          sj_on_mortar_Gauss[is] += n_on_mortar_Gauss[d][is]*n_on_mortar_Gauss[d][is];
        }
        sj_on_mortar_Gauss[is] = sqrt(sj_on_mortar_Gauss[is]);
        if (j_div_sj_mortar_Gauss != NULL)
          j_div_sj_mortar_Gauss[is] = J_on_face_Gauss[i]/sj_on_mortar_Gauss[is];
        for (int d = 0; d < (P4EST_DIM); d++)
          n_on_mortar_Gauss[d][is] /= sj_on_mortar_Gauss[is];
      }  
    }
    else {
      curved_element_data_compute_drst_dxyz(dxyz_drst_on_face_Gauss,
                                          drst_dxyz_on_face_Gauss,
                                          face_mortar_integ_nodes);

    }
    
    face_mortar_integ_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_integ[face_mortar]);
  }

  P4EST_FREE(temp);
  P4EST_FREE(J_on_face_Gauss);

  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      P4EST_FREE(dxyz_drst[i][j]);
      P4EST_FREE(dxyz_drst_on_face_Gauss[i][j]);
    }
}

void
curved_element_data_compute_physical_derivatives_on_face_Gauss_nodes
(
 double* dvec_drst_on_face_Gauss [(P4EST_DIM)], /* should be of mortar length, but not individually rotated */
 curved_element_data_t** e,
 int num_faces_side,
 int num_faces_mortar,
 int* deg_mortar_integ,
 int face_side,
 double* dvec_dxyz_on_face_Gauss [(P4EST_DIM)],
 p4est_geometry_t* geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase
)
{  
  int mortar_face_nodes_integ = 0;
  for (int i = 0; i < num_faces_mortar; i++){
    mortar_face_nodes_integ += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_integ[i]);
  }

  double* drst_dxyz_on_face_Gauss [(P4EST_DIM)][(P4EST_DIM)];
  for (int i = 0; i < (P4EST_DIM); i++){
    for (int j = 0; j < (P4EST_DIM); j++){
      drst_dxyz_on_face_Gauss[i][j] = P4EST_ALLOC(double, mortar_face_nodes_integ);
    }
  }
  
  curved_data_compute_drst_dxyz_Gauss_on_mortar_using_volume_data
    (
     e,
     num_faces_side,
     num_faces_mortar,
     deg_mortar_integ,
     face_side,
     drst_dxyz_on_face_Gauss,
     NULL,
     (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
         , NULL
#endif
         },     
     geom,
     dgmath_jit_dbase,
     NULL
    );


    for (int j = 0; j < (P4EST_DIM); j++){
      for (int k = 0; k < mortar_face_nodes_integ; k++){
        dvec_dxyz_on_face_Gauss[j][k] = 0.;
        for (int i = 0; i < (P4EST_DIM); i++){
          dvec_dxyz_on_face_Gauss[j][k] += drst_dxyz_on_face_Gauss[i][j][k]*dvec_drst_on_face_Gauss[i][k];
          /* printf("drst_dxyz_on_face_Gauss[i][j][k] = %.25f, dvec_drst_on_face_Gauss[i][k] = %.25f\n",drst_dxyz_on_face_Gauss[i][j][k],dvec_drst_on_face_Gauss[i][k]); */
      }
    }    
  }
 

  for (int i = 0; i < (P4EST_DIM); i++){
    for (int j = 0; j < (P4EST_DIM); j++){
      P4EST_FREE(drst_dxyz_on_face_Gauss[i][j]);
    }
  }

  
}

/* void */
/* curved_element_data_map_rst_to_xyz */
/* ( */
/*  d4est_real_t* r [3], */
/*  d4est_real_t* x [3], */
/*  p4est_qcoord_t q [3], */
/*  p4est_qcoord_t dq, */
/*  int deg, */
/*  int which_tree, */
/*  p4est_geometry_t* geom */
/* ) */
/* { */
/*   int nodes = dgmath_get_nodes( (P4EST_DIM) , deg); */
/*   d4est_real_t abc [] = {0,0,0};  /\* [0,1]**DIM *\/ */
/*   d4est_real_t xyz [] = {0,0,0};  /\* curvilinear coordinates *\/ */
/*   for (int i = 0; i < nodes; i++){ */
/*     for (int j = 0; j < (P4EST_DIM); j++){ */
/*       abc[j] = dgmath_rtox(r[j][i], */
/*                            (d4est_real_t)q[j], */
/*                            (d4est_real_t)dq */
/*                           )/(d4est_real_t)P4EST_ROOT_LEN; */

/*     } */
/*     geom->X(p4est_geom, which_tree, abc, xyz); */
/*     for (int j = 0; j < (P4EST_DIM); j++) */
/*       x[j][i] = xyz[j]; */
/*   } */
/* } */

/* void */
/* curved_element_data_compute_xyz */
/* ( */
/*  int deg, */
/*  p4est_qcoord_t qxyz [3], */
/*  p4est_qcoord_t dq, */
/*  d4est_real_t* xyz [3], */
/*  p4est_geometry_t* geom */
/* ) */
/* { */
/*   d4est_real_t* r = dgmath_fetch_xyz( d4est_ops, */
/*                                   (P4EST_DIM), */
/*                                    ed->deg, */
/*                                    0); */
  
/*   d4est_real_t* s = dgmath_fetch_xyz(d4est_ops, */
/*                                   (P4EST_DIM), */
/*                                    ed->deg, */
/*                                    1); */
/* #if (P4EST_DIM) == 3 */
/*   d4est_real_t* t = dgmath_fetch_xyz( d4est_ops, */
/*                                          (P4EST_DIM), */
/*                                          ed->deg, */
/*                                          2); */
/* #else */
/*   d4est_real_t* t = NULL; */
/* #endif */

/*   double* rst [3]; */
/*   rst[0] = r; */
/*   rst[1] = s; */
/*   rst[2] = t; */

/*   d4est_geometry_map_rst_to_xyz(rst, */
/*                                 xyz, */
/*                                 qxyz, */
/*                                 dq, */
/*                                 (P4EST_DIM), */
/*                                 ed->deg, */
/*                                 ed->which_tree, */
/*                                 geom */
/*                                ); */
/* } */



/* void */
/* d4est_geometry_compute_surface_jacobian_and_normal_ver1 */
/* ( */
/*  d4est_real_t* xyz_rst [(P4EST_DIM)][(P4EST_DIM)], */
/*  d4est_real_t* n [(P4EST_DIM)], */
/*  d4est_real_t* sj, */
/*  int face */
/* ) */
/* { */
/*   d4est_face_info_t face_info = d4est_geometry_get_face_info(f_m); */

/*   d4est_real_t* xyz_a_on_f_m [3]; */
/*   d4est_real_t* xyz_b_on_f_m [3]; */

/*   for (int d = 0; d < 3; d++){ */
/*     xyz_a_on_f_m[d] = P4EST_ALLOC_ZERO(d4est_real_t, face_nodes); */
/*     xyz_b_on_f_m[d] = P4EST_ALLOC_ZERO(d4est_real_t, face_nodes); */
/*   } */

/* #if (P4EST_DIM)==2 */
/*   if (face_info.a == 2){ */
/*     for (int i = 0; i < face_nodes; i++) */
/*       xyz_a_on_f_m[2][i] = 1.; */
/*   } */
/*   if (face_info.b == 2){ */
/*     for (int i = 0; i < face_nodes; i++) */
/*       xyz_b_on_f_m[2][i] = 1.; */
/*   } */
/* #endif */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     if (face_info.a < (P4EST_DIM)) */
/*       dgmath_apply_slicer */
/*         ( */
/*          d4est_ops, */
/*          &(xyz_rst[d][face_info.a][0]), */
/*          (P4EST_DIM), */
/*          f_m, */
/*          e_m->deg, */
/*          xyz_a_on_f_m[d] */
/*         ); */

/*     if (face_info.b < (P4EST_DIM)) */
/*       dgmath_apply_slicer */
/*         ( */
/*          d4est_ops, */
/*          &(xyz_rst[d][face_info.b][0]), */
/*          (P4EST_DIM), */
/*          f_m, */
/*          e_m->deg, */
/*          xyz_b_on_f_m[d] */
/*         ); */
/*   } */

/*   for (int i = 0; i < face_nodes; i++){ */
/*     double nx = NAN; */
/*     double ny = NAN; */
/*     double nz = NAN; */
/*     d4est_linalg_cross_prod */
/*       ( */
/*        xyz_a_on_f_m[0][i], */
/*        xyz_a_on_f_m[1][i], */
/*        xyz_a_on_f_m[2][i], */
/*        xyz_b_on_f_m[0][i], */
/*        xyz_b_on_f_m[1][i], */
/*        xyz_b_on_f_m[2][i], */
/*        &nx, */
/*        &ny, */
/*        &nz */
/*       ); */
/*     sj[i] = sqrt((nx*nx) + (ny*ny) + (nz*nz)); */
/*     n[0][i] = nx/sj[i]; */
/*     n[1][i] = ny/sj[i]; */
/*     n[2][i] = nz/sj[i]; */
/*   } */
/*   for (int d = 0; d < 3; d++){ */
/*     P4EST_FREE(xyz_a_on_f_m[d]); */
/*     P4EST_FREE(xyz_b_on_f_m[d]); */
/*   } */
/* } */

/* void */
/* d4est_geometry_compute_surface_jacobian_and_normal_ver2 */
/* ( */
/*  d4est_real_t* rst_xyz [(P4EST_DIM)][(P4EST_DIM)], */
/*  d4est_real_t* J, */
/*  d4est_real_t* n [(P4EST_DIM)], */
/*  d4est_real_t* sj, */
/*  int face, */
/*  int deg */
/* ) */
/* { */
/*   int i0 = -1; */
/*   double sgn = NAN; */
  
/*   if (face == 0 || face == 1){ */
/*     i0 = 0; */
/*     sgn = (face == 0) ? -1. : 1.; */
/*   } */
/*   else if (face == 2 || face == 3){ */
/*     i0 = 1; */
/*     sgn = (face == 2) ? -1. : 1.; */
/*   } */
/*   else if (face == 4 || face == 5){ */
/*     i0 = 2; */
/*     sgn = (face == 4) ? -1. : 1.; */
/*   } */
/*   else { */
/*     P4EST_ABORT("face must be < 6\n"); */
/*   } */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       n[d][i] = sgn*rst_xyz[i0][d][i]; */
/*       sj[i] += n[d][i]*n[d][i]; */
/*     } */
/*     sJ[i] = sqrt(sj[i]); */
/*     for (int d = 0; d < (P4EST_DIM); d++) */
/*       n[d][i] /= sj[i]; */
/*   } */
/* } */

/* void */
/* d4est_geometry_compute_surface_jacobian_and_normal_ver3 */
/* ( */
/*  p4est_qcoord_t qxyz0 [3], */
/*  p4est_qcoord_t dq, */
/*  int face, */
/*  int deg */
/* ) */
/* { */
/*   d4est_face_info_t face_info = d4est_geometry_get_face_info(f_m); */
  
/*   d4est_real_t* a = dgmath_fetch_xyz( d4est_ops, */
/*                                   (P4EST_DIM) - 1, */
/*                                    deg, */
/*                                    0); */
  
/*   d4est_real_t* b = dgmath_fetch_xyz(d4est_ops, */
/*                                   (P4EST_DIM) - 1, */
/*                                    deg, */
/*                                    1); */

/*   d4est_real_t* x = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* y = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* z = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dxda = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dyda = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dzda = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dxdb = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dydb = P4EST_ALLOC(d4est_real_t, face_nodes); */
/*   d4est_real_t* dzdb = P4EST_ALLOC(d4est_real_t, face_nodes); */

/*   for (int i = 0; i < face_nodes; i++){ */
/*     abc[0] = dgmath_rtox(a[i], */
/*                             (d4est_real_t)q[face_info.a], */
/*                             (d4est_real_t)dq)/(d4est_real_t)(P4EST_ROOT_LEN); */
    
/*     abc[1] = dgmath_rtox(b[i], (d4est_real_t)q[face_info.b], */
/*                             (d4est_real_t)dq)/(d4est_real_t)(P4EST_ROOT_LEN); */
/*     abc[2] = (d4est_real_t)q[face_info.c]/(d4est_real_t)(P4EST_ROOT_LEN); */
/*     p4est_geom->X(p4est_geom, which_tree, abc, xyz); */
/*     x[i] = xyz[0]; */
/*     y[i] = xyz[1]; */
/*     z[i] = xyz[2]; */
/*   } */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      x, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      0, */
/*      dxda */
/*     ); */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      y, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      0, */
/*      dyda */
/*     ); */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      z, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      0, */
/*      dzda */
/*     ); */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      x, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      1, */
/*      dxdb */
/*     ); */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      y, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      1, */
/*      dydb */
/*     ); */

/*   dgmath_apply_Dij */
/*     ( */
/*      dgmath_jit_dbase, */
/*      z, */
/*      (P4EST_DIM) - 1, */
/*      deg, */
/*      1, */
/*      dzdb */
/*     ); */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     double nx = 0.0; */
/*     double ny = 0.0; */
/*     double nz = 0.0; */
/*     d4est_linalg_cross_prod */
/*       ( */
/*        dxda[i], */
/*        dyda[i], */
/*        dzda[i], */
/*        dxdb[i], */
/*        dydb[i], */
/*        dzdb[i], */
/*        &nx, */
/*        &ny, */
/*        &nz */
/*       ); */
/*     sj[i] = sqrt((nx*nx) + (ny*ny) + (nz*nz)); */
/*     n[0][i] = nx/sj[i]; */
/*     n[1][i] = ny/sj[i]; */
/*     n[2][i] = nz/sj[i]; */
/*   } */
  
/*   P4EST_FREE(dzdb); */
/*   P4EST_FREE(dydb); */
/*   P4EST_FREE(dxdb); */
/*   P4EST_FREE(dzda); */
/*   P4EST_FREE(dxda); */
/*   P4EST_FREE(dyda); */
/*   P4EST_FREE(z); */
/*   P4EST_FREE(y); */
/*   P4EST_FREE(x); */
/* } */


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
    
    curved_element_data_compute_xyz
      (
       dgmath_jit_dbase,
       d4est_geometry->p4est_geom,
       elem_data->tree,
       deg_Gauss,
       GAUSS,
       elem_data->q,
       elem_data->dq,
       xyz_integ
      );

    curved_element_data_compute_dxyz_drst
      (
       dgmath_jit_dbase,
       elem_data->q,
       elem_data->dq,
       elem_data->tree,
       d4est_geometry->p4est_geom,
       deg_Gauss,
       1,
       xyz_rst_integ,
       (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
           , NULL
#endif
           }

      );

    curved_element_data_compute_J_and_rst_xyz
      (
       xyz_rst_integ,
       J_integ,
       rst_xyz_integ,
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
    
    
    curved_element_data_compute_xyz
      (
       dgmath_jit_dbase,
       d4est_geometry->p4est_geom,
       elem_data->tree,
       deg_Gauss,
       GAUSS,
       elem_data->q,
       elem_data->dq,
       xyz_integ
      );

    curved_element_data_compute_dxyz_drst
      (
       dgmath_jit_dbase,
       elem_data->q,
       elem_data->dq,
       elem_data->tree,
       d4est_geometry->p4est_geom,
       deg_Gauss,
       1,
       xyz_rst_integ,
       (double* [(P4EST_DIM)]){NULL, NULL
#if (P4EST_DIM)==3
           , NULL
#endif
           }

      );

    curved_element_data_compute_J_and_rst_xyz
      (
       xyz_rst_integ,
       J_integ,
       rst_xyz_integ,
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
