#include <pXest.h>
#include <d4est_mesh.h>
#include <d4est_geometry.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_util.h>
#include <sc_reduce.h>
#include <d4est_linalg.h>

d4est_mesh_geometry_storage_t*
d4est_mesh_geometry_storage_init()
{
  d4est_mesh_geometry_storage_t* geometric_factors = P4EST_ALLOC(d4est_mesh_geometry_storage_t, 1);
  geometric_factors->J_quad = NULL; 
  geometric_factors->xyz = NULL;
  geometric_factors->xyz_quad = NULL;
  geometric_factors->xyz_rst_quad = NULL;
  geometric_factors->rst_xyz_quad = NULL;

  return geometric_factors;
}

static void
d4est_mesh_geometry_storage_realloc
(
 p4est_t* p4est,
 d4est_mesh_geometry_storage_t* geometric_factors,
 d4est_local_sizes_t local_sizes
)
{
  int local_nodes = local_sizes.local_nodes;
  int local_nodes_quad = local_sizes.local_nodes_quad;
    
  int vector_nodes = local_nodes*(P4EST_DIM); 
  geometric_factors->xyz = P4EST_REALLOC(geometric_factors->xyz,double,vector_nodes);
  geometric_factors->xyz_quad = P4EST_REALLOC(geometric_factors->xyz_quad, double, (P4EST_DIM)*local_nodes_quad);

  int matrix_nodes_quad = local_nodes_quad*(P4EST_DIM)*(P4EST_DIM);
  geometric_factors->J_quad = P4EST_REALLOC(geometric_factors->J_quad,double,local_nodes_quad);
  geometric_factors->xyz_rst_quad = P4EST_REALLOC(geometric_factors->xyz_rst_quad,double,matrix_nodes_quad);  
  geometric_factors->rst_xyz_quad = P4EST_REALLOC(geometric_factors->rst_xyz_quad,double,matrix_nodes_quad);

}

void
d4est_mesh_geometry_storage_destroy
(
 d4est_mesh_geometry_storage_t* geometric_factors
)
{
  if (geometric_factors != NULL){
    P4EST_FREE(geometric_factors->J_quad);
    P4EST_FREE(geometric_factors->xyz);
    P4EST_FREE(geometric_factors->xyz_quad);
    P4EST_FREE(geometric_factors->xyz_rst_quad);
    P4EST_FREE(geometric_factors->rst_xyz_quad);
    P4EST_FREE(geometric_factors);
  }
}


void
d4est_mesh_get_array_of_degrees
(
 p4est_t* p4est,
 int* deg_array
)
{
  int stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        deg_array[stride] = ed->deg;
        stride++;
      }
    }
}


void
d4est_mesh_get_array_of_estimators
(
 p4est_t* p4est,
 double* eta2_array
)
{
  int stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        eta2_array[stride] = ed->local_estimator;
        stride++;
      }
    }
}

void
d4est_mesh_compute_jacobian_on_lgl_grid
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 double* jacobian_lgl
)
{
  int nodal_stride = 0;
  double* xyz_rst [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(xyz_rst, d4est_lgl_get_nodes((P4EST_DIM),(MAX_DEGREE)));

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg);


        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = elem_data->dq;
        mesh_object.tree = elem_data->tree;
        mesh_object.element_id = elem_data->id;
        mesh_object.q[0] = elem_data->q[0];
        mesh_object.q[1] = elem_data->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = elem_data->q[2];
#endif
      
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg,
                                                                0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL, NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                elem_data->deg, 2);
#endif

        d4est_geometry_compute_dxyz_drst
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           elem_data->tree,
           elem_data->q,
           elem_data->dq,
           elem_data->deg,           
           xyz_rst
          );      


        d4est_geometry_compute_jacobian
          (
           xyz_rst,
           &jacobian_lgl[nodal_stride],
           volume_nodes
          );

        
        nodal_stride += volume_nodes;
      }
    }
                       
  D4EST_FREE_DBYD_MAT(xyz_rst);
}

int d4est_mesh_get_local_matrix_nodes(p4est_t* p4est){

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
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        local_matrix_nodes += volume_nodes*volume_nodes;
      }
    }
  return local_matrix_nodes;
}

void
d4est_mesh_print_number_of_elements_per_tree
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

d4est_element_data_t*
d4est_mesh_get_element_data
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
        d4est_element_data_t* ed = quad->p.user_data;
        if (ed->id == local_element_id)
          return ed;
      }
    }  
  return NULL;
}

int
d4est_mesh_global_node_to_local_node
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
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        if(global_node >= ed->nodal_stride && global_node < (ed->nodal_stride + volume_nodes))
          return global_node - ed->nodal_stride;
      }
    }
  return -1;
}

/* only for serial use */
int
d4est_mesh_debug_find_node
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
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        if(node >= ed->nodal_stride && node < (ed->nodal_stride + volume_nodes))
          return ed->id;
      }
    }
  return -1;
}

void
d4est_mesh_print_element_data_debug
(
 p4est_t* p4est
)
{
/* #ifndef D4EST_DEBUG */
/*   D4EST_ABORT("compile with the debug flag if you want to print curved element data"); */
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
        d4est_element_data_t* ed = quad->p.user_data;

        printf("** Element %d **\n", ed->id);
        printf("deg, deg_quad = %d, %d\n", ed->deg, ed->deg_quad);
        printf("tree, tree_quadid = %d, %d\n", ed->tree, ed->tree_quadid);

        
        int volume_nodes = d4est_lgl_get_nodes( (P4EST_DIM), ed->deg );
       
        
#if (P4EST_DIM)==2 
        printf("q = %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->dq);
        DEBUG_PRINT_2ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes);
        /* DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], volume_nodes); */
#elif (P4EST_DIM)==3
        printf("q = %d, %d, %d, dq = %d\n", ed->q[0], ed->q[1], ed->q[2], ed->dq);
        DEBUG_PRINT_3ARR_DBL(ed->xyz[0], ed->xyz[1], ed->xyz[2], volume_nodes);
#else
        D4EST_ABORT("DIM = 2 or 3");
#endif
        
/* #else */
        /* D4EST_ABORT("DEBUG flag must be set"); */
/* #endif */
      }
    }  
}

double
d4est_mesh_compute_l2_norm_sqr
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 double* nodal_vec,
 int local_nodes,
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
        d4est_element_data_t* ed = quad->p.user_data;
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif

        d4est_quadrature_apply_mass_matrix
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &nodal_vec[ed->nodal_stride],
           ed->deg,
           ed->J_quad,
           ed->deg_quad,
           &Mvec[ed->nodal_stride]
          );
      
        double norm2
          = d4est_linalg_vec_dot(&nodal_vec[ed->nodal_stride],
                                 &Mvec[ed->nodal_stride],
                                 volume_nodes);
        
        if (store_local == STORE_LOCALLY){
          ed->local_estimator = norm2;
        }
        
        l2_norm_sqr += norm2;
      }
    }
  P4EST_FREE(Mvec);
  return l2_norm_sqr;
}



d4est_local_sizes_t
d4est_mesh_init_element_data
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 void(*user_fcn)(d4est_element_data_t*, void*),
 void* user_ctx
)
{
  /* sizes */
  int local_nodes = 0;
  int local_sqr_nodes = 0;
  int local_sqr_mortar_nodes = 0;
  int local_nodes_quad = 0;
  /* int local_sqr_nodes_invM = 0; */

  /* strides */
  int sqr_nodal_stride = 0;
  int sqr_mortar_stride = 0;
  int nodal_stride = 0;
  int quad_stride = 0;
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
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);
        elem_data->tree = tt;
        elem_data->tree_quadid = q;
        elem_data->dq = P4EST_QUADRANT_LEN(quad->level);
        elem_data->q[0] = quad->x;
        elem_data->q[1] = quad->y;
#if (P4EST_DIM)==3  
        elem_data->q[2] = quad->z;
#endif        


        
        elem_data->id = id_stride;
        elem_data->sqr_nodal_stride = sqr_nodal_stride;
        elem_data->sqr_mortar_stride = sqr_mortar_stride;
        elem_data->nodal_stride = nodal_stride;
        elem_data->quad_stride = quad_stride;

        /* user_fcn should set degree, 
           or the degree will be assumed to be set */
        if (user_fcn != NULL){
          /* problem_set_degrees_donald_trump(elem_data, user_ctx); */
          user_fcn(elem_data, user_ctx);
        }

        /* printf("elem_data->deg = %d\n", elem_data->deg); */
        /* printf("elem_data->deg_quad = %d\n", elem_data->deg_quad); */
        
        /* elem_data->deg = 2; */
        /* elem_data->deg_quad = 2; */

        D4EST_ASSERT(elem_data->deg > 0
                   &&
                   elem_data->deg_quad > 0
                   &&
                   elem_data->deg < MAX_DEGREE
                  );

        
        int nodes = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg);
        int nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg_quad);
        int face_nodes = d4est_lgl_get_nodes((P4EST_DIM)-1, elem_data->deg);
        local_nodes += nodes;
        local_sqr_nodes += nodes*nodes;
        local_sqr_mortar_nodes += (P4EST_FACES)*face_nodes*face_nodes;
        local_nodes_quad += d4est_lgl_get_nodes((P4EST_DIM), elem_data->deg_quad);

        sqr_nodal_stride += nodes*nodes;
        sqr_mortar_stride += face_nodes*face_nodes*(P4EST_FACES);
        nodal_stride += nodes;
        quad_stride += nodes_quad;
        id_stride += 1;
      }
    }

  d4est_local_sizes_t local_sizes;
  local_sizes.local_nodes = local_nodes;
  local_sizes.local_sqr_nodes = local_sqr_nodes;
  local_sizes.local_sqr_mortar_nodes = local_sqr_mortar_nodes;
  local_sizes.local_nodes_quad = local_nodes_quad;
  /* local_sizes.local_sqr_nodes_invM = local_sqr_nodes_invM; */
  return local_sizes;
}

static void
d4est_mesh_geometry_storage_initialize_data
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* geometric_factors
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
 
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
      
        d4est_rst_t rst_points_quad;
        rst_points_quad = d4est_quadrature_get_rst_points
                          (
                           d4est_ops,
                           d4est_quad,
                           d4est_geom,
                           &mesh_object,
                           QUAD_OBJECT_VOLUME,
                           QUAD_INTEGRAND_UNKNOWN,
                           ed->deg_quad
                          );


        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ed->deg, 2);
#endif
              
        int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           tt,
           ed->deg,
           ed->q,
           ed->dq,
           ed->xyz
          );

        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_quad,
           tt,
           ed->deg_quad,
           ed->q,
           ed->dq,
           ed->xyz_quad
          );

        d4est_geometry_compute_dxyz_drst
          (
           d4est_ops,
           d4est_geom,
           rst_points_quad,
           ed->tree,
           ed->q,
           ed->dq,
           ed->deg_quad,           
           ed->xyz_rst_quad
          );

    
        d4est_geometry_compute_jacobian
          (
           ed->xyz_rst_quad,
           ed->J_quad,
           volume_nodes_quad
          );

        d4est_geometry_compute_drst_dxyz
          (
           ed->xyz_rst_quad,
           ed->J_quad,
           ed->rst_xyz_quad,
           volume_nodes_quad
          );
      }
    }
}



void
d4est_mesh_geometry_storage_initialize_aliases
(
 p4est_t* p4est,
 d4est_mesh_geometry_storage_t* geometric_factors,
 d4est_local_sizes_t local_sizes
)
{
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int QQ = (p4est_locidx_t) tquadrants->elem_count;
      for (int qq = 0; qq < QQ; ++qq) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);

        elem_data->J_quad = &geometric_factors->J_quad[elem_data->quad_stride];  
        for (int i = 0; i < (P4EST_DIM); i++){
          elem_data->xyz[i] = &geometric_factors->xyz[i*local_sizes.local_nodes + elem_data->nodal_stride];
          elem_data->xyz_quad[i] = &geometric_factors->xyz_quad[i*local_sizes.local_nodes_quad + elem_data->quad_stride];
          for (int j = 0; j < (P4EST_DIM); j++){
            elem_data->xyz_rst_quad[i][j] = &geometric_factors->xyz_rst_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride];
            elem_data->rst_xyz_quad[i][j] = &geometric_factors->rst_xyz_quad[(i*(P4EST_DIM) + j)*local_sizes.local_nodes_quad + elem_data->quad_stride];
          }
        }

      }
    }
}

int
d4est_mesh_update
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_geometry_storage_t* geometric_factors,
 d4est_mesh_quadrature_data_init_option_t quad_init_option,
 d4est_mesh_geometry_data_init_option_t geom_init_option,
 d4est_mesh_geometry_aliases_init_option_t alias_init_option,
 void(*user_fcn)(d4est_element_data_t*, void*),
 void* user_ctx
)
{
  d4est_local_sizes_t local_sizes = d4est_mesh_init_element_data(p4est,
                                                                 d4est_ops,
                                                                 user_fcn,//problem_set_degrees_donald_trump,
                                                                 user_ctx);
  if (quad_init_option == INITIALIZE_QUADRATURE_DATA)
    {
      d4est_quadrature_reinit(
                              p4est,
                              ghost,
                              ghost_data,
                              d4est_ops,
                              d4est_geom,
                              d4est_quad
      );
    }

    
  
  if (geom_init_option == INITIALIZE_GEOMETRY_DATA)
    {
      d4est_mesh_geometry_storage_realloc(
                                          p4est,
                                          geometric_factors,
                                          local_sizes
      );
      d4est_mesh_geometry_storage_initialize_aliases(p4est, geometric_factors, local_sizes);
      d4est_mesh_geometry_storage_initialize_data(
                                                  p4est,
                                                  d4est_ops,
                                                  d4est_geom,
                                                  d4est_quad,
                                                  geometric_factors
      );
    }

  if (alias_init_option == INITIALIZE_GEOMETRY_ALIASES &&
      geom_init_option != INITIALIZE_GEOMETRY_DATA /* if this is false, then this will be a waste of time */
     ){
    d4est_mesh_geometry_storage_initialize_aliases(p4est, geometric_factors, local_sizes);
  }

  return local_sizes.local_nodes; 
}

void
d4est_mesh_init_field
(
 p4est_t* p4est,
 double* node_vec,
 d4est_xyz_fcn_t init_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 void* user
)
{

  double* xyz_temp [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    xyz_temp[d] = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), (MAX_DEGREE)));
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
        d4est_element_data_t* ed = quad->p.user_data;        
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
      
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg,
                                                                0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL, NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 2);
#endif

        
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           tt,
           ed->deg,
           ed->q,
           ed->dq,
           xyz_temp
          );

        
        for (int i = 0; i < volume_nodes; i++){
          node_vec[ed->nodal_stride + i] = init_fcn(xyz_temp[0][i],
                                                    xyz_temp[1][i],
#if (P4EST_DIM)==3
                                                    xyz_temp[2][i],
#endif
                                                    user
                                                   );

        }
      }
    }


  for (int d = 0; d < (P4EST_DIM); d++) {
    P4EST_FREE(xyz_temp[d]);
  }
  
}

void
d4est_mesh_init_field_ext
(
 p4est_t* p4est,
 double* node_vec,
 d4est_xyz_fcn_ext_t fofxyzv,
 double* v,
 double* fofxyzv_user,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{

  double* xyz_temp [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    xyz_temp[d] = P4EST_ALLOC(double, d4est_lgl_get_nodes((P4EST_DIM), (MAX_DEGREE)));
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
        d4est_element_data_t* ed = quad->p.user_data;        
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);


          d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
      
        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg,
                                                                0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL, NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops,
                                                                NULL,
                                                                NULL,
                                                                &mesh_object,
                                                                QUAD_OBJECT_VOLUME,
                                                                QUAD_INTEGRAND_UNKNOWN,
                                                                ed->deg, 2);
#endif

        
        for (int i = 0; i < volume_nodes; i++){
          node_vec[ed->nodal_stride + i] = fofxyzv(xyz_temp[0][i],
                                                    xyz_temp[1][i],
#if (P4EST_DIM)==3
                                                    xyz_temp[2][i],
#endif
                                                    v[i],
                                                    fofxyzv_user
                                                   );

        }
      }
    }


  for (int d = 0; d < (P4EST_DIM); d++) {
    P4EST_FREE(xyz_temp[d]);
  }
  
}

void
d4est_mesh_get_local_nodes_callback
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  p4est_quadrant_t* q = info->quad;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) q->p.user_data;
  int* local_nodes = (int*) user_data;
  *local_nodes += d4est_lgl_get_nodes( (P4EST_DIM),
                                    elem_data->deg);
}

int d4est_mesh_get_local_nodes(p4est_t* p4est)
{
  int local_nodes = 0;
  
  p4est_iterate(p4est,
                NULL,
                (void *) (&local_nodes),
                d4est_mesh_get_local_nodes_callback,
                NULL,
#if (P4EST_DIM)==3
                NULL,       
#endif                
                NULL);
  
  return local_nodes;
}
