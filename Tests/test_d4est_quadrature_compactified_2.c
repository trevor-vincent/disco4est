
#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_cubed_sphere.h>
#include <curved_compute_flux.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <util.h>
#include <limits.h>

#define DEG_LOBATTO 2
#define DEG_QUAD 2

void
problem_set_degrees
(
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  elem_data->deg = DEG_LOBATTO;
  elem_data->deg_quad = DEG_QUAD;
}

typedef struct {
  
  double surface_integral_GL;
  double surface_integral_comp;
  
} test_d4est_quadrature_compactified_surface_integrals_data_t;

static int
uni_refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}


p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}


static void
test_d4est_quadrature_compactified_surface_integrals_bndry
(
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 grid_fcn_t bndry_fcn,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{

  test_d4est_quadrature_compactified_surface_integrals_data_t*  data = params;
  int face_nodes_m_lobatto = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg);
  int face_nodes_m_quad = d4est_operators_get_nodes((P4EST_DIM) - 1, e_m->deg_quad);

  double* sj_on_f_m_quad_GL = P4EST_ALLOC(double, face_nodes_m_quad);
  double* sj_on_f_m_quad_comp = P4EST_ALLOC(double, face_nodes_m_quad);
  
  double* ones_quad = P4EST_ALLOC(double, face_nodes_m_quad);
  d4est_linalg_fill_vec(ones_quad, 1., face_nodes_m_quad);
  

  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_UNKNOWN_INTEGRAND,
     e_m->tree,
     e_m->q,
     e_m->dq,
     mortar_side_id_m,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     NULL,
     sj_on_f_m_quad_GL,
     NULL,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );


  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_JAC_TIMES_POLY_INTEGRAND,
     e_m->tree,
     e_m->q,
     e_m->dq,
     mortar_side_id_m,
     1,
     1,
     &e_m->deg_quad,
     f_m,
     NULL,
     sj_on_f_m_quad_comp,
     NULL,
     NULL,
     NULL,//     J_div_SJ_quad_for_term3,
     COMPUTE_NORMAL_USING_JACOBIAN
    );
  
  
  d4est_quadrature_mortar_t face_object;
  face_object.dq = e_m->dq;
  face_object.tree = e_m->tree;
  face_object.face = f_m;
  face_object.mortar_side_id = mortar_side_id_m;
  face_object.mortar_subface_id = 0;
  
  face_object.q[0] = e_m->q[0];
  face_object.q[1] = e_m->q[1];
#if (P4EST_DIM)==3
  face_object
    .q[2] = e_m->q[2];
#endif



    data->surface_integral_comp +=
    d4est_quadrature_innerproduct
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_MORTAR,
     QUAD_JAC_TIMES_POLY_INTEGRAND,
     ones_quad,
     NULL,
     sj_on_f_m_quad_comp,
     e_m->deg_quad
    );

    data->surface_integral_GL +=
    d4est_quadrature_innerproduct
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_MORTAR,
     QUAD_UNKNOWN_INTEGRAND,
     ones_quad,
     NULL,
     sj_on_f_m_quad_GL,
     e_m->deg_quad
    );
  

  P4EST_FREE(sj_on_f_m_quad_comp);
  P4EST_FREE(sj_on_f_m_quad_GL);
  P4EST_FREE(ones_quad);
}



static void
test_d4est_quadrature_compactified_surface_integrals_interface
(
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* params
)
{
  test_d4est_quadrature_compactified_surface_integrals_data_t*  data = params;
  
  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  /* int deg_p_quad [(P4EST_HALF)]; */
  int face_nodes_p_quad [(P4EST_HALF)];

  int deg_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  /* int deg_m_quad [(P4EST_HALF)]; */
  int face_nodes_m_quad [(P4EST_HALF)];
  
  int nodes_mortar_quad [(P4EST_HALF)];
  int nodes_mortar_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;
  double penalty_mortar [(P4EST_HALF)];


  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
  
  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    /* deg_m_quad[i] = e_m[i]->deg_quad; */
    
    face_nodes_m_lobatto[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg_quad);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_lobatto_porder[i] = e_p[i]->deg;
    /* deg_p_quad[i] = e_p_oriented[i]->deg_quad; */

    face_nodes_p_lobatto[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_operators_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg_quad);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }    

  /* calculate degs and nodes of the mortar faces */
  int total_nodes_mortar_quad = 0;
  int total_nodes_mortar_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = util_max_int( e_m[i]->deg_quad,
                                            e_p_oriented[j]->deg_quad);
      deg_mortar_lobatto[i+j] = util_max_int( e_m[i]->deg,
                                              e_p_oriented[j]->deg );      
      nodes_mortar_quad[i+j] = d4est_operators_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );     
      nodes_mortar_lobatto[i+j] = d4est_operators_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );     
      total_nodes_mortar_quad += nodes_mortar_quad[i+j];
      total_nodes_mortar_lobatto += nodes_mortar_lobatto[i+j];
      
    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int nodes_mortar_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_operators_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    nodes_mortar_quad_porder[inew] = nodes_mortar_quad[i];
  }

   
  double* sjvol_m_on_f_m_mortar_quad_comp = P4EST_ALLOC(double, total_nodes_mortar_quad);
  double* sjvol_m_on_f_m_mortar_quad_GL = P4EST_ALLOC(double, total_nodes_mortar_quad);

  p4est_qcoord_t mortar_q0_forder [(P4EST_HALF)][(P4EST_DIM)];
  p4est_qcoord_t mortar_dq_forder;

  d4est_geometry_compute_qcoords_on_mortar
    (
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     faces_m,
     faces_mortar,
     f_m,
     mortar_q0_forder,
     &mortar_dq_forder
    );


  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_JAC_TIMES_POLY_INTEGRAND,
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     mortar_side_id_m,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     NULL,
     sjvol_m_on_f_m_mortar_quad_comp,
     NULL,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );


  d4est_mortars_compute_geometric_data_on_mortar
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     QUAD_UNKNOWN_INTEGRAND,
     e_m[0]->tree,
     e_m[0]->q,
     e_m[0]->dq,
     mortar_side_id_m,
     faces_m,
     faces_mortar,
     &deg_mortar_quad[0],
     f_m,
     NULL,
     sjvol_m_on_f_m_mortar_quad_GL,
     NULL,
     NULL,
     NULL,
     COMPUTE_NORMAL_USING_JACOBIAN
    );


  double* ones = P4EST_ALLOC(double, total_nodes_mortar_quad);
  for (int i = 0; i < total_nodes_mortar_quad; i++)
    ones[i] = 1.;
  
  for (int i = 0; i < faces_mortar; i++){



    d4est_quadrature_mortar_t face_object;
    face_object.dq = mortar_dq_forder;
    face_object.tree = e_m[0]->tree;
    face_object.mortar_side_id = mortar_side_id_m;
    face_object.mortar_subface_id = i;
    
    face_object.face = f_m;
    face_object.q[0] = mortar_q0_forder[i][0];
    face_object.q[1] = mortar_q0_forder[i][1];
#if (P4EST_DIM)==3
    face_object.q[2] = mortar_q0_forder[i][2];
#endif
    
    data->surface_integral_comp +=
    d4est_quadrature_innerproduct
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_MORTAR,
     QUAD_JAC_TIMES_POLY_INTEGRAND,
     ones,
     NULL,
     sjvol_m_on_f_m_mortar_quad_comp,
     deg_mortar_quad[i]
    );
  }
  for (int i = 0; i < faces_mortar; i++){


    d4est_quadrature_mortar_t face_object;
    face_object.dq = mortar_dq_forder;
    face_object.tree = e_m[0]->tree;
    face_object.mortar_side_id = mortar_side_id_m;
    face_object.mortar_subface_id = i;
    
    face_object.face = f_m;
    face_object.q[0] = mortar_q0_forder[i][0];
    face_object.q[1] = mortar_q0_forder[i][1];
#if (P4EST_DIM)==3
    face_object.q[2] = mortar_q0_forder[i][2];
#endif
    
    data->surface_integral_GL +=
    d4est_quadrature_innerproduct
    (
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &face_object,
     QUAD_MORTAR,
     QUAD_UNKNOWN_INTEGRAND,
     ones,
     NULL,
     sjvol_m_on_f_m_mortar_quad_GL,
     deg_mortar_quad[i]
    );
  }
  
  
  P4EST_FREE(sjvol_m_on_f_m_mortar_quad_GL);
  P4EST_FREE(ones);
  P4EST_FREE(sjvol_m_on_f_m_mortar_quad_comp);
  
}

curved_flux_fcn_ptrs_t
test_d4est_quadrature_compactified_surface_integrals_fetch_fcns
(
 test_d4est_quadrature_compactified_surface_integrals_data_t* data
)
{
  
  curved_flux_fcn_ptrs_t curved_test_mortarjacobianterms_fcns;
  curved_test_mortarjacobianterms_fcns.flux_interface_fcn = test_d4est_quadrature_compactified_surface_integrals_interface;
  curved_test_mortarjacobianterms_fcns.flux_boundary_fcn = test_d4est_quadrature_compactified_surface_integrals_bndry;
  curved_test_mortarjacobianterms_fcns.params = (void*)data;

  return curved_test_mortarjacobianterms_fcns;
}


int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_CUBED_SPHERE_OUTER_SHELL;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = 10;
  sphere_attrs->R2 = 1000;
  sphere_attrs->compactify_outer_shell = 1;
  sphere_attrs->compactify_inner_shell = -1;
  d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geom, sphere_attrs);
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init();
  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG4;
  d4est_quadrature_compactified_new(p4est, d4est_ops, d4est_geom, d4est_quad, "", "");


  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  int num_unifrefs = 1;
  for (int level = 0; level < num_unifrefs; ++level){

      p4est_refine_ext(p4est,
                       0,
                       -1,
                       uni_refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FACE,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);

  }


  d4est_mesh_update
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     geometric_factors,
     INITIALIZE_QUADRATURE_DATA,
     INITIALIZE_GEOMETRY_DATA,
     INITIALIZE_GEOMETRY_ALIASES,
     problem_set_degrees,
     NULL
    );
  
  test_d4est_quadrature_compactified_surface_integrals_data_t test_data;
  test_data.surface_integral_GL = 0.;
  test_data.surface_integral_comp = 0.;
  curved_flux_fcn_ptrs_t ffp = test_d4est_quadrature_compactified_surface_integrals_fetch_fcns(&test_data);
  
  curved_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &ffp,
     EXCHANGE_GHOST_DATA
    );
    

  printf("surface_integral_GL = %.25f\n", test_data.surface_integral_GL);
  printf("surface_integral_comp = %.25f\n", test_data.surface_integral_comp);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_geometry_destroy(d4est_geom);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
