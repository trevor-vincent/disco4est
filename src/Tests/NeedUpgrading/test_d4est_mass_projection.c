#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_brick.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_amr.h>
#include <d4est_poisson.h>
#include <d4est_poisson_flux.h>
#include <d4est_util.h>
#include <limits.h>

#define D4EST_REAL_EPS 100*1e-15
#define TEST_DEG_INIT 2


typedef struct {

  double quad_vs_mass_err;
  
} testd4est_mass_projection_data_t;


static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  elem_data->deg = TEST_DEG_INIT;
  elem_data->deg_quad = TEST_DEG_INIT;
}

double
poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  int deg = TEST_DEG_INIT;
  double poly = pow(x,deg) + pow(y,deg);
#if (P4EST_DIM)==3
  poly += ((P4EST_DIM)==3) ? pow(z,deg) : 0.;
#endif
  return poly;
}


double
laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  int deg = TEST_DEG_INIT;
  double poly = pow(x,deg-2) + pow(y,deg-2);
#if (P4EST_DIM)==3
  poly += ((P4EST_DIM)==3) ? pow(z,deg-2) : 0.;
#endif
  double factor = deg;
  while (factor != 0){
    poly *= factor;
    factor -= 1;
  }
  return poly;
}



static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg_quad = elem_data->deg;

}


static p4est_t*
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
testd4est_mass_projection_interface
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
 d4est_poisson_flux_interface_data_t* mortar_data,
 void* params
)
{
  testd4est_mass_projection_data_t* data = params;
  d4est_quadrature_mortar_t* mortar_face_object = mortar_data->mortar_face_object;
  
  int faces_mortar = mortar_data->faces_mortar;
  int total_side_nodes_m_lobatto = mortar_data->total_side_nodes_m_lobatto;
  int total_nodes_mortar_lobatto = mortar_data->total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad = mortar_data->total_nodes_mortar_quad;
  
  double* u_m_on_f_m_mortar_quad = mortar_data->u_m_on_f_m_mortar_quad;
  double* u_m_on_f_m = mortar_data->u_m_on_f_m;
  double* sj_on_f_m_mortar_quad = mortar_data->sj_on_f_m_mortar_quad;
  double* j_div_sj_on_f_m_mortar_quad = mortar_data->j_div_sj_on_f_m_mortar_quad;
  double* u_p_on_f_p_mortar_quad = mortar_data->u_p_on_f_p_mortar_quad;
  double* j_div_sj_on_f_p_mortar_quad = mortar_data->j_div_sj_on_f_p_mortar_quad;
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_mortar_quad [(P4EST_DIM)];

  D4EST_COPY_DBYD_MAT(mortar_data->drst_dxyz_m_on_mortar_quad, drst_dxyz_m_on_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_m_on_f_m_mortar_quad, dudx_m_on_f_m_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->dudx_p_on_f_p_mortar_quad, dudx_p_on_f_p_mortar_quad);
  D4EST_COPY_DIM_VEC(mortar_data->n_sj_on_f_m_mortar_quad, n_sj_on_f_m_mortar_quad);
  
  int* deg_mortar_quad = mortar_data->deg_mortar_quad;
  int* nodes_mortar_quad = mortar_data->nodes_mortar_quad;
  int* nodes_mortar_lobatto = mortar_data->nodes_mortar_lobatto;
  int* deg_mortar_lobatto = mortar_data->deg_mortar_lobatto;
  int* face_nodes_m_lobatto = mortar_data->deg_mortar_lobatto;
  int* deg_m_lobatto = mortar_data->deg_m_lobatto;
  int* deg_p_lobatto = mortar_data->deg_p_lobatto;

  double* u_m_on_f_m_mortar_lobatto = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* M_u_using_quadrature = P4EST_ALLOC(double, total_nodes_mortar_lobatto);
  double* M_u_using_mass = P4EST_ALLOC(double, total_nodes_mortar_lobatto);

  
  d4est_mortars_project_side_onto_mortar_space
    (
     d4est_ops,
     u_m_on_f_m,
     faces_m,
     deg_m_lobatto,
     u_m_on_f_m_mortar_lobatto,
     faces_mortar,
     deg_mortar_lobatto
    );
  
  int stride = 0;
  int stride_lobatto = 0;
  for (int f = 0; f < faces_mortar; f++){
    for (int k = 0; k < nodes_mortar_quad[f]; k++){
      int ks = k + stride;

      d4est_quadrature_apply_galerkin_integral
        (
         d4est_ops,
         d4est_geom,
         d4est_quad,
         mortar_face_object,
         QUAD_OBJECT_MORTAR,
         QUAD_INTEGRAND_UNKNOWN,
         &u_m_on_f_m_mortar_quad[stride],
         deg_mortar_lobatto[f],
         &sj_on_f_m_mortar_quad[stride],
         deg_mortar_quad[f],
         &M_u_using_quadrature[stride_lobatto]
        );

      d4est_operators_apply_mij(
                                d4est_ops,
                                &u_m_on_f_m_mortar_lobatto[stride_lobatto],
                                (P4EST_DIM)-1,
                                deg_mortar_lobatto[f],
                                &M_u_using_mass[stride_lobatto]
                               );


      d4est_linalg_vec_scale(sj_on_f_m_mortar_quad[stride], &M_u_using_mass[stride_lobatto], nodes_mortar_lobatto[f]);
      
      double biggest_error;
      int biggest_error_id;

      d4est_util_find_biggest_error(&M_u_using_quadrature[stride_lobatto], &M_u_using_mass[stride_lobatto], nodes_mortar_lobatto[f], &biggest_error, &biggest_error_id);


      double* mass = &M_u_using_mass[stride_lobatto];
      double* quad = &M_u_using_quadrature[stride_lobatto];
      
      if (biggest_error > D4EST_REAL_EPS){
        printf ("f = %d\n", f);
        printf ("e_m[0]->id = %d\n", e_m[0]->id);
        printf ("e_p[0]->id = %d\n", e_m[0]->id);
        DEBUG_PRINT_ARR_INT(nodes_mortar_lobatto, faces_mortar);
        DEBUG_PRINT_ARR_INT(nodes_mortar_quad, faces_mortar);
        DEBUG_PRINT_2ARR_DBL(mass, quad, nodes_mortar_lobatto[f]);
      }
      
      data->quad_vs_mass_err += biggest_error;
                         
    }
    stride += nodes_mortar_quad[f];
    stride_lobatto += nodes_mortar_lobatto[f];
  }

double* proj_M_u_using_quadrature = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
double* u_m_on_f_m_proj = P4EST_ALLOC(double, total_side_nodes_m_lobatto);
double* M_u_m_on_f_m_proj = P4EST_ALLOC(double, total_side_nodes_m_lobatto);

  d4est_mortars_project_mass_mortar_onto_side
    (
     d4est_ops,
      M_u_using_quadrature,
     faces_mortar,
     deg_mortar_lobatto,
     proj_M_u_using_quadrature,
     faces_m,
     deg_m_lobatto
    );
  
  d4est_mortars_project_mortar_onto_side
  (
   d4est_ops,
   u_m_on_f_m_mortar_lobatto,
   faces_mortar,
   deg_mortar_lobatto,
   u_m_on_f_m_proj,
   faces_m,
   deg_m_lobatto
  );

  stride_lobatto = 0;
  for (int f = 0; f < faces_m; f++){
    d4est_operators_apply_mij(
                              d4est_ops,
                              &u_m_on_f_m_proj[stride_lobatto],
                              (P4EST_DIM)-1,
                              deg_m_lobatto[f],
                              &M_u_m_on_f_m_proj[stride_lobatto]
    );
    stride_lobatto += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_m_lobatto[f]);
  }

  for (int i = 0; i < total_side_nodes_m_lobatto; i++){
    double ratio = M_u_m_on_f_m_proj[i]/ proj_M_u_using_quadrature[i];
    if (fabs(ratio - 4.0) > D4EST_REAL_EPS && fabs(ratio - 16.0) > D4EST_REAL_EPS && fabs( M_u_m_on_f_m_proj[i] - proj_M_u_using_quadrature[i]) > D4EST_REAL_EPS){
      printf("M_u_m_on_f_m_proj, proj_M_u_using_quadrature, ratio = %.25f, %.25f, %.25f\n",
             M_u_m_on_f_m_proj[i], proj_M_u_using_quadrature[i], M_u_m_on_f_m_proj[i]/ proj_M_u_using_quadrature[i]);

      D4EST_ABORT(" if (fabs(ratio - 4.0) > D4EST_REAL_EPS && fabs(ratio - 16.0) > D4EST_REAL_EPS) failed ");
    }
  }


  /* DEBUG_PRINT_2ARR_DBL(M_u_m_on_f_m_proj, proj_M_u_using_quadrature, total_side_nodes_m_lobatto);  */
  
  P4EST_FREE(u_m_on_f_m_proj);
  P4EST_FREE(M_u_m_on_f_m_proj);
  P4EST_FREE(proj_M_u_using_quadrature);
  P4EST_FREE(u_m_on_f_m_mortar_lobatto);
  P4EST_FREE(M_u_using_quadrature);
  P4EST_FREE(M_u_using_mass);
    
}

void
testd4est_mass_projection_on_interfaces
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_elliptic_data_t* prob_vecs,
 double err_tol
)
{
  d4est_poisson_flux_data_t* d4est_poisson_flux_data = P4EST_ALLOC(d4est_poisson_flux_data_t,1);
  testd4est_mass_projection_data_t* data = P4EST_ALLOC(testd4est_mass_projection_data_t, 1);
  data->quad_vs_mass_err = 0.;
  
  d4est_poisson_flux_data->user = data;
  d4est_poisson_flux_data->interface_fcn = testd4est_mass_projection_interface;
  d4est_poisson_flux_data->boundary_fcn = NULL;
  d4est_poisson_flux_data->boundary_condition = NULL;
  d4est_poisson_flux_data->destroy = NULL;


  d4est_poisson_flux_init_element_data(p4est, d4est_ops, prob_vecs->u, prob_vecs->Au);
  
  d4est_mortars_fcn_ptrs_t flux_fcns = d4est_poisson_flux_fetch_fcns(d4est_poisson_flux_data);
  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     &flux_fcns,
     EXCHANGE_GHOST_DATA
    );

  printf("data->quad_vs_mass_err = %.25f\n", data->quad_vs_mass_err);
  int err = (data->quad_vs_mass_err < err_tol);
  P4EST_FREE(d4est_poisson_flux_data);
  P4EST_FREE(data);
  D4EST_ASSERT(err);
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
  d4est_geom->geom_type = GEOM_BRICK;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_brick_new(proc_rank, "testd4est_poisson_1_brick.input", "geometry", "[Geometry]:", d4est_geom);
    
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    1,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");

  
  d4est_poisson_flux_data_t* flux_data = d4est_poisson_flux_new(p4est, "testd4est_poisson_1_brick.input", zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_with_bc = d4est_poisson_flux_new(p4est, "testd4est_poisson_1_brick.input", poly_vec_fcn, NULL);
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr = d4est_amr_init(
                                          p4est,
                                          "testd4est_poisson_1_brick.input",
                                          NULL
  );


  int local_nodes = d4est_mesh_update
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
                     problem_set_degrees_init,
                     NULL
                    );
  
  double* poly_vec = P4EST_ALLOC(double, local_nodes);
  int same = 1;
  int same2 = 1;
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps; ++level){

    local_nodes = d4est_mesh_update
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
                   problem_set_degrees_amr,
                   NULL
                  );

    
    printf("level = %d, elements = %d, nodes = %d\n", level, p4est->local_num_quadrants, local_nodes);

    if (level == 0){
      d4est_mesh_init_field(p4est, poly_vec, poly_vec_fcn, d4est_ops, d4est_geom, NULL);
    }
    else {

      double* Apoly_vec = P4EST_ALLOC(double, local_nodes);
      d4est_elliptic_data_t elliptic_data;
      elliptic_data.u = poly_vec;
      elliptic_data.Au = Apoly_vec;
      elliptic_data.local_nodes = local_nodes;

      testd4est_mass_projection_on_interfaces
        (
         p4est,
         ghost,
         ghost_data,
         d4est_ops,
         d4est_geom,
         d4est_quad,
         &elliptic_data,
         (D4EST_REAL_EPS)*100
        );
      
      P4EST_FREE(Apoly_vec);
    }

    
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       &poly_vec,
       NULL
      );
    
  }


  P4EST_FREE(poly_vec);

  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_poisson_flux_destroy(flux_data);
  d4est_poisson_flux_destroy(flux_data_with_bc);
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();

}
