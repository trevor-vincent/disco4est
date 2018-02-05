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
#include <d4est_norms.h>
#include <limits.h>

#define D4EST_REAL_EPS 100*1e-15
#if (P4EST_DIM)==2
#define TEST_DEG_INIT 2
#else
#define TEST_DEG_INIT 2
#endif

static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  elem_data->deg = 2;
  elem_data->deg_vol_quad = TEST_DEG_INIT;
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

#if (P4EST_DIM)==3
  return x*x + y*y + z*z;
#else
  return x*x + y*y;
#endif
}


double
boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void *user
)
{
  return poly_vec_fcn(x,
                      y,
#if(P4EST_DIM)==3
                      z,
#endif
                      user);
}

double
neg_laplacian_poly_vec_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
#if (P4EST_DIM)==3
  return -6.;
#else
  return -4.;
#endif
}


/*  */

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

#if (P4EST_DIM)==3
  const char* input_file = "test_d4est_poisson_consistency_curved_3d.input";
#else
  const char* input_file = "test_d4est_poisson_consistency_curved_2d.input";
#endif
  
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                     input_file,
                                                    "geometry",
                                                    "[D4EST_GEOMETRY]");
  /*  */  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");

  
  d4est_poisson_flux_data_t* flux_data_with_homog_bc = d4est_poisson_flux_new(p4est, input_file, zero_fcn, NULL);
  d4est_poisson_flux_data_t* flux_data_with_bc = d4est_poisson_flux_new(p4est, input_file, poly_vec_fcn, NULL);
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);
  d4est_amr_t* d4est_amr = d4est_amr_init(
                                          p4est,
                                          input_file,
                                          "[TEST_D4EST_POISSON_CONSISTENCY_CURVED]:",
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
      /* d4est_mesh_init_field_ext */
      /*   ( */
      /*    p4est, */
      /*    poly_vec, */
      /*    poly_vec_fcn_ext, */
      /*    NULL, */
      /*    d4est_ops, */
      /*    d4est_geom */
      /*   ); */
      
    }
    else {
      /* double* poly_vec_compare = P4EST_ALLOC(double, local_nodes); */
      /* d4est_mesh_init_field(p4est, poly_vec_compare, poly_vec_fcn, d4est_ops, d4est_geom, NULL); */
      /* same = d4est_util_compare_vecs(poly_vec, poly_vec_compare, local_nodes, D4EST_REAL_EPS); */
      /* if (!same){ */
      /*   double biggest_poly_err; */
      /*   int biggest_poly_id; */
      /*   d4est_util_find_biggest_error(poly_vec, poly_vec_compare, local_nodes, &biggest_poly_err, &biggest_poly_id); */
      /*   printf("Apparently poly_vec and poly_vec_compare aren't the same, biggest err = %.25f\n", biggest_poly_err); */
      /*   /\* DEBUG_PRINT_2ARR_DBL(poly_vec, poly_vec_compare, local_nodes); *\/ */
      /* } */
      /* P4EST_FREE(poly_vec_compare); */

      d4est_mesh_init_field(p4est, poly_vec, poly_vec_fcn, d4est_ops, d4est_geom, NULL);
      
      double* Apoly_vec = P4EST_ALLOC(double, local_nodes);
      double* Abc_poly_vec = P4EST_ALLOC(double, local_nodes);
      double* Apoly_vec_compare = P4EST_ALLOC(double, local_nodes);
      double* tmp = P4EST_ALLOC(double, local_nodes);
      d4est_elliptic_data_t elliptic_data;
      elliptic_data.u = poly_vec;
      elliptic_data.Au = Apoly_vec;
      elliptic_data.local_nodes = local_nodes;
  
      d4est_poisson_apply_aij
        (
         p4est,
         ghost,
         ghost_data,
         &elliptic_data,
         flux_data_with_bc,
         d4est_ops,
         d4est_geom,
         d4est_quad
        );

      /* int local_nodes = local_nodes; */
      double* f = P4EST_ALLOC(double, local_nodes);

      d4est_mesh_init_field
        (
         p4est,
         f,
         neg_laplacian_poly_vec_fcn,
         d4est_ops,
         d4est_geom,
         NULL
        );
       
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
               &f[ed->nodal_stride],
               ed->deg,
               ed->J_quad,
               ed->deg_quad,
               &Apoly_vec_compare[ed->nodal_stride]
              );
          }
        }
          P4EST_FREE(f);
      

        double biggest_err;
        int biggest_id;
        d4est_util_find_biggest_error(Apoly_vec, Apoly_vec_compare, local_nodes, &biggest_err, &biggest_id);
        int element = d4est_mesh_debug_find_node(p4est, biggest_id);
        printf("biggest_err, biggest_id, element = %f, %d, %d\n", biggest_err, biggest_id, element);
        printf("Apoly_vec[i], Apoly_vec_compare[i] = %.25f, %.25f\n",
               Apoly_vec[biggest_id], Apoly_vec_compare[biggest_id]);

        same2 = d4est_mesh_compare_two_fields(p4est,
                                      Apoly_vec,
                                      Apoly_vec_compare,
                                      "Apoly_vec and Apoly_vec_compare = ",
                                      DISCARD_NOTHING,
                                      PRINT_ON_ERROR,
                                      D4EST_REAL_EPS
                                     );


      double* error = P4EST_ALLOC(double, local_nodes);
      d4est_mesh_compute_point_error
        (
         Apoly_vec,
         Apoly_vec_compare,
         error,
         local_nodes
        );
      
      d4est_output_vtk
        (p4est,
         d4est_ops,
         d4est_geom,
         Apoly_vec,
         Apoly_vec_compare,
         error,
         input_file,
         "test_d4est_poisson_2_cubed_sphere",
         local_nodes,
         level,0);

      P4EST_FREE(error);
        
      
      P4EST_FREE(Apoly_vec);
      P4EST_FREE(Abc_poly_vec);
      P4EST_FREE(Apoly_vec_compare);
      P4EST_FREE(tmp);
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
    
  d4est_poisson_flux_destroy(flux_data_with_homog_bc);  
  d4est_poisson_flux_destroy(flux_data_with_bc);  
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();
  /* sc_finalize (); */
  if (same && same2)
    return 0;
  else
    return 1;
}
