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
#include <zlog.h>


#define NUM_OF_TRIALS 5
#define D4EST_REAL_EPS 100*1e-15


typedef struct {
  
  int deg;
  int deg_volume_quad;
  int deg_mortar_quad;
  
} test_d4est_stiffness_t;

static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  test_d4est_stiffness_t* data = user_ctx;
  elem_data->deg = data->deg;
  elem_data->deg_vol_quad = data->deg_volume_quad;
}

static int
get_deg_mortar_quad
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  test_d4est_stiffness_t* data = user_ctx;
  return data->deg_mortar_quad;
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
  /* d4est_geometry_t* d4est_geom = user; */
  /* return 1./(x*x + y*y + z*z); */
  return x*x + y*y + z*z;
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
  /* return -2./pow(x*x + y*y + z*z,2); */
  return -6.;
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

static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  test_d4est_stiffness_t* data = user_ctx;
  elem_data->deg = data->deg;
  elem_data->deg_vol_quad = data->deg_volume_quad;
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

  const char* input_file = "test_d4est_stiffness.input";
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                     input_file,
                                                    "geometry",
                                                    c_geom);

  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    1,
                    1
                   );
  
  test_d4est_stiffness_t deg_data;
  deg_data.deg = atoi(argv[1]);
  /* deg_data.deg_mortar_quad = atoi(argv[2]); */
  /* deg_data.deg_quad = atoi(argv[2]); */
  D4EST_ASSERT(argc == 2);
 
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");

  
  d4est_poisson_flux_data_t* flux_data_with_bc = d4est_poisson_flux_new(p4est, input_file, poly_vec_fcn, NULL, get_deg_mortar_quad, &deg_data);
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  d4est_amr_t* d4est_amr = d4est_amr_init(
                                          p4est,
                                          input_file,
                                          NULL
  );


  

  int same = 1;
  int same2 = 1;

  double stiff_last_node_error [NUM_OF_TRIALS];
  double stiff_last_node [NUM_OF_TRIALS];

  double mortar_last_node_error [NUM_OF_TRIALS];
  double mortar_last_node [NUM_OF_TRIALS];

  double rhs_last_node_error [NUM_OF_TRIALS];
  double rhs_last_node [NUM_OF_TRIALS];
  
  for (int level = 0; level < d4est_amr->num_of_amr_steps; ++level){

    
    for (int i = 0; i < NUM_OF_TRIALS; i++){
      deg_data.deg_volume_quad = deg_data.deg + i;
      for (int j = 0; j < NUM_OF_TRIALS; j++){
        deg_data.deg_mortar_quad = deg_data.deg + j;
        D4EST_ASSERT(deg_data.deg > 0 && deg_data.deg_volume_quad >= deg_data.deg && deg_data.deg_mortar_quad >= deg_data.deg);

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
                     problem_set_degrees_amr,
                     &deg_data
                    );

    

      double* poly_vec = P4EST_ALLOC(double, local_nodes);
      d4est_mesh_init_field(p4est, poly_vec, poly_vec_fcn, d4est_ops, d4est_geom, d4est_geom);
      
      double* Apoly_vec = P4EST_ALLOC(double, local_nodes);
      double* Abc_poly_vec = P4EST_ALLOC(double, local_nodes);
      double* Apoly_vec_compare = P4EST_ALLOC(double, local_nodes);
      double* tmp = P4EST_ALLOC(double, local_nodes);
      d4est_elliptic_data_t elliptic_data;
      elliptic_data.u = poly_vec;
      elliptic_data.Au = Apoly_vec;
      elliptic_data.local_nodes = local_nodes;

      d4est_poisson_flux_init_element_data
        (
         p4est,
         d4est_ops,
         elliptic_data.u,
         elliptic_data.Au
        );
      
      d4est_poisson_apply_stiffness_matrix
        (
         p4est,
         d4est_ops,
         d4est_geom,
         d4est_quad
        );
      
      stiff_last_node[i] = Apoly_vec[local_nodes-1];
      if(i == 0){
        stiff_last_node_error[i] = -1.;
      }
      else{
        stiff_last_node_error[i] = fabs(stiff_last_node[i] - stiff_last_node[i-1]);
      }

      d4est_linalg_fill_vec(Apoly_vec, 0., local_nodes);

      d4est_poisson_apply_mortar_matrices
        (
         p4est,
         ghost,
         ghost_data,
         flux_data_with_bc,
         d4est_ops,
         d4est_geom,
         d4est_quad
        );

      mortar_last_node[j] = Apoly_vec[local_nodes-1];
      if(j == 0){
        mortar_last_node_error[j] = -1.;
      }
      else{
        mortar_last_node_error[j] = fabs(mortar_last_node[j] - mortar_last_node[j-1]);
      }
      /* printf("elem %d deg %d deg_vol_quad %d deg_mortar_quad %d stiff_last %.15f mortar_last = %.15f stifferr %.15f morterr %.15f\n", p4est->local_num_quadrants, deg_data.deg, deg_data.deg_volume_quad, deg_data.deg_mortar_quad, stiff_last_node[i],  mortar_last_node[j], stiff_last_node_error[i], mortar_last_node_error[j]); */


      /* int local_nodes = local_nodes; */
      double* f = P4EST_ALLOC(double, local_nodes);

      d4est_mesh_init_field
        (
         p4est,
         f,
         neg_laplacian_poly_vec_fcn,
         d4est_ops,
         d4est_geom,
         d4est_geom
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
               ed->deg_vol_quad,
               &Apoly_vec_compare[ed->nodal_stride]
              );
          }
        }
      P4EST_FREE(f);


      rhs_last_node[i] = Apoly_vec_compare[local_nodes-1];
      if(j == 0){
        rhs_last_node_error[i] = -1.;
      }
      else{
        rhs_last_node_error[i] = fabs(rhs_last_node[i] - rhs_last_node[i-1]);
      }

      d4est_poisson_apply_aij(p4est,ghost, ghost_data, &elliptic_data, flux_data_with_bc, d4est_ops, d4est_geom, d4est_quad);
      double* cons_error = P4EST_ALLOC(double, local_nodes);
      d4est_linalg_vec_axpyeqz(-1., Apoly_vec, Apoly_vec_compare, cons_error, local_nodes);
      double l2_cons_err = d4est_mesh_compute_l2_norm_sqr(p4est,d4est_ops, d4est_geom, d4est_quad, cons_error, local_nodes, DO_NOT_STORE_LOCALLY);
      P4EST_FREE(cons_error);
      printf("elem %d deg %d deg_v %d deg_m %d stiffl %.15f mortl %.15f rhsl %.15f serr %.15f merr %.15f rerr %.15f conserr %.15f\n", p4est->local_num_quadrants, deg_data.deg, deg_data.deg_volume_quad, deg_data.deg_mortar_quad, stiff_last_node[i],  mortar_last_node[j], rhs_last_node[i], stiff_last_node_error[i], mortar_last_node_error[j], rhs_last_node_error[i], l2_cons_err);

      P4EST_FREE(Apoly_vec);
      P4EST_FREE(Abc_poly_vec);
      P4EST_FREE(Apoly_vec_compare);
      P4EST_FREE(tmp);
      P4EST_FREE(poly_vec);
    }
    }

    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr,
       NULL,
       NULL
      );
    
  }




  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
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
