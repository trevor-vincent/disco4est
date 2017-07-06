#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_disk.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <util.h>
#include <limits.h>


#define DEG_LOBATTO 2
#define DEG_QUAD 2

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

static void
problem_set_degrees
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  /* d4est_element_data_t* elem_data = elem_data_tmp; */
  elem_data->deg = DEG_LOBATTO;
  elem_data->deg_quad = DEG_QUAD;
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
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_CUBED_SPHERE_OUTER_SHELL;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_cubed_sphere_attr_t* sphere_attrs
    = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = 10;
  sphere_attrs->R2 = 1000000;
  sphere_attrs->compactify_outer_shell = 1;
  sphere_attrs->compactify_inner_shell = -1;

  d4est_geometry_disk_attr_t* disk_attrs
    = P4EST_ALLOC(d4est_geometry_disk_attr_t, 1);

  disk_attrs->R0 = -1;
  disk_attrs->R1 = 1;
  disk_attrs->R2 = 1000;
  disk_attrs->compactify_outer_wedge = 1;
  
  
  if ((P4EST_DIM)==3){
    d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geom, sphere_attrs);
  }
  else if ((P4EST_DIM)==2){
    d4est_geometry_disk_outer_wedge_new_aux(d4est_geom, disk_attrs);
  }
  else {
    mpi_abort("DIM == 2 or 3 only");
  }
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  /* d4est_quadrature_legendre_new(p4est, d4est_ops, d4est_geom, d4est_quad, "", ""); */
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");


  
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
        printf("elem_data->deg = %d\n", elem_data->deg);
        printf("elem_data->deg = %d\n", elem_data->deg_quad);
      }
    }
  
  

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
