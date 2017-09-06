#include <pXest.h>
#include <d4est_initial_grid_input.h>
#include <problem.h>
#include <d4est_geometry.h>
#include <d4est_element_data.h>
#include <petscsnes.h>


int main(int argc, char *argv[])
{  
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
#ifndef NDEBUG
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE ON\n");
  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE OFF\n");
  p4est_init(NULL, SC_LP_ERROR);
#endif
  
#if (P4EST_DIM)==3
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 3\n");
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
#endif

  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] : "options.input");
 
 
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] : "options.input",
                                                    "geometry",
                                                    "[D4EST_GEOMETRY]");

  d4est_initial_grid_t initial_grid_input = d4est_initial_grid_parse((argc == 2) ? argv[1] : "options.input");
  
  p4est_t* p4est = p4est_new_ext
    (
     mpicomm,
     d4est_geom->p4est_conn,
     initial_grid_input.min_quadrants,
     initial_grid_input.min_level,
     initial_grid_input.fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
  

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);

  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);
   

  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
    printf("[D4EST_INFO]: min_quadrants = %d\n", initial_grid_input.min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", initial_grid_input.min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", initial_grid_input.fill_uniform);
  }
  if(initial_grid_input.print_elements_per_proc){
    sc_MPI_Barrier(mpicomm);
    printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
    sc_MPI_Barrier(mpicomm);
  }
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);  
  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] : "options.input", "quadrature", "[QUADRATURE]");
  
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
                     d4est_initial_grid_degree_init,
                     (void*)&initial_grid_input
                    );
  
  /* Solve Problem */
  problem_init
    (
     p4est,
     &ghost,
     &ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     geometric_factors,
     local_nodes,
     (argc == 2) ? argv[1] : "options.input",
     mpicomm
    );

  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
  
  
  d4est_ops_destroy(d4est_ops);
  
  /* free pXest */
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
