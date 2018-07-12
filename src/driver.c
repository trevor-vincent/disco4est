#include <pXest.h>
#include <problem.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <petscsnes.h>
#include <time.h>
#include <zlog.h>

int main(int argc, char *argv[])
{
  clock_t begin = clock();
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
  // Initialize logging
  if (zlog_init("logging.conf") != 0) 
    printf("Initializing logging failed.\n");
  zlog_category_t *c_default = zlog_get_category("d4est");
  
  if (proc_rank == 0) {
    zlog_info(c_default, "------");
  
#ifdef D4EST_PROBLEM_NAME
    zlog_info(c_default, "# Running problem %s", D4EST_PROBLEM_NAME);
#endif

#ifndef NDEBUG
    zlog_info(c_default, "DEBUG MODE ON");
#else
    zlog_info(c_default, "DEBUG MODE OFF");
#endif
  
#if (P4EST_DIM)==3
    zlog_info(c_default, "DIM = 3");
#else
    zlog_info(c_default, "DIM = 2");
#endif

    zlog_info(c_default, "options file = %s", (argc == 2) ? argv[1] : "options.input");
    zlog_info(c_default, "mpisize = %d", proc_size);
    zlog_info(c_default, "------");
  }

  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] : "options.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] : "options.input", d4est_geom);

  p4est_t* p4est;
  if (initial_grid_input->load_from_checkpoint == 0){
    
    p4est = p4est_new_ext
    (
     mpicomm,
     d4est_geom->p4est_conn,
     initial_grid_input->min_quadrants,
     initial_grid_input->min_level,
     initial_grid_input->fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
  }
  else {
    p4est = d4est_checkpoint_load_mesh
            (
             mpicomm,
             initial_grid_input->checkpoint_prefix,
             &d4est_geom->p4est_conn
            );


    p4est_reset_data(p4est,
                     sizeof(d4est_element_data_t),
                     NULL,
                     NULL
                    );
    
    initial_grid_input->checkpoint_deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_h5_read_dataset(p4est->mpirank, initial_grid_input->checkpoint_prefix, "degree", H5T_NATIVE_INT, initial_grid_input->checkpoint_deg_array);
  }

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  
  /* p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE); */
  /* d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t, */
                                                   /* ghost->ghosts.elem_count); */
   
  d4est_ghost_t* d4est_ghost = NULL;

  if (proc_rank == 0 && initial_grid_input->load_from_checkpoint == 0){
    zlog_debug(c_default, "min_quadrants = %d", initial_grid_input->min_quadrants);
    zlog_debug(c_default, "min_level = %d", initial_grid_input->min_level);
    zlog_debug(c_default, "fill_uniform = %d", initial_grid_input->fill_uniform);
  }
  
  sc_MPI_Barrier(mpicomm);
  zlog_debug(c_default, "elements on proc %d = %d", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(initial_grid_input->max_degree);
  if(proc_rank == 0)
    zlog_debug(c_default, "max_degree = %d", initial_grid_input->max_degree);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] : "options.input", "quadrature");  
  d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update
                                         (
                                          p4est,
                                          &d4est_ghost,
                                          d4est_ops,
                                          d4est_geom,
                                          d4est_quad,
                                          d4est_factors,
                                          initial_grid_input,
                                          INITIALIZE_GHOST,
                                          INITIALIZE_QUADRATURE_DATA,
                                          INITIALIZE_GEOMETRY_DATA,
                                          INITIALIZE_GEOMETRY_ALIASES,
                                          d4est_mesh_set_initial_extents,
                                          (void*)initial_grid_input
                                         );
  
  initial_grid_input->initial_nodes = local_sizes.local_nodes;
  
  /* Solve Problem */
  problem_init
    (
     p4est,
     &d4est_ghost,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     initial_grid_input,
     (argc == 2) ? argv[1] : "options.input",
     mpicomm
    );

  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (d4est_ghost) {
    d4est_ghost_destroy(d4est_ghost);
    /* p4est_ghost_destroy (ghost); */
    /* P4EST_FREE (ghost_data); */
    /* ghost = NULL; */
    /* ghost_data = NULL; */
  }
  
  
  d4est_ops_destroy(d4est_ops);
  
  /* free pXest */
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  PetscFinalize();
  /* sc_finalize (); */

  if (proc_rank == 0) {
    zlog_info(c_default, "Completed garbage collection.");
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    zlog_info(c_default, "Completed problem in %f seconds", time_spent);
  }

  zlog_fini();

  
  return 0;
}
