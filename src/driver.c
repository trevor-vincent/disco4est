#include <version.h>
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
  if (zlog_init("logging.conf") != 0){
    printf("Initializing logging failed.\n");
    D4EST_ABORT("");
  }
  zlog_category_t *c_default = zlog_get_category("d4est");

  char* input_file = P4EST_ALLOC(char, 200);
  /* sprintf(input_file, "%s", (argc == 2) ? argv[1] : "options.input"); */
  sprintf(input_file, "%s", (argc == 2) ? argv[1] : "options.input");

  
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

    zlog_info(c_default, "git commit hash = %s",
              GIT_COMMIT_HASH);
    
#if (P4EST_DIM)==3
    zlog_info(c_default, "DIM = 3");
#else
    zlog_info(c_default, "DIM = 2");
#endif

    zlog_info(c_default, "options file = %s", input_file);
    zlog_info(c_default, "mpisize = %d", proc_size);
    zlog_info(c_default, "------");
  }

  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse(input_file, d4est_geom);

  p4est_t* p4est;
  d4est_operators_t* d4est_ops = d4est_ops_init(initial_grid_input->max_degree);
  
  if (initial_grid_input->load_from_checkpoint == 0 ||
      (
       initial_grid_input->load_from_checkpoint == 1 &&
       initial_grid_input->checkpoint_type == D4EST_CHKPT_HISTORY_H5
      )
     ){
    
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

    if (initial_grid_input->refine_per_region_at_start){
      for (int i = 0; i < initial_grid_input->number_of_regions;
           i++){

        int reg = -1;
        for (int k = 0; k < initial_grid_input->number_of_regions; k++){
          if (initial_grid_input->per_region_order[k] == -1){
            D4EST_ABORT("initial_grid_input->per_region_order[k] == -1");
          }
          if (initial_grid_input->per_region_order[k] == i){
            reg = k;
          }
        }
        
        int num_of_times = initial_grid_input->per_region_number_of_refines[reg];
        if (num_of_times == -1){
          D4EST_ABORT("num_of_times == -1");
        }
        int partition = initial_grid_input->per_region_partition[reg];
        if (partition == -1){
          D4EST_ABORT("partition == -1");
        }
                
        for (int j = 0; j < num_of_times; j++){
          void* tmp = p4est->user_pointer;
          d4est_mesh_refine_in_region_data_t rird;
          rird.d4est_geom = d4est_geom;
          rird.region = reg;        
          p4est->user_pointer = &rird;
          p4est_refine(
                       p4est,
                       0,
                       d4est_mesh_refine_in_region_callback,
                       NULL
                      );

          p4est->user_pointer = tmp;
        }

        if (partition){
          if (initial_grid_input->keep_quad_fams_together){
            p4est_partition(p4est, 1, NULL);
          }
          else {
            p4est_partition(p4est, 0, NULL);
          }          
        }
        
      }
    }
    else {
      if (initial_grid_input->keep_quad_fams_together){
        p4est_partition(p4est, 1, NULL);
      }
      else {
        p4est_partition(p4est, 0, NULL);
      }
    }
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  }

  else if (initial_grid_input->load_from_checkpoint == 1 &&
           initial_grid_input->checkpoint_type == D4EST_CHKPT_P4EST_H5){

    /* D4EST_ABORT("This type of checkpoint may not work anymore, so we abort"); */
    
    p4est = d4est_checkpoint_load_p4est_from_file
            (
             mpicomm,
             initial_grid_input,
             &d4est_geom->p4est_conn
            );

    zlog_info(c_default, "Successfully loaded checkpoint mesh");

    p4est_reset_data(p4est,
                     sizeof(d4est_element_data_t),
                     NULL,
                     NULL
                    );
    
    initial_grid_input->checkpoint_deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_h5_read_dataset(
                          p4est->mpirank,
                          initial_grid_input->checkpoint_prefix,
                          "degree",
                          H5T_NATIVE_INT,
                          initial_grid_input->checkpoint_deg_array
                         );
    zlog_info(c_default, "Successfully read checkpoint degrees");
  }
  if (initial_grid_input->load_from_checkpoint == 1 &&
           initial_grid_input->checkpoint_type == D4EST_CHKPT_HISTORY_H5){
    
    printf("initial_checkpoint_number = %d\n",initial_grid_input->initial_checkpoint_number);
    printf("checkpoint_number = %d\n",initial_grid_input->checkpoint_number);

    /* we need to zero this out for the initial update because 
     * d4est_mesh_set_initial_extents needs to use the load_from_checkpoint == 0 branch
     * in order to set up the initial grid */
    d4est_mesh_update
      (
       p4est,
       NULL,
       d4est_ops,
       d4est_geom,
       NULL,
       NULL,
       initial_grid_input,
       DO_NOT_INITIALIZE_GHOST,
       DO_NOT_INITIALIZE_QUADRATURE_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_DATA,
       DO_NOT_INITIALIZE_GEOMETRY_ALIASES,
       d4est_mesh_set_initial_extents_using_input_file_regions,
       (void*)initial_grid_input
      );
    
    d4est_checkpoint_load_mesh_from_amr_history
      (
       mpicomm,
       p4est,
       d4est_ops,
       d4est_geom,
       initial_grid_input,
       input_file
      );   
    
  }
  else {
    zlog_error(c_default, "Checkpoint parameters not set or something eles is wrong");
    D4EST_ABORT("");
  }
  
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

  if(proc_rank == 0)
    zlog_debug(c_default, "max_degree = %d", initial_grid_input->max_degree);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature");  
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
     input_file,
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

  P4EST_FREE(input_file);
  return 0;
}
