#include <pXest.h>
#include <problem.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_amr.h>
#include <d4est_amr_random.h>
#include <d4est_vtk.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <petscsnes.h>
#include <zlog.h>


static double
init_sinxyz
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "test_d4est_vtk_options.input");
 
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "test_d4est_vtk_options.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "test_d4est_vtk_options.input", d4est_geom);

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
  
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);
   

  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
  }
  if (proc_rank == 0 && initial_grid_input->load_from_checkpoint == 0){
    printf("[D4EST_INFO]: min_quadrants = %d\n", initial_grid_input->min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", initial_grid_input->min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", initial_grid_input->fill_uniform);
  }
  
  sc_MPI_Barrier(mpicomm);
  printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "test_d4est_vtk_options.input", "quadrature", "[QUADRATURE]");
  
  initial_grid_input->initial_nodes = d4est_mesh_update
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
                                d4est_mesh_set_initial_extents,
                                (void*)initial_grid_input
                               );

  d4est_amr_t* d4est_amr_random = d4est_amr_init_random_hp(p4est, 7, 1);

  d4est_amr_step
    (
     p4est,
     &ghost,
     &ghost_data,
     d4est_ops,
     d4est_amr_random,
     NULL,
     NULL
    );

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  int nodes = d4est_mesh_update
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
                   d4est_mesh_set_quadratures_after_amr,
                   (void*)initial_grid_input
                  );

  double* sinvec = P4EST_ALLOC(double, nodes);
  double* element_volume = P4EST_ALLOC(double, p4est->local_num_quadrants);

  int k = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q, ++k) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;
        element_volume[k] = ed->volume;
      }
    }
  
  d4est_mesh_init_field
    (
     p4est,
     sinvec,
     init_sinxyz,
     d4est_ops,
     d4est_geom,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );

  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);
    

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "test_d4est_vtk_options.input",
     "d4est_vtk",
     (const char*[]){"sinvec", NULL},
     (double**)((const double*[]){sinvec, NULL}),
     (const char*[]){"element_vol", NULL},
     (double**)((const double*[]){element_volume, NULL}),
     -1
    );
  /*  */
     
  P4EST_FREE(deg_array);
  P4EST_FREE(element_volume);
  P4EST_FREE(sinvec);
  d4est_amr_destroy(d4est_amr_random);
  

  d4est_mesh_initial_extents_destroy(initial_grid_input);
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
