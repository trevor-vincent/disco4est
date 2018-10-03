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
#include <d4est_ghost_data.h>
#include <d4est_field.h>
#include <d4est_ghost.h>
#include <d4est_util.h>
#include <petscsnes.h>
#include <zlog.h>


static double
poly_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  return x*y*z;
}

static double
sinxyz_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return sin(M_PI*x*y*z);
}


int main(int argc, char *argv[])
{

#ifndef D4EST_TEST
  D4EST_ABORT("D4EST_TEST not defined");
#endif
  
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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_ghost_data.input");
 
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_ghost_data.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_ghost_data.input", d4est_geom);

  p4est_t* p4est;
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


  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  
  d4est_ghost_t* d4est_ghost = d4est_ghost_init(p4est);
     

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
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_ghost_data.input", "quadrature");
  


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

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* d4est_amr_t* d4est_amr_random = d4est_amr_init_uniform_h(p4est, 7, 1); */
  int num_of_amr_steps = 3;
  d4est_amr_t* d4est_amr_random = d4est_amr_init_random_hp(p4est, num_of_amr_steps);

  int nodes = -1;
  for (int i = 0; i < num_of_amr_steps; i++){
  d4est_amr_step
    (
     p4est,
     d4est_ops,
     d4est_amr_random,
     NULL,
     NULL,
     NULL
    );

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
                                          d4est_mesh_set_quadratures_after_amr,
                                          (void*)initial_grid_input
                                         );

  nodes = local_sizes.local_nodes;
  }

  d4est_field_type_t field_type [2] = {NODAL,NODAL};
  d4est_ghost_data_t* d4est_ghost_data = d4est_ghost_data_init(p4est,
                                                               d4est_ghost,
                                                               &field_type[0],
                                                               2);
  
  
  double* vecs = P4EST_ALLOC(double, 2*nodes);
  
  d4est_mesh_init_field
    (
     p4est,
     vecs,
     sinxyz_fcn,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );

  d4est_mesh_init_field
    (
     p4est,
     &vecs[nodes],
     poly_fcn,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
  
  int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
  d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);
     
  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "d4est_test_ghost_data.input",
     "d4est_vtk",
     (const char*[]){"sinvec", "polyvec", NULL},
     (double**)((const double*[]){vecs, &vecs[nodes], NULL}),
     (const char*[]){NULL},
     (double**)((const double*[]){NULL}),
     NULL,
     NULL,
     -1
    );

  /* d4est_element_data_copy_from_vec_to_storage(p4est, sinvec); */

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
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
        d4est_util_copy_1st_to_2nd(&vecs[ed->nodal_stride], &ed->test_vecs[0][0], volume_nodes);
        d4est_util_copy_1st_to_2nd(&vecs[nodes + ed->nodal_stride], &ed->test_vecs[1][0], volume_nodes);
      }
    }
  
  

  d4est_ghost_data_exchange(p4est,d4est_ghost,d4est_ghost_data,vecs);
  /* /\*  *\/ */
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];

    double* ud = d4est_ghost_data_get_field_on_element(ged,0,d4est_ghost_data);
    double* up = &ged->test_vecs[0][0];

    int size = d4est_element_data_get_size_of_field
               (
                ged,
                d4est_ghost_data->transfer_types[0]
               );


    int compare = d4est_util_compare_vecs(ud, up, size, 1e-15);
    
    double* ud1 = d4est_ghost_data_get_field_on_element(ged,1,d4est_ghost_data);
    double* up1 = &ged->test_vecs[1][0];

    int size1 = d4est_element_data_get_size_of_field
               (
                ged,
                d4est_ghost_data->transfer_types[1]
               );
  
    
    int compare1 = d4est_util_compare_vecs(ud1, up1, size1, 1e-15);
    
    if (!compare){
      DEBUG_PRINT_2ARR_DBL(ud, up, size);
    }
    if (!compare1){
      DEBUG_PRINT_2ARR_DBL(ud1, up1, size);
    }
  }
  
  d4est_ghost_destroy(d4est_ghost);
  d4est_ghost_data_destroy(d4est_ghost_data);
  

  P4EST_FREE(deg_array);
  P4EST_FREE(vecs);
  d4est_amr_destroy(d4est_amr_random);
  

  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
 
  
  d4est_ops_destroy(d4est_ops);
  
  /* free pXest */
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
 
