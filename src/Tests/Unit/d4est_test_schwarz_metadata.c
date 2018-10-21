#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_operators.h>
#include <d4est_amr.h>
#include <d4est_amr_random.h>
#include <d4est_vtk.h>
#include <d4est_h5.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <d4est_solver_schwarz.h>
#include <petscsnes.h>
#include <zlog.h>

static void 
d4est_test_iterate_refine_only_one_element
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) info->p4est->user_pointer;
  int* element_to_refine = (int*) (d4est_amr->scheme->amr_scheme_data);
  d4est_element_data_t* elem_data = (d4est_element_data_t*) info->quad->p.user_data;
  if (elem_data->id == *element_to_refine){
    d4est_amr->refinement_log[elem_data->id] = -elem_data->deg;
  }
  else {
    d4est_amr->refinement_log[elem_data->id] = elem_data->deg;
  }
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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_schwarz_metadata.input");
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_schwarz_metadata.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_schwarz_metadata.input", d4est_geom);

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
  

  printf("num_quadrants = %d\n", p4est->local_num_quadrants);

  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  
  d4est_ghost_t* d4est_ghost = NULL;
   

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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_schwarz_metadata.input", "quadrature");
  

  d4est_mesh_local_sizes_t local_sizes= d4est_mesh_update
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



  /* create amr scheme */
  int add_amr_in_one_element = 0;
  if (add_amr_in_one_element){
    int* refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
    int element_to_refine = 0;
    d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
    d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
    scheme->post_balance_callback = NULL;
    scheme->pre_refine_callback = NULL;
    scheme->refine_replace_callback_fcn_ptr = NULL;
    scheme->balance_replace_callback_fcn_ptr = NULL;
    scheme->mark_elements = d4est_test_iterate_refine_only_one_element;
    scheme->amr_scheme_data = &element_to_refine;
    scheme->destroy = NULL;
    
    d4est_amr->mpirank = p4est->mpirank;
    d4est_amr->scheme = scheme;
    d4est_amr->balance_log = NULL;
    d4est_amr->refinement_log = NULL;
    d4est_amr->initial_log = NULL;
    d4est_amr->max_degree = 1000;
    
    d4est_amr_step
      (
       p4est,
       NULL,
       d4est_amr,
       NULL,
       NULL,
       NULL
      );

    P4EST_FREE(refinement_log);
    P4EST_FREE(d4est_amr);
    P4EST_FREE(scheme);

    local_sizes = d4est_mesh_update
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
  
    initial_grid_input->initial_nodes = local_sizes.local_nodes;
  }

  int num_nodes_overlap = 2;
  d4est_solver_schwarz_data_t* schwarz_data
    = d4est_solver_schwarz_init(
                                p4est,
                                d4est_ghost,
                                num_nodes_overlap
                               );
  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "d4est_test_schwarz_metadata.input",
     "d4est_vtk",
     NULL,
     NULL,
     NULL,
     NULL,
     (const char*[]){NULL},
     (int**)((const int*[]){NULL}),
     -1
    );

  d4est_solver_schwarz_destroy
    (
     schwarz_data
    );

  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_mesh_data_destroy(d4est_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  if (d4est_ghost) {
    d4est_ghost_destroy(d4est_ghost);
  }
    
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  PetscFinalize();
  return 0;
}

