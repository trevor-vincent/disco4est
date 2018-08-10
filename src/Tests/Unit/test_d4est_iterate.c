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

typedef struct {

  double* is_corner;
  double* is_edge;
  int element;
  int process;
  int corner;
  int edge;
  
} corner_data_t;


/* static void */
/* iter_edge_callback */
/* ( */
/*  p8est_iter_edge_info_t* info, */
/*  void* user_data */
/* ) */
/* { */
/*   int good_edge = -1; */
/*   p4est_t* p4est = info->p4est; */
/*   corner_data_t* edge_data = user_data; */

/*   for (int side = 0; side < info->sides.elem_count; side++){ */
/*     p8est_iter_edge_side_t* edge_side */
/*       = sc_array_index(&info->sides, */
/*                        side); */
/*     d4est_element_data_t* ed = edge_side->quad->p.user_data; */
/*     printf("element id = %d, edge = %d, side = %d\n", ed->id, edge_side->edge, side); */

/*     if (ed->id == edge_data->element && */
/*         p4est->mpirank == edge_data->process && */
/*         edge_side->edge == edge_data->edge){ */

/*       good_edge = side; */
/*       break; */
/*     } */
/*   } */

/*   for (int side = 0; side < info->sides.elem_count; side++){ */
/*     p8est_iter_edge_side_t* edge_side */
/*       = sc_array_index(&info->sides, */
/*                        side); */
/*     d4est_element_data_t* ed */
/*       = edge_side->quad->p.user_data; */

/*     if (side != good_edge && good_edge >= 0){ */
/*       printf("edge = %d\n", edge_side->edge); */
/*       edge_data->is_edge[ed->id] = edge_side->edge; */
/*     } */
/*   }     */

/* } */
 

static void
iter_corner_callback
(
 p4est_iter_corner_info_t* info,
 void* user_data
)
{
  int good_corner = -1;
  p4est_t* p4est = info->p4est;
  corner_data_t* corner_data = user_data;

  for (int side = 0; side < info->sides.elem_count; side++){
    p4est_iter_corner_side_t* corner_side
      = sc_array_index(&info->sides,
                       side);
    d4est_element_data_t* ed = corner_side->quad->p.user_data;
    printf("element id = %d, corner = %d, side = %d\n", ed->id, corner_side->corner, side);

    if (ed->id == corner_data->element &&
        p4est->mpirank == corner_data->process &&
       corner_side->corner == corner_data->corner){

      printf("\n*** This is the selected corner ****\n");
      printf("ed->id = %d\n", ed->id);
      printf("corner = %d\n", corner_side->corner);
      printf("faces[0] = %d\n", corner_side->faces[0]);
      printf("faces[1] = %d\n", corner_side->faces[1]);
      printf("faces[2] = %d\n", corner_side->faces[2]);
      printf("edges[0] = %d\n", corner_side->edges[0]);
      printf("edges[1] = %d\n", corner_side->edges[1]);
      printf("edges[2] = %d\n\n", corner_side->edges[2]);
      good_corner = side;
      break;
    }
  }

  for (int side = 0; side < info->sides.elem_count; side++){
    p4est_iter_corner_side_t* corner_side
      = sc_array_index(&info->sides,
                       side);
    d4est_element_data_t* ed
      = corner_side->quad->p.user_data;

    if (side != good_corner && good_corner >= 0){

      corner_data->is_corner[ed->id] = corner_side->corner;

      printf("\n**** This is a corner that touches selected ****\n");
      printf("ed->id = %d\n", ed->id);
      printf("corner = %d\n", corner_side->corner);
      printf("faces[0] = %d\n", corner_side->faces[0]);
      printf("faces[1] = %d\n", corner_side->faces[1]);
      printf("faces[2] = %d\n", corner_side->faces[2]);
      printf("edges[0] = %d\n", corner_side->edges[0]);
      printf("edges[1] = %d\n", corner_side->edges[1]);
      printf("edges[2] = %d\n\n", corner_side->edges[2]);
      
    }
  }  
}

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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "test_d4est_iterate_options.input");
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "test_d4est_iterate_options.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "test_d4est_iterate_options.input", d4est_geom);

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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "test_d4est_iterate_options.input", "quadrature");
  

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

  initial_grid_input->initial_nodes = local_sizes.local_nodes;

  /* d4est_amr_t* d4est_amr_random = d4est_amr_init_random_hp(p4est, 1); */

  /* d4est_amr_step */
  /*   ( */
  /*    p4est, */
  /*    NULL, */
  /*    NULL, */
  /*    d4est_ops, */
  /*    d4est_amr_random, */
  /*    NULL, */
  /*    NULL, */
  /*    NULL */
  /*   ); */

  /* p4est_partition(p4est, 1, NULL); */
  /* p4est_balance (p4est, P4EST_CONNECT_FULL, NULL); */

  /* local_sizes = d4est_mesh_update */
  /*                 ( */
  /*                  p4est, */
  /*                  &d4est_ghost, */
  /*                  d4est_ops, */
  /*                  d4est_geom, */
  /*                  d4est_quad, */
  /*                  d4est_factors, */
  /*                  initial_grid_input, */
  /*                  INITIALIZE_GHOST, */
  /*                  INITIALIZE_QUADRATURE_DATA, */
  /*                  INITIALIZE_GEOMETRY_DATA, */
  /*                  INITIALIZE_GEOMETRY_ALIASES, */
  /*                  d4est_mesh_set_quadratures_after_amr, */
  /*                  (void*)initial_grid_input */
  /*                 ); */

  int nodes = local_sizes.local_nodes;
  double* sinvec = P4EST_ALLOC(double, nodes);
  double* element_id = P4EST_ALLOC(double, p4est->local_num_quadrants);

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
        element_id[k] = ed->id;
      }
    }
  
  d4est_mesh_init_field
    (
     p4est,
     sinvec,
     init_sinxyz,
     d4est_ops,
     d4est_geom,
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     NULL
    );
  
  corner_data_t corner_data;
  corner_data.is_corner = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);
  corner_data.is_edge = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);

  for (int i = 0; i < p4est->local_num_quadrants; i++){
    corner_data.is_corner[i] = -1;
  }
  
  corner_data.element = 14;
  corner_data.process = 0;
  corner_data.corner = 6;
  corner_data.edge = 1;

  double* test = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);
  
  p8est_iterate (p4est,
                 d4est_ghost->ghost,
                 &corner_data,
                 NULL,
                 NULL,
                 NULL,
                 iter_corner_callback);


  DEBUG_PRINT_ARR_DBL(corner_data.is_corner, p4est->local_num_quadrants);

  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "test_d4est_iterate_options.input",
     "d4est_vtk",
     (const char*[]){"sinvec", NULL},
     (double**)((const double*[]){sinvec, NULL}),
     (const char*[]){"is_corner", "element_id", NULL},
     (double**)((const double*[]){corner_data.is_corner, element_id, NULL}),
     -1
    );  

  P4EST_FREE(test);
  P4EST_FREE(corner_data.is_corner);
  P4EST_FREE(corner_data.is_edge);
  P4EST_FREE(element_id);
  P4EST_FREE(sinvec);
  /* d4est_amr_destroy(d4est_amr_random); */
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

