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

typedef struct {

  int mpirank;
  int tree;
  int tree_id;
  int id;
  int nodal_stride;
  int faces [3];

  int is_central_element;
  int shares_what; //0 = corner, 1 = edge, 2 = face
  
} schwarz_element_t;


typedef struct{

  schwarz_element_t connections [108];

} schwarz_element_connections_t;


typedef struct {

  /* double* is_corner; */
  /* double* is_edge; */
  int element;
  int process;
  /* int corner; */
  /* int edge; */

  schwarz_element_t connections [108];
  int connections_found;
  
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
    /* printf("element id = %d, corner = %d, side = %d\n", ed->id, corner_side->corner, side); */

    if (ed->id == corner_data->element &&
        p4est->mpirank == corner_data->process &&
        corner_side->corner == corner_data->corner){

      printf("\n*** This is the selected corner ****\n");
      printf("ed->id = %d\n", ed->id);
      printf("corner = %d\n", corner_side->corner);
      printf("faces[0] = %d\n", corner_side->faces[0]);
      printf("faces[1] = %d\n", corner_side->faces[1]);
#if (P4EST_DIM)==3
      printf("faces[2] = %d\n", corner_side->faces[2]);
      printf("edges[0] = %d\n", corner_side->edges[0]);
      printf("edges[1] = %d\n", corner_side->edges[1]);
      printf("edges[2] = %d\n\n", corner_side->edges[2]);
#endif
      good_corner = side;
      break;
    }
  }

  for (int side = 0; side < info->sides.elem_count; side++){

    /* printf("info->sides.elem_count = %d\n", info->sides.elem_count); */
    p4est_iter_corner_side_t* corner_side
      = sc_array_index(&info->sides,
                       side);
    d4est_element_data_t* ed
      = corner_side->quad->p.user_data;

    if (side != good_corner && good_corner >= 0){

      corner_data->is_corner[ed->id] = corner_side->corner;

      p4est_iter_corner_side_t* good_side
        = sc_array_index(&info->sides,
                         good_corner);

      printf("\n**** This is a corner that touches selected ****\n");
      
      int found_face [] = {0,0,0};
      int found_edge [] = {0,0,0};
      for (int f = 0; f < (P4EST_DIM); f++){
        if ((good_side->faces[f] == corner_side->faces[f])){
          found_face[f] = 1;
        }
      }

#if (P4EST_DIM)==3         
      for (int e = 0; e < (P4EST_DIM); e++){
        if ((good_side->edges[e] == corner_side->edges[e])){
          found_edge[e] = 1;          
        }
      }
#endif
      
      int shares_face = 0;
      int face_id;
      int shares_edge = 0;
      int edge_id;
      int shares_corner = 1;
      int face_edge_sum = 0;
      
      for (int ef = 0; ef < (P4EST_DIM); ef++){
        if (found_face[ef] == 1){
#if (P4EST_DIM)==3
          face_id = p8est_corner_faces[corner_side->corner][ef];
#else
          face_id = p4est_corner_faces[corner_side->corner][ef];
#endif
          printf("This side shares a face that touches the corner = %d\n", face_id);
          shares_face++;
          face_edge_sum++;
          break;
        }
#if (P4EST_DIM)==3       
        if (found_edge[ef] == 1){
          printf("This side shares an edge that touches the corner = %d\n", edge_id);
   
          edge_id = p8est_corner_edges[corner_side->corner][ef];
          printf("The two faces that touch this edge are %d, %d\n",
                 p8est_edge_faces[edge_id][0],p8est_edge_faces[edge_id][1]);          



          /* corner_data->connections[connections_found].id; */
          shares_edge++;
          face_edge_sum++;
        }
#endif
      }

      
      d4est_element_data_t* ed = corner_side->quad->p.user_data; 
      int do_not_save = 0;
      for (int i = 0; i < corner_data->connections_found; i++){
        if(corner_data->connections[i].id == ed->id){
          do_not_save = 1;
        }
      }

 
      if (face_edge_sum == 0){
        printf("This side only shares a corner\n");
#if (P4EST_DIM)==3   
        printf("The three faces that share this corner are %d, %d, %d\n",
               p8est_corner_faces[corner_side->corner][0],
               p8est_corner_faces[corner_side->corner][1],
               p8est_corner_faces[corner_side->corner][2]);
#else
        printf("The two faces that share this corner are %d, %d\n",
               p4est_corner_faces[corner_side->corner][0],
               p4est_corner_faces[corner_side->corner][1]);
        /* p8est_corner_faces[corner_side->corner][2]);         */
#endif
      }

      /* if(do_not_save == 0){ */
        if (shares_face){
          corner_data->connections[corner_data->connections_found].faces[0] = face_id;
          corner_data->connections[corner_data->connections_found].faces[1] = -1;
          corner_data->connections[corner_data->connections_found].faces[2] = -1;
        }
#if (P4EST_DIM)==3
        else if (shares_edge){    
          corner_data->connections[corner_data->connections_found].faces[0] = p8est_edge_faces[edge_id][0];
          corner_data->connections[corner_data->connections_found].faces[1] = p8est_edge_faces[edge_id][1];
          corner_data->connections[corner_data->connections_found].faces[2] = -1;
        }
#endif
        else {

#if (P4EST_DIM)==3
          corner_data->connections[corner_data->connections_found].faces[0] = p8est_corner_faces[corner_side->corner][0];
          corner_data->connections[corner_data->connections_found].faces[1] = p8est_corner_faces[corner_side->corner][1];
          corner_data->connections[corner_data->connections_found].faces[2] = p8est_corner_faces[corner_side->corner][2];
#else
          corner_data->connections[corner_data->connections_found].faces[0] = p4est_corner_faces[corner_side->corner][0];
          corner_data->connections[corner_data->connections_found].faces[1] = p4est_corner_faces[corner_side->corner][1];        
#endif
        }

        corner_data->connections[corner_data->connections_found].id = ed->id;
        corner_data->connections[corner_data->connections_found].is_central_element = 0;
        corner_data->connections[corner_data->connections_found].mpirank = ed->mpi_rank;
        corner_data->connections[corner_data->connections_found].tree = ed->tree;
        corner_data->connections[corner_data->connections_found].nodal_stride = ed->nodal_stride;
        corner_data->connections[corner_data->connections_found].shares_what = (shares_face > 0) ? 2 :
                                                                               ((shares_edge > 0) ? 1 : 0);
 
      
      /* } */
      
      printf("ed->id = %d\n", corner_data->connections[corner_data->connections_found].id);
      /* printf("corner = %d\n", corner_side->corner); */
      printf("faces[0] = %d\n", corner_data->connections[corner_data->connections_found].faces[0]);
      printf("faces[1] = %d\n", corner_data->connections[corner_data->connections_found].faces[1]);
#if (P4EST_DIM)==3
      printf("faces[2] = %d\n", corner_data->connections[corner_data->connections_found].faces[2]);
      /* printf("edges[0] = %d\n", corner_data->connections[corner_data->connections_found].edges[0]); */
      /* printf("edges[1] = %d\n", corner_data->connections[corner_data->connections_found].edges[1]); */
      /* printf("edges[2] = %d\n\n", corner_data->connections[corner_data->connections_found].edges[2]); */
#endif


      corner_data->connections_found++;

      /* if(corner_data->connections_found > 26){ */
        /* D4EST_ABORT("Connections should not be > 26"); */
      /* } */
      
    }
  }  
}

/* static double */
/* init_sinxyz */
/* ( */
/*  double x, */
/*  double y, */
/* #if (P4EST_DIM)==3 */
/*  double z, */
/* #endif */
/*  void* user */
/* ) */
/* { */
/*   return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z); */
/* } */

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
    printf("[D4EST_INFO]: options file = %s\n", (argc == 2) ? argv[1] :      "d4est_test_iterate_options.input");
 
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    (argc == 2) ? argv[1] :      "d4est_test_iterate_options.input",
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse((argc == 2) ? argv[1] :      "d4est_test_iterate_options.input", d4est_geom);

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
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, (argc == 2) ? argv[1] :      "d4est_test_iterate_options.input", "quadrature");
  

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


  local_sizes= d4est_mesh_update
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

  int nodes = local_sizes.local_nodes;
  double* sinvec = P4EST_ALLOC(double, nodes);
  double* element_id = P4EST_ALLOC(double, p4est->local_num_quadrants);

  for (int i = 0; i < p4est->local_num_quadrants; i++){
    d4est_element_data_t* ed = d4est_factors->element_data[i];
    element_id[i] = ed->id;
    /* printf("ed->deg from factors = %d\n", ed->deg); */
  }
    
  corner_data_t corner_data;
  corner_data.is_corner = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);
  corner_data.is_edge = P4EST_ALLOC_ZERO(double, p4est->local_num_quadrants);

  for (int i = 0; i < p4est->local_num_quadrants; i++){
    corner_data.is_corner[i] = -1;
  }

  for (int i = 0; i < 27*4; i++){
    corner_data.connections[i].id = -1;
    corner_data.connections[i].faces[0] = -1;
    corner_data.connections[i].faces[1] = -1;
    corner_data.connections[i].faces[2] = -1;
  }
  
  corner_data.connections_found = 0;

  int schwarz_center = 10;
  /* for (int i = 0; i < (P4EST_CHILDREN); i++){ */
    
    corner_data.element = schwarz_center;
    corner_data.process = 0;
    /* corner_data.corner = i; */
    /* corner_data.edge = -1; */
#if (P4EST_DIM)==3
    p8est_iterate (p4est,
                   d4est_ghost->ghost,
                   &corner_data,
                   NULL,
                   NULL,
                   NULL,
                   iter_corner_callback);

#else
    p4est_iterate (p4est,
                   d4est_ghost->ghost,
                   &corner_data,
                   NULL,
                   /* NULL, */
                   NULL,
                   iter_corner_callback);
#endif

  /* } */
  
  for (int i = 0; i < corner_data.connections_found; i++){
    if(corner_data.connections[i].id != -1){
      printf("***********\n");
      printf("Connection %d\n",i);
      printf("Element id = %d\n", corner_data.connections[i].id);
      printf("Element tree = %d\n", corner_data.connections[i].tree);
      printf("Faces 0 = %d\n", corner_data.connections[i].faces[0]);
      printf("Faces 1 = %d\n", corner_data.connections[i].faces[1]);
      printf("Faces 2 = %d\n", corner_data.connections[i].faces[2]);
    }
  }

 
  
  int local_nodes = local_sizes.local_nodes;
  double* ones = P4EST_ALLOC_ZERO(double, local_nodes);
  for (int i = 0; i < local_nodes; i++){
    ones[i] = 1.;
  }
  double* restricted_ones = P4EST_ALLOC_ZERO(double, local_nodes);
  double* restricted_ones_transpose = P4EST_ALLOC_ZERO(double, local_nodes);
  
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

        int connection = -1;
        int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);

        for (int i = 0; i < corner_data.connections_found; i++){
          if(ed->id == corner_data.connections[i].id){
            connection = i;
          }
        }

        if (ed->id == schwarz_center){
          d4est_util_copy_1st_to_2nd(&ones[ed->nodal_stride],&restricted_ones[ed->nodal_stride],volume_nodes);
          d4est_util_copy_1st_to_2nd(&ones[ed->nodal_stride],&restricted_ones_transpose[ed->nodal_stride],volume_nodes);
        }
        else if (connection != -1){
          double* res_rst_lobatto [3];
          double* res_rst_lobatto_trans [3];
          double* rst_lobatto [3];
          for (int i = 0; i < (P4EST_DIM); i++){
            rst_lobatto[i] = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM), ed->deg, i);
            res_rst_lobatto[i] = P4EST_ALLOC(double, volume_nodes);
            res_rst_lobatto_trans[i] = P4EST_ALLOC(double, volume_nodes);
          }
          
          printf("***********\n");
          printf("Connection %d\n",connection);
          printf("Element id = %d\n", corner_data.connections[connection].id);
          printf("Element tree = %d\n", corner_data.connections[connection].tree);
          printf("Faces 0 = %d\n", corner_data.connections[connection].faces[0]);
          printf("Faces 1 = %d\n", corner_data.connections[connection].faces[1]);
          printf("Faces 2 = %d\n", corner_data.connections[connection].faces[2]);
          /* printf("P4EST_DIM = %d\n", (P4EST_DIM)); */

          /* double *onesptr = &ones[ed->nodal_stride]; */
          /* DEBUG_PRINT_ARR_DBL(onesptr, volume_nodes); */
          
          int faces [3];
          faces[0] = corner_data.connections[connection].faces[0];
          faces[1] = corner_data.connections[connection].faces[1];
          faces[2] = corner_data.connections[connection].faces[2];

          double* ones_vol = P4EST_ALLOC(double, volume_nodes);
          for (int i = 0 ; i < volume_nodes; i++)
            ones_vol[i] = 1;
          
          d4est_operators_apply_schwarz_restrictor(d4est_ops,
                                                   ones_vol,
                                                   (P4EST_DIM),
                                                   &(faces[0]),
                                                   ed->deg,
                                                   ed->deg,/* (ed->deg == 1) ? ed->deg + 1 : ed->deg, */
                                                   D4OPS_NO_TRANSPOSE,
                                                   &restricted_ones[ed->nodal_stride]);


          
          double *restrictedonesptr = &restricted_ones[ed->nodal_stride];
          int res_volume_nodes = d4est_lgl_get_nodes((P4EST_DIM),ed->deg);
          DEBUG_PRINT_ARR_DBL(restrictedonesptr, res_volume_nodes);

          
          d4est_operators_apply_schwarz_restrictor(d4est_ops,
                                                   &restricted_ones[ed->nodal_stride],
                                                   (P4EST_DIM),
                                                   &(faces[0]),
                                                   ed->deg,
                                                   ed->deg,
                                                   D4OPS_TRANSPOSE,
                                                   &restricted_ones_transpose[ed->nodal_stride]);



          double *restrictedones_transposeptr = &restricted_ones_transpose[ed->nodal_stride];
          DEBUG_PRINT_ARR_DBL(restrictedones_transposeptr, volume_nodes);
          
          for (int i = 0; i < (P4EST_DIM); i++){
            

            d4est_operators_apply_schwarz_restrictor
              (
               d4est_ops,
               rst_lobatto[i],
               (P4EST_DIM),
               &faces[0],
               ed->deg,
               ed->deg,
               D4OPS_NO_TRANSPOSE,
               res_rst_lobatto[i]
              );
          
            d4est_operators_apply_schwarz_restrictor(d4est_ops,
                                                     res_rst_lobatto[i],
                                                     (P4EST_DIM),
                                                     &(faces[0]),
                                                     ed->deg,
                                                     ed->deg,
                                                     D4OPS_TRANSPOSE,
                                                     res_rst_lobatto_trans[i]);
          
          }

#if (P4EST_DIM)==3
          DEBUG_PRINT_3ARR_DBL(res_rst_lobatto[0], res_rst_lobatto[1], res_rst_lobatto[2], volume_nodes);
          DEBUG_PRINT_3ARR_DBL(res_rst_lobatto[0], res_rst_lobatto[1], res_rst_lobatto[2], res_volume_nodes);
#else
          DEBUG_PRINT_2ARR_DBL(res_rst_lobatto[0], res_rst_lobatto[1], volume_nodes);
          DEBUG_PRINT_2ARR_DBL(res_rst_lobatto[0], res_rst_lobatto[1], res_volume_nodes);
#endif
          for (int i = 0; i < (P4EST_DIM); i++){
            P4EST_FREE(res_rst_lobatto[i]);
            P4EST_FREE(res_rst_lobatto_trans[i]);
          }
        }
        else{
          for (int i = 0; i< volume_nodes; i++){
            restricted_ones[ed->nodal_stride + i] = 0.;
            restricted_ones_transpose[ed->nodal_stride + i] = 0.;
          }
        }
      }
      
    }

  printf("This element is on bndry? = %d\n", d4est_factors->element_touches_boundary[corner_data.element]);
  
  d4est_vtk_save
    (
     p4est,
     d4est_ops,
     "d4est_test_iterate_options.input",
     "d4est_vtk",
     (const char*[]){"restricted_ones_transpose", NULL},
     (double**)((const double*[]){restricted_ones_transpose, NULL}),
     (const char*[]){"is_corner", "element_id", NULL},
     (double**)((const double*[]){corner_data.is_corner, element_id, NULL}),
     -1
    );

  P4EST_FREE(restricted_ones);
  P4EST_FREE(ones);
  P4EST_FREE(restricted_ones_transpose);
  P4EST_FREE(corner_data.is_corner);
  P4EST_FREE(corner_data.is_edge);
  P4EST_FREE(element_id);
  P4EST_FREE(sinvec);
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

