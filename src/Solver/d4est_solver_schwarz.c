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
d4est_solver_schwarz_corner_callback
(
 p4est_iter_corner_info_t* info,
 void* user_data
)
{
  int good_corner = -1;
  p4est_t* p4est = info->p4est;
  d4est_solver_schwarz_center_data_t* center_data = user_data;

  for (int side = 0; side < info->sides.elem_count; side++){
    p4est_iter_corner_side_t* corner_side
      = sc_array_index(&info->sides,
                       side);
    d4est_element_data_t* ed = corner_side->quad->p.user_data;

    if (ed->id == center_data->element &&
        /* p4est->mpirank == center_data->process */
       )
      {
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

      center_data->is_corner[ed->id] = corner_side->corner;

      p4est_iter_corner_side_t* good_side
        = sc_array_index(&info->sides,
                         good_corner);

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
          shares_face++;
          face_edge_sum++;
          break;
        }
#if (P4EST_DIM)==3       
        if (found_edge[ef] == 1){
          edge_id = p8est_corner_edges[corner_side->corner][ef];
          shares_edge++;
          face_edge_sum++;
        }
#endif
      }

      d4est_element_data_t* ed = corner_side->quad->p.user_data; 
      int do_not_save = 0;
      for (int i = 0; i < center_data->connections_found; i++){
        if(center_data->connections[i].id == ed->id){
          do_not_save = 1;
        }
      }

      if (shares_face){
        center_data->connections[center_data->connections_found].faces[0] = face_id;
        center_data->connections[center_data->connections_found].faces[1] = -1;
        center_data->connections[center_data->connections_found].faces[2] = -1;
      }
#if (P4EST_DIM)==3
      else if (shares_edge){    
        center_data->connections[center_data->connections_found].faces[0] = p8est_edge_faces[edge_id][0];
        center_data->connections[center_data->connections_found].faces[1] = p8est_edge_faces[edge_id][1];
        center_data->connections[center_data->connections_found].faces[2] = -1;
      }
#endif
      else {

#if (P4EST_DIM)==3
        center_data->connections[center_data->connections_found].faces[0] = p8est_corner_faces[corner_side->corner][0];
        center_data->connections[center_data->connections_found].faces[1] = p8est_corner_faces[corner_side->corner][1];
        center_data->connections[center_data->connections_found].faces[2] = p8est_corner_faces[corner_side->corner][2];
#else
        center_data->connections[center_data->connections_found].faces[0] = p4est_corner_faces[corner_side->corner][0];
        center_data->connections[center_data->connections_found].faces[1] = p4est_corner_faces[corner_side->corner][1];        
#endif
      }

      center_data->connections[center_data->connections_found].id = ed->id;
      center_data->connections[center_data->connections_found].mpirank = ed->mpi_rank;
      center_data->connections[center_data->connections_found].tree = ed->tree; 
      center_data->connections_found++;      
    }
  }  
}


void
d4est_solver_schwarz_get_center_data
(
 p4est_t* p4est,
 int element,
 d4est_solver_schwarz_subdomain_data_t* subdomain_data
)
{
#if (P4EST_DIM)==3
  int max_connections = 57;
#else
  int max_connections = 13;
#endif

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
  
}

void
d4est_solver_schwarz_restrict_field
(
 p4est_t* p4est,
 d4est_solver_schwarz_center_data_t* center_data, /* restricted field will have same sorting as center_data */
 d4est_factors_t* factors,
 double* field,
 double* restricted_field
){
  D4EST_ASSERT(p4est->mpisize == 1);
  for (int i = 0; i < subdomain_size; i++){

  }

  
}

void
d4est_solver_schwarz_restrict_transpose_field_single_core
(
 p4est_t* p4est,
 d4est_solver_schwarz_center_data_t* center_data, /* restricted field needs to be sorted in same manner as center data */
 double* field,
 double* restricted_field
)
{
  D4EST_ASSERT(p4est->mpisize == 1);
}

void
d4est_solver_schwarz_apply_lhs_single_core
(


){
  /* transpose restrict */
  /* apply op */
  /* restrict */
}

/* void */
/* d4est_solver_schwarz_add_correction */
/* ( */

/* ) */
/* { */
/*   for (int i = 0; i < schwarz_data->num_subdomain; i++){ */
/*     d4est_solver_schwarz_get_subdomain_data(p4est, i, &schwarz_data->subdomain_data[i]); */
/*   } */
/* } */


d4est_solver_schwarz_data_t*
d4est_solver_schwarz_init
(
 p4est_t* p4est,
 int restricted_size,
 d4est_solver_schwarz_weight_type_t weight_type
)
{
  d4est_solver_schwarz_data_t* schwarz_data = P4EST_ALLOC(d4est_solver_schwarz_data_t, 1);
  schwarz_data->num_subdomains = p4est->local_num_quadrants;
  schwarz_data->subdomain_data = P4EST_ALLOC(d4est_solver_subdomain_data_t, schwarz_data->num_subdomains);
  schwarz_data->weight_type = weight_type;

  /* fill subdomains */
  for (int i = 0; i < schwarz_data->num_subdomain; i++){
    d4est_solver_schwarz_get_subdomain_data(p4est, i, &schwarz_data->subdomain_data[i]);
  }
  
  return schwarz_data;
}

void
d4est_solver_schwarz_destroy
(
 d4est_solver_schwarz_data_t* schwarz_data
)
{
  
}

