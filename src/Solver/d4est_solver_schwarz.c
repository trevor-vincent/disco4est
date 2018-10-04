// TODO: implement sort over tree, quadid
// TODO: add face info to element datas from mesh_data
// TODO: test

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
#include <d4est_solver_schwarz.h>
#include <zlog.h>


static void
d4est_solver_schwarz_corner_callback
(
 p4est_iter_corner_info_t* info,
 void* user_data
)
{
  int core_side = -1;
  p4est_t* p4est = info->p4est;
  d4est_solver_schwarz_data_t* schwarz_data = user_data;

  for (int core = 0; core < info->sides.elem_count; core++){
    
    p4est_iter_corner_side_t* core_info
      = sc_array_index(&info->sides,
                       core);

    
    d4est_element_data_t* ed_core = core_info->quad->p.user_data;
    if (ed_core->mpirank != p4est->mpirank){ /* ghost elements cannot be core elements */
      continue; 
    }

    d4est_solver_schwarz_subdomain_data_t sub_data = schwarz_data->subdomain_data[ed_core->id];
    /* add core element to subdomain */
    sub_data.element_data[0].id = ed_core->id;
    sub_data.element_data[0].mpirank = ed_core->mpirank;
    sub_data.element_data[0].tree = ed_core->tree;
    sub_data.element_data[0].tree_quadid = ed_core->tree_quadid;
    sub_data.element_data[0].is_core = 1;
    sub_data.element_data[0].deg = ed_core->deg;
    sub_data.element_data[0].faces[0] = -1;
    sub_data.element_data[0].faces[1] = -1;
    sub_data.element_data[0].faces[2] = -1;
    sub_data.element_data[0].core_faces[0] = -1;
    sub_data.element_data[0].core_faces[1] = -1;
    sub_data.element_data[0].core_faces[2] = -1;    
    sub_data.num_elements++;
    
    for (int side = 0; side < info->sides.elem_count; side++){

      p4est_iter_corner_side_t* side_info
        = sc_array_index(&info->sides,
                         side);
      d4est_element_data_t* ed
        = side_info->quad->p.user_data;

      if (side != core){

        int found_side_face [] = {0,0,0};
        int found_side_edge [] = {0,0,0};
        int found_core_face [] = {0,0,0};
        int found_core_edge [] = {0,0,0};
        
        for (int f_side = 0; f_side < (P4EST_DIM); f_side++){
          for (int f_core = 0; f_core < (P4EST_DIM); f_core++){
            if ((core_info->faces[f_core] == side_info->faces[f_side])){
              found_side_face[f_side] = 1;
              found_core_face[f_core] = 1;
            }
          }
        }
#if (P4EST_DIM)==3         
          for (int e_side = 0; e_side < (P4EST_DIM); e_side++){
            for (int e_core = 0; e_core < (P4EST_DIM); e_core++){
              if ((core_info->edges[e_core] == side_info->edges[e_side])){
                found_side_edge[e_side] = 1;          
                found_core_edge[e_core] = 1;          
              }
            }
          }
#endif

            int side_face_id = -1;
            int side_edge_id = -1;
            int core_face_id = -1;
            int core_edge_id = -1;
        
            int shares_face = 0;
            int shares_edge = 0;
            int shares_corner = 1;
            int shares_face_or_edge = 0;
      
            for (int ef = 0; ef < (P4EST_DIM); ef++){
              if (found_side_face[ef] == 1){
#if (P4EST_DIM)==3
                side_face_id = p8est_corner_faces[side_info->corner][ef];
#else
                side_face_id = p4est_corner_faces[side_info->corner][ef];
#endif
                shares_face++;
                shares_face_or_edge++;
              }

              if (found_core_face[ef] == 1){
#if (P4EST_DIM)==3
                core_face_id = p8est_corner_faces[core_info->corner][ef];
#else
                core_face_id = p4est_corner_faces[core_info->corner][ef];
#endif
              }
          
#if (P4EST_DIM)==3       
              if (found_side_edge[ef] == 1){
                side_edge_id = p8est_corner_edges[side_info->corner][ef];
                shares_edge++;
                shares_face_or_edge++;
              }

              if (found_core_edge[ef] == 1){
                core_edge_id = p8est_corner_edges[core_info->corner][ef];
              }
#endif
            }

            d4est_element_data_t* ed_side = side_info->quad->p.user_data; 

            if (shares_face){
              sub_data.element_data[sub_data.num_elements].faces[0] = side_face_id;
              sub_data.element_data[sub_data.num_elements].faces[1] = -1;
              sub_data.element_data[sub_data.num_elements].faces[2] = -1;

              sub_data.element_data[sub_data.num_elements].core_faces[0] = core_face_id;
              sub_data.element_data[sub_data.num_elements].core_faces[1] = -1;
              sub_data.element_data[sub_data.num_elements].core_faces[2] = -1;
          
            }
#if (P4EST_DIM)==3
            else if (shares_edge){    
              sub_data.element_data[sub_data.num_elements].faces[0] = p8est_edge_faces[side_edge_id][0];
              sub_data.element_data[sub_data.num_elements].faces[1] = p8est_edge_faces[side_edge_id][1];
              sub_data.element_data[sub_data.num_elements].faces[2] = -1;


              sub_data.element_data[sub_data.num_elements].core_faces[0] = p8est_edge_faces[core_edge_id][0];
              sub_data.element_data[sub_data.num_elements].core_faces[1] = p8est_edge_faces[core_edge_id][1];
              sub_data.element_data[sub_data.num_elements].core_faces[2] = -1;
          
            }
#endif
            else {

#if (P4EST_DIM)==3
              sub_data.element_data[sub_data.num_elements].faces[0] = p8est_corner_faces[side_info->corner][0];
              sub_data.element_data[sub_data.num_elements].faces[1] = p8est_corner_faces[side_info->corner][1];
              sub_data.element_data[sub_data.num_elements].faces[2] = p8est_corner_faces[side_info->corner][2];

              sub_data.element_data[sub_data.num_elements].core_faces[0] = p8est_corner_faces[core_info->corner][0];
              sub_data.element_data[sub_data.num_elements].core_faces[1] = p8est_corner_faces[core_info->corner][1];
              sub_data.element_data[sub_data.num_elements].core_faces[2] = p8est_corner_faces[core_info->corner][2];

#else
              sub_data.element_data[sub_data.num_elements].faces[0] = p4est_corner_faces[side_info->corner][0];
              sub_data.element_data[sub_data.num_elements].faces[1] = p4est_corner_faces[side_info->corner][1];

              sub_data.element_data[sub_data.num_elements].core_faces[0] = p4est_corner_faces[core_info->corner][0];
              sub_data.element_data[sub_data.num_elements].core_faces[1] = p4est_corner_faces[core_info->corner][1];
#endif
            }
            /* add side element to subdomain */
            sub_data.element_data[sub_data.num_elements].id = ed_side->id;
            sub_data.element_data[sub_data.num_elements].mpirank = ed_side->mpirank;
            sub_data.element_data[sub_data.num_elements].tree = ed_side->tree;
            sub_data.element_data[sub_data.num_elements].tree_quadid = ed_side->tree_quadid;
            sub_data.element_data[sub_data.num_elements].is_core = 0;
            sub_data.element_data[sub_data.num_elements].deg = ed_side->deg;
            sub_data.num_elements++;      
          }
        }
      }
    }

void
d4est_solver_schwarz_init_subdomain_metadata
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_data_t* schwarz_data
)
{
    
#if (P4EST_DIM)==3
  p8est_iterate (p4est,
                 d4est_ghost->ghost,
                 &schwarz_data,
                 NULL,
                 NULL,
                 NULL,
                 d4est_solver_schwarz_corner_callback);

#else
  p4est_iterate (p4est,
                 d4est_ghost->ghost,
                 &schwarz_data,
                 NULL,
                 /* NULL, */
                 NULL,
                 d4est_solver_schwarz_corner_callback);
#endif

  /* sort */
  /* build strides and sizes */
}

d4est_solver_schwarz_data_t*
d4est_solver_schwarz_init
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 int num_nodes_overlap
)
{
  D4EST_ASSERT(num_nodes_overlap > 0);
  d4est_solver_schwarz_data_t* schwarz_data = P4EST_ALLOC(d4est_solver_schwarz_data_t, 1);
  schwarz_data->num_subdomains = p4est->local_num_quadrants;
  schwarz_data->subdomain_data = P4EST_ALLOC(d4est_solver_schwarz_subdomain_data_t, schwarz_data->num_subdomains);
  schwarz_data->num_nodes_overlap = num_nodes_overlap;
  int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13; 
  /* fill subdomains */
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    schwarz_data->subdomain_data->element_data = P4EST_ALLOC(d4est_solver_schwarz_element_data_t, max_connections);
  }

  d4est_solver_schwarz_init_subdomain_metadata(p4est, d4est_ghost, schwarz_data);
  
  return schwarz_data;
}

void
d4est_solver_schwarz_destroy
(
 d4est_solver_schwarz_data_t* schwarz_data
)
{
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    P4EST_FREE(schwarz_data->subdomain_data[i].element_data);
  }
  P4EST_FREE(schwarz_data->subdomain_data);
  P4EST_FREE(schwarz_data);
}



