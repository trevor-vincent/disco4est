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
#include <d4est_util.h>
#include <d4est_checkpoint.h>
#include <d4est_element_data.h>
#include <petscsnes.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_ghost_data_ext.h>
#include <zlog.h>
#include <ini.h>

static int schwarz_subdomain_field_size_of_ghost_fcn
(
 d4est_ghost_t* d4est_ghost,
 int gid,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_solver_schwarz_subdomain_metadata_t);
}

static int schwarz_subdomain_field_size_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_solver_schwarz_subdomain_metadata_t);
}

static int schwarz_subdomain_field_stride_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  return d4est_ghost->mirror_elements[mirr].id*sizeof(d4est_solver_schwarz_subdomain_metadata_t);
}



static int schwarz_element_field_size_of_ghost_fcn
(
 d4est_ghost_t* d4est_ghost,
 int gid,
 int tn,
 void* user_ctx
)
{
  /* return sizeof(d4est_solver_schwarz_element_metadata_t); */


  /* d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx; */
  /* d4est_solver_schwarz_subdomain_metadata_t* sub_data = */
    /* (d4est_solver_schwarz_subdomain_metadata_t*) */
    /* d4est_ghost_data_ext_get_field_on_element */
    /* ( */
     /* &d4est_ghost->ghost_elements[gid], */
     /* 0, */
     /* schwarz_metadata->subdomain_ghostdata */
    /* ); */
        
  /* return sizeof(d4est_solver_schwarz_element_metadata_t)*sub_data->num_elements; */

  int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13;
  return sizeof(d4est_solver_schwarz_element_metadata_t)*max_connections;
}

static int schwarz_element_field_size_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  /* d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx; */
  /* int id = d4est_ghost->mirror_elements[mirr].id; */
  /* return sizeof(d4est_solver_schwarz_element_metadata_t)*schwarz_metadata->subdomain_metadata[id].num_elements; */

  int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13;
  return sizeof(d4est_solver_schwarz_element_metadata_t)*max_connections;
  
}

static int schwarz_element_field_stride_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  /* d4est_solver_schwarz_metadata_t* schwarz_metadata = user_ctx; */
  int id = d4est_ghost->mirror_elements[mirr].id;
  /* return sizeof(d4est_solver_schwarz_element_metadata_t)*schwarz_metadata->subdomain_metadata[id].element_stride; */

  int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13;
  return sizeof(d4est_solver_schwarz_element_metadata_t)*max_connections*id;
  
}

static
int d4est_solver_schwarz_metadata_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_schwarz_metadata_t* pconfig = (d4est_solver_schwarz_metadata_t*)user;
  if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"num_nodes_overlap")) {
    D4EST_ASSERT(pconfig->num_nodes_overlap == -1);
    pconfig->num_nodes_overlap = atoi(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"subdomain_atol")) {
    D4EST_ASSERT(pconfig->subdomain_atol == -1);
    pconfig->subdomain_atol = atof(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"subdomain_rtol")) {
    D4EST_ASSERT(pconfig->subdomain_rtol == -1);
    pconfig->subdomain_rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"subdomain_iter")) {
    D4EST_ASSERT(pconfig->subdomain_iter == -1);
    pconfig->subdomain_iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"schwarz_iter")) {
    D4EST_ASSERT(pconfig->schwarz_iter == -1);
    pconfig->schwarz_iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,"d4est_solver_schwarz",name,"print_info")) {
    D4EST_ASSERT(pconfig->print_info == -1);
    pconfig->print_info = atoi(value);
  } 
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_schwarz_metadata_input
(
 p4est_t* p4est,
 const char* input_file,
 d4est_solver_schwarz_metadata_t* schwarz_data
)
{
  zlog_category_t* c_default = zlog_get_category("d4est_solver_schwarz_metadata");
  schwarz_data->num_nodes_overlap = -1;
  schwarz_data->subdomain_iter = -1;
  schwarz_data->subdomain_rtol = -1;
  schwarz_data->subdomain_atol = -1;
  schwarz_data->schwarz_iter = -1;
  schwarz_data->print_info = -1;
  
  if (
      ini_parse(input_file,
                d4est_solver_schwarz_metadata_input_handler,
                schwarz_data) < 0
  ) {
    D4EST_ABORT("Can't load input file");
  }
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->num_nodes_overlap, -1);
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->schwarz_iter, -1);
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->subdomain_iter, -1);
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->subdomain_rtol, -1);
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->subdomain_atol, -1);
  D4EST_CHECK_INPUT("d4est_solver_schwarz", schwarz_data->print_info, -1);
  
  if (schwarz_data->num_nodes_overlap <= 0){
    D4EST_ABORT("num_nodes_overlap <= 0");
  }
  
  zlog_info(c_default, "num_nodes_overlap = %d\n", schwarz_data->num_nodes_overlap);
}

static void
d4est_solver_schwarz_metadata_corner_callback
(
 p4est_iter_corner_info_t* info,
 void* user_data
)
{
  int core_side = -1;
  p4est_t* p4est = info->p4est;
  d4est_solver_schwarz_metadata_t* schwarz_data = user_data;
  d4est_element_data_t* ghost_data = schwarz_data->d4est_ghost->ghost_elements;

  
  for (int core = 0; core < info->sides.elem_count; core++){
    
    p4est_iter_corner_side_t* core_info
      = sc_array_index(&info->sides,
                       core);
    
    if (!core_info->is_ghost){ /* ghost elements cannot be core elements */

      d4est_element_data_t* ed_core = core_info->quad->p.user_data;
    
      d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[ed_core->id];

      /* add core element to subdomain */
      if(sub_data->num_elements == 0){
        sub_data->num_elements++;
        sub_data->core_deg = ed_core->deg;
        sub_data->core_id = ed_core->id;
        sub_data->element_metadata[0].id = ed_core->id;
        sub_data->element_metadata[0].mpirank = ed_core->mpirank;
        sub_data->element_metadata[0].tree = ed_core->tree;
        sub_data->element_metadata[0].tree_quadid = ed_core->tree_quadid;
        sub_data->element_metadata[0].is_core = 1;
        sub_data->element_metadata[0].deg = ed_core->deg;
        sub_data->element_metadata[0].faces[0] = -1;
        sub_data->element_metadata[0].faces[1] = -1;
        sub_data->element_metadata[0].faces[2] = -1;
        sub_data->element_metadata[0].core_faces[0] = -1;
        sub_data->element_metadata[0].core_faces[1] = -1;
        sub_data->element_metadata[0].core_faces[2] = -1;    
      }

      /* printf("ed_core->id = %d\n", ed_core->id); */
      /* printf("ed_core->deg = %d\n", ed_core->deg); */
    
      for (int side = 0; side < info->sides.elem_count; side++){

        if (side != core){

          p4est_iter_corner_side_t* side_info
            = sc_array_index(&info->sides,
                             side);

          d4est_element_data_t* ed_side = NULL;
          if (side_info->is_ghost){
            ed_side = &ghost_data[side_info->quadid];
          }
          else {
            ed_side = side_info->quad->p.user_data; 
          }
          int found = 0;
          for (int i = 0; i < sub_data->num_elements; i++){
            /* already added this element */          
            if (sub_data->element_metadata[i].tree == ed_side->tree &&
                sub_data->element_metadata[i].tree_quadid == ed_side->tree_quadid &&
                sub_data->element_metadata[i].mpirank == ed_side->mpirank){
              found = 1;
              break;
            }
          }
          /* if (side_info->is_ghost && p4est->mpirank == 0){ */
            /* printf("local mpirank = %d, ghost mpirank = %d, core_id = %d, ghost_quadid = %d, tree, tree_quadid = %d,%d, found = %d\n", p4est->mpirank, */
                   /* ed_side->mpirank, ed_core->id, side_info->quadid, ed_side->tree, ed_side->tree_quadid, found); */
          /* } */
          
          if (found == 1){
            continue;
          }
        
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

          /* d4est_element_data_t* ed_side = side_info->quad->p.user_data;  */
          sub_data->num_elements++;      
          if (shares_face){
            sub_data->element_metadata[sub_data->num_elements - 1].faces[0] = side_face_id;
            sub_data->element_metadata[sub_data->num_elements - 1].faces[1] = -1;
            sub_data->element_metadata[sub_data->num_elements - 1].faces[2] = -1;

            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[0] = core_face_id; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[1] = -1; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[2] = -1; */
          
          }
#if (P4EST_DIM)==3
          else if (shares_edge){    
            sub_data->element_metadata[sub_data->num_elements - 1].faces[0] = p8est_edge_faces[side_edge_id][0];
            sub_data->element_metadata[sub_data->num_elements - 1].faces[1] = p8est_edge_faces[side_edge_id][1];
            sub_data->element_metadata[sub_data->num_elements - 1].faces[2] = -1;


            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[0] = p8est_edge_faces[core_edge_id][0]; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[1] = p8est_edge_faces[core_edge_id][1]; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[2] = -1; */
          
          }
#endif
          else {

#if (P4EST_DIM)==3
            sub_data->element_metadata[sub_data->num_elements - 1].faces[0] = p8est_corner_faces[side_info->corner][0];
            sub_data->element_metadata[sub_data->num_elements - 1].faces[1] = p8est_corner_faces[side_info->corner][1];
            sub_data->element_metadata[sub_data->num_elements - 1].faces[2] = p8est_corner_faces[side_info->corner][2];

            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[0] = p8est_corner_faces[core_info->corner][0]; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[1] = p8est_corner_faces[core_info->corner][1]; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[2] = p8est_corner_faces[core_info->corner][2]; */

#else
            sub_data->element_metadata[sub_data->num_elements - 1].faces[0] = p4est_corner_faces[side_info->corner][0];
            sub_data->element_metadata[sub_data->num_elements - 1].faces[1] = p4est_corner_faces[side_info->corner][1];

            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[0] = p4est_corner_faces[core_info->corner][0]; */
            /* sub_data->element_data[sub_data->num_elements - 1].core_faces[1] = p4est_corner_faces[core_info->corner][1]; */
#endif
          }

          D4EST_ASSERT(sub_data->num_elements <= (((P4EST_DIM) ==3) ? 57 : 13));

          for (int face_it = 0; face_it < (P4EST_DIM); face_it++){
            int test_face = sub_data->element_metadata[sub_data->num_elements - 1].faces[face_it];
            sub_data->element_metadata[sub_data->num_elements - 1].core_faces[face_it] = -1;
            if(test_face != -1){
              sub_data->element_metadata[sub_data->num_elements - 1].core_faces[face_it] =
                d4est_reference_get_mirrored_face(test_face);
            }
          }
        
          /* add side element to subdomain */
          sub_data->element_metadata[sub_data->num_elements - 1].id = ed_side->id;
          sub_data->element_metadata[sub_data->num_elements - 1].mpirank = ed_side->mpirank;
          sub_data->element_metadata[sub_data->num_elements - 1].tree = ed_side->tree;
          sub_data->element_metadata[sub_data->num_elements - 1].tree_quadid = ed_side->tree_quadid;
          sub_data->element_metadata[sub_data->num_elements - 1].is_core = 0;
          sub_data->element_metadata[sub_data->num_elements - 1].deg = ed_side->deg;
        }
      }
    }
  }
}

static void
d4est_solver_schwarz_metadata_init_subdomain_metadata
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_metadata_t* schwarz_data
)
{
  schwarz_data->d4est_ghost = d4est_ghost;
  
#if (P4EST_DIM)==3
  p8est_iterate (p4est,
                 (void*)d4est_ghost->ghost,
                 schwarz_data,
                 NULL,
                 NULL,
                 NULL,
                 d4est_solver_schwarz_metadata_corner_callback);

#else
  p4est_iterate (p4est,
                 (void*)d4est_ghost->ghost,
                 schwarz_data,
                 NULL,
                 /* NULL, */
                 NULL,
                 d4est_solver_schwarz_metadata_corner_callback);
#endif

  /* sort */
  /* build strides and sizes */
}

static int
d4est_solver_schwarz_metadata_sort_callback(const void *p,
                                   const void *q)

{
  int ret;
  d4est_solver_schwarz_element_metadata_t x = *(const d4est_solver_schwarz_element_metadata_t *)p;
  d4est_solver_schwarz_element_metadata_t y = *(const d4est_solver_schwarz_element_metadata_t *)q;

  /* Avoid return x - y, which can cause undefined behaviour
       because of signed integer overflow. */

  if (x.tree < y.tree)
    ret = -1;
  else if (x.tree == y.tree){
    if (x.tree_quadid == y.tree_quadid)
      return 0;
    else if (x.tree_quadid < y.tree_quadid)
      return -1;
    else
      return 1;
  }
  else
    ret = 1;

  return ret;
}


static void
d4est_solver_schwarz_metadata_sort_elements
(
 d4est_solver_schwarz_metadata_t* schwarz_data
)
{  
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[i];
    qsort(sub_data.element_metadata, sub_data.num_elements,
          sizeof(d4est_solver_schwarz_element_metadata_t),d4est_solver_schwarz_metadata_sort_callback);
  }
}




static void
d4est_solver_schwarz_metadata_compute_strides_and_sizes
(
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_ghost_t* d4est_ghost //set to NULL if setting strides for local data, non-NULL for ghost data strides
)
{
  int sub_nodal_stride = 0;
  int sub_restricted_nodal_stride = 0;
  int sub_element_stride = 0;

  int num_subdomains = schwarz_data->num_subdomains;
  if (d4est_ghost != NULL){
    num_subdomains = schwarz_data->subdomain_ghostdata->num_ghosts;
  }
  
  for (int i = 0; i < num_subdomains; i++){

    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[i];
    if (d4est_ghost != NULL){
      sub_data =
        (d4est_solver_schwarz_subdomain_metadata_t*)
        d4est_ghost_data_ext_get_field_on_element
        (
         &d4est_ghost->ghost_elements[i],
         0,
         schwarz_data->subdomain_ghostdata
        );
    }

    sub_data->restricted_nodal_size = 0;
    sub_data->nodal_size = 0;
    sub_data->nodal_stride = sub_nodal_stride;
    sub_data->restricted_nodal_stride = sub_restricted_nodal_stride;
    sub_data->element_stride = sub_element_stride;
    
    int element_nodal_stride = 0;
    int element_restricted_nodal_stride = 0;

    /* printf("sub_data->num_elements = %d\n", sub_data->num_elements); */
    
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_metadata_t* ed = &sub_data->element_metadata[j];
      int res_volume_nodes = 1;
      for (int k = 0; k < (P4EST_DIM); k++){
        if (ed->faces[k] != -1){
          res_volume_nodes *= schwarz_data->num_nodes_overlap;
        }
        else {
          res_volume_nodes *= (ed->deg + 1);
        }
      }
      ed->nodal_size = d4est_lgl_get_nodes((P4EST_DIM), ed->deg);
      ed->nodal_stride = element_nodal_stride;
      ed->restricted_nodal_size = res_volume_nodes;
      ed->restricted_nodal_stride = element_restricted_nodal_stride;

      element_nodal_stride += ed->nodal_size;
      element_restricted_nodal_stride += ed->restricted_nodal_size;
      sub_data->nodal_size += ed->nodal_size;
      sub_data->restricted_nodal_size += ed->restricted_nodal_size;
    }

    sub_nodal_stride += sub_data->nodal_size;
    sub_restricted_nodal_stride += sub_data->restricted_nodal_size;
    /* sub_element_stride += sub_data->num_elements; */
    int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13;
    sub_element_stride += max_connections;
    
    schwarz_data->restricted_nodal_size += sub_data->restricted_nodal_size;
    schwarz_data->nodal_size += sub_data->nodal_size;
    schwarz_data->num_elements += sub_data->num_elements;
  }
}

d4est_solver_schwarz_metadata_t*
d4est_solver_schwarz_metadata_init
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 const char* input_file
 /* int num_nodes_overlap */
)
{
  /* TODO: Get rid of max_connections and actually just use proper number of elements */
  int max_connections = ((P4EST_DIM)==3 ) ? 57 : 13;

  d4est_solver_schwarz_metadata_t* schwarz_data = P4EST_ALLOC(d4est_solver_schwarz_metadata_t, 1);
  schwarz_data->restricted_nodal_size = 0;
  schwarz_data->nodal_size = 0;
  schwarz_data->num_elements = 0;
  schwarz_data->num_subdomains = p4est->local_num_quadrants;
  schwarz_data->subdomain_metadata = P4EST_ALLOC(d4est_solver_schwarz_subdomain_metadata_t, schwarz_data->num_subdomains);
  schwarz_data->element_metadata = P4EST_ALLOC(d4est_solver_schwarz_element_metadata_t, schwarz_data->num_subdomains*max_connections);

  d4est_solver_schwarz_metadata_input
    (
     p4est,
     input_file,
     schwarz_data
    );
  
  /* schwarz_data->num_nodes_overlap = num_nodes_overlap; */

  
  /* fill subdomains */
  int stride = 0;
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    schwarz_data->subdomain_metadata[i].num_elements = 0;
    schwarz_data->subdomain_metadata[i].element_metadata = &schwarz_data->element_metadata[stride];
    stride += max_connections;
  }

  d4est_solver_schwarz_metadata_init_subdomain_metadata(p4est, d4est_ghost, schwarz_data);
  d4est_solver_schwarz_metadata_sort_elements(schwarz_data);
  d4est_solver_schwarz_metadata_compute_strides_and_sizes(schwarz_data,NULL);

  
  schwarz_data->subdomain_ghostdata =
    d4est_ghost_data_ext_init
    (
     p4est,
     d4est_ghost,
     1,
     schwarz_subdomain_field_size_of_ghost_fcn,
     schwarz_subdomain_field_size_of_mirror_fcn,
     schwarz_subdomain_field_stride_of_mirror_fcn,
     NULL
    );
  
  d4est_ghost_data_ext_exchange
    (
     p4est,
     d4est_ghost,
     schwarz_data->subdomain_ghostdata,
     (char**)&schwarz_data->subdomain_metadata
    );

  schwarz_data->element_ghostdata =
    d4est_ghost_data_ext_init
    (
     p4est,
     d4est_ghost,
     1,
     schwarz_element_field_size_of_ghost_fcn,
     schwarz_element_field_size_of_mirror_fcn,
     schwarz_element_field_stride_of_mirror_fcn,
     NULL
    );
  
  d4est_ghost_data_ext_exchange
    (
     p4est,
     d4est_ghost,
     schwarz_data->element_ghostdata,
     (char**)&schwarz_data->element_metadata
    );

  if(d4est_ghost != NULL && p4est->mpisize > 1){
    int stride = 0;
    /* printf("Ghost subdomains on proc %d\n", p4est->mpirank); */
    /* printf("***************************\n", p4est->mpirank); */
    for (int gid = 0; gid < schwarz_data->subdomain_ghostdata->num_ghosts; gid++){
      d4est_solver_schwarz_subdomain_metadata_t* sub_data =
        (d4est_solver_schwarz_subdomain_metadata_t*)
        d4est_ghost_data_ext_get_field_on_element
        (
         &d4est_ghost->ghost_elements[gid],
         0,
         schwarz_data->subdomain_ghostdata
        );

       sub_data->element_metadata = d4est_ghost_data_ext_get_field_on_element
        (
         &d4est_ghost->ghost_elements[gid],
         0,
         schwarz_data->element_ghostdata
        );

       sub_data->element_stride = stride;
       stride += max_connections;      
    }
    printf("***************************\n", p4est->mpirank);
  }

  d4est_solver_schwarz_metadata_compute_strides_and_sizes(schwarz_data,d4est_ghost);
  
  
  return schwarz_data;
}

/* void */
/* d4est_solver_schwarz_metadata_get_face_neighbours */
/* ( */

/* ){ */
/*   for (int i = 0; i < schwarz_data->num_subdomains; i++){ */
/*     d4est_solver_schwarz_subdomain_metadata_t sub_data = schwarz_data->subdomain_metadata[subdomain]; */
/*     for (int j = 0; j < sub_data->num_elements; j++){ */
/*       d4est_solver_schwarz_element_metadata_t schwarz_ed = sub_data.element_metadata[j]; */

/*       d4est_element_data_t* ed_metadata = d4est_factors->element_data[schwarz_ed->id]; */

/*       for (int f = 0; f < (P4EST_FACES); f++){ */
/*         for (int s_p = 0; s_p < (P4EST_HALF): s_p++){ */
          
/*           p_tree_that_touch_face [f][s_p]; */
/*           p_tree_quadid_that_touch_face[f][s_p]; */
/*           p_face_that_touch_face [ */
          
/*         } */
/*       } */
/*     } */
/*   } */
/* } */

void
d4est_solver_schwarz_metadata_destroy
(
 d4est_solver_schwarz_metadata_t* schwarz_data
)
{
  P4EST_FREE(schwarz_data->subdomain_metadata);
  P4EST_FREE(schwarz_data->element_metadata);

  d4est_ghost_data_ext_destroy
    (
     schwarz_data->subdomain_ghostdata
    );


  d4est_ghost_data_ext_destroy
    (
     schwarz_data->element_ghostdata
    );

  
  P4EST_FREE(schwarz_data);
}

static void
d4est_solver_schwarz_metadata_print_subdomain
(
 d4est_solver_schwarz_subdomain_metadata_t* sub_data,
 int i,
 int ghost
)
{
    printf("\n");
    if (ghost){
      printf("Ghost Subdomain %d info\n", i);
    }
    else {
      printf("Subdomain %d info\n", i);
    }
    printf("********************************\n");
    printf("mpirank = %d\n", sub_data->mpirank);
    printf("num_elements = %d\n",sub_data->num_elements);
    printf("restricted_nodal_size = %d\n", sub_data->restricted_nodal_size);
    printf("nodal_size = %d\n", sub_data->nodal_size);
    printf("restricted_nodal_stride = %d\n", sub_data->restricted_nodal_stride);
    printf("nodal_stride = %d\n", sub_data->nodal_stride);
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_metadata_t* ed = &sub_data->element_metadata[j];
      printf("Element %d info\n", j);
      printf("*****************\n");
      printf("mpirank = %d\n", ed->mpirank);
      printf("local mesh id = %d\n", ed->id);
      printf("tree = %d\n", ed->tree);
      printf("tree_quadid = %d\n", ed->tree_quadid);
      printf("deg = %d\n", ed->deg);
      printf("nodal_size = %d\n", ed->nodal_size);
      printf("nodal_stride = %d\n", ed->nodal_stride);
      printf("restricted_nodal_stride = %d\n", ed->restricted_nodal_stride);
      printf("restricted_nodal_size = %d\n", ed->restricted_nodal_size);
      printf("is_core = %d\n", ed->is_core);
      printf("faces[0] = %d\n", ed->faces[0]);
      printf("faces[1] = %d\n", ed->faces[1]);
#if (P4EST_DIM)==3
      printf("faces[2] = %d\n", ed->faces[2]);
#endif
      printf("core_faces[0] = %d\n", ed->core_faces[0]);
      printf("core_faces[1] = %d\n", ed->core_faces[1]);
#if (P4EST_DIM)==3
      printf("core_faces[2] = %d\n", ed->core_faces[2]);
#endif
    }
}

void
d4est_solver_schwarz_metadata_print
(
 p4est_t* p4est,
 d4est_solver_schwarz_metadata_t* schwarz_data,
 d4est_ghost_t* d4est_ghost
)
{
  printf("Schwarz data on proc %d\n", p4est->mpirank);
  printf("***********************\n");
  printf("num subdomains = %d\n", schwarz_data->num_subdomains);
  printf("num_nodes_overlap = %d\n", schwarz_data->num_nodes_overlap);
  printf("nodal_size = %d\n", schwarz_data->nodal_size);
  printf("restricted_nodal_size = %d\n", schwarz_data->restricted_nodal_size);
  for (int i = 0; i < schwarz_data->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_data->subdomain_metadata[i];
    d4est_solver_schwarz_metadata_print_subdomain(sub_data,i,0);
  }
  
  if(d4est_ghost != NULL){
    printf("Ghost subdomains on proc %d\n", p4est->mpirank);
    printf("***************************\n", p4est->mpirank);
    for (int gid = 0; gid < schwarz_data->subdomain_ghostdata->num_ghosts; gid++){
      d4est_solver_schwarz_subdomain_metadata_t* sub_data =
        (d4est_solver_schwarz_subdomain_metadata_t*)
        d4est_ghost_data_ext_get_field_on_element
        (
         &d4est_ghost->ghost_elements[gid],
         0,
         schwarz_data->subdomain_ghostdata
        );

      printf("sub_data->num_elements = %d\n", sub_data->num_elements);
      printf("sub_data->num_elements = %d\n", sub_data->restricted_nodal_size);
      d4est_solver_schwarz_metadata_print_subdomain
        (
         sub_data,
         gid,
         1
        );
    }
    printf("***************************\n", p4est->mpirank);
  }
}
