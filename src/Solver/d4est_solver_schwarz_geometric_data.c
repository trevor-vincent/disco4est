#include <pXest.h>
#include <d4est_ghost.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_operators.h>
#include <d4est_mortars.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_helpers.h>
#include <ini.h>

#include <sc_reduce.h>
static int field_size_of_ghost_fcn
(
 d4est_ghost_t* d4est_ghost,
 int gid,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES);
}

static int field_size_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirror,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES);
}

static int field_stride_of_mirror_fcn
(
 d4est_ghost_t* d4est_ghost,
 int mirr,
 int tn,
 void* user_ctx
)
{
  return sizeof(d4est_mortar_side_data_t)*(P4EST_FACES)*mirr;
}


static
int d4est_solver_schwarz_geometric_data_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_schwarz_geometric_data_t* pconfig = (d4est_solver_schwarz_geometric_data_t*)user;
  const char* input_section = pconfig->input_section;
  if (d4est_util_match_couple(section,"mesh_parameters",name,"face_h_type")) {
    D4EST_ASSERT(pconfig->face_h_type == FACE_H_EQ_NOT_SET);
    if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_QUAD")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_QUAD;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_MIN_LOBATTO;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO")){
      pconfig->face_h_type = FACE_H_EQ_J_DIV_SJ_MEAN_LOBATTO;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_VOLUME_DIV_AREA")){
      pconfig->face_h_type = FACE_H_EQ_VOLUME_DIV_AREA;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_FACE_DIAM")){
      pconfig->face_h_type = FACE_H_EQ_FACE_DIAM;
    }
    else if(d4est_util_match(value, "FACE_H_EQ_TREE_H")){
      pconfig->face_h_type = FACE_H_EQ_TREE_H;
    }
    else {
      printf("face_h_type = %s\n", value);
      D4EST_ABORT("face_h_type is not set to a supported value\n");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}



void
d4est_solver_schwarz_geometric_data_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
)
{  
  schwarz_geometric_data->input_section = input_section;
  schwarz_geometric_data->face_h_type = FACE_H_EQ_NOT_SET;

  if(
     ini_parse(input_file,
               d4est_solver_schwarz_geometric_data_input_handler,
               schwarz_geometric_data) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("mesh_parameters", schwarz_geometric_data->face_h_type, FACE_H_EQ_NOT_SET);
  
}




static int
get_mesh_id_if_local_or_first_layer_ghost_id
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_element_metadata_t* ed
)
{
  /* local */
  if (ed->mpirank == p4est->mpirank){
    return ed->id;
  }
  else
    /* might be 1st layer ghost if this is true */
    if (ed->id < d4est_ghost->ghost->ghosts.elem_count){
      d4est_element_data_t* ged = &d4est_ghost->ghost_elements[ed->id];
      /* it is a first layer ghost if this is true */
      if (
          ged->mpirank == ed->mpirank &&
          ged->tree == ed->tree &&
          ged->tree_quadid == ed->tree_quadid
      ){
        return p4est->local_num_quadrants + ged->id;
      }
      /* Not in local mesh or local ghost mesh */
      else {
        return -1;
      }
    }
  /* Not in local mesh or local ghost mesh */
    else {
      return -1;
    }
}


void
d4est_solver_schwarz_geometric_data_boundary_callback
(
 p4est_t* p4est,
 d4est_element_data_t* e_m,
 int f_m,
 int mortar_side_id_m,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params
)
{
  d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data = params;
  int total_mortar_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg_quad);
  int total_mortar_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM)-1, e_m->deg);

  d4est_mortar_side_data_t* mortar_side_data = &schwarz_geometric_data->mortar_side_data[schwarz_geometric_data->mortar_side_stride];

  mortar_side_data->faces_m = 1;
  mortar_side_data->faces_p = 0;
  mortar_side_data->f_p = -1;
  mortar_side_data->f_m = f_m;
  mortar_side_data->tree_p = -1;
  mortar_side_data->tree_m = e_m->tree;
  mortar_side_data->orientation = -1;
  mortar_side_data->e_m[0] = *e_m;
    
  mortar_side_data->mortar_side_id = schwarz_geometric_data->mortar_side_stride;
  mortar_side_data->mortar_quad_stride = schwarz_geometric_data->mortar_quad_stride;
  mortar_side_data->boundary_quad_stride = schwarz_geometric_data->boundary_quad_stride;
  mortar_side_data->total_mortar_nodes_quad = total_mortar_nodes_quad;
  mortar_side_data->total_mortar_nodes_lobatto = total_mortar_nodes_lobatto;
  mortar_side_data->is_ghost = 0;
  
  schwarz_geometric_data->mortar_which_touches_face[e_m->id*(P4EST_FACES) + f_m] = mortar_side_data->mortar_side_id;
  schwarz_geometric_data->mortar_quad_stride += total_mortar_nodes_quad;
  schwarz_geometric_data->boundary_quad_stride += total_mortar_nodes_quad;
  schwarz_geometric_data->mortar_side_stride += 1;
}

void
d4est_solver_schwarz_geometric_data_interface_callback
(
 p4est_t* p4est,
 d4est_element_data_t** e_m,
 int faces_m,
 int f_m,
 int mortar_side_id_m,
 d4est_element_data_t** e_p,
 int faces_p,
 int f_p,
 int mortar_side_id_p,
 int* e_m_is_ghost,
 int orientation,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* params
){

  int stride;
  int deg_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_lobatto [(P4EST_HALF)];
  int face_nodes_p_quad [(P4EST_HALF)];
  int deg_m_lobatto [(P4EST_HALF)];
  int deg_m_quad [(P4EST_HALF)];
  int deg_p_quad [(P4EST_HALF)];
  int face_nodes_m_lobatto [(P4EST_HALF)];
  int face_nodes_m_quad [(P4EST_HALF)];
  int mortar_nodes_quad [(P4EST_HALF)];
  int mortar_nodes_lobatto [(P4EST_HALF)];
  int deg_mortar_quad [(P4EST_HALF)];
  int deg_mortar_lobatto [(P4EST_HALF)];
  int deg_p_lobatto_porder [(P4EST_HALF)];
  int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

  d4est_element_data_t* e_p_oriented [(P4EST_HALF)];
  d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);

  /* calculate degs and nodes of each face of (-) side */
  int total_side_nodes_m_lobatto = 0;
  int total_side_nodes_m_quad = 0;
  for (int i = 0; i < faces_m; i++){
    deg_m_lobatto[i] = e_m[i]->deg;
    deg_m_quad[i] = e_m[i]->deg_quad;
    
    face_nodes_m_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_m[i]->deg);
    face_nodes_m_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_m_quad[i]);
    
    total_side_nodes_m_lobatto += face_nodes_m_lobatto[i];
    total_side_nodes_m_quad += face_nodes_m_quad[i];
  }

  /* calculate degs and nodes of each face of (+) side  */
  int total_side_nodes_p_lobatto = 0;
  int total_side_nodes_p_quad = 0;
  for (int i = 0; i < faces_p; i++){
    deg_p_lobatto[i] = e_p_oriented[i]->deg;
    deg_p_quad[i] = e_p_oriented[i]->deg_quad;
    deg_p_lobatto_porder[i] = e_p[i]->deg;

    face_nodes_p_lobatto[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, e_p_oriented[i]->deg );
    face_nodes_p_quad[i] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_p_quad[i]);
    
    total_side_nodes_p_lobatto += face_nodes_p_lobatto[i];
    total_side_nodes_p_quad += face_nodes_p_quad[i];
  }

  /* calculate degs and nodes of the mortar faces */
  int total_mortar_nodes_quad = 0;
  int total_mortar_nodes_lobatto = 0;
  for (int i = 0; i < faces_m; i++)
    for (int j = 0; j < faces_p; j++){
      /* find max degree for each face pair of the two sides*/
      deg_mortar_quad[i+j] = d4est_util_max_int(deg_m_quad[i],
                                                deg_p_quad[j]);
      deg_mortar_lobatto[i+j] = d4est_util_max_int( e_m[i]->deg,
                                                    e_p_oriented[j]->deg );
      mortar_nodes_quad[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_quad[i+j] );
      mortar_nodes_lobatto[i+j] = d4est_lgl_get_nodes( (P4EST_DIM) - 1, deg_mortar_lobatto[i+j] );
      total_mortar_nodes_quad += mortar_nodes_quad[i+j];
      total_mortar_nodes_lobatto += mortar_nodes_lobatto[i+j];

    }

  int deg_mortar_quad_porder [(P4EST_HALF)];
  int deg_mortar_lobatto_porder [(P4EST_HALF)];
  int mortar_nodes_quad_porder [(P4EST_HALF)];
  
  for(int i = 0; i < faces_mortar; i++){
    int inew = i;
    if (faces_mortar == (P4EST_HALF)){
      inew = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, i);
    }
    deg_mortar_quad_porder[inew] = deg_mortar_quad[i];
    deg_mortar_lobatto_porder[inew] = deg_mortar_lobatto[i];
    mortar_nodes_quad_porder[inew] = mortar_nodes_quad[i];
  }
  d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data = params;
  d4est_mortar_side_data_t* mortar_side_data = &schwarz_geometric_data->mortar_side_data[schwarz_geometric_data->mortar_side_stride];

  mortar_side_data->faces_m = faces_m;
  mortar_side_data->faces_p = faces_p;
  mortar_side_data->f_p = f_p;
  mortar_side_data->f_m = f_m;
  mortar_side_data->tree_p = e_p[0]->tree;
  mortar_side_data->tree_m = e_m[0]->tree;
  mortar_side_data->orientation = orientation;
  mortar_side_data->is_ghost = 0;
  mortar_side_data->mortar_side_id = schwarz_geometric_data->mortar_side_stride;
  mortar_side_data->mortar_quad_stride = schwarz_geometric_data->mortar_quad_stride;
  mortar_side_data->boundary_quad_stride = -1;
  mortar_side_data->total_mortar_nodes_quad = total_mortar_nodes_quad;
  mortar_side_data->total_mortar_nodes_lobatto = total_mortar_nodes_lobatto;
  
  for(int f = 0; f < faces_m; f++){
    mortar_side_data->e_m[f] = *e_m[f];
    schwarz_geometric_data->mortar_which_touches_face[e_m[f]->id*(P4EST_FACES) + f_m] = mortar_side_data->mortar_side_id;
  }
  for(int f = 0; f < faces_p; f++){
    mortar_side_data->e_p[f] = *e_p[f];

  }
  
  schwarz_geometric_data->mortar_quad_stride += total_mortar_nodes_quad;
  schwarz_geometric_data->mortar_side_stride += 1;

}

void
d4est_solver_schwarz_geometric_data_destroy
(
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
)
{
  P4EST_FREE(schwarz_geometric_data->mortar_which_touches_face);
  d4est_ghost_data_ext_destroy(schwarz_geometric_data->mortar_side_ghost_data);
  
  P4EST_FREE(schwarz_geometric_data->num_of_mortars_per_subdomain);
  P4EST_FREE(schwarz_geometric_data->mortar_strides_per_subdomain);
  P4EST_FREE(schwarz_geometric_data->subdomain_mortars);
  P4EST_FREE(schwarz_geometric_data->zero_and_skip_m);
  P4EST_FREE(schwarz_geometric_data->nodal_stride_m);
  P4EST_FREE(schwarz_geometric_data->nodal_stride_p);
  P4EST_FREE(schwarz_geometric_data->zero_and_skip_p);
  P4EST_FREE(schwarz_geometric_data->J_quad_ghost);
  P4EST_FREE(schwarz_geometric_data->rst_xyz_quad_ghost);
  P4EST_FREE(schwarz_geometric_data->volume_quad_strides_per_ghost);
  /* P4EST_FREE(schwarz_geometric_data->skip_p_sum); */
  
  P4EST_FREE(schwarz_geometric_data->drst_dxyz_m_mortar_quad);
  P4EST_FREE(schwarz_geometric_data->drst_dxyz_p_mortar_quad_porder);
  P4EST_FREE(schwarz_geometric_data->sj_m_mortar_quad);
  P4EST_FREE(schwarz_geometric_data->n_m_mortar_quad);
  P4EST_FREE(schwarz_geometric_data->hm_mortar_quad);
  P4EST_FREE(schwarz_geometric_data->hp_mortar_quad);
  P4EST_FREE(schwarz_geometric_data->xyz_on_f_m_quad);
  P4EST_FREE(schwarz_geometric_data->xyz_on_f_m_lobatto);
  
  P4EST_FREE(schwarz_geometric_data->mortar_side_data);
  P4EST_FREE(schwarz_geometric_data);
}

void
d4est_solver_schwarz_geometric_data_reduce_to_minimal_set
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
)
{
  /* printf("schwarz_metadata->num_elements = %d\n", */
         /* schwarz_metadata->num_elements); */
  int* sides_done = P4EST_ALLOC_ZERO(int, schwarz_metadata->num_elements*(P4EST_FACES));

  int stride = 0;
  for (int i = 0; i < schwarz_metadata->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[i];
    int num_of_mortars = 0;
  
    for (int j = 0; j < sub_data->num_elements; j++){
      d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[j];

      for (int f = 0; f < (P4EST_FACES); f++){
        int mortar_face_stride = (stride + num_of_mortars)*(P4EST_HALF);
        /* already computed for this subdomain */
        /* printf("j = %d\n", j); */
        /* printf("sub_data->element_stride = %d\n", sub_data->element_stride); */
        /* printf("(sub_data->element_stride + j)*(P4EST_FACES) + f = %d\n", */
               /* (sub_data->element_stride + j)*(P4EST_FACES) + f); */
        if (sides_done[(sub_data->element_stride + j)*(P4EST_FACES) + f] == 1){
          continue;
        }

        d4est_mortar_side_data_t* mortar_side_data = NULL;
        if (schwarz_ed->mpirank == p4est->mpirank){          
          int mortar_side_id = schwarz_geometric_data->mortar_which_touches_face[schwarz_ed->id*(P4EST_FACES) + f];
          mortar_side_data = &schwarz_geometric_data->mortar_side_data[mortar_side_id];
        }
        else {
          d4est_mortar_side_data_t* ghost_mortar_side_data =
            d4est_ghost_data_ext_get_field_on_element
            (
             &d4est_ghost->ghost_elements[schwarz_ed->id],
             0,
             schwarz_geometric_data->mortar_side_ghost_data
            );
          mortar_side_data = &ghost_mortar_side_data[f];
        }
          
          
        for (int fm = 0; fm < mortar_side_data->faces_m; fm++){

          int is_in_subdomain =
            d4est_solver_schwarz_is_element_in_subdomain(sub_data,
                                                         mortar_side_data->e_m[fm].mpirank,
                                                         mortar_side_data->e_m[fm].tree,
                                                         mortar_side_data->e_m[fm].tree_quadid
                                                        );
          if (is_in_subdomain >= 0){
            schwarz_geometric_data->nodal_stride_m[mortar_face_stride + fm]
              = sub_data->element_metadata[is_in_subdomain].nodal_stride;
            sides_done[(sub_data->element_stride + is_in_subdomain)*(P4EST_FACES) + mortar_side_data->f_m] == 1;

            D4EST_ASSERT(schwarz_geometric_data->nodal_stride_m[mortar_face_stride + fm] < sub_data->nodal_size);
          }
          else {
            schwarz_geometric_data->zero_and_skip_m[mortar_face_stride + fm] = 1;
          }
        }

        int skip_p_sum = 0;
        for (int fp = 0; fp < mortar_side_data->faces_p; fp++){

          int is_in_subdomain =
            d4est_solver_schwarz_is_element_in_subdomain(sub_data,
                                                         mortar_side_data->e_p[fp].mpirank,
                                                         mortar_side_data->e_p[fp].tree,
                                                         mortar_side_data->e_p[fp].tree_quadid
                                                        );

          skip_p_sum += (is_in_subdomain < 0);
          if (is_in_subdomain >= 0){
            schwarz_geometric_data->nodal_stride_p[mortar_face_stride + fp]
              = sub_data->element_metadata[is_in_subdomain].nodal_stride;
            D4EST_ASSERT(schwarz_geometric_data->nodal_stride_p[mortar_face_stride + fp] < sub_data->nodal_size);
              
            sides_done[(sub_data->element_stride + is_in_subdomain)*(P4EST_FACES) + mortar_side_data->f_p] == 1;
          }
          else {
            schwarz_geometric_data->zero_and_skip_p[mortar_face_stride + fp] = 1;
          }           
        }
        /* potential speed up */
        /* if (skip_p_sum == mortar_side_data->faces_p */
            /* && */
            /* overlap_m_sum == mortar_side_data->faces_m){ */
          /* continue; */
        /* } */
        
        /* save mortar then move stride one forward */
        schwarz_geometric_data->subdomain_mortars[stride + num_of_mortars] = mortar_side_data;
        num_of_mortars++;
      }
    }
    
    schwarz_geometric_data->num_of_mortars_per_subdomain[i] = num_of_mortars;
    schwarz_geometric_data->mortar_strides_per_subdomain[i] = stride;
    stride += num_of_mortars;
  }

  P4EST_FREE(sides_done);

  
}

d4est_solver_schwarz_geometric_data_t*
d4est_solver_schwarz_geometric_data_init
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 const char* input_file,
 const char* input_section
)
{
  d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data =
    P4EST_ALLOC(d4est_solver_schwarz_geometric_data_t, 1);

d4est_solver_schwarz_geometric_data_input
(
 p4est,
 input_file,
 input_section,
 schwarz_geometric_data
);
  
  /* printf(" d4est_factors->local_sizes.local_mortar_sides = %d\n",  d4est_factors->local_sizes.local_mortar_sides); */
  schwarz_geometric_data->mortar_side_data = P4EST_ALLOC(d4est_mortar_side_data_t, d4est_factors->local_sizes.local_mortar_sides);
  schwarz_geometric_data->mortar_which_touches_face = P4EST_ALLOC(int, p4est->local_num_quadrants*(P4EST_FACES));


  schwarz_geometric_data->mortar_side_stride = 0;
  schwarz_geometric_data->mortar_quad_stride = 0;
  schwarz_geometric_data->boundary_quad_stride = 0;
  
  d4est_mortars_fcn_ptrs_t flux_fcns;
  flux_fcns.flux_interface_fcn = d4est_solver_schwarz_geometric_data_interface_callback;
  flux_fcns.flux_boundary_fcn = d4est_solver_schwarz_geometric_data_boundary_callback;
  flux_fcns.user_ctx = (void*)schwarz_geometric_data;

  d4est_mortars_compute_flux_on_local_elements
    (
     p4est,
     d4est_ghost,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     &flux_fcns
    );

  d4est_mortar_side_data_t* mirror_data = P4EST_ALLOC(d4est_mortar_side_data_t,
                                                      (P4EST_FACES)*d4est_ghost->ghost->mirrors.elem_count);

  for (int i = 0; i < d4est_ghost->ghost->mirrors.elem_count; i++){
    /* get mirror data */
    d4est_element_data_t* mirror_ed_ref = &d4est_ghost->mirror_elements[i];

    d4est_element_data_t* mirror_ed
      = d4est_element_data_get_ptr
      (
       p4est,
       mirror_ed_ref->tree,
       mirror_ed_ref->tree_quadid
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){
      int mortar_side_id = schwarz_geometric_data->mortar_which_touches_face[mirror_ed->id*(P4EST_FACES) + f];
      mirror_data[i*(P4EST_FACES) + f] = schwarz_geometric_data->mortar_side_data[mortar_side_id];
    }
  }
  
  schwarz_geometric_data->mortar_side_ghost_data =
    d4est_ghost_data_ext_init
    (
     p4est,
     d4est_ghost,
     1,
     field_size_of_ghost_fcn,
     field_size_of_mirror_fcn,
     field_stride_of_mirror_fcn,
     NULL
    );

  d4est_ghost_data_ext_exchange
    (
     p4est,
     d4est_ghost,
     schwarz_geometric_data->mortar_side_ghost_data,
     (char**)&mirror_data
    );

  P4EST_FREE(mirror_data);
  
  /* loop through data and set strides */
  
  /* int mortar_quad_ghost_size = 0; */
  /* int boundary_quad_ghost_size = 0; */
  
  int boundary_quad_stride = 0;
  int mortar_quad_stride = 0;

  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    
    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       schwarz_geometric_data->mortar_side_ghost_data
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){
      ghost_mortar_side_data[f].is_ghost = 1;
      int total_mortar_nodes_quad = 0;
 
      if (ghost_mortar_side_data[f].faces_p == 0){
        ghost_mortar_side_data[f].boundary_quad_stride = boundary_quad_stride;
        ghost_mortar_side_data[f].mortar_quad_stride = mortar_quad_stride;
        boundary_quad_stride += ghost_mortar_side_data[f].total_mortar_nodes_quad;
      }
      else {
        ghost_mortar_side_data[f].mortar_quad_stride = mortar_quad_stride;
        ghost_mortar_side_data[f].boundary_quad_stride = -1;
      }      
      mortar_quad_stride += ghost_mortar_side_data[f].total_mortar_nodes_quad;
    }
  }

  int mortar_quad_ghost_size = mortar_quad_stride;
  int boundary_quad_ghost_size = boundary_quad_stride;

  /* printf("mortar_quad_ghost_size = %d\n", mortar_quad_ghost_size); */
  
  schwarz_geometric_data->sj_m_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_geometric_data->n_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_geometric_data->hp_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_geometric_data->hm_mortar_quad = P4EST_ALLOC_ZERO(double, mortar_quad_ghost_size);
  schwarz_geometric_data->drst_dxyz_m_mortar_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_geometric_data->drst_dxyz_p_mortar_quad_porder = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*(P4EST_DIM)*mortar_quad_ghost_size);
  schwarz_geometric_data->xyz_on_f_m_quad = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);
  schwarz_geometric_data->xyz_on_f_m_lobatto = P4EST_ALLOC_ZERO(double, (P4EST_DIM)*boundary_quad_ghost_size);

  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){

    d4est_mortar_side_data_t* ghost_mortar_side_data =
      d4est_ghost_data_ext_get_field_on_element
      (
       &d4est_ghost->ghost_elements[gid],
       0,
       schwarz_geometric_data->mortar_side_ghost_data
      );
    
    for (int f = 0; f < (P4EST_FACES); f++){

      int scalar_stride = ghost_mortar_side_data[f].mortar_quad_stride;
      int vector_stride = ghost_mortar_side_data[f].mortar_quad_stride*(P4EST_DIM);
      int matrix_stride = ghost_mortar_side_data[f].mortar_quad_stride*(P4EST_DIM)*(P4EST_DIM);
      int boundary_vector_stride = ghost_mortar_side_data[f].boundary_quad_stride*(P4EST_DIM);

      /* printf("boundary_vector_stride = %d\n", boundary_vector_stride); */
      double* xyz_on_f_m_quad [(P4EST_DIM)];
      double* xyz_on_f_m_lobatto [(P4EST_DIM)];
      double* n_m_mortar_quad [(P4EST_DIM)];
      double* drst_dxyz_m_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
      double* drst_dxyz_p_mortar_quad_porder [(P4EST_DIM)][(P4EST_DIM)];

      /* printf("scalar stride = %d\n", scalar_stride); */
      double* sj_m_mortar_quad = &schwarz_geometric_data->sj_m_mortar_quad[scalar_stride];
      double* hm_mortar_quad = &schwarz_geometric_data->hm_mortar_quad[scalar_stride];
      double* hp_mortar_quad = &schwarz_geometric_data->hp_mortar_quad[scalar_stride];
  
      for (int i = 0; i < (P4EST_DIM); i++){
        n_m_mortar_quad[i] =
          &schwarz_geometric_data->n_m_mortar_quad
          [vector_stride + i*ghost_mortar_side_data[f].total_mortar_nodes_quad];

        for (int j = 0; j < (P4EST_DIM); j++){
          drst_dxyz_m_mortar_quad[i][j] =
            &schwarz_geometric_data->drst_dxyz_m_mortar_quad
            [matrix_stride + (i + j*(P4EST_DIM))*ghost_mortar_side_data[f].total_mortar_nodes_quad];
          
          drst_dxyz_p_mortar_quad_porder[i][j] =
            &schwarz_geometric_data->drst_dxyz_p_mortar_quad_porder
            [matrix_stride + (i + j*(P4EST_DIM))*ghost_mortar_side_data[f].total_mortar_nodes_quad];
        }
      }

      if (ghost_mortar_side_data[f].faces_p == 0){

        d4est_element_data_t* e_m = &d4est_ghost->ghost_elements[gid];
        /* D4EST_ASSERT(ghost_mortar_side_data[f].faces_m == 1); */
        /* D4EST_ASSERT(ghost_mortar_side_data[f].tree_m == e_m->tree); */
        /* D4EST_ASSERT(ghost_mortar_side_data[f].tree_quadid_m[0] == e_m->tree_quadid); */

        for (int i = 0; i < (P4EST_DIM); i++){
          xyz_on_f_m_quad [i] = &schwarz_geometric_data->xyz_on_f_m_quad[boundary_vector_stride + i*ghost_mortar_side_data[f].total_mortar_nodes_quad];
          xyz_on_f_m_lobatto [i] = &schwarz_geometric_data->xyz_on_f_m_lobatto[boundary_vector_stride + i*ghost_mortar_side_data[f].total_mortar_nodes_quad];
        }

        double* xyz_lobatto [(P4EST_DIM)];
        D4EST_ALLOC_DIM_VEC(xyz_lobatto, d4est_lgl_get_nodes((P4EST_DIM),ghost_mortar_side_data[f].e_m[0].deg));


        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq =  e_m->dq;
        mesh_object.tree = e_m->tree;
        
        mesh_object.q[0] = e_m->q[0];
        mesh_object.q[1] = e_m->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = e_m->q[2];
#endif

        

        d4est_rst_t rst_points_lobatto;
        rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 0);
        rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 1);
        rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
        rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst(d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, e_m->deg, 2);
#endif
        
        
        d4est_geometry_compute_xyz
          (
           d4est_ops,
           d4est_geom,
           rst_points_lobatto,
           ghost_mortar_side_data[f].tree_m,
           ghost_mortar_side_data[f].e_m[0].deg,
           ghost_mortar_side_data[f].e_m[0].q,
           ghost_mortar_side_data[f].e_m[0].dq,
           xyz_lobatto
          );
  
        for (int d = 0; d < (P4EST_DIM); d++){

          d4est_operators_apply_slicer(d4est_ops,
                                       xyz_lobatto[d],
                                       (P4EST_DIM),
                                       ghost_mortar_side_data[f].f_m,
                                       ghost_mortar_side_data[f].e_m[0].deg,
                                       xyz_on_f_m_lobatto[d]);

        }

        D4EST_FREE_DIM_VEC(xyz_lobatto);
 
        double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           e_m->tree,
           e_m->q,
           e_m->dq,
           ghost_mortar_side_data[f].mortar_side_id,
           1,
           1,
           &e_m->deg,
           &e_m->deg_quad,
           ghost_mortar_side_data[f].f_m,
           xyz_on_f_m_quad,
           drst_dxyz_m_mortar_quad,
           sj_m_mortar_quad,
           n_m_mortar_quad,
           NULL,
           j_div_sj_m_mortar_quad,
           COMPUTE_NORMAL_USING_JACOBIAN
          );


        /* DEBUG_PRINT_ARR_DBL(j_div_sj_m_mortar_quad, ghost_mortar_side_data[f].total_mortar_nodes_quad); */
        d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
        
        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &ged,
           ghost_mortar_side_data[f].f_m,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           schwarz_geometric_data->face_h_type,
           j_div_sj_m_mortar_quad,
           hm_mortar_quad,
           1,
           1,
           &ghost_mortar_side_data[f].total_mortar_nodes_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );

        P4EST_FREE(j_div_sj_m_mortar_quad);
        
      }
      else {

        int faces_m = ghost_mortar_side_data[f].faces_m;
        int faces_p = ghost_mortar_side_data[f].faces_p;
        int orientation = ghost_mortar_side_data[f].orientation;
        int f_m = ghost_mortar_side_data[f].f_m;
        int f_p = ghost_mortar_side_data[f].f_p;
        int faces_mortar = (faces_m > faces_p) ? faces_m : faces_p;

        int deg_mortar_quad_porder [P4EST_HALF];
        int deg_mortar_lobatto_porder [P4EST_HALF];
        int deg_mortar_quad_morder [P4EST_HALF];
        int deg_mortar_lobatto_morder [P4EST_HALF];
        
        
        for (int i = 0; i < faces_m; i++){
          for (int j = 0; j < faces_p; j++){
            /* find max degree for each face pair of the two sides*/
            
            int jnew = j;
            if (faces_p == (P4EST_HALF)){
              jnew = d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                         f_m,
                                                         f_p,
                                                         orientation,
                                                         j);
            }
            deg_mortar_quad_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data[f].e_m[i].deg_quad,
                                                               ghost_mortar_side_data[f].e_p[jnew].deg_quad);
            deg_mortar_lobatto_morder[i + j] = d4est_util_max_int(ghost_mortar_side_data[f].e_m[i].deg,
                                                                  ghost_mortar_side_data[f].e_p[jnew].deg);
            
          }
        }

        for (int i = 0; i < faces_mortar; i++){
          int inew = i;
          if (faces_mortar == P4EST_HALF){
            d4est_reference_reorient_face_order((P4EST_DIM)-1,
                                                f_m,
                                                f_p,
                                                orientation, i);
          }
          deg_mortar_quad_porder[inew] = deg_mortar_quad_morder[i];
          deg_mortar_lobatto_porder[inew] = deg_mortar_lobatto_morder[i];
        }
        double* j_div_sj_m_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        double* j_div_sj_p_mortar_quad = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        double* j_div_sj_p_mortar_quad_porder = P4EST_ALLOC(double, ghost_mortar_side_data[f].total_mortar_nodes_quad);
        
        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           ghost_mortar_side_data[f].tree_m,
           ghost_mortar_side_data[f].e_m[0].q,
           ghost_mortar_side_data[f].e_m[0].dq,
           ghost_mortar_side_data[f].mortar_side_id,
           ghost_mortar_side_data[f].faces_m,
           faces_mortar,
           &deg_mortar_lobatto_morder[0],
           &deg_mortar_quad_morder[0],
           f_m,
           NULL,
           drst_dxyz_m_mortar_quad,
           sj_m_mortar_quad,
           n_m_mortar_quad,
           NULL,
           j_div_sj_m_mortar_quad,
           COMPUTE_NORMAL_USING_JACOBIAN
          );

        d4est_mortars_compute_geometric_data_on_mortar
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           QUAD_INTEGRAND_UNKNOWN,
           ghost_mortar_side_data[f].tree_p,
           ghost_mortar_side_data[f].e_p[0].q,
           ghost_mortar_side_data[f].e_p[0].dq,
           ghost_mortar_side_data[f].mortar_side_id,
           ghost_mortar_side_data[f].faces_p,
           faces_mortar,
           &deg_mortar_lobatto_porder[0],
           &deg_mortar_quad_porder[0],
           ghost_mortar_side_data[f].f_p,
           NULL,
           drst_dxyz_p_mortar_quad_porder,
           NULL,
           NULL,
           NULL,
           j_div_sj_p_mortar_quad_porder,
           COMPUTE_NORMAL_USING_JACOBIAN
          );
  
        int face_mortar_stride = 0;
        for (int face = 0; face < faces_mortar; face++){
          int face_p = face;
          if (faces_mortar == (P4EST_HALF))
            face_p = d4est_reference_reorient_face_order((P4EST_DIM)-1, f_m, f_p, orientation, face);

          int oriented_face_mortar_stride = 0;
          for (int b = 0; b < face_p; b++){
            oriented_face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_porder[b]);
          }

          d4est_operators_reorient_face_data
            (
             d4est_ops,
             &j_div_sj_p_mortar_quad_porder[oriented_face_mortar_stride],
             (P4EST_DIM)-1,
             deg_mortar_quad_morder[face],
             ghost_mortar_side_data[f].orientation,
             ghost_mortar_side_data[f].f_m,
             ghost_mortar_side_data[f].f_p,
             &j_div_sj_p_mortar_quad[face_mortar_stride]
            );
    
          face_mortar_stride += d4est_lgl_get_nodes((P4EST_DIM)-1, deg_mortar_quad_morder[face]);
        }

        int nodes_mortar_quad [(P4EST_HALF)];

        for (int i = 0; i < faces_mortar; i++){
          nodes_mortar_quad[i]
            = d4est_lgl_get_nodes((P4EST_DIM)-1,
                                  deg_mortar_quad_morder[i]);
        }

        d4est_element_data_t* e_m [P4EST_HALF];
        d4est_element_data_t* e_p [P4EST_HALF];

        for (int i = 0; i < faces_m; i++){
          e_m[f] = &ghost_mortar_side_data[f].e_m[i];
          D4EST_ASSERT(e_m[f]->deg > 0);
        }

        for (int i = 0; i < faces_p; i++){
          e_p[f] = &ghost_mortar_side_data[f].e_p[i];
          D4EST_ASSERT(e_p[f]->deg > 0);
        }
        
        d4est_element_data_t* e_p_oriented [P4EST_HALF];
        d4est_element_data_reorient_f_p_elements_to_f_m_order(e_p, (P4EST_DIM)-1, f_m, f_p, orientation, faces_p, e_p_oriented);
        
        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &e_m[0],
           f_m,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           schwarz_geometric_data->face_h_type,
           j_div_sj_m_mortar_quad,
           hm_mortar_quad,
           faces_mortar,
           faces_m,
           nodes_mortar_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );
  
        d4est_mesh_calculate_mortar_h
          (
           p4est,
           &e_p_oriented[0],
           f_p,
           d4est_ops,
           d4est_geom,
           d4est_quad,
           schwarz_geometric_data->face_h_type,
           j_div_sj_p_mortar_quad,
           hp_mortar_quad,
           faces_mortar,
           faces_p,
           nodes_mortar_quad,
           d4est_mesh_get_size_parameters(d4est_factors)
          );

        for (int b = 0; b < ghost_mortar_side_data->total_mortar_nodes_quad; b++){
          D4EST_ASSERT(hp_mortar_quad[b] > 0);
        }
        
        /* DEBUG_PRINT_ARR_DBL(j_div_sj_p_mortar_quad, ghost_mortar_side_data->total_mortar_nodes_quad); */
        
        P4EST_FREE(j_div_sj_m_mortar_quad);
        P4EST_FREE(j_div_sj_p_mortar_quad);
        P4EST_FREE(j_div_sj_p_mortar_quad_porder);
      }
    }
  }

  schwarz_geometric_data->subdomain_mortars = P4EST_ALLOC(d4est_mortar_side_data_t*, schwarz_metadata->num_elements*(P4EST_FACES));
  schwarz_geometric_data->zero_and_skip_m = P4EST_ALLOC_ZERO(int, schwarz_metadata->num_elements*(P4EST_FACES)*(P4EST_HALF));
  schwarz_geometric_data->nodal_stride_m = P4EST_ALLOC_ZERO(int, schwarz_metadata->num_elements*(P4EST_FACES)*(P4EST_HALF));
  schwarz_geometric_data->nodal_stride_p = P4EST_ALLOC_ZERO(int, schwarz_metadata->num_elements*(P4EST_FACES)*(P4EST_HALF));
  schwarz_geometric_data->zero_and_skip_p = P4EST_ALLOC_ZERO(int, schwarz_metadata->num_elements*(P4EST_FACES)*(P4EST_HALF));

  for (int i = 0; i < schwarz_metadata->num_elements*(P4EST_FACES)*(P4EST_HALF); i++){
    schwarz_geometric_data->nodal_stride_m[i] = -1;
    schwarz_geometric_data->nodal_stride_p[i] = -1;
  }
    
  for (int i = 0; i < schwarz_metadata->num_elements*(P4EST_FACES); i++){
    schwarz_geometric_data->subdomain_mortars[i] = NULL;
  }

  schwarz_geometric_data->num_of_mortars_per_subdomain
    = P4EST_ALLOC(int,
                  schwarz_metadata->num_subdomains);

  schwarz_geometric_data->mortar_strides_per_subdomain
    = P4EST_ALLOC(int,
                  schwarz_metadata->num_subdomains);

  d4est_solver_schwarz_geometric_data_reduce_to_minimal_set
    (
     p4est,
     d4est_ghost,
     schwarz_metadata,
     schwarz_geometric_data
    );


  int ghost_volume_quad_size = 0;
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
    ghost_volume_quad_size += d4est_lgl_get_nodes((P4EST_DIM), ged->deg_quad);
  }

  schwarz_geometric_data->total_ghost_volume_quad_size
    = ghost_volume_quad_size;
  
  schwarz_geometric_data->J_quad_ghost = P4EST_ALLOC(double, ghost_volume_quad_size);
  schwarz_geometric_data->rst_xyz_quad_ghost = P4EST_ALLOC(double, (P4EST_DIM)*(P4EST_DIM)*ghost_volume_quad_size);
  schwarz_geometric_data->volume_quad_strides_per_ghost
    = P4EST_ALLOC(
     int,
     d4est_ghost->ghost->ghosts.elem_count
    );

  int volume_ghost_quad_stride = 0;
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
    int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ged->deg_quad);
    int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM), ged->deg);


    double* xyz_quad_ghost [P4EST_DIM];
    D4EST_ALLOC_DIM_VEC(xyz_quad_ghost, volume_nodes_quad);
    double* xyz_lobatto_ghost [P4EST_DIM];
    D4EST_ALLOC_DIM_VEC(xyz_lobatto_ghost, volume_nodes);
    double* xyz_rst_quad_ghost [P4EST_DIM][P4EST_DIM];
    D4EST_ALLOC_DBYD_MAT(xyz_rst_quad_ghost, volume_nodes_quad);

    double* J_quad_ghost
      = &schwarz_geometric_data->J_quad_ghost[volume_ghost_quad_stride];

    double* rst_xyz_quad_ghost [(P4EST_DIM)][(P4EST_DIM)];

    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
          rst_xyz_quad_ghost[i][j] =
            &schwarz_geometric_data->rst_xyz_quad_ghost
            [volume_ghost_quad_stride + (i*(P4EST_DIM) + j)*ghost_volume_quad_size];
      }
    }

    
    d4est_quadrature_volume_t mesh_object;
    mesh_object.dq =  ged->dq;
    mesh_object.tree = ged->tree;
    mesh_object.element_id = ged->id;
        
    mesh_object.q[0] = ged->q[0];
    mesh_object.q[1] = ged->q[1];
#if (P4EST_DIM)==3
    mesh_object.q[2] = ged->q[2];
#endif

    d4est_rst_t rst_points_quad;
    rst_points_quad = d4est_quadrature_get_rst_points
                      (
                       d4est_ops,
                       d4est_quad,
                       d4est_geom,
                       &mesh_object,
                       QUAD_OBJECT_VOLUME,
                       QUAD_INTEGRAND_UNKNOWN,
                       ged->deg_quad
                      );


    d4est_rst_t rst_points_lobatto;
    rst_points_lobatto.r = d4est_quadrature_lobatto_get_rst
                           (d4est_ops, NULL, NULL, &mesh_object,
                            QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ged->deg, 0);
    rst_points_lobatto.s = d4est_quadrature_lobatto_get_rst
                           (d4est_ops, NULL, NULL, &mesh_object,
                            QUAD_OBJECT_VOLUME, QUAD_INTEGRAND_UNKNOWN, ged->deg, 1);
    rst_points_lobatto.t = NULL;
#if (P4EST_DIM)==3
    rst_points_lobatto.t = d4est_quadrature_lobatto_get_rst
                           (d4est_ops, NULL, NULL, &mesh_object, QUAD_OBJECT_VOLUME,
                            QUAD_INTEGRAND_UNKNOWN, ged->deg, 2);
#endif
              



    d4est_geometry_compute_xyz
      (
       d4est_ops,
       d4est_geom,
       rst_points_lobatto,
       ged->tree,
       ged->deg,
       ged->q,
       ged->dq,
       xyz_lobatto_ghost
      );

    d4est_geometry_compute_xyz
      (
       d4est_ops,
       d4est_geom,
       rst_points_quad,
       ged->tree,
       ged->deg_quad,
       ged->q,
       ged->dq,
       xyz_quad_ghost
      );

        
    if (d4est_geom->DX_compute_method == GEOM_COMPUTE_NUMERICAL){
      double* tmp = P4EST_ALLOC(double, volume_nodes);
      for (int d = 0; d < (P4EST_DIM); d++){
        for (int d1 = 0; d1 < (P4EST_DIM); d1++){
          d4est_operators_apply_dij(d4est_ops, &(xyz_lobatto_ghost[d][0]), (P4EST_DIM), ged->deg, d1, tmp);
          d4est_quadrature_interpolate(
                                       d4est_ops,
                                       d4est_quad,
                                       d4est_geom,
                                       &mesh_object,
                                       QUAD_OBJECT_VOLUME,
                                       QUAD_INTEGRAND_UNKNOWN,
                                       tmp,
                                       ged->deg,
                                       &(xyz_rst_quad_ghost[d][d1][0]),
                                       ged->deg_quad
                                      );
        }
      }
      P4EST_FREE(tmp);
    }
    else if (d4est_geom->DX_compute_method == GEOM_COMPUTE_ANALYTIC){
      d4est_geometry_compute_dxyz_drst_analytic
        (
         d4est_ops,
         d4est_geom,
         rst_points_quad,
         ged->tree,
         ged->q,
         ged->dq,
         ged->deg_quad,
         xyz_rst_quad_ghost
        );
    }
    else {
      D4EST_ABORT("Not a supported compute method for DX");
    }

    d4est_geometry_compute_jacobian
      (
       xyz_rst_quad_ghost,
       J_quad_ghost,
       volume_nodes_quad
      );

    d4est_geometry_compute_drst_dxyz
      (
       xyz_rst_quad_ghost,
       J_quad_ghost,
       rst_xyz_quad_ghost,
       volume_nodes_quad
      );

    /* printf("rank %d, gid %d, stride %d, rst[1][0][0] = %.15f\n",  p4est->mpirank, gid, volume_ghost_quad_stride, rst_xyz_quad_ghost[1][0][0]); */

    D4EST_FREE_DIM_VEC(xyz_quad_ghost);
    D4EST_FREE_DIM_VEC(xyz_lobatto_ghost);
    D4EST_FREE_DBYD_MAT(xyz_rst_quad_ghost);
    schwarz_geometric_data->volume_quad_strides_per_ghost[gid] = volume_ghost_quad_stride;

   /*  double* rst_xyz_quad_ghost_test [(P4EST_DIM)][(P4EST_DIM)]; */
   /* for (int i = 0; i < (P4EST_DIM); i++){ */
   /*    for (int j = 0; j < (P4EST_DIM); j++){ */
   /*        rst_xyz_quad_ghost_test[i][j] = */
   /*          &schwarz_geometric_data->rst_xyz_quad_ghost */
   /*          [schwarz_geometric_data->volume_quad_strides_per_ghost[gid] + (i + j*(P4EST_DIM))*ghost_volume_quad_size]; */
   /*    } */
   /*  } */

   
   /* printf("rank %d, gid %d, rst[1][0][0] = %.15f, test = %.15f\n", p4est->mpirank, gid, rst_xyz_quad_ghost[1][0][0], rst_xyz_quad_ghost_test[1][0][0]); */
    
    volume_ghost_quad_stride += volume_nodes_quad;    
  }
  
  return schwarz_geometric_data;
}


void
d4est_solver_schwarz_geometric_data_check_hp
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
)  
{   
  for (int i = 0; i < schwarz_metadata->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[i];

    int num_mortars = schwarz_geometric_data->num_of_mortars_per_subdomain[i];
    int mortar_stride = schwarz_geometric_data->mortar_strides_per_subdomain[i];
    /* total_num_of_mortars += num_mortars; */
    /* printf("**************Subdomain %d**************\n", i); */
    for (int m = 0; m < num_mortars; m++){

      d4est_mortar_side_data_t* mortar_side_data =
        schwarz_geometric_data->subdomain_mortars[mortar_stride + m];

      int total_mortar_nodes_quad = mortar_side_data->total_mortar_nodes_quad;
      int total_mortar_nodes_lobatto = mortar_side_data->total_mortar_nodes_lobatto;
      int mortar_quad_stride = mortar_side_data->mortar_quad_stride;
      int boundary_vector_stride = (P4EST_DIM)*mortar_side_data->boundary_quad_stride;
      int boundary_quad_stride = mortar_side_data->boundary_quad_stride;
      /* printf("**************Mortar %d**************\n", m); */
      /* printf("total mortar nodes quad = %d\n", */
             /* mortar_side_data->total_mortar_nodes_quad); */
      /* printf("mpirank = %d\n", p4est->mpirank); */
      /* printf("subdomain = %d\n", i); */
      /* printf("is ghost = %d\n", mortar_side_data->is_ghost); */
      /* printf("faces_m = %d\n", mortar_side_data->faces_m); */
      /* printf("faces_p = %d\n", mortar_side_data->faces_p); */

      double* hp = NULL;
      if (mortar_side_data->faces_p != 0){
        if (mortar_side_data->is_ghost){
          hp = &schwarz_geometric_data->hp_mortar_quad[mortar_quad_stride];
        }
        else {
          hp = &d4est_factors->hm_mortar_quad[mortar_quad_stride];
        }
        /* for ( */
             /* int n = 0; */
             /* n < mortar_side_data->total_mortar_nodes_quad; */
             /* n++ */
        /* ){ */
          /* printf("hp[%d] = %.15f\n", n, hp[n]); */
        /* } */
      }
    }
  }
}


void
d4est_solver_schwarz_geometric_data_volume_sum_test
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
){
  double sum = 0;
  for (int subdomain = 0; subdomain < schwarz_metadata->num_subdomains; subdomain++){
  d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[subdomain];
 for (int el = 0; el < sub_data->num_elements; el++){

   d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[el];

   d4est_element_data_t* mesh_ed = NULL;
   double* J_quad = NULL;
   double* rst_xyz_quad [P4EST_DIM][P4EST_DIM];
   if (schwarz_ed->mpirank == p4est->mpirank){
     mesh_ed = d4est_element_data_get_ptr
               (
                p4est,
                schwarz_ed->tree,
                schwarz_ed->tree_quadid
               );

     J_quad = d4est_mesh_get_jacobian_on_quadrature_points
              (
               d4est_factors,
               mesh_ed
              );

     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad
                              [(i*(P4EST_DIM) + j)*d4est_factors->local_sizes.local_nodes_quad
                               + mesh_ed->quad_stride];
       }
     }     
   }
   else {
     mesh_ed = &d4est_ghost->ghost_elements[schwarz_ed->id];
     int ghost_quad_stride = schwarz_geometric_data->volume_quad_strides_per_ghost[schwarz_ed->id];
     int total_ghost_size = schwarz_geometric_data->total_ghost_volume_quad_size;
     J_quad = &schwarz_geometric_data->J_quad_ghost[ghost_quad_stride];
     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &schwarz_geometric_data->rst_xyz_quad_ghost
                              [(i*(P4EST_DIM) + j)*total_ghost_size
                               + ghost_quad_stride];
       }
     }     
     
   }
   int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM),
                                               schwarz_ed->deg);


   for (int k = 0; k < volume_nodes_quad; k++){
     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         sum += rst_xyz_quad[i][j][k];
       }
     }
     sum += J_quad[k];
   }

   /* d4est_laplacian_apply_stiffness_matrix_on_element */
   /*   ( */
   /*    d4est_ops, */
   /*    d4est_geom, */
   /*    d4est_quad, */
   /*    d4est_factors, */
   /*    mesh_ed, */
   /*    &restrict_transpose_u_over_subdomain[stride], */
   /*    &Au[stride] */
   /*   ); */
   /* stride += volume_nodes_lobatto; */
 }
  }

   double global_sum;
   sc_reduce(
             &sum,
             &global_sum,
             4,
             sc_MPI_DOUBLE,
             sc_MPI_SUM,
             0,
             sc_MPI_COMM_WORLD
   );

  if (p4est->mpirank == 0){
    printf("global_volume_data_sum = %.15f\n", global_sum);
  }
   
}


void
d4est_solver_schwarz_geometric_data_volume_sum_test_2
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
){
  /* double sum = 0; */
  for (int subdomain = 0; subdomain < schwarz_metadata->num_subdomains; subdomain++){
  d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[subdomain];
 for (int el = 0; el < sub_data->num_elements; el++){

   d4est_solver_schwarz_element_metadata_t* schwarz_ed = &sub_data->element_metadata[el];

   d4est_element_data_t* mesh_ed = NULL;
   double* J_quad = NULL;
   double* rst_xyz_quad [P4EST_DIM][P4EST_DIM];
   if (schwarz_ed->mpirank == p4est->mpirank){
     mesh_ed = d4est_element_data_get_ptr
               (
                p4est,
                schwarz_ed->tree,
                schwarz_ed->tree_quadid
               );

     J_quad = d4est_mesh_get_jacobian_on_quadrature_points
              (
               d4est_factors,
               mesh_ed
              );

     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &d4est_factors->rst_xyz_quad
                              [(i*(P4EST_DIM) + j)*d4est_factors->local_sizes.local_nodes_quad
                               + mesh_ed->quad_stride];
       }
     }

     /* printf("core_tree, elem, schwarz_ed_id, rst[1][0][0] = %d %d %d %.15f\n", */
            /* sub_data->core_tree, el, schwarz_ed->id, rst_xyz_quad[1][0][0]); */
   }
   else {
     mesh_ed = &d4est_ghost->ghost_elements[schwarz_ed->id];
     int ghost_quad_stride = schwarz_geometric_data->volume_quad_strides_per_ghost[schwarz_ed->id];
     int total_ghost_size = schwarz_geometric_data->total_ghost_volume_quad_size;

     
     J_quad = &schwarz_geometric_data->J_quad_ghost[ghost_quad_stride];
     
     for (int i = 0; i < (P4EST_DIM); i++){
       for (int j = 0; j < (P4EST_DIM); j++){
         rst_xyz_quad[i][j] = &schwarz_geometric_data->rst_xyz_quad_ghost
                              [(i*(P4EST_DIM) + j)*total_ghost_size
                               + ghost_quad_stride];
       }
     }     
     /* printf("rank %d core_tree, elem-isg, schwarz_ed_id, ghost_stride, rst[1][0][0] = %d %d %d %d %.15f\n",p4est->mpirank, */
            /* sub_data->core_tree, el, schwarz_ed->id, (1*(P4EST_DIM) + 0)*total_ghost_size + ghost_quad_stride,rst_xyz_quad[1][0][0]); */
   }
   int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM),
                                               schwarz_ed->deg);

   double sum = 0;
   for (int k = 0; k < volume_nodes_quad; k++){
     sum += J_quad[k];
   }

   /* D4EST_ASSERT(sub_data->element_metadata[0].is_core == 1); */
   if (schwarz_ed->is_core == 1){
   d4est_element_data_t*  core_ed = d4est_element_data_get_ptr
               (
                p4est,
                schwarz_ed->tree,
                schwarz_ed->tree_quadid
               );

   
/*    int ctree = core_ed->tree; */
/*    int q0 = core_ed->q[0]; */
/*    int q1 = core_ed->q[1]; */
/* #if (P4EST_DIM)==3 */
/*    int q2 = core_ed->q[2]; */
/* #endif */

/* #if (P4EST_DIM)==3 */
/*    printf("mpirank %d, subdomain %d, core tree %d, core q %d %d %d\n", p4est->mpirank, subdomain, ctree, q0, q1, q2); */
/* #else */
/*    printf("mpirank %d, subdomain %d, core tree %d, core q %d %d\n", p4est->mpirank, subdomain, ctree, q0, q1); */
/* #endif */
   }
   /* printf("Mpirank %d, Subdomain %d, Element %d, jac sum %.15f\n", p4est->mpirank, subdomain,  i, sum); */
   
   /* d4est_laplacian_apply_stiffness_matrix_on_element */
   /*   ( */
   /*    d4est_ops, */
   /*    d4est_geom, */
   /*    d4est_quad, */
   /*    d4est_factors, */
   /*    mesh_ed, */
   /*    &restrict_transpose_u_over_subdomain[stride], */
   /*    &Au[stride] */
   /*   ); */
   /* stride += volume_nodes_lobatto; */
 }
  }

   /* double global_sum; */
   /* sc_reduce( */
             /* &eesum, */
             /* &global_sum, */
             /* 4, */
             /* sc_MPI_DOUBLE, */
             /* sc_MPI_SUM, */
             /* 0, */
             /* sc_MPI_COMM_WORLD */
   /* ); */

  /* if (p4est->mpirank == 0){ */
    /* printf("global_volume_data_sum = %.15f\n", global_sum); */
  /* } */
   
}


void
d4est_solver_schwarz_geometric_data_sum_test
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost,
 d4est_mesh_data_t* d4est_factors,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data
)
{
  double sj_sum = 0;
  double xyz_sum = 0;
  double xyz_lobatto_sum = 0;
  double hm_sum = 0;
  double hp_sum = 0;
  int nodal_m_sum = 0;
  int nodal_p_sum = 0;
  
  int total_num_of_mortars = 0;
  int num_of_bndry_mortars = 0;
  
  for (int i = 0; i < schwarz_metadata->num_subdomains; i++){
    d4est_solver_schwarz_subdomain_metadata_t* sub_data = &schwarz_metadata->subdomain_metadata[i];

    int num_mortars = schwarz_geometric_data->num_of_mortars_per_subdomain[i];
    int mortar_stride = schwarz_geometric_data->mortar_strides_per_subdomain[i];
    total_num_of_mortars += num_mortars;
    
    for (int m = 0; m < num_mortars; m++){

      d4est_mortar_side_data_t* mortar_side_data =
        schwarz_geometric_data->subdomain_mortars[mortar_stride + m];

      int total_mortar_nodes_quad = mortar_side_data->total_mortar_nodes_quad;
      int total_mortar_nodes_lobatto = mortar_side_data->total_mortar_nodes_lobatto;
      int mortar_quad_stride = mortar_side_data->mortar_quad_stride;
      int boundary_vector_stride = (P4EST_DIM)*mortar_side_data->boundary_quad_stride;
      int boundary_quad_stride = mortar_side_data->boundary_quad_stride;
      
      double* sj = NULL;
      double* hm = NULL;
      double* hp = NULL;
      double* xyz = NULL;
      double* xyz_lobatto = NULL;

      for (int fm = 0; fm < mortar_side_data->faces_m; fm++){
        nodal_m_sum += schwarz_geometric_data->nodal_stride_m[(m + mortar_stride)*(P4EST_HALF) + fm];
      }
      for (int fp = 0; fp < mortar_side_data->faces_p; fp++){
        nodal_p_sum += schwarz_geometric_data->nodal_stride_p[(m + mortar_stride)*(P4EST_HALF) + fp];
      }
      
      if(mortar_side_data->is_ghost){
        sj = &schwarz_geometric_data->sj_m_mortar_quad[mortar_quad_stride];
        hm = &schwarz_geometric_data->hm_mortar_quad[mortar_quad_stride];
        hp = &schwarz_geometric_data->hp_mortar_quad[mortar_quad_stride];

      }
      else {
        d4est_mesh_local_strides_t* strides_check =
          &d4est_factors->local_strides[mortar_side_data->e_m[0].id];
        D4EST_ASSERT(strides_check->mortar_quad_stride[mortar_side_data->f_m] == mortar_quad_stride);
        D4EST_ASSERT(strides_check->boundary_quad_stride[mortar_side_data->f_m] == boundary_quad_stride);
            
        sj = &d4est_factors->sj_m_mortar_quad[mortar_quad_stride];
        hm = &d4est_factors->hm_mortar_quad[mortar_quad_stride];
        hp = &d4est_factors->hp_mortar_quad[mortar_quad_stride];
      }

      /* if (mortar_side_data->tree_m != mortar_side_data->tree_p && mortar_side_data->tree_p != -1){ */
        /* DEBUG_PRINT_3ARR_DBL(sj,hm,hp, mortar_side_data->total_mortar_nodes_quad); */
      /* } */

      
      if (mortar_side_data->faces_p == 0){
        num_of_bndry_mortars++;
        if(mortar_side_data->is_ghost){
          xyz = &schwarz_geometric_data->xyz_on_f_m_quad[boundary_vector_stride];
          xyz_lobatto = &schwarz_geometric_data->xyz_on_f_m_lobatto[boundary_vector_stride];
        }
        else {
          xyz = &d4est_factors->xyz_m_mortar_quad[boundary_vector_stride];
          xyz_lobatto = &d4est_factors->xyz_m_mortar_lobatto[boundary_vector_stride];
        }
      }

      if (xyz != NULL){
        for (int k = 0; k < (P4EST_DIM)*total_mortar_nodes_quad; k++){
          xyz_sum += xyz[k];

        }
        for (int d = 0; d < (P4EST_DIM); d++){
          for (int k = 0; k < total_mortar_nodes_lobatto; k++){
            xyz_lobatto_sum += xyz_lobatto[d*(total_mortar_nodes_quad) + k];
          }
        }
      }

      for (int k = 0; k < total_mortar_nodes_quad; k++){
        sj_sum += sj[k];
        hm_sum += hm[k];
        if (mortar_side_data->faces_p!=0){
        hp_sum += hp[k];
        if (hp[k] <= 0){
          printf("mpirank %d, subdomain %d, mortar %d, mortar_side_data->is_ghost = %d, hp[%d] > 0 = %.15f\n",
                 p4est->mpirank, i, m,
                 mortar_side_data->is_ghost, k, hp[k]);
        }
        D4EST_ASSERT(hp[k] > 0);
        }

      }      
    }
  }

  int local_num_mortars [] = {total_num_of_mortars, num_of_bndry_mortars, nodal_m_sum, nodal_p_sum};
  int global_num_mortars [4];

  double global_sums [5];
  double local_sums [] = {sj_sum, xyz_sum, xyz_lobatto_sum, hm_sum, hp_sum};
  sc_reduce(
            &local_sums[0],
            &global_sums,
            5,
            sc_MPI_DOUBLE,
            sc_MPI_SUM,
            0,
            sc_MPI_COMM_WORLD
  );

  sc_reduce(
            &local_num_mortars,
            &global_num_mortars,
            4,
            sc_MPI_INT,
            sc_MPI_SUM,
            0,
            sc_MPI_COMM_WORLD
  );


  
  if (p4est->mpirank == 0){
    printf("xyz sum = %.15f, sj_sum = %.15f, xyz_lobatto_sum = %.15f, hm_sum = %.15f, hp_sum = %.15f\n",
           global_sums[1], global_sums[0], global_sums[2], global_sums[3], global_sums[4]);
    printf("num of mortars = %d, bndyr = %d, nodal_m_sum = %d, nodal_p_sum = %d\n", global_num_mortars[0], global_num_mortars[1],
           global_num_mortars[2], global_num_mortars[3]);
  }
}

 
