#include <pXest.h>
#include <d4est_amr.h>
#include <d4est_xyz_functions.h>
#include <d4est_util.h>
#include <d4est_mesh.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
#include <d4est_amr_smooth_pred.h>
#include <d4est_amr_uniform.h>
#include <d4est_amr_random.h>
#include <ini.h>
#include <zlog.h>


static int
d4est_amr_refine_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{
  d4est_amr_t* d4est_amr = (d4est_amr_t*) p4est->user_pointer;
  d4est_element_data_t* elem_data = (d4est_element_data_t *) quadrant->p.user_data;
  
  int* refinement_log = d4est_amr->refinement_log;
  int* initial_log = d4est_amr->initial_log;
  initial_log[elem_data->id] = elem_data->deg;
  
  /* h-refine */
  if (refinement_log[elem_data->id] < 0){
    return 1;
  }
  /* p-refine, p-coarsen or do nothing */
  else {
    if (refinement_log[elem_data->id] > d4est_amr->max_degree){
      refinement_log[elem_data->id] = d4est_amr->max_degree;
    }
    elem_data->deg = refinement_log[elem_data->id];
    return 0;
  }
}


static void
d4est_amr_refine_replace_callback (
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
  D4EST_ASSERT(num_outgoing == 1 && num_incoming == (P4EST_CHILDREN));
  d4est_amr_t* d4est_amr = (d4est_amr_t*) p4est->user_pointer;
  d4est_operators_t* d4est_ops = d4est_amr->d4est_ops;
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  
  for (int i = 0; i < P4EST_CHILDREN; i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    child_data->deg = abs(d4est_amr->refinement_log[parent_data->id]);
  }

  if(d4est_amr->scheme->refine_replace_callback_fcn_ptr != NULL)
    d4est_amr->scheme->refine_replace_callback_fcn_ptr(p4est, which_tree, num_outgoing, outgoing, num_incoming, incoming);
}

static void
d4est_amr_balance_replace_callback (
			     p4est_t * p4est,
			     p4est_topidx_t which_tree,
			     int num_outgoing,
			     p4est_quadrant_t * outgoing[],
			     int num_incoming,
			     p4est_quadrant_t * incoming[]
			     )
{
  D4EST_ASSERT(num_outgoing == 1 && num_incoming == (P4EST_CHILDREN));
  d4est_amr_t* d4est_amr = (d4est_amr_t*) p4est->user_pointer;
  d4est_operators_t* d4est_ops = d4est_amr->d4est_ops;
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  
  int i;
  d4est_amr->balance_log[parent_data->id] = -parent_data->deg;
  for (i = 0; i < P4EST_CHILDREN; i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
  }

  if(d4est_amr->scheme->balance_replace_callback_fcn_ptr != NULL)
    d4est_amr->scheme->balance_replace_callback_fcn_ptr(p4est, which_tree, num_outgoing, outgoing, num_incoming, incoming);
}

static void
d4est_amr_balance_elements
(
 p4est_t* p4est
)
{
  d4est_amr_t* d4est_amr = p4est->user_pointer;
  d4est_amr->balance_log_size = p4est->local_num_quadrants;
  d4est_amr->balance_log = D4EST_REALLOC(d4est_amr->balance_log,
                                          int,
                                         d4est_amr->balance_log_size);
  
  /* we need the quad-ids to store the balance info, so update them */
  int id_stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* elem_data = (d4est_element_data_t*)(quad->p.user_data);
        elem_data->id = id_stride;
        d4est_amr->balance_log[id_stride] = elem_data->deg;
        id_stride++;
      }
    }
  
  p4est_balance_ext
    (
     p4est,
     P4EST_CONNECT_FULL,
     NULL,
     d4est_amr_balance_replace_callback
    );
}

static void
d4est_amr_refine_elements
(
 p4est_t* p4est
)
{
  p4est_refine_ext
    (
     p4est,
     0,
     -1,
     d4est_amr_refine_callback,
     NULL,
     d4est_amr_refine_replace_callback
    );
}

static void
d4est_amr_mark_elements
(
 p4est_t* p4est
)
{
  d4est_amr_t* d4est_amr = p4est->user_pointer;
  
  d4est_amr->initial_log_size = p4est->local_num_quadrants;

  d4est_amr->initial_log = D4EST_REALLOC
                     (
                      d4est_amr->initial_log,
                      int,
                      d4est_amr->initial_log_size
                     );

  d4est_amr->refinement_log = D4EST_REALLOC
                        (
                         d4est_amr->refinement_log,
                         int,
                         d4est_amr->initial_log_size
                        );

  p4est_iterate(p4est,
                NULL,
                NULL,
                d4est_amr->scheme->mark_elements,
		NULL,
#if P4EST_DIM==3
                NULL,
#endif
		NULL);
}


static void
d4est_amr_interpolate_field_on_element
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 double* in,
 int in_deg,
 double* out,
 int interp_code
)
{
  int* out_deg = D4EST_ALLOC(int, (P4EST_CHILDREN));
  if (interp_code < 0){
    /* we only support same-degree h-refining atm */
    D4EST_ASSERT(-interp_code == in_deg);
    for (int i = 0; i < (P4EST_CHILDREN); i++)
      out_deg[i] = -interp_code;
  }
  else {
    out_deg[0] = interp_code;
  }

  /* h-refine interpolation */
  if (interp_code < 0){
    d4est_operators_apply_hp_prolong
      (
       d4est_ops,
       in,
       in_deg,
       (P4EST_DIM),
       &(out_deg[0]),
       out
      );
  }

  /* p-refine interpolation */
  else if (interp_code >= in_deg){
    d4est_operators_apply_p_prolong
      (
       d4est_ops,
       in,
       in_deg,
       (P4EST_DIM),
       out_deg[0],
       out
      );
  }

  else {
    zlog_category_t *c_default = zlog_get_category("d4est_amr");
    zlog_info(c_default, "interp_code = %d", interp_code);
    zlog_info(c_default, "indeg = %d", in_deg);
    zlog_info(c_default, "out_deg[0] = %d", out_deg[0]);
    D4EST_ABORT("hp amr code should be >= deg or -deg, coarsening is currently not supported in amr");
  }
  
  D4EST_FREE(out_deg);
}

static void
d4est_amr_interpolate_field
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_amr_t* d4est_amr,
 double** field,
 int post_balance_nodes
)
{
  double *aux = D4EST_ALLOC(double, post_balance_nodes);

  /* interpolate from initial to auxiliary (unbalanced) grid */
  int pre_refine_stride = 0;
  int post_refine_stride = 0;
  for (int i = 0; i < d4est_amr->initial_log_size; i++){
    int pre_volume_nodes =  d4est_lgl_get_nodes
                            (
                             (P4EST_DIM),
                             d4est_amr->initial_log[i]
                            );
        
    int post_volume_nodes = d4est_lgl_get_nodes
                            (
                             (P4EST_DIM),
                             abs(d4est_amr->refinement_log[i])
                            );

        
    if(d4est_amr->refinement_log[i] < 0)
      post_volume_nodes *= (P4EST_CHILDREN);


    d4est_amr_interpolate_field_on_element
      (
       p4est,
       d4est_ops,
       &((*field)[pre_refine_stride]),
       d4est_amr->initial_log[i],
       &aux[post_refine_stride],
       d4est_amr->refinement_log[i]
      );

    pre_refine_stride += pre_volume_nodes;
    post_refine_stride += post_volume_nodes;
  }

  *field = D4EST_REALLOC(*field, double, post_balance_nodes);
  
  /* interpolate from auxiliary to balanced grid */
  int pre_balance_stride = 0;
  int post_balance_stride = 0;
  for (int i = 0; i < d4est_amr->balance_log_size; i++){
    int pre_balance_volume_nodes =  d4est_lgl_get_nodes
                                    (
                                     (P4EST_DIM),
                                     abs(d4est_amr->balance_log[i])
                                    );
        
    int post_balance_volume_nodes = d4est_lgl_get_nodes
                                    (
                                     (P4EST_DIM),
                                     abs(d4est_amr->balance_log[i])
                                    );

        
    if(d4est_amr->balance_log[i] < 0)
      post_balance_volume_nodes *= (P4EST_CHILDREN);


    d4est_amr_interpolate_field_on_element
      (
       p4est,
       d4est_ops,
       &aux[pre_balance_stride],
       abs(d4est_amr->balance_log[i]),
       &((*field)[post_balance_stride]),
       d4est_amr->balance_log[i]
      );

    pre_balance_stride += pre_balance_volume_nodes;
    post_balance_stride += post_balance_volume_nodes;
  }
    
  D4EST_FREE(aux);
}
 


static
int d4est_amr_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_amr");

  d4est_amr_t* pconfig = (d4est_amr_t*)user;

  if (d4est_util_match_couple(section,"amr",name,"scheme")) {
    if(d4est_util_match(value, "smooth_pred")){
      pconfig->scheme->amr_scheme_type = AMR_SMOOTH_PRED;
      if(pconfig->mpirank == 0)
        zlog_info(c_default, "loading scheme = %s", "smooth_pred");
    }
    else if (d4est_util_match(value, "uniform_h")){
      pconfig->scheme->amr_scheme_type = AMR_UNIFORM_H;
      if(pconfig->mpirank == 0)
        zlog_info(c_default, "loading scheme = %s", "uniform_h");
    }
    else if (d4est_util_match(value, "uniform_p")){
      pconfig->scheme->amr_scheme_type = AMR_UNIFORM_P;
      if(pconfig->mpirank == 0)
        zlog_info(c_default, "loading scheme = %s", "uniform_p");
    }
    else if (d4est_util_match(value, "random_h")){
      pconfig->scheme->amr_scheme_type = AMR_RANDOM_H;
      if(pconfig->mpirank == 0)
        zlog_info(c_default, "loading scheme = %s", "random_h");
    }
    else if (d4est_util_match(value, "random_hp")){
      pconfig->scheme->amr_scheme_type = AMR_RANDOM_HP;
      if(pconfig->mpirank == 0)
        zlog_info(c_default, "loading scheme = %s", "random_hp");
    }
    else {
      zlog_error(c_default, "You tried to use %s as a mapping type.", value);
      D4EST_ABORT("This mapping is not supported.");
    }
  }
  else if (d4est_util_match_couple(section,"amr",name,"num_of_amr_steps")) {
    D4EST_ASSERT(pconfig->num_of_amr_steps == -1);
    pconfig->num_of_amr_steps = atoi(value);
    D4EST_ASSERT(pconfig->num_of_amr_steps >= 0);
  }
  else if (d4est_util_match_couple(section,"amr",name,"max_degree")) {
    D4EST_ASSERT(pconfig->max_degree == -1);
    pconfig->max_degree = atoi(value);
    D4EST_ASSERT(atoi(value) > 0);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_amr_input
(
 const char* input_file,
 d4est_amr_t* d4est_amr
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_amr");

  d4est_amr->max_degree = -1;
  d4est_amr->num_of_amr_steps = -1;
  d4est_amr->scheme->amr_scheme_type = AMR_NOT_SET;

  if (ini_parse(input_file, d4est_amr_input_handler, d4est_amr) < 0) {
    D4EST_ABORT("Can't load input file in d4est_amr_input.");
  }

  D4EST_CHECK_INPUT("amr", d4est_amr->scheme->amr_scheme_type, AMR_NOT_SET);
  D4EST_CHECK_INPUT("amr", d4est_amr->max_degree, -1);
  D4EST_CHECK_INPUT("amr", d4est_amr->num_of_amr_steps, -1);

  if(d4est_amr->mpirank == 0){
    zlog_info(c_default, "num_of_amr_steps = %d",  d4est_amr->num_of_amr_steps);
    zlog_info(c_default, "max_degree = %d",  d4est_amr->max_degree);
  }
}

d4est_amr_t*
d4est_amr_init
(
 p4est_t* p4est,
 const char* input_file,
 void* scheme_data
)
{
  d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
  d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);

  d4est_amr->mpirank = p4est->mpirank;
  d4est_amr->scheme = scheme;
  d4est_amr->balance_log = NULL;
  d4est_amr->refinement_log = NULL;
  d4est_amr->initial_log = NULL;
  
  d4est_amr_input(input_file, d4est_amr);
  
  if (scheme->amr_scheme_type == AMR_SMOOTH_PRED) {
    d4est_amr_smooth_pred_init(p4est, input_file, scheme, scheme_data);
  }
  else if (scheme->amr_scheme_type == AMR_UNIFORM_H ||
          scheme->amr_scheme_type == AMR_UNIFORM_P){
    d4est_amr_uniform_init(p4est, input_file, scheme, scheme_data);
  }
  else if (scheme->amr_scheme_type == AMR_RANDOM_H ||
          scheme->amr_scheme_type == AMR_RANDOM_HP){
    d4est_amr_random_init(p4est, input_file, scheme, scheme_data);
  }
  else {
    D4EST_ABORT("This amr scheme is currently not supported.");
  }

  return d4est_amr;
}

d4est_amr_t*
d4est_amr_init_uniform_h
(
 p4est_t* p4est,
 int max_degree,
 int num_of_amr_steps
)
{
  d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
  d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
  
  d4est_amr->scheme = scheme;
  d4est_amr->balance_log = NULL;
  d4est_amr->refinement_log = NULL;
  d4est_amr->initial_log = NULL;
  d4est_amr->max_degree = max_degree;
  d4est_amr->num_of_amr_steps = num_of_amr_steps;
  scheme->amr_scheme_type = AMR_UNIFORM_H;
  d4est_amr_uniform_init(p4est, NULL, scheme, NULL);

  return d4est_amr;

}

d4est_amr_t*
d4est_amr_init_uniform_p
(
 p4est_t* p4est,
 int max_degree,
 int num_of_amr_steps
)
{
  d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
  d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
  
  d4est_amr->scheme = scheme;
  d4est_amr->balance_log = NULL;
  d4est_amr->refinement_log = NULL;
  d4est_amr->initial_log = NULL;
  d4est_amr->max_degree = max_degree;
  d4est_amr->num_of_amr_steps = num_of_amr_steps;
  scheme->amr_scheme_type = AMR_UNIFORM_P;
  d4est_amr_uniform_init(p4est, NULL, scheme, NULL);

  return d4est_amr;

}


d4est_amr_t*
d4est_amr_init_random_hp
(
 p4est_t* p4est,
 int max_degree,
 int num_of_amr_steps
)
{
  d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
  d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
  
  d4est_amr->scheme = scheme;
  d4est_amr->balance_log = NULL;
  d4est_amr->refinement_log = NULL;
  d4est_amr->initial_log = NULL;
  d4est_amr->max_degree = max_degree;
  d4est_amr->num_of_amr_steps = num_of_amr_steps;
  scheme->amr_scheme_type = AMR_RANDOM_HP;
  d4est_amr_random_init(p4est, NULL, scheme, NULL);

  return d4est_amr;
}



static void
d4est_amr_custom_destroy(d4est_amr_scheme_t* scheme)
{
  P4EST_FREE(scheme);
}

d4est_amr_t*
d4est_amr_custom_init
(
 p4est_t* p4est,
 int max_degree,
 int num_of_amr_steps,
 void(*d4est_amr_custom_mark_elements)(p4est_iter_volume_info_t*,void*),
 void* user
)
{
  d4est_amr_t* d4est_amr = P4EST_ALLOC(d4est_amr_t, 1);
  d4est_amr_scheme_t* scheme = P4EST_ALLOC(d4est_amr_scheme_t, 1);
  
  d4est_amr->scheme = scheme;
  d4est_amr->balance_log = NULL;
  d4est_amr->refinement_log = NULL;
  d4est_amr->initial_log = NULL;
  d4est_amr->max_degree = max_degree;
  d4est_amr->num_of_amr_steps = num_of_amr_steps;

  scheme->amr_scheme_type = AMR_CUSTOM;
  scheme->pre_refine_callback
    = NULL;
  
  scheme->balance_replace_callback_fcn_ptr
    = NULL;

  scheme->refine_replace_callback_fcn_ptr
    = NULL;

  scheme->amr_scheme_data
    = user;

  scheme->post_balance_callback
    = NULL;

  scheme->mark_elements
    = d4est_amr_custom_mark_elements;
  
  scheme->destroy
    = d4est_amr_custom_destroy;
  
  return d4est_amr;
}
 
 


void
d4est_amr_step
(
 p4est_t* p4est,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_amr_t* d4est_amr,
 double** field,
 double* d4est_estimator,
 d4est_estimator_stats_t** stats
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_amr");

  d4est_amr->d4est_estimator_stats = stats;
  d4est_amr->d4est_estimator = d4est_estimator;
  
  void* backup = p4est->user_pointer;
  p4est->user_pointer = d4est_amr;
    
  if(d4est_amr->scheme->pre_refine_callback != NULL){
    d4est_amr->scheme->pre_refine_callback(p4est, d4est_amr->scheme->amr_scheme_data);
  }
  if (p4est->mpirank == 0)
    zlog_info(c_default, "Starting to mark elements");

  d4est_amr_mark_elements(p4est);

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Starting to refine elements");

  d4est_amr_refine_elements(p4est);

  if (p4est->mpirank == 0)
    zlog_info(c_default, "Starting to balance elements");

  d4est_amr_balance_elements(p4est);

  if(d4est_amr->scheme->post_balance_callback != NULL){
    d4est_amr->scheme->post_balance_callback(p4est, d4est_amr->scheme->amr_scheme_data);
  }

  if(field != NULL)
    d4est_amr_interpolate_field
      (
       p4est,
       d4est_ops,
       d4est_amr,
       field,
       d4est_mesh_get_local_nodes(p4est)
      );
  
  p4est->user_pointer = backup;

  p4est_ghost_destroy(*ghost);
  P4EST_FREE(*ghost_data);
  *ghost = NULL;
  *ghost_data = NULL;
  *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
  *ghost_data = P4EST_ALLOC(d4est_element_data_t, (*ghost)->ghosts.elem_count);
}

void
d4est_amr_destroy
(
 d4est_amr_t* d4est_amr
)
{
  D4EST_FREE(d4est_amr->balance_log);
  D4EST_FREE(d4est_amr->refinement_log);
  D4EST_FREE(d4est_amr->initial_log);
  d4est_amr->scheme->destroy(d4est_amr->scheme);
  D4EST_FREE(d4est_amr);
}
