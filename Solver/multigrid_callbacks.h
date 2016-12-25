#ifndef MULTIGRID_CALLBACKS_H
#define MULTIGRID_CALLBACKS_H 

static int
multigrid_coarsen (p4est_t * p4est,
			 p4est_topidx_t which_tree,
			 p4est_quadrant_t * children[])
{
  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = (multigrid_refine_data_t*) mg_data->coarse_grid_refinement;
  int* stride = &mg_data->stride;

  if (children[1] != NULL){
    coarse_grid_refinement[*stride].hrefine = 1;
    
    int i;
    for (i = 0; i < (P4EST_CHILDREN); i++){
      coarse_grid_refinement[*stride].degh[i] = ((element_data_t*)children[i]->p.user_data)->deg;
    }
  }
  else{
    coarse_grid_refinement[*stride].hrefine = 0;
    coarse_grid_refinement[*stride].degh[0] = ((element_data_t*)children[0]->p.user_data)->deg;
    coarse_grid_refinement[*stride].degH = ((element_data_t*)children[0]->p.user_data)->deg;
  }
  (*stride)++;
  return 1;
}

static void
multigrid_coarsen_init
(
 p4est_t *p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  element_data_t* element_data = (element_data_t*)quadrant->p.user_data;
  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;

  int* stride = &mg_data->stride;
  int children = (coarse_grid_refinement[*stride - 1].hrefine == 1) ? (P4EST_CHILDREN) : 1;
  int min_deg = coarse_grid_refinement[*stride - 1].degh[0];

  int i;
  if (children == (P4EST_CHILDREN))
    for(i = 0; i < children; i++){
      if (coarse_grid_refinement[*stride - 1].degh[i] < min_deg)
        min_deg = coarse_grid_refinement[*stride - 1].degh[i];
    }

  element_data->deg = min_deg;
  coarse_grid_refinement[*stride - 1].degH = min_deg;
}

static void
multigrid_store_balance_changes
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
)
{
  if (num_outgoing != 1)
    {
      printf("num_outgoing != 1 in store_balance_changes");
      exit(1);
    }
  
  element_data_t* parent_data = (element_data_t*) outgoing[0]->p.user_data;
  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;
  int stride = mg_data->stride;
  
  coarse_grid_refinement[stride + parent_data->id].hrefine = 2;

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    element_data_t* child_data = (element_data_t*) incoming[i]->p.user_data;
    child_data->deg = coarse_grid_refinement[stride + parent_data->id].degh[i];
  }
}

static void
multigrid_apply_restriction
(
 p4est_iter_volume_info_t* info,
 void* user_data
)
{
  multigrid_data_t* mg_data = (multigrid_data_t*) info->p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;

  if(mg_data->user_defined_fields && mg_data->mg_restrict_user_callback != NULL){
    mg_data->mg_restrict_user_callback(info, user_data);
  }
  
  double* x = mg_data->intergrid_ptr;
  int* stride = &mg_data->stride;
  int* fine_stride = &mg_data->fine_stride;
  int* coarse_stride = &mg_data->coarse_stride;
  int* temp_stride = &mg_data->temp_stride;
  int fine_nodes = mg_data->fine_nodes;

  double* x_children = &x[*fine_stride];
  double* x_parent = &x[fine_nodes + *coarse_stride];

  if (coarse_grid_refinement[(*stride)].hrefine == 0){
    int degh = coarse_grid_refinement[(*stride)].degh[0];
    int degH = coarse_grid_refinement[(*stride)].degH;

    dgmath_apply_p_prolong_transpose
      (
       dgmath_jit_dbase,
       x_children,
       degh,
       (P4EST_DIM),
       degH,
       x_parent
      );

    (*stride)++;
    (*fine_stride)+=dgmath_get_nodes((P4EST_DIM), degh);
    (*coarse_stride)+=dgmath_get_nodes((P4EST_DIM), degH);
  }
  
  else if (coarse_grid_refinement[(*stride)].hrefine == 1){
    int* degh = &(coarse_grid_refinement[(*stride)].degh[0]);
    int degH = coarse_grid_refinement[(*stride)].degH;

    dgmath_apply_hp_prolong_transpose
      (
       dgmath_jit_dbase,
       x_children,
       degh,
       (P4EST_DIM),
       degH,
       x_parent
      );

    (*stride)++;
    for (int i = 0; i < (P4EST_CHILDREN); i++){
      (*fine_stride) += dgmath_get_nodes( (P4EST_DIM) , degh[i] );
    }
    (*coarse_stride) += dgmath_get_nodes( (P4EST_DIM) , degH);
  }
  else {
    int deg = coarse_grid_refinement[(*stride)].degh[*temp_stride];
    int nodes = dgmath_get_nodes( (P4EST_DIM) , deg );
    linalg_copy_1st_to_2nd(x_children,
                           x_parent,
                           nodes
                          );
    if (*temp_stride == (P4EST_CHILDREN)-1){
      (*fine_stride) += nodes;
      (*coarse_stride) += nodes;
      (*stride)++;
      *temp_stride = 0;
    }
    else{
      (*fine_stride) += nodes;
      (*coarse_stride) += nodes;
      *temp_stride = *temp_stride + 1;
    }
  }
}

static void
multigrid_refine_and_apply_prolongation_replace
(
 p4est_t *p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t *outgoing[],
 int num_incoming,
 p4est_quadrant_t *incoming[]
)
{
  /* printf("multigrid refine and apply prolongation replace function\n"); */
  if (num_outgoing != 1)
    {
      printf("num_outgoing != 1 in store_balance_changes");
      exit(1);
    }
  
  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;

  /* the refine function increments the stride, so the stride of
     the parent element we want is -1 */
  int stride = mg_data->stride - 1;
  /* printf("mg_data->stride = %d\n", mg_data->stride); */
  /* printf("stride = %d\n", stride); */
  
  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++){
    element_data_t* child_data = (element_data_t*) incoming[i]->p.user_data;
    child_data->deg = coarse_grid_refinement[stride].degh[i];
  }
}

static int
multigrid_refine_and_apply_prolongation
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{  
  multigrid_data_t* mg_data = (multigrid_data_t*) p4est->user_pointer;
  multigrid_refine_data_t* coarse_grid_refinement = mg_data->coarse_grid_refinement;
  dgmath_jit_dbase_t* dgmath_jit_dbase = mg_data->dgmath_jit_dbase;

  if(mg_data->user_defined_fields && mg_data->mg_prolong_user_callback != NULL)
    mg_data->mg_prolong_user_callback(p4est,which_tree,quadrant);
  
  double* x = mg_data->intergrid_ptr;

  int* stride = &mg_data->stride;
  int* fine_stride = &mg_data->fine_stride;
  int* coarse_stride = &mg_data->coarse_stride;
  int* temp_stride = &mg_data->temp_stride;
  int fine_nodes = mg_data->fine_nodes;

  double* x_children = &x[*fine_stride];
  double* x_parent = &x[fine_nodes + *coarse_stride];
  
  if (coarse_grid_refinement[(*stride)].hrefine == 0){
    int degh = coarse_grid_refinement[(*stride)].degh[0];
    int degH = coarse_grid_refinement[(*stride)].degH;
    dgmath_apply_p_prolong
      (
       dgmath_jit_dbase,
       x_parent,
       degH,
       (P4EST_DIM),
       degh,
       x_children
      );
    
    (*stride)++;
    (*fine_stride)+=dgmath_get_nodes((P4EST_DIM), degh);
    (*coarse_stride)+=dgmath_get_nodes((P4EST_DIM), degH);
    
    return 0;
  }
  else if (coarse_grid_refinement[(*stride)].hrefine == 1){
 
   int* degh = &(coarse_grid_refinement[(*stride)].degh[0]);
   int degH = coarse_grid_refinement[(*stride)].degH;

    dgmath_apply_hp_prolong
      (
       dgmath_jit_dbase,
       x_parent,
       degH,
       (P4EST_DIM),
       degh,
       x_children
      );   
    
    (*stride)++;
    int i;
    for (i = 0; i < (P4EST_CHILDREN); i++)
      (*fine_stride) += dgmath_get_nodes( (P4EST_DIM) , degh[i] );
    (*coarse_stride) += dgmath_get_nodes( (P4EST_DIM) , degH);

    return 1;
  }
  else {
    int deg = coarse_grid_refinement[(*stride)].degh[*temp_stride];
    int nodes = dgmath_get_nodes( (P4EST_DIM) , deg );
    
    linalg_copy_1st_to_2nd(x_parent,
                           x_children,
                           nodes
                          );    
    if (*temp_stride == (P4EST_CHILDREN)-1){
      (*fine_stride) += nodes;
      (*coarse_stride) += nodes;
      (*stride)++;
      *temp_stride = 0;
    }
    else{
      (*fine_stride) += nodes;
      (*coarse_stride) += nodes;
      *temp_stride = *temp_stride + 1;
    }

    return 0;
  }
}

#endif
