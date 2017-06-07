#define SAFETY
#define NDEBUG

#include "../hpAMR/curved_hp_amr.h"
#include "../GridFunctions/grid_functions.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../ElementData/d4est_element_data.h"
#include "../Support/Visualization/curved_hacked_p4est_vtk.h"

/**
 * @file   curved_hp_amr.c
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Mon Nov  2 01:08:28 2015
 * 
 * @brief  Here are the rules for the refinement_log
 * array :
 *
 * This is the interface for running hp-amr loops.
 * 
 * Suppose the element currently has degree p, then:
 *
 * if refinement_log[i] = p' with p' < 0
 * then we h-refine in element i with the absolute
 * value of p as the degree of the parent
 * and children octants
 *
 * if refinement_log[i] = p' with p' > 0 then we
 * p-refine if p' > p and we p-coarsen if p' < p
 * We do nothing if p' = p.
 *
 */
void
curved_hp_amr_replace_callback_default
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 int num_outgoing,
 p4est_quadrant_t * outgoing[],
 int num_incoming,
 p4est_quadrant_t * incoming[]
)
{
#ifdef SAFETY  
  mpi_assert(num_outgoing == 1 && num_incoming == (P4EST_CHILDREN));
#endif

  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*)p4est->user_pointer;
  d4est_operators_t* d4est_ops = curved_hp_amr_data->d4est_ops;
  d4est_element_data_t* parent_data = (d4est_element_data_t*) outgoing[0]->p.user_data;
  d4est_element_data_t* child_data;
  
  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;
      
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degH);
  
  /* to store h-prolongation data */
  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  
  d4est_operators_apply_hp_prolong
    (
     d4est_ops,
     &(parent_data->u_elem[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );

  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (d4est_element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    d4est_linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_elem[0], volume_nodes);
  }
    
  /* we're done with this refined parent, move on to next */
  P4EST_FREE(temp_data);
}

static int
curved_hp_amr_refine_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{
  curved_hp_amr_data_t* curved_hp_amr_data = (curved_hp_amr_data_t*) p4est->user_pointer;
  d4est_element_data_t* elem_data = (d4est_element_data_t*) quadrant->p.user_data;
  d4est_operators_t* d4est_ops = curved_hp_amr_data->d4est_ops;
  
  int* refinement_log = curved_hp_amr_data->refinement_log;
  int* refinement_log_stride = &curved_hp_amr_data->refinement_log_stride;

  if (elem_data->debug_flag == 1){
    printf("refinement_log[*refinement_log_stride] = %d\n",refinement_log[*refinement_log_stride]);
  }
  
  /* h-refine */
  if (refinement_log[*refinement_log_stride] < 0){
    (*refinement_log_stride)++;
    return 1;
  }

  /* p-refine, p-coarsen or do nothing */
  else {
    int degH = elem_data->deg;
    int degh = refinement_log[*refinement_log_stride];
    int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), degh);
    
    double* temp_data = P4EST_ALLOC(double, volume_nodes);
    
    if (elem_data->deg < degh){
      d4est_operators_apply_p_prolong
        (
         d4est_ops,
         &elem_data->u_elem[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data        
        );
      d4est_linalg_copy_1st_to_2nd(temp_data, &elem_data->u_elem[0], volume_nodes);
    }
    else if (elem_data->deg > degh){
      d4est_operators_apply_p_restrict
        (
         d4est_ops,
         &elem_data->u_elem[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data
        );
      d4est_linalg_copy_1st_to_2nd(temp_data, &elem_data->u_elem[0], volume_nodes);
    }     

    elem_data->deg = degh;
    P4EST_FREE(temp_data);
    (*refinement_log_stride)++;
    return 0;
  }  
}


/* void */
/* curved_hp_amr_save_to_vtk */
/* ( */
/*  p4est_t* p4est, */
/*  p4est_geometry_t* geom, */
/*  char* filename, */
/*  int save_local_est */
/* ) */
/* { */
/*   if (save_local_est == 1){ */
/*     /\* mpi_abort("NOT SUPPORTED YET"); *\/ */
/*     /\* double* eta_corner = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN)); *\/ */
/*     /\* d4est_element_data_store_local_estimator_in_corner_array *\/ */
/*     /\*   ( *\/ */
/*     /\*    p4est, *\/ */
/*     /\*    eta_corner *\/ */
/*     /\*   ); *\/ */

/*     /\* /\\* util_print_matrix( eta_corner, p4est->local_num_quadrants*(P4EST_CHILDREN), 1, "eta_corner = ", 0); *\\/ *\/ */
    
/*     /\* curved_hacked_p4est_vtk_write_all *\/ */
/*     /\*   (p4est, geom,     /\\* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer *\\/ *\/ */
/*     /\*    0.99,    /\\* draw each quadrant at almost full scale *\\/ *\/ */
/*     /\*    0,       /\\* do not write the tree id's of each quadrant (there is only one tree in this example) *\\/ *\/ */
/*     /\*    1,       /\\* do write the refinement level of each quadrant *\\/ *\/ */
/*     /\*    1,       /\\* do write the mpi process id of each quadrant *\\/ *\/ */
/*     /\*    0,       /\\* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) *\\/ *\/ */
/*     /\*    1,       /\\* write one scalar field: the solution value *\\/ *\/ */
/*     /\*    0,       /\\* write no vector fields *\\/ *\/ */
/*     /\*    filename, "eta2", eta_corner *\/ */
/*     /\*   ); *\/ */
  
/*     /\* P4EST_FREE(eta_corner); *\/ */
/*   } */

/*   else { */
/*     curved_hacked_p4est_vtk_write_all */
/*       (p4est, geom,     /\* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer *\/ */
/*        0.99,    /\* draw each quadrant at almost full scale *\/ */
/*        0,       /\* do not write the tree id's of each quadrant (there is only one tree in this example) *\/ */
/*        1,       /\* do write the refinement level of each quadrant *\/ */
/*        1,       /\* do write the mpi process id of each quadrant *\/ */
/*        0,       /\* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) *\/ */
/*        0,       /\* write one scalar field: the solution value *\/ */
/*        0,       /\* write no vector fields *\/ */
/*        filename */
/*       ); */
/*   }   */
/* } */

/** 
 * Assumes iter_volume will impose a refinement criterion
 * 
 * @param p4est 
 * @param data_to_hp_refine 
 * @param iter_volume 
 * @param curved_hp_amr_refine_store_fcn_ptr if NULL, we use default
 * @param curved_hp_amr_balance_store_fcn_ptr if NULL, we use default
 */
void
curved_hp_amr
(
 p4est_t* p4est,
 double** data_to_hp_refine,
 p4est_iter_volume_t iter_volume,
 p4est_replace_t refine_replace_callback_fcn_ptr, 
 p4est_replace_t balance_replace_callback_fcn_ptr, 
 void* curved_hp_amr_scheme_data,
 estimator_stats_t* stats,
 d4est_operators_t* d4est_ops
)
{
  curved_hp_amr_data_t curved_hp_amr_data;
  curved_hp_amr_data.refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
  curved_hp_amr_data.refinement_log_stride = 0;
  curved_hp_amr_data.data = *data_to_hp_refine;
  curved_hp_amr_data.curved_hp_amr_scheme_data = curved_hp_amr_scheme_data;
  curved_hp_amr_data.d4est_ops = d4est_ops;
  curved_hp_amr_data.estimator_stats = stats;
  /* set defaults if NULL, the extra data stored is only the deg */
  /* util_print_matrix(*data_to_hp_refine, d4est_element_data_get_local_nodes(p4est), 1, "data_to_hp_refine = ", 0); */
  
  if (refine_replace_callback_fcn_ptr == NULL){
    curved_hp_amr_data.refine_replace_callback_fcn_ptr
      = curved_hp_amr_replace_callback_default;
    printf("[DISCO4EST_WARNING]: Using default hp_amr callback for replacing quads after refining\n");
  }
  else{
    curved_hp_amr_data.refine_replace_callback_fcn_ptr
      = refine_replace_callback_fcn_ptr;
  }

  if (balance_replace_callback_fcn_ptr == NULL){
    curved_hp_amr_data.balance_replace_callback_fcn_ptr
      = curved_hp_amr_replace_callback_default;
    printf("[DISCO4EST_WARNING]: Using default hp_amr callback for replacing quads after balancing\n");
  }
  else{
    curved_hp_amr_data.balance_replace_callback_fcn_ptr
      = balance_replace_callback_fcn_ptr;
  }
  
  /* store vector */
  d4est_element_data_copy_from_vec_to_storage
    (
     p4est,
     *data_to_hp_refine
    );

  /* store p4est user pointer */
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = &curved_hp_amr_data;

  /* iterate and store refinement boolean */
  p4est_iterate(p4est,
                NULL,
                NULL,
                iter_volume,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);


  /* iter_volume will increment it */
  curved_hp_amr_data.refinement_log_stride = 0;
  
  /* int i; */
  /* for(i = 0; i < p4est->local_num_quadrants; i++){ */
  /* printf("refinement_log[%d] = %d\n", i, curved_hp_amr_data.refinement_log[i]); */
  /* } */
  /* interpolate solution and refine mesh */
  p4est_refine_ext
    (
     p4est,
     0,
     -1,
     curved_hp_amr_refine_callback,
     NULL,
     curved_hp_amr_data.refine_replace_callback_fcn_ptr
    );
    
  /* 2-1 balance mesh */
  p4est_balance_ext
    (
     p4est,
     P4EST_CONNECT_FACE,
     NULL,
     curved_hp_amr_data.balance_replace_callback_fcn_ptr
    );
  
  /* restore pointer back to original state */
  p4est->user_pointer = tmp;

  int new_nodes = d4est_element_data_get_local_nodes(p4est);
  /* realloc room for new refined mesh data */
  *data_to_hp_refine = P4EST_REALLOC(*data_to_hp_refine,
                                     double,
                                     new_nodes
                                    );
  
  d4est_element_data_copy_from_storage_to_vec
    (
     p4est,
     *data_to_hp_refine
    );

  
  /* int world_rank; */
  /* sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank); */
  /* if (world_rank == 0) */
  /*   util_print_matrix(*data_to_hp_refine, new_nodes, 1, "*data_to_hp_refine = ", 0); */
  P4EST_FREE(curved_hp_amr_data.refinement_log);
}

