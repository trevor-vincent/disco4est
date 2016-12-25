#define SAFETY
#define NDEBUG

#include "../hpAMR/hp_amr.h"
#include "../GridFunctions/grid_functions.h"
#include "../Utilities/util.h"
#include "../LinearAlgebra/linalg.h"
#include "../ElementData/element_data.h"
#include <hacked_p4est_vtk.h>


void
hp_amr_save_to_vtk
(
 p4est_t* p4est,
 char* filename,
 int save_local_est
)
{
  if (save_local_est == 1){
    double* eta_corner = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN));
    element_data_store_local_estimator_in_corner_array
      (
       p4est,
       eta_corner
      );

    /* util_print_matrix( eta_corner, p4est->local_pfnum_quadrants*(P4EST_CHILDREN), 1, "eta_corner = ", 0); */
    
    hacked_p4est_vtk_write_all
      (p4est, NULL,     /* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer */
       0.99,    /* draw each quadrant at almost full scale */
       0,       /* do not write the tree id's of each quadrant (there is only one tree in this example) */
       1,       /* do write the refinement level of each quadrant */
       1,       /* do write the mpi process id of each quadrant */
       0,       /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
       1,       /* write one scalar field: the solution value */
       0,       /* write no vector fields */
       filename, "eta2", eta_corner
      );
  
    P4EST_FREE(eta_corner);
  }

  else {
    hacked_p4est_vtk_write_all
      (p4est, NULL,     /* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer */
       0.99,    /* draw each quadrant at almost full scale */
       0,       /* do not write the tree id's of each quadrant (there is only one tree in this example) */
       1,       /* do write the refinement level of each quadrant */
       1,       /* do write the mpi process id of each quadrant */
       0,       /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
       0,       /* write one scalar field: the solution value */
       0,       /* write no vector fields */
       filename
      );
  }  
}


/**
 * @file   hp_amr.c
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
hp_amr_replace_callback_default
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

  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*)p4est->user_pointer;
  dgmath_jit_dbase_t* dgmath_jit_dbase = hp_amr_data->dgmath_jit_dbase;
  element_data_t* parent_data = (element_data_t*) outgoing[0]->p.user_data;
  element_data_t* child_data;
  
  int degh [(P4EST_CHILDREN)];
  int degH = parent_data->deg;

  int i;
  for (i = 0; i < (P4EST_CHILDREN); i++)
    degh[i] = degH;
      
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), degH);
  
  /* to store h-prolongation data */
  double* temp_data = P4EST_ALLOC(double, volume_nodes*(P4EST_CHILDREN));
  
  dgmath_apply_hp_prolong
    (
     dgmath_jit_dbase,
     &(parent_data->u_storage[0]),
     degH,
     (P4EST_DIM),
     &degh[0],
     temp_data
    );

  /* printf("element */
  /* util_print_matrix( &(parent_data->u_storage[0]), volume_nodes, 1, "u_storage = ", 0); */
  /* util_print_matrix( temp_data, volume_nodes*(P4EST_CHILDREN), 1, "u_storage new = ", 0); */
  
  for (i = 0; i < (P4EST_CHILDREN); i++){
    child_data = (element_data_t*) incoming[i]->p.user_data;
    child_data->deg = parent_data->deg;
    linalg_copy_1st_to_2nd(&temp_data[volume_nodes*i], &child_data->u_storage[0], volume_nodes);
  }
    
  /* we're done with this refined parent, move on to next */
  P4EST_FREE(temp_data);
}

static int
hp_amr_refine_callback
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t * quadrant
)
{
  hp_amr_data_t* hp_amr_data = (hp_amr_data_t*) p4est->user_pointer;
  element_data_t* elem_data = (element_data_t*) quadrant->p.user_data;
  dgmath_jit_dbase_t* dgmath_jit_dbase = hp_amr_data->dgmath_jit_dbase;
  
  int* refinement_log = hp_amr_data->refinement_log;
  int* refinement_log_stride = &hp_amr_data->refinement_log_stride;
  
  /* h-refine */
  if (refinement_log[*refinement_log_stride] < 0){
    (*refinement_log_stride)++;
    hp_amr_data->elements_marked_for_hrefine += 1;
    /* printf("Element will be refined, rank=%d, id=%d, estimator=%.25f\n", p4est->mpirank, elem_data->id, elem_data->local_estimator); */
    return 1;
  }

  /* p-refine, p-coarsen or do nothing */
  else {
    int degH = elem_data->deg;
    int degh = refinement_log[*refinement_log_stride];
    int volume_nodes = dgmath_get_nodes((P4EST_DIM), degh);
    
    double* temp_data = P4EST_ALLOC(double, volume_nodes);
    
    if (elem_data->deg < degh){
      hp_amr_data->elements_marked_for_prefine += 1;
      dgmath_apply_p_prolong
        (
         dgmath_jit_dbase,
         &elem_data->u_storage[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data        
        );
      linalg_copy_1st_to_2nd(temp_data, &elem_data->u_storage[0], volume_nodes);
    }
    else if (elem_data->deg > degh){
      dgmath_apply_p_restrict
        (
         dgmath_jit_dbase,
         &elem_data->u_storage[0],
         degH,
         (P4EST_DIM),
         degh,
         temp_data
        );
      linalg_copy_1st_to_2nd(temp_data, &elem_data->u_storage[0], volume_nodes);
    }     

    elem_data->deg = degh;
    P4EST_FREE(temp_data);
    (*refinement_log_stride)++;
    return 0;
  }  
}

/* void */
/* hp_amr_save_to_vtk */
/* ( */
/*  p4est_t* p4est, */
/*  char* filename, */
/*  int save_local_est */
/* ) */
/* { */
/*   if (save_local_est == 1){ */
/*     double* eta_corner = P4EST_ALLOC(double, p4est->local_num_quadrants*(P4EST_CHILDREN)); */
/*     element_data_store_local_estimator_in_corner_array */
/*       ( */
/*        p4est, */
/*        eta_corner */
/*       ); */

/*     /\* util_print_matrix( eta_corner, p4est->local_pfnum_quadrants*(P4EST_CHILDREN), 1, "eta_corner = ", 0); *\/ */
    
/*     hacked_p4est_vtk_write_all */
/*       (p4est, NULL,     /\* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer *\/ */
/*        0.99,    /\* draw each quadrant at almost full scale *\/ */
/*        0,       /\* do not write the tree id's of each quadrant (there is only one tree in this example) *\/ */
/*        1,       /\* do write the refinement level of each quadrant *\/ */
/*        1,       /\* do write the mpi process id of each quadrant *\/ */
/*        0,       /\* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) *\/ */
/*        1,       /\* write one scalar field: the solution value *\/ */
/*        0,       /\* write no vector fields *\/ */
/*        filename, "eta2", eta_corner */
/*       ); */
  
/*     P4EST_FREE(eta_corner); */
/*   } */

/*   else { */
/*     hacked_p4est_vtk_write_all */
/*       (p4est, NULL,     /\* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer *\/ */
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
 * @param hp_amr_refine_store_fcn_ptr if NULL, we use default
 * @param hp_amr_balance_store_fcn_ptr if NULL, we use default
 */
void
hp_amr
(
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double** data_to_hp_refine,
 estimator_stats_t* stats,
 hp_amr_scheme_t* scheme
)
{
  hp_amr_data_t hp_amr_data;
  hp_amr_data.refinement_log = P4EST_ALLOC(int, p4est->local_num_quadrants);
  hp_amr_data.refinement_log_stride = 0;
  hp_amr_data.data = *data_to_hp_refine;
  hp_amr_data.hp_amr_scheme_data = scheme->hp_amr_scheme_data;
  hp_amr_data.dgmath_jit_dbase = dgmath_jit_dbase;
  hp_amr_data.estimator_stats = stats;
  hp_amr_data.elements_marked_for_hrefine = 0;
  hp_amr_data.elements_marked_for_prefine = 0;
  
  
  if (scheme->refine_replace_callback_fcn_ptr == NULL){
    hp_amr_data.refine_replace_callback_fcn_ptr
      = hp_amr_replace_callback_default;
    printf("[WARNING]: Using default hp amr refine replace callback. We suggest setting it explicitly to the default, to avoid this warning\n ");
  }
  else{
    hp_amr_data.refine_replace_callback_fcn_ptr
      = scheme->refine_replace_callback_fcn_ptr;
  }

  if (scheme->balance_replace_callback_fcn_ptr == NULL){
    hp_amr_data.balance_replace_callback_fcn_ptr
      = hp_amr_replace_callback_default;
    printf("[WARNING]: Using default hp amr balance replace callback. We suggest setting it explicitly to the default, to avoid this warning\n ");
  }
  else{
    hp_amr_data.balance_replace_callback_fcn_ptr
      = scheme->balance_replace_callback_fcn_ptr;
  }
  
  if(scheme->pre_refine_callback != NULL){
    scheme->pre_refine_callback(p4est, scheme->hp_amr_scheme_data);
  }

  
  /* iterate and store refinement boolean */
  element_data_copy_from_vec_to_storage
    (
     p4est,
     *data_to_hp_refine
    );

  /* store p4est user pointer */
  void* tmp = p4est->user_pointer;
  p4est->user_pointer = &hp_amr_data;
  
  p4est_iterate(p4est,
                NULL,
                NULL,
                scheme->iter_volume,
		NULL,
#if (P4EST_DIM)==3
                 NULL,       
#endif                
		NULL);

  /* printf("Elements before refine = %d\n", p4est->local_num_quadrants); */
  
  /* interpolate solution and refine mesh */
  p4est_refine_ext
    (
     p4est,
     0,
     -1,
     hp_amr_refine_callback,
     NULL,
     hp_amr_data.refine_replace_callback_fcn_ptr
    );

  /* printf("Elements marked for h-refine = %d\n", hp_amr_data.elements_marked_for_hrefine); */
  /* printf("Elements marked for p-refine = %d\n", hp_amr_data.elements_marked_for_prefine); */
  /* printf("Elements after refine = %d\n", p4est->local_num_quadrants); */
  
  /* 2-1 balance mesh */
  p4est_balance_ext
    (
     p4est,
     P4EST_CONNECT_FULL,
     NULL,
     hp_amr_data.balance_replace_callback_fcn_ptr
    );


  /* printf("Elements after balance = %d\n", p4est->local_num_quadrants); */
  
  /* restore pointer back to original state */
  p4est->user_pointer = tmp;

  int new_nodes = element_data_get_local_nodes(p4est);
  /* realloc room for new refined mesh data */
  *data_to_hp_refine = P4EST_REALLOC(*data_to_hp_refine,
                                     double,
                                     new_nodes
                                    );
  
  element_data_copy_from_storage_to_vec
    (
     p4est,
     *data_to_hp_refine
    );

  if(scheme->post_balance_callback != NULL){
    scheme->post_balance_callback(p4est, scheme->hp_amr_scheme_data);
  }
  
  /* int world_rank; */
  /* sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank); */
  /* if (world_rank == 0) */
  /*   util_print_matrix(*data_to_hp_refine, new_nodes, 1, "*data_to_hp_refine = ", 0); */
  P4EST_FREE(hp_amr_data.refinement_log);
}

