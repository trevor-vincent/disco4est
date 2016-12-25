/**
 * @file   compute_flux.c
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Fri Sep 25 01:53:31 2015
 * 
 * @brief  
 * 
 * 
 */

#include "../Utilities/util.h"
#include "../ElementData/element_data.h"
#include "../Flux/compute_flux.h"

void compute_flux_on_local_elements(p4est_iter_face_info_t * info, void *user_data){
  int i;
  int s_p, s_m;
 
  compute_flux_user_data_t* compute_flux_user_data = (compute_flux_user_data_t*) info->p4est->user_pointer;
  flux_fcn_ptrs_t* flux_fcn_ptrs = (flux_fcn_ptrs_t*) compute_flux_user_data->flux_fcn_ptrs;
  dgmath_jit_dbase_t * dgmath_jit_dbase = (dgmath_jit_dbase_t*) compute_flux_user_data->dgmath_jit_dbase;

 
  element_data_t *ghost_data = (element_data_t *) user_data;
  element_data_t * e_p [(P4EST_HALF)];
  element_data_t * e_m [(P4EST_HALF)];
  int e_m_is_ghost [(P4EST_HALF)];
  
  sc_array_t *sides = &(info->sides);

  /* check if it's an interface boundary, otherwise it's a physical boundary */
  if (sides->elem_count == 2){
    
    p4est_iter_face_side_t *side[2];
 
    side[0] = p4est_iter_fside_array_index_int(sides,0);
    side[1] = p4est_iter_fside_array_index_int(sides,1);

    /* iterate through both sides of the interface */
    for (s_m = 0; s_m < 2; s_m++){

      /* initialize everything to non-ghost */
      for (i = 0; i < (P4EST_HALF); i++) e_m_is_ghost[i] = 0;
      
      s_p = (s_m == 0) ? 1 : 0;    

      /* the minus side is hanging */
      if (side[s_m]->is_hanging) {

        for (i = 0; i < (P4EST_HALF); i++) {
          /* we skip ghosts, b.c those are handled by other processes */
	  if (!side[s_m]->is.hanging.is_ghost[i]) {  
	    e_m[i] = (element_data_t *) side[s_m]->is.hanging.quad[i]->p.user_data;
	    e_m_is_ghost[i] = 0;
	  }
	  else{
	    /* e_m[i] = NULL; */
            e_m[i] = (element_data_t *) &ghost_data[side[s_m]->is.hanging.quadid[i]];
            /* printf("ghost_id = %d, ghost e_m[i]->deg = %d\n", side[s_m]->is.hanging.quadid[i],  e_m[i]->deg); */
	    e_m_is_ghost[i] = 1;
	  }
        }	    
        
        /* The + side must be full if the - side is hanging */
        /* We allow the + side to be a ghost, b.c we just copy data for this case */
        if (side[s_p]->is.full.is_ghost)
          e_p[0] = (element_data_t *) &ghost_data[side[s_p]->is.full.quadid];
        else
          e_p[0] = (element_data_t *) side[s_p]->is.full.quad->p.user_data;

        int sum_ghost_array = util_sum_array_int(e_m_is_ghost, (P4EST_HALF));

        if(sum_ghost_array != 0 && sum_ghost_array != (P4EST_HALF)){
          int dd;
          printf("sum_ghost_array != 0 || sum_ghost_array != %d, it = %d\n", (P4EST_HALF), sum_ghost_array);
          for (dd = 0; dd < (P4EST_HALF); dd++){
            printf(" ghost_id = %d\n", side[s_m]->is.hanging.quadid[dd]); 
          }
          /* mpi_abort("sum_ghost_array != 0 || sum_ghost_array != (P4EST_HALF)"); */
        }
        
        /* unless every hanging face is a ghost, we calculate the flux */
        if (sum_ghost_array < (P4EST_HALF)){
          flux_fcn_ptrs->flux_interface_fcn
            (
             e_m,
             (P4EST_HALF),
             side[s_m]->face,
             e_p,
             1,
             side[s_p]->face,
             e_m_is_ghost,
             dgmath_jit_dbase,
             flux_fcn_ptrs->params
            );
        }
      }

      /* If the - side is full */
      else {
	/* We only do calculations if the side is not a ghost */
	if (!side[s_m]->is.full.is_ghost) {
	  e_m[0] = (element_data_t *) side[s_m]->is.full.quad->p.user_data;
	  //if + side is hanging
	  if (side[s_p]->is_hanging){
            for (i = 0; i < (P4EST_HALF); i++) {
	      if (side[s_p]->is.hanging.is_ghost[i])
		e_p[i] = (element_data_t *) &ghost_data[side[s_p]->is.hanging.quadid[i]];	      
	      else
		e_p[i] = (element_data_t *) side[s_p]->is.hanging.quad[i]->p.user_data;
            }
            flux_fcn_ptrs->flux_interface_fcn
              (
               e_m,
               1,
               side[s_m]->face,
               e_p,
               (P4EST_HALF),
               side[s_p]->face,
               e_m_is_ghost,
               dgmath_jit_dbase,
               flux_fcn_ptrs->params
              );    
            
	  }
	  
	  //if + side is full
	  else {
	    if (side[s_p]->is.full.is_ghost){
	      e_p[0] = (element_data_t *) &ghost_data[side[s_p]->is.full.quadid];
            }
	    else{
	      e_p[0] = (element_data_t *) side[s_p]->is.full.quad->p.user_data;
            }
            flux_fcn_ptrs->flux_interface_fcn
              (
               e_m,
               1,
               side[s_m]->face,
               e_p,
               1,
               side[s_p]->face,
               e_m_is_ghost,
               dgmath_jit_dbase,
               flux_fcn_ptrs->params
              );  
	  }
	}        
      }	
    }
  }
  
  /* Weak Enforcement of Boundary Conditions */
  else {
    p4est_iter_face_side_t *side;
    side = p4est_iter_fside_array_index_int(sides,0);
    e_m[0] = (element_data_t *) side->is.full.quad->p.user_data;
/* #ifndef NDEBUG */
    /* printf("e_m[0]->deg = %d\n", e_m[0]->deg); */
/* #endif */
    flux_fcn_ptrs->flux_boundary_fcn(e_m[0],
                                     side->face,
                                     flux_fcn_ptrs->bndry_fcn,
                                     dgmath_jit_dbase,
                                     flux_fcn_ptrs->params
                                    );
  }
  
}
