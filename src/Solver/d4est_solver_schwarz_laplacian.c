#include <d4est_solver_schwarz_laplacian.h>

void
d4est_solver_schwarz_laplacian_apply_aij_over_subdomain
(
 d4est_schwarz_metadata_t* schwarz_data,
 d4est_schwarz_operators_t* schwarz_ops,
 double* Au_over_subdomain,
 double* u_over_subdomain,
 int subdomain
)
{
  /* restrict transpose subdomain */

  for (int i = 0; i < sub_data->num_elements; i++){

    /* apply sitffness */
    
    for (int f = 0; f < (P4EST_FACES); f++){

      /* d4est laplacian flux boundary is the exact same */
      /* d4est interface is exact same except set all elements inside subdomain to non-ghost and set local nodes to
 restricted_transpose_nodes and set ed->nodal_stride to the restricted_transpose stride in the subdomain, then set all
the rest of the stuff to the same as the element for outside the element set u to big zero vector and all elements to non-ghost and zero vector otherwise rewrite laplacian boundary flux for them <---- I think this is a better idea*/
      
      /* get element data together or NULL if outside subdomain */
      /* check if boundary */
      /* send to interface or boundary routine */
      /* have different interface routine for outside */
      
    }

  }

}
 
