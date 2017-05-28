#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H 

/* #include "../Flux/compute_flux.h" */
#include "../Flux/curved_compute_flux.h"

typedef struct {

  int mpi_rank;
  
  /* total nodes for this CPU */
  int local_nodes;
  
  /* used to store Aij*uj (pointer alias)*/
  double* Au;

  /* primary node vector to 
     hold potential solution (pointer alias)*/
  double* u;
    
  /* auxiliary node vector to
   * hold previous iterate or
   * initial iterate for Newton
   * Raphson (pointer alias) 
   */
  double* u0;

  /* used to store rhs of weak eqns (pointer alias)*/
  double* rhs;
  
  /* function pointers to calculate flux for scalar variable */
  /* flux_fcn_ptrs_t scalar_flux_fcn_data; */
  curved_flux_fcn_ptrs_t curved_scalar_flux_fcn_data;

  /* function pointers to calculate flux for vector variable, will be deprecated eventually */
  /* flux_fcn_ptrs_t vector_flux_fcn_data; */
  curved_flux_fcn_ptrs_t curved_vector_flux_fcn_data;

  /* convenience pointer for the user */
  void* user;
  
} problem_data_t;

void
problem_data_copy_ptrs
(
 problem_data_t* pd1,
 problem_data_t* pd2
);

#endif
