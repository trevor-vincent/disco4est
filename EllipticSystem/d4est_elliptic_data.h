#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H 


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
  
  /* convenience pointer for the user */
  void* user;
  
} d4est_elliptic_problem_data_t;

void
problem_data_copy_ptrs
(
 d4est_elliptic_problem_data_t* pd1,
 d4est_elliptic_problem_data_t* pd2
);

#endif
