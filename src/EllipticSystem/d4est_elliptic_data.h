#ifndef D4EST_ELLIPTIC_DATA_H
#define D4EST_ELLIPTIC_DATA_H 

#include <d4est_field.h>

typedef struct {

  int mpirank;
    
  /* total nodes for this CPU */
  int local_nodes;

  int num_of_fields;
  d4est_field_type_t* field_types;
  
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
  
} d4est_elliptic_data_t;

void
d4est_elliptic_data_copy_ptrs
(
 d4est_elliptic_data_t* pd1,
 d4est_elliptic_data_t* pd2
);

#endif
