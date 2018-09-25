#ifndef D4EST_SOLVER_SCHWARZ_H
#define D4EST_SOLVER_SCHWARZ_H 

#include <pXest.h>

typedef struct {

  int mpirank;
  int tree;
  int tree_quadid;
  int id;

  /* the faces that touch either a face, edge, or corner of the schwarz center
   * these are in the orientation of the tree that owns the connection element */
  int faces [3];

  /* gives local or ghost element that touches a face or -1 if it's a 2nd-layer ghost or no face, gives same id as element if boundary 
   * the id is between 0 and local_num_quadrants + ghost_num_quadrants - 1*/
  int elements_that_touch_face[(P4EST_FACES)][(P4EST_HALF)];
  
} d4est_solver_schwarz_subdomain_element_t;

typedef struct {

  int center_id;
  d4est_solver_schwarz_subdomain_element_t* subdomain_elements;
  int subdomain_size;
  
} d4est_solver_schwarz_subdomain_data_t;


typedef struct {

  int overlap_size;   /* Quantifies how much the schwarz subdomain 
                       * overlaps the central element. In LGL node units */
  int num_subdomains;
  d4est_solver_schwarz_subdomain_data_t* subdomain_data;

  d4est_solver_schwarz_weight_function_t weight_function_type;
  
} d4est_solver_schwarz_data_t;

#endif
