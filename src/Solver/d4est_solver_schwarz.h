
#ifndef D4EST_SOLVER_SCHWARZ_H
#define D4EST_SOLVER_SCHWARZ_H 

#include <pXest.h>

typedef struct {

  int process_p;
  int id_p; /* either local id on process that holds m side or ghost id on process that holds m side otherwise -1 */
  int face_p;
  int tree_p;
  int quadid_p;
  int is_hanging;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;

  
} d4est_schwarz_plus_face_info_t;

typedef struct {

  int mpirank;
  int tree;
  int tree_quadid;
  int id;
  int deg;
  /* the faces of this subdomain element (w.r.t it's tree coordinates) that touch either a face, edge, or corner of the schwarz center, faces = {-1,-1,-1} if subdomain element is the core element  */
  int faces [3];

  /* these are the faces of the core element (w.r.t it's tree coordinates) that touch a face, edge or corner of this subdomain element
   * this is useful for applying the weighting function when we need to know if a subdomain element is left or right of the core, core_faces = {-1,-1,-1} if subdomain element is the core element */
  int core_faces [3];

  /* 1 if core element, 0 if side element */
  int is_core;
  
  /* gives local or ghost element that touches a face or -1 if it's a 2nd-layer ghost or no face, gives same id as element if boundary 
   * the id is between 0 and local_num_quadrants + ghost_num_quadrants - 1*/
  /* d4est_mortar_data_t elements_that_touch_face[(P4EST_FACES)][(P4EST_HALF)]; */

  int nodal_size; 
  int nodal_stride; /* stride into just the nodal field on this subdomain, it is not a stride into a local nodal field and is not a stride into a field over all subdomains, i.e. it is zero on the first node of the first element of the containing subdomain*/
  int restricted_nodal_size; 
  int restricted_nodal_stride; /* stride into just the restricted nodal field on this subdomain, it is not a stride into a local nodal field and is not a stride into a restricted field over all subdomains, i.e. it is zero on the first node of the first element of the containing subdomain*/
  
} d4est_solver_schwarz_element_metadata_t;

typedef struct {

  int mpirank;
  int subdomain_id;
  int core_id; /* id stride in element_data for core */
  d4est_solver_schwarz_element_metadata_t*  element_metadata;
  int core_deg;
  int num_elements;
  int restricted_nodal_size;
  int restricted_nodal_stride;
  int nodal_size;
  int nodal_stride;
  
} d4est_solver_schwarz_subdomain_metadata_t;

typedef struct {

  int num_nodes_overlap;   /* Quantifies how many 1-D nodes are in the overlap of the the schwarz subdomain 
                            * this is a number between 1 and (minimum mesh degree + 1) */

  int restricted_nodal_size; /* restricted nodal size of all subdomains combined */
  int nodal_size; /* nodal size of all subdomains combined */
  int num_subdomains; /* Equivalent to number of local quadrants */
  d4est_solver_schwarz_subdomain_metadata_t* subdomain_metadata; /* The elements and their connections in the subdomain */

  
} d4est_solver_schwarz_metadata_t;

/* This file was automatically generated.  Do not edit! */
void d4est_solver_schwarz_metadata_print(p4est_t *p4est,d4est_solver_schwarz_metadata_t *schwarz_data);
void d4est_solver_schwarz_metadata_destroy(d4est_solver_schwarz_metadata_t *schwarz_data);
d4est_solver_schwarz_metadata_t *d4est_solver_schwarz_metadata_init(p4est_t *p4est,d4est_ghost_t *d4est_ghost,const char *input_file);
void d4est_solver_schwarz_metadata_input(p4est_t *p4est,const char *input_file,d4est_solver_schwarz_metadata_t *schwarz_data);

#endif
