#ifndef D4EST_ELEMENT_DATA_H
#define D4EST_ELEMENT_DATA_H 

#include <pXest.h>
#include <d4est_xyz_functions.h>
#include <d4est_operators.h>
#include <d4est_field.h>

#ifdef D4EST_TEST
#define MAX_NODES 729
#endif

typedef struct {

  /* processor information for element */
  int id;
  int mpirank;
  int tree_quadid; /* local id of quadrant in tree */
  
  /* TODO: Convenience strides, will be removed 
   * and places in d4est_mesh_data eventually*/
  int sqr_nodal_stride;
  int nodal_stride;
  int quad_stride;

  /* topological information for element */
  int region;
  int tree;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;
  
  /* nodal dG information for element */
  int deg; /* nodal degree */
  int deg_quad; /* deg for quadrature */

  /* only used in tests */
#ifdef D4EST_TEST
  double test_vecs[3][MAX_NODES];
  /* double area [P4EST_FACES]; */
  /* double volume; */
  /* double j_div_sj_min_lobatto [P4EST_FACES]; */
  /* double j_div_sj_avg_lobatto [P4EST_FACES]; */
  /* double j_div_sj_max_lobatto [P4EST_FACES]; */
#endif
  
} d4est_element_data_t;

/* This file was automatically generated.  Do not edit! */
d4est_element_data_t *d4est_element_data_get_ptr(p4est_t *p4est,int tree,int tree_quadid);
int d4est_element_data_get_stride_for_field(d4est_element_data_t *ed,d4est_field_type_t type);
int d4est_element_data_get_size_of_field(d4est_element_data_t *ed,d4est_field_type_t type);
void d4est_element_data_reorient_f_p_elements_to_f_m_order(d4est_element_data_t **e_p,int face_dim,int f_m,int f_p,int o,int faces_p,d4est_element_data_t *e_p_oriented[(P4EST_HALF)]);
void d4est_element_data_store_nodal_vec_in_vertex_array(p4est_t *p4est,double *nodal_vec,double *corner_vec);
void d4est_element_data_store_element_scalar_in_vertex_array(p4est_t *p4est,double *vertex_array,double(*get_local_scalar_fcn)(d4est_element_data_t *));

#endif
