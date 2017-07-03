#ifndef CURVED_ELEMENT_DATA_H
#define CURVED_ELEMENT_DATA_H 

#define MAX_DEGREE 20
#if (P4EST_DIM) == 3
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#else
#define MAX_NODES (MAX_DEGREE + 1) * (MAX_DEGREE + 1)
#endif

#include <d4est_xyz_functions.h>
#include <d4est_operators.h>
#include <d4est_quadrature.h>


typedef struct {

  /* identification */
  int id;
  int mpi_rank;
  
  int sqr_nodal_stride;
  int sqr_mortar_stride;
  int nodal_stride;
  int quad_stride;
  
  int tree;
  int tree_quadid;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;
  
  /* geometric factors */
  double* J_quad; /* Jacobian */
  double* xyz [(P4EST_DIM)]; /* points on lobatto grid */
  double* xyz_quad [(P4EST_DIM)]; /* points on quadrature grid */
  double* xyz_rst_quad[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  double* rst_xyz_quad[(P4EST_DIM)][(P4EST_DIM)]; /* inverse mapping derivatives */
 
  /* estimator variables */
  double local_estimator;
  double local_predictor;

  /* local nodal fields */
  double u_elem[MAX_NODES];   /* storage for MPI transfers */
  double dudr_elem[(P4EST_DIM)][MAX_NODES];   /* storage for MPI transfers */
  double* Au_elem;  /* alias for Au */
  
  int deg; /* nodal degree */
  int deg_quad; /* deg for quadrature */
    
} d4est_element_data_t;

/* This file was automatically generated.  Do not edit! */
void d4est_element_data_print_local_estimator(p4est_t *p4est);
void d4est_element_data_reorient_f_p_elements_to_f_m_order(d4est_element_data_t **e_p,int face_dim,int f_m,int f_p,int o,int faces_p,d4est_element_data_t *e_p_oriented[(P4EST_HALF)]);
void d4est_element_data_store_nodal_vec_in_vertex_array(p4est_t *p4est,double *nodal_vec,double *corner_vec);
void d4est_element_data_store_element_scalar_in_vertex_array(p4est_t *p4est,double *vertex_array,double(*get_local_scalar_fcn)(d4est_element_data_t *));
void d4est_element_data_copy_from_storage_to_vec(p4est_t *p4est,double *vec);
void d4est_element_data_copy_from_vec_to_storage(p4est_t *p4est,double *vec);


#endif
