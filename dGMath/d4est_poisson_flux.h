#ifndef D4EST_POISSON_FLUX_H
#define D4EST_POISSON_FLUX_H 

#include <d4est_mortars.h>

typedef struct {

  d4est_quadrature_mortar_t* mortar_face_object;
  
  int total_side_nodes_m_lobatto;
  int total_side_nodes_p_lobatto;
  int total_nodes_mortar_lobatto;
  int total_nodes_mortar_quad;
  int faces_mortar;
  
  double* u_m_on_f_m_mortar_quad;
  double* u_m_on_f_m;
  double* u_p_on_f_p;
  double* sj_on_f_m_mortar_quad;
  double* j_div_sj_on_f_m_mortar_quad;
  double* u_p_on_f_p_mortar_quad;
  double* j_div_sj_on_f_p_mortar_quad;
  
  double* drst_dxyz_m_on_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_mortar_quad [(P4EST_DIM)];
  double* dudx_p_on_f_p_mortar_quad [(P4EST_DIM)];
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_mortar_quad [(P4EST_DIM)];

  int* face_nodes_m_lobatto;
  int* face_nodes_p_lobatto;
  int* deg_mortar_quad;
  int* nodes_mortar_quad;
  int* deg_mortar_lobatto;
  int* nodes_mortar_lobatto;
  int* deg_m_lobatto;
  int* deg_p_lobatto;

} d4est_poisson_flux_interface_data_t;


typedef struct {

  d4est_quadrature_mortar_t* face_object;

  int face_nodes_m_lobatto;
  int face_nodes_m_quad;

  double* xyz_on_f_m_lobatto [(P4EST_DIM)];
  double* u_m_on_f_m_quad;  
  double* u_m_on_f_m;  
  double* drst_dxyz_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* dudx_m_on_f_m_quad [(P4EST_DIM)];
  double* n_on_f_m_quad [(P4EST_DIM)];
  double* n_sj_on_f_m_quad [(P4EST_DIM)];
  double* sj_on_f_m_quad;
  double* j_div_sj_quad;
  
} d4est_poisson_flux_boundary_data_t;


typedef
void (*d4est_poisson_flux_interface_fcn_t)
(
 d4est_element_data_t**,
 int,
 int,
 int,
 d4est_element_data_t**,
 int,
 int,
 int,
 int*,
 int,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_poisson_flux_interface_data_t*,
 void*
);

typedef
void (*d4est_poisson_flux_boundary_fcn_t)
(
 d4est_element_data_t*,
 int,
 int,
 d4est_xyz_fcn_t,
 d4est_operators_t*,
 d4est_geometry_t*,
 d4est_quadrature_t*,
 d4est_poisson_flux_boundary_data_t*,
 void*
);

typedef enum {FLUX_SIPG, FLUX_NIPG, FLUX_IIPG, FLUX_NOT_SET} d4est_poisson_flux_type_t;
typedef struct d4est_poisson_flux_data d4est_poisson_flux_data_t;

struct d4est_poisson_flux_data{

  d4est_poisson_flux_type_t flux_type; 
  d4est_poisson_flux_interface_fcn_t interface_fcn;
  d4est_poisson_flux_boundary_fcn_t boundary_fcn;
  d4est_xyz_fcn_t boundary_condition;
  void* user;
  void (*destroy)(d4est_poisson_flux_data_t*);
  
};

/* This file was automatically generated.  Do not edit! */
void d4est_poisson_flux_destroy(d4est_poisson_flux_data_t *data);
d4est_poisson_flux_data_t *d4est_poisson_flux_new(p4est_t *p4est,const char *input_file,d4est_xyz_fcn_t boundary_condition,void *user);
d4est_mortar_fcn_ptrs_t d4est_poisson_flux_fetch_fcns(d4est_poisson_flux_data_t *data);
void d4est_poisson_flux_init_element_data(p4est_t *p4est,d4est_operators_t *d4est_ops,double *u,double *Au);


#endif
