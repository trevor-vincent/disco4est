/**
 * @file   element_data.h
 * @author Trevor Vincent <tvincent@cita.utoronto.ca>
 * @date   Tue Nov 10 22:49:00 2015
 *
 * @brief  element data for mixed-variable 2nd order
 * dG elliptic problems on cubic domains.
 *
 *
 */

#ifndef ELEMENTDATA_H
#define ELEMENTDATA_H

#include "../GridFunctions/grid_functions.h"
#include "../Flux/ip_flux.h"
#include "../dGMath/dgmath.h"

typedef struct {

  /* identification */
  int stride;
  int id;
  int mpirank;
  
  /* front corner */
  double xyz_corner[3];
  double h;

  /* primary node vector alias */
  double *u_elem;

  /* primary trace vector alias */
  double *ustar_min_u;

  /* auxiliary node vector aliases */
  double *q_elem[(P4EST_DIM)];
  double *du_elem[(P4EST_DIM)];
  double *Au_elem;

  /* auxiliary trace vector alias */
  double *qstar_min_q[(P4EST_DIM)];

  /* jacobians */
  double surface_jacobian;
  double jacobian;

  /* aposteriori/apriori error indicator for hp_amr or h_amr */
  double local_estimator;
  double local_predictor;

  /* For Bi amr algorithm */
  /* double* J1 [(P4EST_DIM)]; */
  /* double* J2 [(P4EST_DIM)]; */
  /* int p_JN; /\* max{|pmax - pneighbour|} over all neighbours *\/ */
  /* int p_J; /\* max{|p - pneighbour|} over all neighbours *\/ */

  /* storage for MPI transfers */
  double u_storage[MAX_NODES];

#ifdef D4EST_DEBUG
  /* only for the central flux */
  double q_storage[(P4EST_DIM)][(MAX_NODES)];
#endif

  
  /* storage for varying coefficient problems */
  /* double mu_storage [MAX_NODES]; */

  /* nodal degree */
  int deg;

  /* degree for Gauss integration if needed for non-linear terms */
  int deg_integ;
  
} element_data_t;

typedef int (*slice_fcn_t)(double, double
#if (P4EST_DIM) == 3
                           ,
                           double
#endif // (P4EST_DIM)==3
                           );

typedef struct {

  /* slicing condition */
  slice_fcn_t scf;

  /* pointer to vector to slice */
  double *node_vec;

  /* slice and attributes */
  double *slice;
  double **xyz;
  int slice_nodes;

  /* convenience stride */
  int stride;

} slice_data_t;

void element_data_get_slice_data(p4est_t *p4est, slice_data_t *slice_data,
                                 dgmath_jit_dbase_t *dgmath_jit_dbase);

void element_data_get_slice_data_count(p4est_t *p4est,
                                       slice_data_t *slice_data);

void element_data_init_node_vec(p4est_t *p4est, double *nodal_vec,
                                grid_fcn_t init_fcn,
                                dgmath_jit_dbase_t *dgmath_jit_dgbase);

double element_data_compute_l2_norm_sqr(p4est_t *p4est, double *nodal_vec,
                                        dgmath_jit_dbase_t *dgmath_jit_dbase);

void element_data_init(p4est_t *p4est, int deg);

int element_data_get_local_nodes(p4est_t *p4est);

void
element_data_integrate_auv_andaddto
(
 p4est_t* p4est,
 double a,
 double* u,
 double* v,
 double* a_uv,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void element_data_integrate_au_andaddto(p4est_t *p4est, double a, double *u,
                                        double *a_u,
                                        dgmath_jit_dbase_t *dgmath_jit_dbase);

void element_data_integrate_fofuv_andaddto(
    p4est_t *p4est, double *u, double *v, double *out,
    /* int nonlin_deg, */
    grid_fcn_ext_t f, dgmath_jit_dbase_t *dgmath_jit_dbase);

void
element_data_print_node_vec
(
 p4est_t* p4est,
 double* nodal_vec,
 char* name_string,
 int mpi_rank,
 int save_to_file,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void element_data_init_quadid(p4est_t *p4est);

double
element_data_compute_l2_norm_sqr_no_local(p4est_t *p4est, double *nodal_vec,
                                          dgmath_jit_dbase_t *dgmath_jit_dbase);

void element_data_get_xyz(p4est_t *p4est, double *xyz[(P4EST_DIM)]);

double element_data_compute_H1_norm_sqr(p4est_t *p4est, double *nodal_vec);

double
element_data_compute_DG_norm_sqr
(
 p4est_t* p4est,
 double* nodal_vec,
 grid_fcn_t bndry_fcn,
 ip_flux_params_t* ip_flux_params,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_ghost_t* ghost,
 element_data_t* ghost_data
);

void element_data_store_local_estimator_in_corner_array(p4est_t *p4est,
                                                        double *est_corner);
void element_data_print_local_estimator(p4est_t *p4est);

void element_data_print(p4est_t *p4est);

void element_data_copy_from_vec_to_storage(
 p4est_t* p4est,
 double* vec
);

void element_data_copy_from_storage_to_vec
(
 p4est_t* p4est,
 double* vec
);


void
element_data_compute_f_of_uxyz
(
 p4est_t* p4est,
 double* f,
 double* u,
 grid_fcn_ext_t f_fcn,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
element_data_apply_Mij_on_vec
(
 p4est_t* p4est,
 double* u,
 double* Mu,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);


typedef double
(*boundary_integral_fcn_t)
(
 element_data_t*,
 int,
 void*,
 dgmath_jit_dbase_t*
);

typedef struct {

  boundary_integral_fcn_t boundary_integral_fcn;
  void* ctx;
  double* boundary_integral;

} boundary_integral_data_t;


double
element_data_compute_boundary_integral
(
 p4est_t* p4est,
 boundary_integral_fcn_t boundary_integral_fcn,
 void* ctx,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

void
element_data_store_nodal_vec_in_vertex_array
(
 p4est_t* p4est,
 double* nodal_vec,
 double* corner_vec
);

int
element_data_which_quadrant_of_root
(
 element_data_t* elem_data
);

void
element_data_apply_Mij_on_f_of_vec1_x_vec2
(
 p4est_t* p4est,
 double* vec1,
 double* vec2,
 double* M_f_of_vec1_x_vec2,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 grid_fcn_ext_t f_fcn,
 int proj_deltap 
);

void
element_data_apply_Mij_on_f_of_vec
(
 p4est_t* p4est,
 double* u,
 double* Mu,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 grid_fcn_ext_t f_fcn,
 int proj_deltap 
);

double
element_data_get_local_h_min
(
 p4est_t* p4est
);


element_data_t*
element_data_get_element_data
(
 p4est_t* p4est,
 int local_element_id
);


void
element_data_get_level_range
(
 p4est_t* p4est,
 int* min_level,
 int* max_level
);

double
element_data_compute_l2_norm_error_no_local
(
 p4est_t* p4est,
 double* vec,
 int nodes,
 grid_fcn_t analytical_solution,
 dgmath_jit_dbase_t* dgmath_jit_dbase
);

int element_data_get_local_matrix_nodes
(
 p4est_t* p4est
);

void element_data_init_ext
(
 p4est_t* p4est,
 int set_deg,
 int set_deg_integ
);

int
element_data_init_new
(
 p4est_t* p4est,
 void(*user_fcn)(element_data_t*,void*),
 /* element_data_user_fcn_t user_fcn, */
 void* user_ctx
);

#endif
