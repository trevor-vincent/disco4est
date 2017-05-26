#ifndef CURVED_ELEMENT_DATA_H
#define CURVED_ELEMENT_DATA_H 

#include <ip_flux_params.h>
#include "../dGMath/d4est_operators.h"
#include "../GridFunctions/grid_functions.h"
#include <d4est_geometry.h>



typedef enum { DO_NOT_STORE_LOCALLY, STORE_LOCALLY } norm_storage_option_t;
typedef enum { DIAM_APPROX, NO_DIAM_APPROX, DIAM_APPROX_CUBE} diam_compute_option_t;

typedef struct {

  /* double* J; */

  double* xyz;
  double* xyz_quad;
  double* xyz_rst_quad;
  double* J_quad;
  double* rst_xyz_quad;

  /* double* custom_abscissas; */
  /* double* custom_weights; */
  /* double* custom_interp; */
  
} geometric_factors_t;


typedef struct {

  /* identification */
  int id;
  int mpi_rank;
  
  int sqr_nodal_stride;
  int sqr_trace_stride;
  int nodal_stride;
  int integ_stride;
  
  int tree;
  int tree_quadid;
  p4est_qcoord_t q [(P4EST_DIM)];
  p4est_qcoord_t dq;
  
  /* auxiliary node vector aliases */
  double *Au_elem;

  /* Storage aliases for trace space calculations */
  double *M_ustar_min_u_n [(P4EST_DIM)];
  double *M_qstar_min_q_dot_n;

  /* geometric factors */
  /* double* J; /\* Jacobian *\/ */
  double* J_quad; /* Jacobian */
  /* double* invM; */
  /* double* invMface [P4EST_FACES]; */
  double* xyz [(P4EST_DIM)]; /* collocation points on physical grid */
  double* xyz_quad [(P4EST_DIM)]; /* collocation points on physical grid */
  /* double* xyz_rst[(P4EST_DIM)][(P4EST_DIM)]; /\* mapping derivatives *\/ */
  double* xyz_rst_quad[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  double* xyz_rst_Lobatto_quad[(P4EST_DIM)][(P4EST_DIM)]; /* mapping derivatives */
  /* double* rst_xyz[(P4EST_DIM)][(P4EST_DIM)]; /\* inverse mapping derivatives *\/ */
  double* rst_xyz_quad[(P4EST_DIM)][(P4EST_DIM)]; /* inverse mapping derivatives */
  double diam; /* approximate value of element diameter*/

  double volume;
  double surface_area [(P4EST_FACES)];
  
  /* aposteriori/apriori error indicator for hp_amr or h_amr */
  double local_estimator;
  double local_predictor;
  
  /* storage for MPI transfers */
  double u_storage[MAX_NODES];
  /* double du_elem[(P4EST_DIM)][MAX_NODES]; */
  double dudr_elem[(P4EST_DIM)][MAX_NODES];
  /* double q_elem[(P4EST_DIM)][MAX_NODES]; */
  
  /* nodal degree */
  int deg;
  int deg_quad; /* for integration other than stiffness matrix */
  int deg_stiffness; /* for integration of the stiffness matrix */
  
#ifndef NDEBUG
  /* useful flag for debugging */
  int debug_flag;
  int on_bdry;
#endif
  
} curved_element_data_t;

typedef void
(*curved_element_data_user_fcn_t)
(
 void*,
 void*
);

typedef struct {
  
  int local_nodes;
  int local_sqr_nodes;
  int local_sqr_trace_nodes;
  int local_nodes_quad;
  /* int local_sqr_nodes_invM; */
  
} curved_element_data_local_sizes_t;



#endif
