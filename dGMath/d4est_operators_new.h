#ifndef D4EST_OPERATORS_H
#define D4EST_OPERATORS_H 

typedef struct {
  int max_degree;
  
  double ** lobatto_nodes_1d_table;
  double ** lobatto_weights_1d_table;
  double ** gauss_nodes_1d_table;
  double ** gauss_weights_1d_table;

  double ** mij_1d_table; /*  */
  double ** invmij_1d_table; /*  */
  double ** dij_1d_table; /*  */
  double ** lift_1d_table; /*  */
  double ** flip_1d_table; /*  */

  double ** lobatto_rst_2d_table;
  double ** lobatto_rst_3d_table;
  double ** gauss_rst_2d_table;
  double ** gauss_rst_3d_table;
  double ** vtk_rst_2d_table; /*  */
  double ** vtk_rst_3d_table;/*  */
  double ** vtk_interp_1d_table; /*  */
  
  double *** lobatto_to_gauss_interp_1d_table;
  double *** lobatto_to_gauss_interp_trans_1d_table;
  double *** hp_prolong_1d_table;
  double *** hp_prolong_transpose_1d_table;
  double *** hp_restrict_1d_table;
  double *** hp_restrict_interp_1d_table;
  double *** p_prolong_1d_table;
  double *** p_prolong_transpose_1d_table;
  double *** p_restrict_1d_table;
  double *** p_restrict_interp_1d_table;


} d4est_operators_t;


typedef void (*d4est_operators_hp_restrict_fcn_t)(d4est_operators_t *, double *, int *,
                                         int, int, double *);

typedef void (*d4est_operators_p_restrict_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int, double *);

typedef void (*d4est_operators_hp_prolong_fcn_t)(d4est_operators_t *, double *, int,
                                        int, int *, double *);

typedef void (*d4est_operators_p_prolong_fcn_t)(d4est_operators_t *, double *, int, int,
                                       int, double *);

typedef struct {

  double* r;
  double* s;
  double* t;
  
} d4est_rst_t;

typedef enum {QUADRATURE_GAUSS, QUADRATURE_LOBATTO} quadrature_type_t; 

#endif
