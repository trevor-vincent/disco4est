#ifndef D4EST_OPERATORS_H
#define D4EST_OPERATORS_H 

typedef struct {
  int max_degree;
  
  double ** Lobatto_nodes_1d_table;
  double ** Lobatto_weights_1d_table;
  double ** Gauss_nodes_1d_table;
  double ** Gauss_weights_1d_table;

  double ** Mij_1d_table; /*  */
  double ** invMij_1d_table; /*  */
  double ** Dij_1d_table; /*  */
  double ** lift_1d_table; /*  */
  double ** flip_1d_table; /*  */

  double ** Lobatto_rst_2d_table;
  double ** Lobatto_rst_3d_table;
  double ** Gauss_rst_2d_table;
  double ** Gauss_rst_3d_table;
  double ** vtk_rst_2d_table; /*  */
  double ** vtk_rst_3d_table;/*  */
  double ** vtk_interp_1d_table; /*  */
  
  double *** GLL_to_GL_interp_1d_table;
  double *** GLL_to_GL_interp_trans_1d_table;
  double *** hp_prolong_1d_table;
  double *** p_prolong_1d_table;
  double *** hp_restrict_1d_table;
  double *** p_restrict_1d_table;
  double *** hp_prolong_transpose_1d_table;
  double *** p_prolong_transpose_1d_table;
  double *** hp_restrict_interp_1d_table;
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
  
} d4est_operators_rst_t;

typedef enum {QUADRATURE_GAUSS, QUADRATURE_LOBATTO} quadrature_type_t; 

#endif
