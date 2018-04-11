
#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_util.h>
#include <limits.h>

#define D4EST_TEST_EPS DBL_EPSILON
#define DEG_SUB_TEST_1 2
#define DEG_SUB_TEST_2 2
#define DEG_SUB_TEST_3 2
#define DEG_SUB_TEST_4 2

void
test_d4est_quadrature_compactified_weights_and_abscissas
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom
)
{
  char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas";
  printf("\nSTARTING %s sub-test\n", test_name);
  
  int deg = DEG_SUB_TEST_1;
  double* abscissas = P4EST_ALLOC(double, deg+1);
  double* weights = P4EST_ALLOC(double, deg+1);

  d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = cubed_sphere_attrs->R2;
  double R1 = cubed_sphere_attrs->R1;

 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);

        p4est_qcoord_t q [(P4EST_DIM)];
        q[0] = quad->x;
        q[1] = quad->y;
        q[2] = quad->z;
        p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level);
        
        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           d4est_quad,
           abscissas,
           weights,
           tt,
           q,
           dq,
           deg
          );

        double amin = q[0];
        double amax = q[0] + dq;
        double bmin = q[1];
        double bmax = q[1] + dq;
        double cmin = q[2];
        double cmax = q[2] + dq;

        /* transform element corners to [0,1]^3 topological space */
        amin /= (double)P4EST_ROOT_LEN;
        amax /= (double)P4EST_ROOT_LEN;
        bmin /= (double)P4EST_ROOT_LEN;
        bmax /= (double)P4EST_ROOT_LEN;
        cmin /= (double)P4EST_ROOT_LEN;
        cmax /= (double)P4EST_ROOT_LEN;

        /* transform element corners to [-1,1]^2 x [1,2] topological space */
        amin = 2.*amin - 1.;
        amax = 2.*amax - 1.;
        bmin = 2.*bmin - 1.;
        bmax = 2.*bmax - 1.;
        cmin = cmin + 1.;
        cmax = cmax + 1.;
        
        double a = (R2-R1)*(cmax-cmin);
        double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;

        
        double integral = 0.;
        for (int i = 0; i < deg + 1; i++){
          integral += weights[i]*abscissas[i]*abscissas[i]*(1./pow((a*abscissas[i]+b),4));
        }
        /* double integral_exact = (-2*(pow(a,2) + 3*pow(b,2)))/(3.*pow(a - b,3)*pow(a + b,3)); */
        /* double integral_exact = (8*a*b)/(3.*pow(a - b,3)*pow(a + b,3)); */
        double integral_exact = (-2*(3*pow(a,2) + pow(b,2)))/(3.*pow(a - b,3)*pow(a + b,3));

        if (fabs(integral_exact - integral) > D4EST_TEST_EPS){
          printf("fabs(integral_exact - integral) > D4EST_TEST_EPS\n");
          printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral));
          printf("D4EST_TEST_EPS = %.25f\n",(D4EST_TEST_EPS));
          printf("integral = %.25f\n",integral);
          printf("integral_exact = %.25f\n",integral_exact);
          printf("SUBTEST FAILED\n");
          P4EST_FREE(abscissas);
          P4EST_FREE(weights); 
          exit(1);
        }
        else {
          printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral));
          printf("SUBTEST PASSED\n");
        }

        /* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); */
     }
    }

 P4EST_FREE(abscissas);
 P4EST_FREE(weights); 
}


/* void */
/* test_d4est_quadrature_compactified_weights_and_abscissas */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_quadrature_t* d4est_quad, */
/*  d4est_geometry_t* d4est_geom */
/* ) */
/* { */
/*   char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas"; */
/*   printf("\nSTARTING %s sub-test\n", test_name); */
  
/*   int deg = DEG_SUB_TEST_1; */
/*   double* abscissas = P4EST_ALLOC(double, deg+1); */
/*   double* weights = P4EST_ALLOC(double, deg+1); */

/*   d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user); */
/*   double R2 = cubed_sphere_attrs->R2; */
/*   double R1 = cubed_sphere_attrs->R1; */

/*  for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */

/*         p4est_qcoord_t q [(P4EST_DIM)]; */
/*         q[0] = quad->x; */
/*         q[1] = quad->y; */
/*         q[2] = quad->z; */
/*         p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level); */
        
/*         d4est_quadrature_compactified_compute_abscissas_and_weights */
/*           ( */
/*            d4est_geom, */
/*            d4est_quad, */
/*            abscissas, */
/*            weights, */
/*            tt, */
/*            q, */
/*            dq, */
/*            deg */
/*           ); */

/*         double amin = q[0]; */
/*         double amax = q[0] + dq; */
/*         double bmin = q[1]; */
/*         double bmax = q[1] + dq; */
/*         double cmin = q[2]; */
/*         double cmax = q[2] + dq; */

/*         /\* transform element corners to [0,1]^3 topological space *\/ */
/*         amin /= (double)P4EST_ROOT_LEN; */
/*         amax /= (double)P4EST_ROOT_LEN; */
/*         bmin /= (double)P4EST_ROOT_LEN; */
/*         bmax /= (double)P4EST_ROOT_LEN; */
/*         cmin /= (double)P4EST_ROOT_LEN; */
/*         cmax /= (double)P4EST_ROOT_LEN; */

/*         /\* transform element corners to [-1,1]^2 x [1,2] topological space *\/ */
/*         amin = 2.*amin - 1.; */
/*         amax = 2.*amax - 1.; */
/*         bmin = 2.*bmin - 1.; */
/*         bmax = 2.*bmax - 1.; */
/*         cmin = cmin + 1.; */
/*         cmax = cmax + 1.; */
        
/*         double a = (R2-R1)*(cmax-cmin); */
/*         double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1; */

        
/*         double integral = 0.; */
/*         for (int i = 0; i < deg + 1; i++){ */
/*           integral += weights[i]*abscissas[i]*abscissas[i]*(1./pow((a*abscissas[i]+b),4)); */
/*         } */
/*         /\* double integral_exact = (-2*(pow(a,2) + 3*pow(b,2)))/(3.*pow(a - b,3)*pow(a + b,3)); *\/ */
/*         /\* double integral_exact = (8*a*b)/(3.*pow(a - b,3)*pow(a + b,3)); *\/ */
/*         double integral_exact = (-2*(3*pow(a,2) + pow(b,2)))/(3.*pow(a - b,3)*pow(a + b,3)); */

/*         if (fabs(integral_exact - integral) > D4EST_TEST_EPS){ */
/*           printf("fabs(integral_exact - integral) > D4EST_TEST_EPS\n"); */
/*           printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral)); */
/*           printf("D4EST_TEST_EPS = %.25f\n",(D4EST_TEST_EPS)); */
/*           printf("integral = %.25f\n",integral); */
/*           printf("integral_exact = %.25f\n",integral_exact); */
/*           printf("SUBTEST FAILED\n"); */
/*           P4EST_FREE(abscissas); */
/*           P4EST_FREE(weights);  */
/*           exit(1); */
/*         } */
/*         else { */
/*           printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral)); */
/*           printf("SUBTEST PASSED\n"); */
/*         } */

/*         /\* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); *\/ */
/*      } */
/*     } */

/*  P4EST_FREE(abscissas); */
/*  P4EST_FREE(weights);  */
/* } */

static double
poly_fcn
(double r,
 double s,
 double t,
 double a,
 double b,
 int deg
){
  return  (1./pow((a*t+b),4))*(pow(r, deg) + pow(s, deg) + pow(t, deg));
}

static double
poly_fcn_no1ot
(double r,
 double s,
 double t,
 double a,
 double b,
 int deg
){
  return (pow(r, deg) + pow(s, deg) + pow(t, deg));
}


void
test_d4est_quadrature_compactified_weights_and_abscissas_volume
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom
)
{
  char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas_volume";
  printf("\nSTARTING %s sub-test\n", test_name);
  
  int deg = DEG_SUB_TEST_2;
  double* abscissas = P4EST_ALLOC(double, deg+1);
  double* weights = P4EST_ALLOC(double, deg+1);

  d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = cubed_sphere_attrs->R2;
  double R1 = cubed_sphere_attrs->R1;

  
 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);

        p4est_qcoord_t q [(P4EST_DIM)];
        q[0] = quad->x;
        q[1] = quad->y;
        q[2] = quad->z;
        p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level);
        
        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           d4est_quad,
           abscissas,
           weights,
           tt,
           q,
           dq,
           deg
          );


        double amin = q[0];
        double amax = q[0] + dq;
        double bmin = q[1];
        double bmax = q[1] + dq;
        double cmin = q[2];
        double cmax = q[2] + dq;

        /* transform element corners to [0,1]^3 topological space */
        amin /= (double)P4EST_ROOT_LEN;
        amax /= (double)P4EST_ROOT_LEN;
        bmin /= (double)P4EST_ROOT_LEN;
        bmax /= (double)P4EST_ROOT_LEN;
        cmin /= (double)P4EST_ROOT_LEN;
        cmax /= (double)P4EST_ROOT_LEN;

        /* transform element corners to [-1,1]^2 x [1,2] topological space */
        amin = 2.*amin - 1.;
        amax = 2.*amax - 1.;
        bmin = 2.*bmin - 1.;
        bmax = 2.*bmax - 1.;
        cmin = cmin + 1.;
        cmax = cmax + 1.;
        
        double a = (R2-R1)*(cmax-cmin);
        double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;

        double* rst [(P4EST_DIM)];
        double volume_nodes_quad = (deg+1)*(deg+1)*(deg+1);
        rst[2] = P4EST_ALLOC(double, volume_nodes_quad);
        rst[1] = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg, 1);
        rst[0] = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM), deg, 0);
        d4est_quadrature_compactified_compute_rst_volume(abscissas, deg, rst[2], 2);

        double* weights_volume = P4EST_ALLOC(double, volume_nodes_quad);
        double* gauss_weights = d4est_operators_fetch_gauss_weights_1d(d4est_ops, deg);
        
        d4est_kron_AoBoC(weights,
                                gauss_weights,
                                gauss_weights,
                                weights_volume,
                                1, deg+1,
                                1, deg+1,
                                1, deg+1);
        
        double integral = 0.;
        for (int i = 0; i < volume_nodes_quad; i++){
          double poly = poly_fcn(rst[0][i], rst[1][i], rst[2][i], a, b, deg);
          integral += weights_volume[i]*poly;
        }
        /* double integral_exact = (-2*(pow(a,2) + 3*pow(b,2)))/(3.*pow(a - b,3)*pow(a + b,3)); */
        /* double integral_exact = (8*a*b)/(3.*pow(a - b,3)*pow(a + b,3)); */
        double integral_exact = (-8*(11*pow(a,2) + 9*pow(b,2)))/(9.*pow(a - b,3)*pow(a + b,3));

        if (fabs(integral_exact - integral) > D4EST_TEST_EPS){
          printf("fabs(integral_exact - integral) > D4EST_TEST_EPS\n");
          printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral));
          printf("D4EST_TEST_EPS = %.25f\n",(D4EST_TEST_EPS));
          printf("integral = %.25f\n",integral);
          printf("integral_exact = %.25f\n",integral_exact);
          printf("SUBTEST FAILED\n");
          P4EST_FREE(abscissas);
          P4EST_FREE(weights); 
          exit(1);
        }
        else {
          printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral));
          printf("SUBTEST PASSED\n");
        }

        /* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); */
     }
    }

 P4EST_FREE(abscissas);
 P4EST_FREE(weights); 
}

void
test_d4est_quadrature_compactified_weights_and_abscissas_face
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom
)
{
  char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas_face";
  printf("\nSTARTING %s sub-test\n", test_name);
  
  int deg = DEG_SUB_TEST_3;
  int deg_GL = deg + 5;
  double* abscissas = P4EST_ALLOC(double, deg+1);
  double* weights = P4EST_ALLOC(double, deg+1);

  d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = cubed_sphere_attrs->R2;
  double R1 = cubed_sphere_attrs->R1;

  
 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);

        p4est_qcoord_t q [(P4EST_DIM)];
        q[0] = quad->x;
        q[1] = quad->y;
        q[2] = quad->z;
        p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level);
        
        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           d4est_quad,
           abscissas,
           weights,
           tt,
           q,
           dq,
           deg
          );


        double amin = q[0];
        double amax = q[0] + dq;
        double bmin = q[1];
        double bmax = q[1] + dq;
        double cmin = q[2];
        double cmax = q[2] + dq;

        /* transform element corners to [0,1]^3 topological space */
        amin /= (double)P4EST_ROOT_LEN;
        amax /= (double)P4EST_ROOT_LEN;
        bmin /= (double)P4EST_ROOT_LEN;
        bmax /= (double)P4EST_ROOT_LEN;
        cmin /= (double)P4EST_ROOT_LEN;
        cmax /= (double)P4EST_ROOT_LEN;

        /* transform element corners to [-1,1]^2 x [1,2] topological space */
        amin = 2.*amin - 1.;
        amax = 2.*amax - 1.;
        bmin = 2.*bmin - 1.;
        bmax = 2.*bmax - 1.;
        cmin = cmin + 1.;
        cmax = cmax + 1.;
        
        double a = (R2-R1)*(cmax-cmin);
        double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;

        int face_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, deg);
        int face_nodes_quad_GL = d4est_lgl_get_nodes((P4EST_DIM)-1, deg_GL);
        double* sj_face = P4EST_ALLOC(double, face_nodes_quad);
        double* j_div_sj_face = P4EST_ALLOC(double, face_nodes_quad);
        double* weights_face = P4EST_ALLOC(double, face_nodes_quad);
        double* weights_face_GL = P4EST_ALLOC(double, face_nodes_quad_GL);

        double* rst_face [(P4EST_DIM)-1];
        for (int i = 0; i < (P4EST_DIM)-1; i++){
          rst_face[i] = P4EST_ALLOC(double, face_nodes_quad);
        }
        
        
        d4est_quadrature_compactified_compute_rst_face
          (
           abscissas,
           deg,
           rst_face
          );

        double face_integral [(P4EST_FACES)];
        double face_integral_GL [(P4EST_FACES)];
        for (int f = 0; f < (P4EST_FACES); f++){

          d4est_rst_t rst_quad_mortar;
          d4est_rst_t rst_quad_mortar_GL;
          d4est_geometry_face_info_t face_info;
          d4est_geometry_get_face_info(f, &face_info);
          
          if (f == 0 || f == 1 || f == 2 || f == 3){
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0);
            rst_quad_mortar.s = rst_face[1];
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;
          
            double* gauss_weights = d4est_operators_fetch_gauss_weights_1d(d4est_ops, deg);
            double* gauss_weights_GL = d4est_operators_fetch_gauss_weights_1d(d4est_ops, deg_GL);
        
            d4est_kron_AoB(weights,
                                  gauss_weights,
                                  weights_face,
                                  1, deg+1,
                                  1, deg+1);
          

            d4est_kron_AoB(gauss_weights_GL,
                                  gauss_weights_GL,
                                  weights_face_GL,
                                  1, deg_GL+1,
                                  1, deg_GL+1);
          
            
          }
          else{
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0); 
            rst_quad_mortar.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 1); 
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;

            double* gauss_weights = d4est_operators_fetch_gauss_weights_1d(d4est_ops, deg);
            double* gauss_weights_GL = d4est_operators_fetch_gauss_weights_1d(d4est_ops, deg_GL);
        
            d4est_kron_AoB(gauss_weights,
                                  gauss_weights,
                                  weights_face,
                                  1, deg+1,
                                  1, deg+1);
          
           

            d4est_kron_AoB(gauss_weights_GL,
                                  gauss_weights_GL,
                                  weights_face_GL,
                                  1, deg_GL+1,
                                  1, deg_GL+1);
            
          }
          
          d4est_mortars_compute_geometric_data_on_mortar_aux
            (
             d4est_ops,
             d4est_geom,
             tt,
             q,
             dq,
             &rst_quad_mortar,
             1,
             1,
             &deg,
             f,
             NULL,
             sj_face,
             NULL,
             NULL,
             j_div_sj_face,
             COMPUTE_NORMAL_USING_JACOBIAN
            );




          double rst_3d_comp [(P4EST_DIM)];
          double rst_3d_GL [(P4EST_DIM)];
          
        face_integral[f] = 0.;
        face_integral_GL[f] = 0.;
        for (int i = 0; i < face_nodes_quad; i++){
          rst_3d_comp[face_info.c] = face_info.sgn;
          rst_3d_comp[face_info.a] = rst_quad_mortar.r[i];
          rst_3d_comp[face_info.b] = rst_quad_mortar.s[i];

          double poly = poly_fcn(rst_3d_comp[0], rst_3d_comp[1], rst_3d_comp[2], a, b, deg);
          face_integral[f] += weights_face[i]*poly;

        }

        printf("a,b = %f,%f\n", a,b);

        for (int i = 0; i < face_nodes_quad_GL; i++){
          rst_3d_GL[face_info.c] = face_info.sgn;
          rst_3d_GL[face_info.a] = rst_quad_mortar_GL.r[i];
          rst_3d_GL[face_info.b] = rst_quad_mortar_GL.s[i];
          double poly_GL = poly_fcn(rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], a, b, deg);
          printf("poly_GL, weights_face_GL = %.25f, %.25f, %.25f, %.25f, %.25f\n",rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], poly_GL, weights_face_GL[i]);
          face_integral_GL[f] += weights_face_GL[i]*poly_GL;
        }
        
        }
        
        double face_integral_analytic [(P4EST_FACES)];
          
        face_integral_analytic[0] = (-4*(13*pow(a,2) + 15*pow(b,2)))/(9.*pow(pow(a,2) - pow(b,2),3));
        face_integral_analytic[1] = (-4*(13*pow(a,2) + 15*pow(b,2)))/(9.*pow(pow(a,2) - pow(b,2),3));
        face_integral_analytic[2] = (-4*(13*pow(a,2) + 15*pow(b,2)))/(9.*pow(a - b,3)*pow(a + b,3));
        face_integral_analytic[3] = (-4*(13*pow(a,2) + 15*pow(b,2)))/(9.*pow(a - b,3)*pow(a + b,3));
        face_integral_analytic[4] = 20/(3.*pow(a - b,4));
        face_integral_analytic[5] = 20/(3.*pow(a + b,4));

        for (int f = 0; f < (P4EST_FACES); f++){
          printf("face num integral, face_gl_integral, face_ana integral = %.16f, %.16f, %.16f\n", face_integral[f], face_integral_GL[f], face_integral_analytic[f]);
        }                                                            
        /* if (fabs(integral_exact - integral) > D4EST_TEST_EPS){ */
        /*   printf("fabs(integral_exact - integral) > D4EST_TEST_EPS\n"); */
        /*   printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral)); */
        /*   printf("D4EST_TEST_EPS = %.25f\n",(D4EST_TEST_EPS)); */
        /*   printf("integral = %.25f\n",integral); */
        /*   printf("integral_exact = %.25f\n",integral_exact); */
        /*   printf("SUBTEST FAILED\n"); */
        /*   P4EST_FREE(abscissas); */
        /*   P4EST_FREE(weights);  */
        /*   exit(1); */
        /* } */
        /* else { */
        /*   printf("fabs(integral_exact - integral) = %.25f\n",fabs(integral_exact - integral)); */
        /*   printf("SUBTEST PASSED\n"); */
        /* } */


        for (int i = 0; i < (P4EST_DIM)-1; i++){
          P4EST_FREE(rst_face[i]);
        }

        P4EST_FREE(sj_face);
        P4EST_FREE(j_div_sj_face);
        P4EST_FREE(weights_face);
        
        
        /* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); */
     }
    }

 P4EST_FREE(abscissas);
 P4EST_FREE(weights); 
}

void
test_d4est_quadrature_compactified_weights_and_abscissas_interp_face
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom
)
{
  char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas_face";
  printf("\nSTARTING %s sub-test\n", test_name);
  
  int deg = DEG_SUB_TEST_4;
  int deg_GL = DEG_SUB_TEST_4;
  int deg_GLL = DEG_SUB_TEST_4;
  double* abscissas = P4EST_ALLOC(double, deg+1);
  double* weights = P4EST_ALLOC(double, deg+1);

  d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = cubed_sphere_attrs->R2;
  double R1 = cubed_sphere_attrs->R1;

  
 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);

        p4est_qcoord_t q [(P4EST_DIM)];
        q[0] = quad->x;
        q[1] = quad->y;
        q[2] = quad->z;
        p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level);
        
        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           d4est_quad,
           abscissas,
           weights,
           tt,
           q,
           dq,
           deg
          );


        double amin = q[0];
        double amax = q[0] + dq;
        double bmin = q[1];
        double bmax = q[1] + dq;
        double cmin = q[2];
        double cmax = q[2] + dq;

        /* transform element corners to [0,1]^3 topological space */
        amin /= (double)P4EST_ROOT_LEN;
        amax /= (double)P4EST_ROOT_LEN;
        bmin /= (double)P4EST_ROOT_LEN;
        bmax /= (double)P4EST_ROOT_LEN;
        cmin /= (double)P4EST_ROOT_LEN;
        cmax /= (double)P4EST_ROOT_LEN;

        /* transform element corners to [-1,1]^2 x [1,2] topological space */
        amin = 2.*amin - 1.;
        amax = 2.*amax - 1.;
        bmin = 2.*bmin - 1.;
        bmax = 2.*bmax - 1.;
        cmin = cmin + 1.;
        cmax = cmax + 1.;
        
        double a = (R2-R1)*(cmax-cmin);
        double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;

        int face_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, deg);
        int face_nodes_quad_GL = d4est_lgl_get_nodes((P4EST_DIM)-1, deg_GL);
        double* poly_GL = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_GLL = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_comp = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_lobatto_to_comp = P4EST_ALLOC(double, face_nodes_quad);
        double* lobatto_to_comp_interp = P4EST_ALLOC(double, (deg+1)*(deg+1));

        double* rst_face [(P4EST_DIM)-1];
        for (int i = 0; i < (P4EST_DIM)-1; i++){
          rst_face[i] = P4EST_ALLOC(double, face_nodes_quad);
        }
        
        
        d4est_quadrature_compactified_compute_rst_face
          (
           abscissas,
           deg,
           rst_face
          );


        d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
          (
           d4est_ops,
           abscissas,
           lobatto_to_comp_interp,
           NULL,
           deg,
           deg
          );

        double* lobatto_to_gauss_interp = d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg, deg);
        

        for (int f = 0; f < (P4EST_FACES); f++){

          d4est_rst_t rst_quad_mortar;
          d4est_rst_t rst_quad_mortar_GL;
          d4est_rst_t rst_quad_mortar_GLL;
          d4est_geometry_face_info_t face_info;
          d4est_geometry_get_face_info(f, &face_info);
          
          if (f < 4){
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0);
            rst_quad_mortar.s = rst_face[1];
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;


            rst_quad_mortar_GLL.r = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 0);
            rst_quad_mortar_GLL.s = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 1);
            rst_quad_mortar_GLL.t = NULL;
                       
            
          }
          else{
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0); 
            rst_quad_mortar.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 1); 
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;


            rst_quad_mortar_GLL.r = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 0);
            rst_quad_mortar_GLL.s = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 1);
            rst_quad_mortar_GLL.t = NULL;
                        
          }
          
          double rst_3d_comp [(P4EST_DIM)];
          double rst_3d_GL [(P4EST_DIM)];
          double rst_3d_GLL [(P4EST_DIM)];
          for (int i = 0; i < face_nodes_quad; i++){
            rst_3d_comp[face_info.c] = face_info.sgn;
            rst_3d_comp[face_info.a] = rst_quad_mortar.r[i];
            rst_3d_comp[face_info.b] = rst_quad_mortar.s[i];
            poly_comp[i] = poly_fcn_no1ot(rst_3d_comp[0], rst_3d_comp[1], rst_3d_comp[2], a, b, deg);
          }

          for (int i = 0; i < face_nodes_quad_GL; i++){
            rst_3d_GL[face_info.c] = face_info.sgn;
            rst_3d_GL[face_info.a] = rst_quad_mortar_GL.r[i];
            rst_3d_GL[face_info.b] = rst_quad_mortar_GL.s[i];
            poly_GL[i] = poly_fcn_no1ot(rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], a, b, deg);
            printf(" r_GL, s_GL, t_GL, poly_GL = %f, %f, %f, %f\n", rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], poly_GL[i]);
          }

          for (int i = 0; i < face_nodes_quad_GL; i++){
            rst_3d_GLL[face_info.c] = face_info.sgn;
            rst_3d_GLL[face_info.a] = rst_quad_mortar_GLL.r[i];
            rst_3d_GLL[face_info.b] = rst_quad_mortar_GLL.s[i];
            poly_GLL[i] = poly_fcn_no1ot(rst_3d_GLL[0], rst_3d_GLL[1], rst_3d_GLL[2], a, b, deg);
          }

          if (f < 4){
            d4est_kron_A1A2x_nonsqr(poly_lobatto_to_comp,
                                           lobatto_to_comp_interp,
                                           lobatto_to_gauss_interp,
                                           poly_GLL,
                                           deg + 1,
                                           deg + 1,
                             deg + 1, deg + 1);
          }
          else {
            d4est_kron_A1A2x_nonsqr(poly_lobatto_to_comp, lobatto_to_gauss_interp, lobatto_to_gauss_interp, poly_GLL, deg + 1, deg + 1,
                             deg + 1, deg + 1);
          }

          DEBUG_PRINT_4ARR_DBL(poly_comp, poly_GL, poly_GLL, poly_lobatto_to_comp, face_nodes_quad_GL);
          
        
        }


        for (int i = 0; i < (P4EST_DIM)-1; i++){
          P4EST_FREE(rst_face[i]);
        }


        P4EST_FREE(poly_comp);
        P4EST_FREE(poly_GL);
        P4EST_FREE(poly_GLL);
        P4EST_FREE(poly_lobatto_to_comp);
        P4EST_FREE(lobatto_to_comp_interp);
        
        
        /* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); */
     }
    }

 P4EST_FREE(abscissas);
 P4EST_FREE(weights); 
}


void
test_d4est_quadrature_compactified_weights_and_abscissas_galerkin_integral_face
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quad,
 d4est_geometry_t* d4est_geom
)
{
  char test_name [] = "test_d4est_quadrature_compactified_weights_and_abscissas_face";
  printf("\nSTARTING %s sub-test\n", test_name);
  
  int deg = DEG_SUB_TEST_4;
  int deg_GL = DEG_SUB_TEST_4;
  int deg_GLL = DEG_SUB_TEST_4;
  double* abscissas = P4EST_ALLOC(double, deg+1);
  double* weights = P4EST_ALLOC(double, deg+1);

  d4est_geometry_cubed_sphere_attr_t* cubed_sphere_attrs = ((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = cubed_sphere_attrs->R2;
  double R1 = cubed_sphere_attrs->R1;

  
 for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);

        p4est_qcoord_t q [(P4EST_DIM)];
        q[0] = quad->x;
        q[1] = quad->y;
        q[2] = quad->z;
        p4est_qcoord_t dq = P4EST_QUADRANT_LEN(quad->level);
        
        d4est_quadrature_compactified_compute_abscissas_and_weights
          (
           d4est_geom,
           d4est_quad,
           abscissas,
           weights,
           tt,
           q,
           dq,
           deg
          );


        double amin = q[0];
        double amax = q[0] + dq;
        double bmin = q[1];
        double bmax = q[1] + dq;
        double cmin = q[2];
        double cmax = q[2] + dq;

        /* transform element corners to [0,1]^3 topological space */
        amin /= (double)P4EST_ROOT_LEN;
        amax /= (double)P4EST_ROOT_LEN;
        bmin /= (double)P4EST_ROOT_LEN;
        bmax /= (double)P4EST_ROOT_LEN;
        cmin /= (double)P4EST_ROOT_LEN;
        cmax /= (double)P4EST_ROOT_LEN;

        /* transform element corners to [-1,1]^2 x [1,2] topological space */
        amin = 2.*amin - 1.;
        amax = 2.*amax - 1.;
        bmin = 2.*bmin - 1.;
        bmax = 2.*bmax - 1.;
        cmin = cmin + 1.;
        cmax = cmax + 1.;
        
        double a = (R2-R1)*(cmax-cmin);
        double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;

        int face_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, deg);
        int face_nodes_quad_GL = d4est_lgl_get_nodes((P4EST_DIM)-1, deg_GL);
        double* poly_GL = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_GLL = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_comp = P4EST_ALLOC(double, face_nodes_quad);
        double* poly_lobatto_to_comp = P4EST_ALLOC(double, face_nodes_quad);
        double* lobatto_to_comp_interp = P4EST_ALLOC(double, (deg+1)*(deg+1));

        double* rst_face [(P4EST_DIM)-1];
        for (int i = 0; i < (P4EST_DIM)-1; i++){
          rst_face[i] = P4EST_ALLOC(double, face_nodes_quad);
        }
        
        
        d4est_quadrature_compactified_compute_rst_face
          (
           abscissas,
           deg,
           rst_face
          );


        d4est_quadrature_compactified_compute_interpolation_matrix_and_transpose
          (
           d4est_ops,
           abscissas,
           lobatto_to_comp_interp,
           NULL,
           deg,
           deg
          );

        double* lobatto_to_gauss_interp = d4est_operators_fetch_lobatto_to_gauss_interp_1d(d4est_ops, deg, deg);
        

        for (int f = 0; f < (P4EST_FACES); f++){

          d4est_rst_t rst_quad_mortar;
          d4est_rst_t rst_quad_mortar_GL;
          d4est_rst_t rst_quad_mortar_GLL;
          d4est_geometry_face_info_t face_info;
          d4est_geometry_get_face_info(f, &face_info);
          
          if (f < 4){
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0);
            rst_quad_mortar.s = rst_face[1];
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;


            rst_quad_mortar_GLL.r = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 0);
            rst_quad_mortar_GLL.s = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 1);
            rst_quad_mortar_GLL.t = NULL;
                       
            
          }
          else{
            rst_quad_mortar.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 0); 
            rst_quad_mortar.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg, 1); 
            rst_quad_mortar.t = NULL;

            rst_quad_mortar_GL.r = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 0);
            rst_quad_mortar_GL.s = d4est_operators_fetch_gauss_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GL, 1);
            rst_quad_mortar_GL.t = NULL;


            rst_quad_mortar_GLL.r = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 0);
            rst_quad_mortar_GLL.s = d4est_operators_fetch_lobatto_rst_nd(d4est_ops, (P4EST_DIM)-1, deg_GLL, 1);
            rst_quad_mortar_GLL.t = NULL;
                        
          }
          
          double rst_3d_comp [(P4EST_DIM)];
          double rst_3d_GL [(P4EST_DIM)];
          double rst_3d_GLL [(P4EST_DIM)];
          for (int i = 0; i < face_nodes_quad; i++){
            rst_3d_comp[face_info.c] = face_info.sgn;
            rst_3d_comp[face_info.a] = rst_quad_mortar.r[i];
            rst_3d_comp[face_info.b] = rst_quad_mortar.s[i];
            poly_comp[i] = poly_fcn_no1ot(rst_3d_comp[0], rst_3d_comp[1], rst_3d_comp[2], a, b, deg);
          }

          for (int i = 0; i < face_nodes_quad_GL; i++){
            rst_3d_GL[face_info.c] = face_info.sgn;
            rst_3d_GL[face_info.a] = rst_quad_mortar_GL.r[i];
            rst_3d_GL[face_info.b] = rst_quad_mortar_GL.s[i];
            poly_GL[i] = poly_fcn_no1ot(rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], a, b, deg);
            printf(" r_GL, s_GL, t_GL, poly_GL = %f, %f, %f, %f\n", rst_3d_GL[0], rst_3d_GL[1], rst_3d_GL[2], poly_GL[i]);
          }

          for (int i = 0; i < face_nodes_quad_GL; i++){
            rst_3d_GLL[face_info.c] = face_info.sgn;
            rst_3d_GLL[face_info.a] = rst_quad_mortar_GLL.r[i];
            rst_3d_GLL[face_info.b] = rst_quad_mortar_GLL.s[i];
            poly_GLL[i] = poly_fcn_no1ot(rst_3d_GLL[0], rst_3d_GLL[1], rst_3d_GLL[2], a, b, deg);
          }

          if (f < 4){
            d4est_kron_A1A2x_nonsqr(poly_lobatto_to_comp,
                                           lobatto_to_comp_interp,
                                           lobatto_to_gauss_interp,
                                           poly_GLL,
                                           deg + 1,
                                           deg + 1,
                             deg + 1, deg + 1);
          }
          else {
            d4est_kron_A1A2x_nonsqr(poly_lobatto_to_comp, lobatto_to_gauss_interp, lobatto_to_gauss_interp, poly_GLL, deg + 1, deg + 1,
                             deg + 1, deg + 1);
          }

          DEBUG_PRINT_4ARR_DBL(poly_comp, poly_GL, poly_GLL, poly_lobatto_to_comp, face_nodes_quad_GL);
          
        
        }


        for (int i = 0; i < (P4EST_DIM)-1; i++){
          P4EST_FREE(rst_face[i]);
        }


        P4EST_FREE(poly_comp);
        P4EST_FREE(poly_GL);
        P4EST_FREE(poly_GLL);
        P4EST_FREE(poly_lobatto_to_comp);
        P4EST_FREE(lobatto_to_comp_interp);
        
        
        /* DEBUG_PRINT_2ARR_DBL(abscissas, weights, deg+1); */
     }
    }

 P4EST_FREE(abscissas);
 P4EST_FREE(weights); 
}


p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     0,
     NULL,
     NULL
    );
}

int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_CUBED_SPHERE_OUTER_SHELL;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;

  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = 10;
  sphere_attrs->R2 = 1000000;
  sphere_attrs->compactify_outer_shell = 1;
  sphere_attrs->compactify_inner_shell = -1;
  d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geom, sphere_attrs);
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    0,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init();  
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE_COMPACTIFIED_C1PC2T_NEG4;
  d4est_quadrature_compactified_new(p4est, d4est_ops, d4est_geom, d4est_quad, "", "");

  
  test_d4est_quadrature_compactified_weights_and_abscissas(p4est, d4est_ops, d4est_quad, d4est_geom);
  test_d4est_quadrature_compactified_weights_and_abscissas_volume(p4est, d4est_ops, d4est_quad, d4est_geom);
  test_d4est_quadrature_compactified_weights_and_abscissas_face(p4est, d4est_ops, d4est_quad, d4est_geom);
  test_d4est_quadrature_compactified_weights_and_abscissas_interp_face(p4est, d4est_ops, d4est_quad, d4est_geom);

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_geometry_destroy(d4est_geom);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
