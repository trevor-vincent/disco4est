
#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_quadrature_compactified.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <util.h>
#include <limits.h>

#define D4EST_TEST_EPS DBL_EPSILON
#define DEG_SUB_TEST_1 2
#define DEG_SUB_TEST_2 2
#define DEG_SUB_TEST_3 2

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
        rst[1] = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), deg, 1);
        rst[0] = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), deg, 0);
        d4est_quadrature_compactified_compute_rst_volume(abscissas, deg, rst[2], 2);

        double* weights_volume = P4EST_ALLOC(double, volume_nodes_quad);
        double* GL_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg);
        
        d4est_linalg_kron_AoBoC(weights,
                                GL_weights,
                                GL_weights,
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
        rst[1] = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), deg, 1);
        rst[0] = d4est_operators_fetch_Gauss_xyz_nd(d4est_ops, (P4EST_DIM), deg, 0);
        d4est_quadrature_compactified_compute_rst_volume(abscissas, deg, rst[2], 2);

        double* weights_volume = P4EST_ALLOC(double, volume_nodes_quad);
        double* GL_weights = d4est_operators_fetch_GL_weights_1d(d4est_ops, deg);
        
        d4est_linalg_kron_AoBoC(weights,
                                GL_weights,
                                GL_weights,
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


        for (int f = 0; f < (P4EST_FACES); f++){
          
          
        }


        
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
  d4est_quad->quad_type = QUAD_GAUSS_LEGENDRE_COMPACTIFIED;
  d4est_quadrature_compactified_new(p4est, d4est_ops, d4est_geom, d4est_quad, "");

  
  test_d4est_quadrature_compactified_weights_and_abscissas(p4est, d4est_ops, d4est_quad, d4est_geom);
  test_d4est_quadrature_compactified_weights_and_abscissas_volume(p4est, d4est_ops, d4est_quad, d4est_geom);

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_geometry_destroy(d4est_geom);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  
  PetscFinalize();
  /* sc_finalize (); */
  return 0;
}
