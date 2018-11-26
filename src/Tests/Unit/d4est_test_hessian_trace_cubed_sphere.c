#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_brick.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_amr.h>
#include <d4est_laplacian.h>
#include <d4est_hessian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_solver_matrix_symmetry.h>
#include <d4est_util.h>
#include <d4est_norms.h>
#include <d4est_vtk.h>
#include <sc_reduce.h>
#include <limits.h>
#include <zlog.h>
#include <ini.h>

typedef struct {

  int use_sinx;
  int use_r2;
  int use_lorentzian;
  int test_dg_laplacian;

} d4est_test_hessian_trace_params_t;

static
int d4est_test_hessian_trace_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_test_hessian_trace_params_t* pconfig = (d4est_test_hessian_trace_params_t*)user;
  if (d4est_util_match_couple(section,"test_params",name,"test_dg_laplacian")) {
    D4EST_ASSERT(pconfig->test_dg_laplacian == -1);
    pconfig->test_dg_laplacian = atoi(value);
  }
  else if (d4est_util_match_couple(section,"test_params",name,"use_sinx")) {
    D4EST_ASSERT(pconfig->use_sinx == -1);
    pconfig->use_sinx = atoi(value);
  }
  else if (d4est_util_match_couple(section,"test_params",name,"use_r2")) {
    D4EST_ASSERT(pconfig->use_r2 == -1);
    pconfig->use_r2 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"test_params",name,"use_lorentzian")) {
    D4EST_ASSERT(pconfig->use_lorentzian == -1);
    pconfig->use_lorentzian = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
d4est_test_hessian_trace_params_t
d4est_test_hessian_trace_params_input
(
 const char* input_file
)
{
  d4est_test_hessian_trace_params_t input;
  input.test_dg_laplacian = -1;
  input.use_sinx = -1;
  input.use_r2 = -1;
  input.use_lorentzian = -1;
  if (ini_parse(input_file,
                d4est_test_hessian_trace_params_handler,
                &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }
  D4EST_CHECK_INPUT("test_params", input.use_sinx, -1);
  D4EST_CHECK_INPUT("test_params", input.use_r2, -1);
  D4EST_CHECK_INPUT("test_params", input.use_lorentzian, -1);
  D4EST_CHECK_INPUT("test_params", input.test_dg_laplacian, -1);

  int sum = input.use_sinx + input.use_r2 + input.use_lorentzian;
  D4EST_ASSERT(sum == 1);
  
  return input;
}


double
u_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  d4est_test_hessian_trace_params_t* params = user;
  if (params->use_sinx){
#if (P4EST_DIM)==3
    double pi = 3.1415926535897932384626433832795;
    return sin(pi*x)*sin(pi*y)*sin(pi*z);
#else
    double pi = 3.1415926535897932384626433832795;
    return sin(pi*x)*sin(pi*y);
#endif
  }
  if (params->use_r2){
#if (P4EST_DIM)==3
    return x*x + y*y + z*z;
#else
    return x*x + y*y;   
#endif
  }
  if (params->use_lorentzian){
#if (P4EST_DIM)==3
    double r2 = x*x + y*y + z*z;
#else
    double r2 = x*x + y*y;
#endif
    return 1./sqrt(1.+r2);
  }
}

double
laplacian_u_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
){
  d4est_test_hessian_trace_params_t* params = user;
  if (params->use_sinx){
    double PIE = 3.1415926535897932384626433832795;
    return -3*PIE*PIE*sin(PIE*x)*sin(PIE*y)*sin(PIE*z);
  }
  if (params->use_r2){
#if (P4EST_DIM)==3
    return 6.;
#else
    return 4.;
#endif
  }
  if (params->use_lorentzian){
#if (P4EST_DIM)==3
    return -3./pow((1. + x*x + y*y + z*z),2.5);
#else
    return (-2 + x*x + y*y)/pow(1 + x*x + y*y,2.5);
#endif
  }
}


int main(int argc, char *argv[])
{

#ifndef D4EST_TEST
  D4EST_ABORT("D4EST_TEST not defined");
#endif
  
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  
#ifndef NDEBUG
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE ON\n");
  p4est_init(NULL, SC_LP_ERROR);
  /* p4est_init(NULL, SC_LP_ALWAYS); */
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DEBUG MODE OFF\n");
  p4est_init(NULL, SC_LP_ERROR);
#endif
  
#if (P4EST_DIM)==3
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 3\n");
#else
  if(proc_rank == 0)
    printf("[D4EST_INFO]: DIM = 2\n");
#endif

  char* input_file = P4EST_ALLOC(char, 100);
  sprintf(input_file, "%s", (argc == 2) ? argv[1] : "d4est_test_hessian_trace_cubed_sphere.input");
  
  if (proc_rank == 0)
    printf("[D4EST_INFO]: options file = %s\n", input_file);
    
  zlog_category_t *c_geom = zlog_get_category("d4est_geometry");
  d4est_geometry_t* d4est_geom = d4est_geometry_new(proc_rank,
                                                    input_file,
                                                    "geometry",
                                                    c_geom);

  d4est_mesh_initial_extents_t* initial_grid_input = d4est_mesh_initial_extents_parse(input_file, d4est_geom);

  p4est_t* p4est;
  p4est = p4est_new_ext
          (
           mpicomm,
           d4est_geom->p4est_conn,
           initial_grid_input->min_quadrants,
           initial_grid_input->min_level,
           initial_grid_input->fill_uniform,
           sizeof(d4est_element_data_t),
           NULL,
           NULL
          );


  p4est_partition(p4est, 1, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  

     
  if (proc_rank == 0){
    printf("[D4EST_INFO]: mpisize = %d\n", proc_size);
  }
  if (proc_rank == 0 && initial_grid_input->load_from_checkpoint == 0){
    printf("[D4EST_INFO]: min_quadrants = %d\n", initial_grid_input->min_quadrants);
    printf("[D4EST_INFO]: min_level = %d\n", initial_grid_input->min_level);
    printf("[D4EST_INFO]: fill_uniform = %d\n", initial_grid_input->fill_uniform);
  }
  
  sc_MPI_Barrier(mpicomm);
  printf("[D4EST_INFO]: elements on proc %d = %d\n", proc_rank, p4est->local_num_quadrants);
  sc_MPI_Barrier(mpicomm);
  
  /* start just-in-time dg-math */
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* d4est_factors = d4est_mesh_data_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature");
  


  d4est_test_hessian_trace_params_t params = d4est_test_hessian_trace_params_input
                                             (input_file);
  
  d4est_ghost_t* d4est_ghost = NULL;
  
  d4est_mesh_local_sizes_t local_sizes = d4est_mesh_update
                                         (
                                          p4est,
                                          &d4est_ghost,
                                          d4est_ops,
                                          d4est_geom,
                                          d4est_quad,
                                          d4est_factors,
                                          initial_grid_input,
                                          INITIALIZE_GHOST,
                                          INITIALIZE_QUADRATURE_DATA,
                                          INITIALIZE_GEOMETRY_DATA,
                                          INITIALIZE_GEOMETRY_ALIASES,
                                          d4est_mesh_set_initial_extents,
                                          (void*)initial_grid_input
                                         );


  
  int local_nodes = local_sizes.local_nodes;    
  double* u = P4EST_ALLOC(double, local_nodes);
  d4est_mesh_init_field
    (
     p4est,
     u,
     u_fcn,
     d4est_ops, // unnecessary?
     d4est_geom, // unnecessary?
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     &params
    );

  double* f_quad = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
  double* f = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes);
  double* hessian_trace_u_ana = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
  double* hessian_trace_u_num = P4EST_ALLOC(double, d4est_factors->local_sizes.local_nodes_quad);
  d4est_mesh_init_field
    (
     p4est,
     f,
     laplacian_u_fcn,
     d4est_ops, // unnecessary?
     d4est_geom, // unnecessary?
     d4est_factors,
     INIT_FIELD_ON_LOBATTO,
     &params
    );

  for (int i = 0; i < p4est->local_num_quadrants; i++){

    d4est_element_data_t* ed = d4est_factors->element_data[i];
    
    d4est_quadrature_volume_t mesh_object;
    mesh_object.dq =  ed->dq;
    mesh_object.tree = ed->tree;
    mesh_object.element_id = ed->id;
        
    mesh_object.q[0] = ed->q[0];
    mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
    mesh_object.q[2] = ed->q[2];
#endif
    d4est_quadrature_interpolate
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &mesh_object,
       QUAD_OBJECT_VOLUME,
       QUAD_INTEGRAND_UNKNOWN,
       &f[ed->nodal_stride],
       ed->deg,
       &f_quad[ed->quad_stride],
       ed->deg_quad
    );

  }

  d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     HESSIAN_ANALYTICAL,
     u,
     hessian_trace_u_ana
    );

  d4est_hessian_compute_hessian_trace_of_field_on_quadrature_points
    (
     p4est,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors,
     HESSIAN_NUMERICAL,
     u,
     hessian_trace_u_num
    );

  int number_of_regions = d4est_geom->get_number_of_regions(d4est_geom);
  double* max_error_local = P4EST_ALLOC_ZERO(double, number_of_regions);
  double* max_error_global = P4EST_ALLOC_ZERO(double, number_of_regions);
  
  int num = 0;
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    d4est_element_data_t* ed = d4est_factors->element_data[i];
    int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
    for (int i = 0; i < volume_nodes_quad; i++){
      double f = f_quad[ed->quad_stride + i];
      double f_ana = hessian_trace_u_ana[ed->quad_stride + i];
      double f_num = hessian_trace_u_num[ed->quad_stride + i];
      /* printf("%.15f %.15f %.15f\n", f, f_ana, f_num); */
      double error = fabs(f - f_ana);
      max_error_local [ed->region] = (error > max_error_local[ed->region]) ? error : max_error_local[ed->region];
    }
  }

    /* for (int i = 0; i < number_of_regions; i++){ */
    /*   printf("quad, nodes, error_region_%d = %d %d %.15f\n", i, p4est->local_num_quadrants, local_nodes, max_error_local[i]); */
    /* } */
  
  sc_reduce
    (
     max_error_local,
     max_error_global,
     number_of_regions,
     sc_MPI_DOUBLE,
     sc_MPI_MAX,
     0,
     mpicomm
    );

  int global_nodes = 0;
  sc_reduce
    (
     &local_nodes,
     &global_nodes,
     1,
     sc_MPI_INT,
     sc_MPI_SUM,
     0,
     mpicomm
    );

  
  if (p4est->mpirank == 0){
    printf("quad, nodes, errors = %d %d ", p4est->global_num_quadrants, global_nodes);
    double max_error = max_error_global[0];
    for (int i = 0; i < number_of_regions; i++){
      max_error = (max_error < max_error_global[i]) ? max_error_global[i] : max_error;
      /* printf("quad, nodes, error_region_%d = %.15f\n", i,  max_error_global[i]); */
      printf("%.15f ", max_error_global[i]);
    }
    printf("\n", max_error);
  }

  P4EST_FREE(max_error_local);
  P4EST_FREE(max_error_global);
  
  if(d4est_quad->quad_type == QUAD_TYPE_GAUSS_LEGENDRE_LOBATTO){
    d4est_vtk_save
    (
     p4est,
     d4est_ops,
     input_file,
     "d4est_vtk",
     (const char*[]){"f_quad, u, hessian_trace_ana, hessian_trace_num", NULL},
     (double**)((const double*[]){f_quad, u, hessian_trace_u_ana, hessian_trace_u_num, NULL}),
     NULL,
     NULL,
     NULL, // (const char*[]){"element_id", NULL},
     NULL, //(int**)((const int*[]){element_id, NULL}),
     -1
    );
  }

  P4EST_FREE(f_quad);
  P4EST_FREE(f);
  P4EST_FREE(u);
  P4EST_FREE(hessian_trace_u_num);
  P4EST_FREE(hessian_trace_u_ana);


  /* APPLY DG LAPLACIAN */
  /* APPLY DG LAPLACIAN */
  /* APPLY DG LAPLACIAN */
  /* APPLY DG LAPLACIAN */
  /* APPLY DG LAPLACIAN */
  /* APPLY DG LAPLACIAN */

  int test_dg_laplacian = 0;
  if (test_dg_laplacian){
    
/*     double* f = P4EST_ALLOC(double, local_nodes); */
/*     double* Au = P4EST_ALLOC(double, local_nodes); */
/*     double* Au_compare = P4EST_ALLOC(double, local_nodes); */
/*     dirichlet_bndry_eval_method_t eval_method = EVAL_BNDRY_FCN_ON_LOBATTO; */

/*     /\* / Setup boundary conditions *\/ */
/*     d4est_laplacian_dirichlet_bc_t bc_data_for_lhs; */
/*     bc_data_for_lhs.dirichlet_fcn = u_fcn; */
/*     bc_data_for_lhs.eval_method = eval_method;   */

/*     d4est_laplacian_flux_data_t* flux_data_for_apply_lhs = d4est_laplacian_flux_new(p4est, input_file, BC_DIRICHLET, &bc_data_for_lhs); */


/*     d4est_elliptic_eqns_t prob_fcns; */
/*     prob_fcns.build_residual = NULL; */
/*     prob_fcns.apply_lhs = d4est_test_apply_lhs; */
/*     prob_fcns.user = flux_data_for_apply_lhs; */

/*     d4est_elliptic_data_t elliptic_data; */
/*     elliptic_data.u = u; */
/*     elliptic_data.Au = Au; */
/*     elliptic_data.local_nodes = local_nodes; */
/*     elliptic_data.field_types = &field_type; */
/*     elliptic_data.num_of_fields = 1; */
    
/*     d4est_elliptic_eqns_apply_lhs */
/*       ( */
/*        p4est, */
/*        d4est_ghost, */
/*        d4est_ghost_data, */
/*        &prob_fcns, */
/*        &elliptic_data, */
/*        d4est_ops, */
/*        d4est_geom, */
/*        d4est_quad, */
/*        d4est_factors */
/*       ); */

/*     d4est_mesh_init_field */
/*       ( */
/*        p4est, */
/*        f, */
/*        laplacian_poly_vec_fcn, */
/*        d4est_ops, // unnecessary? */
/*        d4est_geom, // unnecessary? */
/*        d4est_factors, */
/*        INIT_FIELD_ON_LOBATTO, */
/*        NULL */
/*       ); */
  
/*     for (int e = 0; e < p4est->local_num_quadrants; e++){ */
/*       d4est_element_data_t* ed = d4est_factors->element_data[e]; */
/*       d4est_quadrature_volume_t mesh_object; */
/*       mesh_object.dq = ed->dq; */
/*       mesh_object.tree = ed->tree; */
/*       mesh_object.element_id = ed->id; */
/*       mesh_object.q[0] = ed->q[0]; */
/*       mesh_object.q[1] = ed->q[1]; */
/* #if (P4EST_DIM)==3 */
/*       mesh_object.q[2] = ed->q[2]; */
/* #endif */
/*       element_id[ed->id] = ed->id; */
          
/*       double* J_quad = d4est_mesh_get_jacobian_on_quadrature_points(d4est_factors,ed); */

/*       d4est_quadrature_apply_mass_matrix */
/*         ( */
/*          d4est_ops, */
/*          d4est_geom, */
/*          d4est_quad, */
/*          &mesh_object, */
/*          QUAD_OBJECT_VOLUME, */
/*          QUAD_INTEGRAND_UNKNOWN, */
/*          &f[ed->nodal_stride], */
/*          ed->deg, */
/*          J_quad, */
/*          ed->deg_quad, */
/*          &Au_compare[ed->nodal_stride] */
/*         );      */
/*     } */
  }
  /* END APPLY DG LAPLACIAN */
  /* END APPLY DG LAPLACIAN */
  /* END APPLY DG LAPLACIAN */
  /* END APPLY DG LAPLACIAN */
  /* END APPLY DG LAPLACIAN */
  /* END APPLY DG LAPLACIAN */
    
  if (d4est_ghost != NULL)
    d4est_ghost_destroy(d4est_ghost);
  
  d4est_mesh_initial_extents_destroy(initial_grid_input);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);

  P4EST_FREE(input_file);
  PetscFinalize();
  return 0;
}
