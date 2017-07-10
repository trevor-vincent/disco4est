#include <sc_reduce.h>
#include <pXest.h>
#include <util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
#include <sipg_flux_vector_fcns.h>
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
#include <central_flux_params.h>
#include <curved_bi_estimator.h>
#include <krylov_petsc.h>
#include <matrix_sym_tester.h>
#include <dg_norm.h>
#include <hp_amr.h>
#include <hp_amr_curved_smooth_pred.h>
#include <hp_amr_curved_uniform.h>
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_disk.h>
#include <curved_poisson_debug_vecs.h>
#include <d4est_vtk.h>
#include <bi_estimator_flux_fcns.h>
#include <newton_petsc.h>
#include <ini.h>
#include <d4est_poisson.h>
#include <curved_Gauss_primal_sipg_kronbichler_flux_fcns.h>
#include <multigrid_matrix_operator.h>
#include <multigrid_smoother_cheby_d4est.h>
#include <multigrid_smoother_krylov_petsc.h>
#include <multigrid_bottom_solver_cg_d4est.h>
#include <multigrid_bottom_solver_cheby_d4est.h>
#include <multigrid_bottom_solver_krylov_petsc.h>
#include <multigrid_logger_residual.h>
#include <multigrid_element_data_updater_curved.h>
#include <krylov_pc_multigrid.h>
#include <krylov_pc_invstiff.h>
#include <ip_flux_params.h>
#include <ip_flux.h>
#include <d4est_cg.h>
#include "time.h"
#include "util.h"

d4est_geometry_cubed_sphere_attr_t global_cubed_sphere_attrs;
d4est_poisson_flux_sipg_params_t global_ip_flux_params;

typedef struct {

  double* u;
  double* u_analytic;
  double* error;
  double* jacobian;

} vtk_nodal_vecs_t;


typedef struct {

  double* u;
  double* Au;
  double* Au_error;
  double* jacobian;

} vtk_nodal_vecs_debug_t;

void
vtk_field_plotter
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
  vtk_nodal_vecs_t* vecs = user;
  vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx,
                                           4,
                                           0,
                                           "u",
                                           vecs->u,
                                           "u_analytic",
                                           vecs->u_analytic,
                                           "error",
                                           vecs->error,
                                           "jacobian",
                                           vecs->jacobian,
                                           vtk_ctx
                                          );


   vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 0,
                 0,
                 vtk_ctx
                );
  
}


void
vtk_field_plotter_debug
(
 d4est_vtk_context_t* vtk_ctx,
 void* user
)
{
  vtk_nodal_vecs_debug_t* vecs = user;
  vtk_ctx = d4est_vtk_write_dg_point_dataf(vtk_ctx,
                                           4,
                                           0,
                                           "u",
                                           vecs->u,
                                           "Au",
                                           vecs->Au,
                                           "error",
                                           vecs->Au_error,
                                           "jacobian",
                                           vecs->jacobian,
                                           vtk_ctx
                                          );


   vtk_ctx = d4est_vtk_write_dg_cell_dataf
                (
                 vtk_ctx,
                 1,
                 1,
                 1,
                 0,
                 1,
                 0,
                 0,
                 vtk_ctx
                );
  
}


static int
uni_refine_function
(
 p4est_t * p4est,
 p4est_topidx_t which_tree,
 p4est_quadrant_t *quadrant
)
{
  return 1;
}

typedef struct {

  int num_unifrefs;
  int deg_R0;
  int deg_quad_R0;
  int deg_stiffness_R0;
  int deg_R1;
  int deg_quad_R1;
  int deg_stiffness_R1;
  int deg_R2;
  int deg_quad_R2;
  int deg_stiffness_R2;
  char rhs_compute_method [50];
  
} problem_input_t;

static
int problem_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  problem_input_t* pconfig = (problem_input_t*)user;
  if (util_match_couple(section,"problem",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"rhs_compute_method")) {
    D4EST_ASSERT(pconfig->rhs_compute_method[0] == '*');
    snprintf (pconfig->rhs_compute_method, sizeof(pconfig->rhs_compute_method), "%s", value);
    D4EST_ASSERT(util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_GAUSS") ||
               util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_LOBATTO") );
  }
  else if (util_match_couple(section,"problem",name,"deg_R0")) {
    D4EST_ASSERT(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R0")) {
    D4EST_ASSERT(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_stiffness_R0")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R0 == -1);
    pconfig->deg_stiffness_R0 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_R1")) {
    D4EST_ASSERT(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R1")) {
    D4EST_ASSERT(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_stiffness_R1")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R1 == -1);
    pconfig->deg_stiffness_R1 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_R2")) {
    D4EST_ASSERT(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_quad_R2")) {
    D4EST_ASSERT(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
  }
  else if (util_match_couple(section,"problem",name,"deg_stiffness_R2")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R2 == -1);
    pconfig->deg_stiffness_R2 = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static
problem_input_t
problem_input
(
 const char* input_file
)
{
  
  problem_input_t input;
  input.num_unifrefs = -1;
  input.deg_R0 = -1;
  input.deg_quad_R0 = -1;
  input.deg_stiffness_R0 = -1;
  input.deg_R1 = -1;
  input.deg_quad_R1 = -1;
  input.deg_stiffness_R1 = -1;
  input.deg_R2 = -1;
  input.deg_quad_R2 = -1; 
  input.deg_stiffness_R2 = -1; 
  input.rhs_compute_method[0] = '*'; 

  if (ini_parse(input_file, problem_input_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.num_unifrefs, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R0, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R1, -1);
  D4EST_CHECK_INPUT("problem", input.deg_R2, -1);
  D4EST_CHECK_INPUT("problem", input.deg_quad_R2, -1);
  D4EST_CHECK_INPUT("problem", input.deg_stiffness_R2, -1);  
  D4EST_CHECK_INPUT("problem", input.rhs_compute_method[0], '*');  

  printf("[PROBLEM]: test_uniform_cubedsphere\n");
  printf("[PROBLEM]: rhs_compute_method = %s\n", input.rhs_compute_method);
  
  
  return input;
}



/* d4est_rst_t */
/* problem_get_rst_points_custom */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_element_data_t* ed, */
/*  int deg, */
/*  int dim, */
/*  void* user */
/* ){ */
/*   d4est_rst_t rst; */
/*   double* rsttmp [3] = {NULL}; */


/*   for (int d = 0; d < 2; d++){ */
/*     rsttmp[d] = d4est_operators_fetch_gauss_rst_nd */
/*                 ( */
/*                  d4est_ops, */
/*                  dim, */
/*                  deg, */
/*                  d */
/*                 ); */
/*   } */

/*   rsttmp[2] =  */
  
/*   rst.r = rsttmp[0]; */
/*   rst.s = rsttmp[1]; */
/*   rst.t = rsttmp[2];   */
/*   return rst; */
/* } */

/* d4est_rst_t */
/* problem_get_weights_custom */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_element_data_t* ed, */
/*  int deg, */
/*  int dim, */
/*  void* user */
/* ){ */
/*   d4est_rst_t rst; */
/*   double* rsttmp [3] = {NULL}; */

/*   for (int d = 0; d < 2; d++){ */
/*     rsttmp[d] = d4est_operators_fetch_gauss_rst_nd */
/*                 ( */
/*                  d4est_ops, */
/*                  dim, */
/*                  deg, */
/*                  d */
/*                 ); */
/*   } */

/*   rsttmp[2] =  */
  
/*   rst.r = rsttmp[0]; */
/*   rst.s = rsttmp[1]; */
/*   rst.t = rsttmp[2];   */
/*   return rst; */
/* } */

void
problem_build_custom_weights_and_abscissas
(
 double* custom_weights,
 double* custom_abscissas,
 double* custom_interp
)
{


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
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}


void
problem_set_degrees
(
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg = input->deg_R2;
    elem_data->deg_quad = input->deg_quad_R2;
    elem_data->deg_stiffness = input->deg_stiffness_R2;
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_quad = input->deg_quad_R1;
    elem_data->deg_stiffness = input->deg_stiffness_R1;
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_quad = input->deg_quad_R0;
    elem_data->deg_stiffness = input->deg_stiffness_R0;
  } 
}




void
problem_init
(
 const char* input_file,
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_operators_t* d4est_ops,
 int proc_size,
 sc_MPI_Comm mpicomm
)
{
  problem_input_t input = problem_input(input_file);
  global_cubed_sphere_attrs = *((d4est_geometry_cubed_sphere_attr_t*)d4est_geom->user);
  double R2 = global_cubed_sphere_attrs.R2;
  double R1 = global_cubed_sphere_attrs.R1;

  
  int level;
  
  D4EST_ASSERT((P4EST_DIM) == 2 || (P4EST_DIM) == 3);
  int world_rank, world_size;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &world_rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &world_size);

  double* Au = P4EST_ALLOC_ZERO(double, 1);
  double* rhs = P4EST_ALLOC_ZERO(double, 1);
  double* u = P4EST_ALLOC_ZERO(double, 1);
  double* u_analytic = P4EST_ALLOC_ZERO(double, 1);
  double* error = P4EST_ALLOC_ZERO(double, 1);
  int local_nodes = 1;

  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);

  p4est_partition(p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
  
  d4est_elliptic_problem_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  ip_flux_t* ip_flux = ip_flux_dirichlet_new(p4est, "[IP_FLUX]", input_file, zero_fcn);
  /* ip_flux_t* ip_flux = ip_flux_dirichlet_new(p4est, "[IP_FLUX]", input_file, boundary_fcn); */
  prob_vecs.flux_fcn_data = ip_flux->mortar_fcn_ptrs;
  global_ip_flux_params = *(ip_flux->ip_flux_params);
  prob_vecs.user = &input;

  
  d4est_mesh_geometry_storage_t* geometric_factors = geometric_factors_init(p4est);


  for (level = 0; level < input.num_unifrefs; ++level){

    if (world_rank == 0)
      printf("[D4EST_INFO]: UNIFORM REFINEMENT LEVEL %d\n", level);
    
      p4est_refine_ext(p4est,
                       0,
                       -1,
                       uni_refine_function,
                       NULL,
                       NULL
                      );

      p4est_partition(p4est, 0, NULL);
      p4est_balance_ext
        (
         p4est,
         P4EST_CONNECT_FACE,
         NULL,
         NULL
        );

      p4est_ghost_destroy(ghost);
      P4EST_FREE(ghost_data);

      ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
      ghost_data = P4EST_ALLOC(d4est_element_data_t, ghost->ghosts.elem_count);




      d4est_element_data_init_new(p4est,
                                   geometric_factors,
                                   d4est_ops,
                                   d4est_geom,
                                   problem_set_degrees,
                                   (void*)&input, 1, 1);
    
      local_nodes = d4est_element_data_get_local_nodes(p4est);

      Au = P4EST_REALLOC(Au, double, local_nodes);
      u = P4EST_REALLOC(u, double, local_nodes);
      rhs = P4EST_REALLOC(rhs, double, local_nodes);
      prob_vecs.Au = Au;
      prob_vecs.u = u;
      prob_vecs.rhs = rhs;
      prob_vecs.local_nodes = local_nodes;
      
  }


  d4est_element_data_init_new(p4est,
                               geometric_factors,
                               d4est_ops,
                               d4est_geom,
                               problem_set_degrees,
                               (void*)&input,
                               1,
                               1);


    
    local_nodes = d4est_element_data_get_local_nodes(p4est);

    for (p4est_topidx_t tt = p4est->first_local_tree;
         tt <= p4est->last_local_tree;
         ++tt)
      {
        p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
        sc_array_t* tquadrants = &tree->quadrants;
        int QQ = (p4est_locidx_t) tquadrants->elem_count;
        for (int qq = 0; qq < QQ; ++qq) {
          p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq);
          d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);
          int volume_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM), ed->deg_quad);
          
          double amin = ed->q[0];
          double amax = ed->q[0] + ed->dq;
          double bmin = ed->q[1];
          double bmax = ed->q[1] + ed->dq;
          double cmin = ed->q[2];
          double cmax = ed->q[2] + ed->dq;

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

          long double c = ((cmax + cmin)/(cmax - cmin)) - ((4*R2 - 
     2*R1)/((R2 - R1)*(cmax - cmin)));

          printf("c = %Le\n",c);
          printf("R1 - cmin R1 + (-2 + cmin) R2 = %f\n", R1 - cmin*R1 + (-2 + cmin)*R2);
          printf("R1 - cmax R1 + (-2 + cmax) R2 = %f\n", R1 - cmax*R1 + (-2 + cmax)*R2);

          d4est_rst_t rst_points = d4est_operators_get_rst_points(d4est_ops,
                                                          ed->deg_quad,
                                                          (P4EST_DIM),
                                                          QUAD_GAUSS);


          

          printf("R1 = %.15f\n",R1);
          printf("R2 = %.15f\n",R2);
          printf("cmax = %.15f\n",cmax);
          printf("cmin = %.15f\n",cmin);


          double* gauss_weights = d4est_operators_fetch_gauss_weights_1d(d4est_ops, ed->deg_quad);
          double* gauss_nodes = d4est_operators_fetch_gauss_nodes_1d(d4est_ops, ed->deg_quad);
          long double* jacmap_abscissas = P4EST_ALLOC(long double, ed->deg_quad + 1);
          long double* jacmap_weights = P4EST_ALLOC(long double, ed->deg_quad + 1);
          double* jacmap_abscissas_dbl = P4EST_ALLOC(double, ed->deg_quad + 1);
          double* jacmap_weights_dbl = P4EST_ALLOC(double, ed->deg_quad + 1);
          double* jacmap_rst [(P4EST_DIM)];
          D4EST_ALLOC_DIM_VEC(jacmap_rst, volume_nodes_quad);


          d4est_geometry_cubed_sphere_outer_shell_block_get_custom_quadrature
            (
             d4est_geom,
             tt,
             ed->q,
             ed->dq,
             ed->deg_quad,
             jacmap_abscissas,
             jacmap_weights,
             1
            );

          int nodes =  ed->deg_quad + 1;
          for (int i = 0; i<nodes; i++){
            jacmap_abscissas_dbl[i] = (double)jacmap_abscissas[i];
            jacmap_weights_dbl[i] = (double)jacmap_weights[i];
            printf("%d jacmap abscissas = %.25Lf\n", i, jacmap_abscissas[i]);
            /* printf("%d jacmap abscissas = %e\n", i,  jacmap_abscissas_dbl[i]); */
            /* printf("%d jacmap weights = %.25Lf\n", i,  jacmap_weights[i]); */
            /* printf("%d jacmap weights = %e\n", i, jacmap_weights_dbl[i]); */
          }


          /* int nodes =  ed->deg_quad + 1; */
          for (int i = 0; i<nodes; i++){
            /* jacmap_abscissas_dbl[i] = (double)jacmap_abscissas[i]; */
            /* jacmap_weights_dbl[i] = (double)jacmap_weights[i]; */
            /* printf("%d jacmap abscissas = %.25Lf\n", i, jacmap_abscissas[i]); */
            /* printf("%d jacmap abscissas = %e\n", i,  jacmap_abscissas_dbl[i]); */
            printf("%d jacmap weights = %.25Lf\n", i,  jacmap_weights[i]);
            /* printf("%d jacmap weights = %e\n", i, jacmap_weights_dbl[i]); */
          }


          for (int i = 0; i<nodes; i++){
            /* jacmap_abscissas_dbl[i] = (double)jacmap_abscissas[i]; */
            /* jacmap_weights_dbl[i] = (double)jacmap_weights[i]; */
            /* printf("%d jacmap abscissas = %.25Lf\n", i, jacmap_abscissas[i]); */
            /* printf("%d jacmap abscissas = %e\n", i,  jacmap_abscissas_dbl[i]); */
            printf("%d gauss weights = %.25f\n", i,  gauss_weights[i]);
            /* printf("%d jacmap weights = %e\n", i, jacmap_weights_dbl[i]); */
          }

          

          double* eye = P4EST_ALLOC_ZERO(double, nodes);
          for (int i = 0; i < nodes; i++) eye[i] = 1.;
          d4est_linalg_kron_AoBoC(eye, eye, gauss_nodes, jacmap_rst[0], nodes, 1, nodes, 1, nodes,
                            1);
          d4est_linalg_kron_AoBoC(eye, gauss_nodes, eye, jacmap_rst[1], nodes, 1, nodes, 1,
                            nodes, 1);
          d4est_linalg_kron_AoBoC(jacmap_abscissas_dbl, eye, eye, jacmap_rst[2], nodes, 1,
                            nodes, 1, nodes, 1);


          double a = (R2-R1)*(cmax-cmin);
          double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1;
          
          double* jac = P4EST_ALLOC(double, volume_nodes_quad);
          double* jac_times_apbtto4 = P4EST_ALLOC(double, volume_nodes_quad);
          double* wgau_wgau_wjac_jac = P4EST_ALLOC(double, volume_nodes_quad);
          double* wgau_wgau_wgau_jac = P4EST_ALLOC(double, volume_nodes_quad);
          double rst [(P4EST_DIM)];
          
          for (int i = 0; i < volume_nodes_quad; i++){
            rst[0] = jacmap_rst[0][i];
            rst[1] = jacmap_rst[1][i];
            rst[2] = jacmap_rst[2][i];
            d4est_geom->JAC(d4est_geom, tt, ed->q, ed->dq, rst, &jac[i]);
          /*   jac_times_apbtto4[i] = jac[i]*powl(a*rst[2] + b,4); */
          /*   printf("r,s,t = %f,%f,%f jac = %f, jac*(at+b)^4 = %f\n", rst[0], rst[1], rst[2], jac[i], jac_times_apbtto4[i]); */
            printf("%d,r,s,t,jac,jac_gauss = %.16f, %.16f, %.16f, %.16f, %.16f\n", i, rst[0], rst[1], rst[2], jac[i], ed->J_quad[i]);
          }

    
          d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_x(jacmap_weights_dbl,
                                               gauss_weights,
                                               gauss_weights,
                                               jac,
                                               ed->deg_quad + 1,
                                               ed->deg_quad + 1,
                                               ed->deg_quad + 1,
                                               wgau_wgau_wjac_jac
                                              );

          d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_x(gauss_weights,
                                               gauss_weights,
                                               gauss_weights,
                                               ed->J_quad,
                                               ed->deg_quad + 1,
                                               ed->deg_quad + 1,
                                               ed->deg_quad + 1,
                                               wgau_wgau_wgau_jac
                                              );

          /* NEXT: figure out how to make points and weights for each of the faces */
          

          DEBUG_PRINT_2ARR_DBL(wgau_wgau_wjac_jac, wgau_wgau_wgau_jac, volume_nodes_quad);
          
          double jacmap_sum = 0.;
          double gauss_sum = 0.;
          for (int i = 0; i < volume_nodes_quad; i++){
            jacmap_sum += wgau_wgau_wjac_jac[i];
            gauss_sum += wgau_wgau_wgau_jac[i];
          }


          int face_nodes_quad = d4est_lgl_get_nodes((P4EST_DIM)-1, ed->deg_quad);
          d4est_rst_t rst_points_face
            = d4est_operators_get_rst_points(d4est_ops, ed->deg_quad, (P4EST_DIM) - 1, QUAD_GAUSS);

          double* s_points_inverse_map = P4EST_ALLOC(double, face_nodes_quad);
          
          d4est_linalg_kron_AoB(jacmap_abscissas_dbl,
                          eye,
                          s_points_inverse_map,
                          nodes,
                          1,
                          nodes,
                          1
                         );
          
          for (int f = 0; f < (P4EST_FACES); f++){

            d4est_geometry_face_info_t face_info;
            d4est_geometry_get_face_info(f, &face_info);

            printf("Face %d\n", f);
            printf("face_info.a = %d\n", face_info.a);
            printf("face_info.b = %d\n", face_info.b);
            printf("face_info.c = %d\n", face_info.c);
            printf("face_info.sgn = %f\n", face_info.sgn);
            
            for (int i = 0; i < face_nodes_quad; i++){
              /* x face, r coord is gauss points */
              if ( (f == 0 || f == 1) || (f == 2 || f == 3) ){
                rst[face_info.a] = rst_points_face.r[i];
                rst[face_info.b] = s_points_inverse_map[i]; /* compactified coord */
                rst[face_info.c] = face_info.sgn;
              }
              /* z face, only gauss points */
              else{
                rst[face_info.a] = rst_points.r[i];
                rst[face_info.b] = rst_points.s[i]; /* compactified coord */
                rst[face_info.c] = face_info.sgn;
              }
              printf("face r, s, t = %.15f, %15f, %15f\n", rst[0], rst[1], rst[2]); 
            }

            
          }
          

          printf("jacmap_sum = %.25f\n", jacmap_sum);
          printf("gauss_sum = %.25f\n", gauss_sum);
          printf("theoretical = %.25f\n", (1./6.)*(4./3.)*M_PI*(R2*R2*R2 - R1*R1*R1));

          P4EST_FREE(s_points_inverse_map);
          P4EST_FREE(eye);
          D4EST_FREE_DIM_VEC(jacmap_rst);
          P4EST_FREE(jacmap_abscissas);
          P4EST_FREE(jacmap_abscissas_dbl);
          P4EST_FREE(jacmap_weights);
          P4EST_FREE(jacmap_weights_dbl);
          P4EST_FREE(wgau_wgau_wjac_jac);
          P4EST_FREE(wgau_wgau_wgau_jac);
          
        }               
      }
    

       
          /* for (int i = 0; i < volume_nodes_quad; i++){ */

          /*   double fac = pow(R2*(-4 + cmax + cmin + cmax*rst_points.t[i] - cmin*rst_points.t[i]) - R1*(-2 + cmax */
          /*                                                                                 + cmin + cmax*rst_points.t[i] - cmin*rst_points.t[i]),4); */

          /*   double a = (R2-R1)*(cmax-cmin); */
          /*   double b = (R2-R1)*(cmax+cmin) - 4*R2 + 2*R1; */
            
          /*   printf(" i, x, y, z, t, J, J*(t+c)^4 J*fac= %d %.15f %.15f %.15f %.15f %.15f %.15f %f %f\n", */
          /*          i, */
          /*          ed->xyz_quad[0][i], */
          /*          ed->xyz_quad[1][i], */
          /*          ed->xyz_quad[2][i], */
          /*          rst_points.t[i], */
          /*          ed->J_quad[i], */
          /*          (double)(ed->J_quad[i]*powl(rst_points.t[i] + c, 4)), */
          /*          ed->J_quad[i]*fac, */
          /*          (double)(ed->J_quad[i]*powl(a*rst_points.t[i] + b, 4)) */
          /*         ); */
          /* } */

 
    
    /* clock_t begin = 0; */
    /* clock_t end = -1; */

    /* if (world_rank == 0){ */
    /*   begin = clock(); */
    /* } */

    
    /* double* jacobian_lgl = P4EST_ALLOC(double, local_nodes); */
    /* d4est_element_data_compute_jacobian_on_lgl_grid(p4est,d4est_geom, d4est_ops, jacobian_lgl); */


    /* DEBUG_PRINT_ARR_DBL(jacobian_lgl, local_nodes); */
    
    /* d4est_linalg_vec_fabs(error, local_nodes); /\* not needed, except for vtk output *\/ */
    /* vtk_nodal_vecs_t vtk_nodal_vecs; */
    /* vtk_nodal_vecs.u = u; */
    /* vtk_nodal_vecs.u_analytic = u_analytic; */
    /* vtk_nodal_vecs.error = error; */
    /* vtk_nodal_vecs.jacobian = jacobian_lgl; */

    /* int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants); */
    /* d4est_element_data_get_array_of_degrees(p4est, deg_array); */

    /* char save_as [500]; */
    /* sprintf(save_as, "%s_level_%d", "test_uniform_cubedsphere", input.num_unifrefs); */
    
    /* /\* vtk output *\/ */
    /* d4est_vtk_save_geometry_and_dg_fields */
    /*   ( */
    /*    save_as, */
    /*    p4est, */
    /*    d4est_ops, */
    /*    deg_array, */
    /*    input_file, */
    /*    "d4est_vtk_geometry", */
    /*    vtk_field_plotter, */
    /*    (void*)&vtk_nodal_vecs */
    /*   ); */

    /* P4EST_FREE(jacobian_lgl); */
    /* P4EST_FREE(deg_array); */

  geometric_factors_destroy(geometric_factors);
  ip_flux_dirichlet_destroy(ip_flux);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  P4EST_FREE(u_analytic);
  P4EST_FREE(error);
  P4EST_FREE(u);
}
