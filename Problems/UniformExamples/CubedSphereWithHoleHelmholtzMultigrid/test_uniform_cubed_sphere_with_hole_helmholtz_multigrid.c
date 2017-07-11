#include <sc_reduce.h>
#include <pXest.h>
#include <d4est_util.h>
#include <d4est_linalg.h>
#include <d4est_element_data.h>
/* #include <sipg_flux_vector_fcns.h> */
#include <problem.h>
#include <problem_data.h>
#include <d4est_elliptic_eqns.h>
/* #include <central_flux_params.h> */
/* #include <curved_bi_estimator.h> */
#include <krylov_petsc.h>
#include <d4est_cg.h>
/* #include <matrix_sym_tester.h> */
/* #include <dg_norm.h> */
/* #include <hp_amr.h> */
/* #include <hp_amr_curved_smooth_pred.h> */
/* #include <hp_amr_curved_uniform.h> */
#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_geometry_disk.h>
/* #include <curved_poisson_debug_vecs.h> */
#include <d4est_vtk.h>
/* #include <bi_estimator_flux_fcns.h> */
#include <newton_petsc.h>
#include <ini.h>
#include <d4est_poisson.h>
#include <curved_Gauss_primal_sipg_kronbichler_flux_fcns.h>
/* #include <multigrid_matrix_operator.h> */
/* #include <multigrid_smoother_cheby_d4est.h> */
/* #include <multigrid_smoother_krylov_petsc.h> */
/* #include <multigrid_bottom_solver_cg_d4est.h> */
/* #include <multigrid_bottom_solver_cheby_d4est.h> */
/* #include <multigrid_bottom_solver_krylov_petsc.h> */
/* #include <multigrid_logger_residual.h> */
/* #include <multigrid_element_data_updater_curved.h> */
/* #include <krylov_pc_multigrid.h> */
/* #include <krylov_pc_invstiff.h> */
#include <d4est_mesh.h>
#include <ip_flux_params.h>
#include <ip_flux.h>
#include "time.h"
#include "d4est_util.h"

/* #define NASTY_DEBUG */

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
  if (which_tree == 0)
    return 1;
  else
    return 0;
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
  int solve_with_multigrid;
  int use_mg_as_pc_for_ksp;
  
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
    if (d4est_util_match_couple(section,"problem",name,"solve_with_multigrid")) {
    D4EST_ASSERT(pconfig->solve_with_multigrid == -1);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    pconfig->solve_with_multigrid = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"use_mg_as_pc_for_ksp")) {
    D4EST_ASSERT(pconfig->use_mg_as_pc_for_ksp == -1);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    pconfig->use_mg_as_pc_for_ksp = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"num_unifrefs")) {
    D4EST_ASSERT(pconfig->num_unifrefs == -1);
    pconfig->num_unifrefs = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"rhs_compute_method")) {
    D4EST_ASSERT(pconfig->rhs_compute_method[0] == '*');
    snprintf (pconfig->rhs_compute_method, sizeof(pconfig->rhs_compute_method), "%s", value);
    D4EST_ASSERT(d4est_util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_QUADRATURE") ||
               d4est_util_match(pconfig->rhs_compute_method, "COMPUTE_RHS_ON_LOBATTO") );
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R0")) {
    D4EST_ASSERT(pconfig->deg_R0 == -1);
    pconfig->deg_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R0")) {
    D4EST_ASSERT(pconfig->deg_quad_R0 == -1);
    pconfig->deg_quad_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R0")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R0 == -1);
    pconfig->deg_stiffness_R0 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R1")) {
    D4EST_ASSERT(pconfig->deg_R1 == -1);
    pconfig->deg_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R1")) {
    D4EST_ASSERT(pconfig->deg_quad_R1 == -1);
    pconfig->deg_quad_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R1")) {
    D4EST_ASSERT(pconfig->deg_stiffness_R1 == -1);
    pconfig->deg_stiffness_R1 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_R2")) {
    D4EST_ASSERT(pconfig->deg_R2 == -1);
    pconfig->deg_R2 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_quad_R2")) {
    D4EST_ASSERT(pconfig->deg_quad_R2 == -1);
    pconfig->deg_quad_R2 = atoi(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"deg_stiffness_R2")) {
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
  input.solve_with_multigrid = -1;
  input.use_mg_as_pc_for_ksp = -1;
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
  D4EST_CHECK_INPUT("problem", input.solve_with_multigrid, -1);  
  D4EST_CHECK_INPUT("problem", input.use_mg_as_pc_for_ksp, -1);  
  
  printf("[PROBLEM]: test_uniform_cubedsphere\n");
  printf("[PROBLEM]: rhs_compute_method = %s\n", input.rhs_compute_method);
  
  
  return input;
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
    /* elem_data->deg_stiffness = input->deg_stiffness_R2; */
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg = input->deg_R1;
    elem_data->deg_quad = input->deg_quad_R1;
    /* elem_data->deg_stiffness = input->deg_stiffness_R1; */
  }
  /* center cube */
  else {
    elem_data->deg = input->deg_R0;
    elem_data->deg_quad = input->deg_quad_R0;
    /* elem_data->deg_stiffness = input->deg_stiffness_R0; */
  } 
}

void
set_deg_quad
(
 void* elem_data_tmp,
 void* user_ctx
)
{
  d4est_element_data_t* elem_data = elem_data_tmp;
  problem_input_t* input = user_ctx;
  /* outer shell */
  if (elem_data->tree < 6){
    elem_data->deg_quad = input->deg_quad_R2;
    /* elem_data->deg_stiffness = input->deg_stiffness_R2; */
  }
  /* inner shell */
  else if(elem_data->tree < 12){
    elem_data->deg_quad = input->deg_quad_R1;
    /* elem_data->deg_stiffness = input->deg_stiffness_R1; */
  }
  /* center cube */
  else {
    elem_data->deg_quad = input->deg_quad_R0;
    /* elem_data->deg_stiffness = input->deg_stiffness_R0; */
  } 
}



static
double helmholtz_fcn
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  /* return 2.*x*y*z*x*y; */
  return 0.;
}

static
double analytic_solution_fcn
(
 double x,
 double y,
 double z
)
{
  double R2 = global_cubed_sphere_attrs.R2;
  return 1./sqrt(x*x + y*y + z*z) - 1./R2;
}

static
double u_poly
(
 double x,
 double y,
 double z
)
{
  return x*y*z*x*y*z;
}


static
double boundary_fcn
(
 double x,
 double y,
 double z
)
{
  /* printf("x, y, z, boundary_fcn = %f, %f, %f, %f\n", x,y,z,analytic_solution_fcn(x,y,z)); */
  return analytic_solution_fcn(x,y,z);
}

static
double f_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,
 double z
#endif
)
{
  return 0. + helmholtz_fcn(x,y,z,0,NULL)*analytic_solution_fcn(x,y,z);
}

static
double f_fcn_ext
(
 double x,
 double y,
 double z,
 double u0,
 void* user
)
{
  return 0. + helmholtz_fcn(x,y,z,0,NULL)*analytic_solution_fcn(x,y,z);
}

static
void problem_apply_lhs
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  /* prob_vecs->flux_fcn_data = curved_gauss_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns */
                                             /* (zero_fcn, &global_ip_flux_params); */

  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad);

 /* double* M_helmf_u = P4EST_ALLOC(double, prob_vecs->local_nodes); */
 

 /* int nodal_stride = 0; */
 /* for (p4est_topidx_t tt = p4est->first_local_tree; */
 /*      tt <= p4est->last_local_tree; */
 /*      ++tt) */
 /*   { */
 /*     p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
 /*     sc_array_t* tquadrants = &tree->quadrants; */
 /*     int Q = (p4est_locidx_t) tquadrants->elem_count; */
 /*     for (int q = 0; q < Q; ++q) { */
 /*       p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
 /*       d4est_element_data_t* ed = quad->p.user_data; */
 /*       int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM), ed->deg); */

 /*       d4est_element_data_apply_fofufofvlilj */
 /*         ( */
 /*          d4est_ops, */
 /*          d4est_geom, */
 /*          &prob_vecs->u[nodal_stride], */
 /*          NULL, */
 /*          NULL, */
 /*          ed, */
 /*          ed->deg_stiffness, // + params->deg_offset_for_nonlinear_quad, */
 /*          d4est_geom->geom_quad_type, */
 /*          (P4EST_DIM), */
 /*          &M_helmf_u[nodal_stride],            */
 /*          helmholtz_fcn, */
 /*          prob_vecs->user, */
 /*          NULL, */
 /*          NULL */
 /*         ); */

 /*       nodal_stride += volume_nodes_lobatto; */
 /*     } */
 /*   } */
  
 /*  d4est_linalg_vec_axpy(1.0, M_helmf_u, prob_vecs->Au, prob_vecs->local_nodes); */
 /*  P4EST_FREE(M_helmf_u); */
}

static
void problem_build_rhs
(
 p4est_t* p4est,
 d4est_elliptic_data_t* prob_vecs,
 d4est_elliptic_eqns_t* prob_fcns,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 problem_input_t* input,
 void* user
)
{
  double* f = P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_mesh_init_field
    (
     p4est,
     f,
     f_fcn,
     d4est_ops,
     d4est_geom
    );

  prob_vecs->flux_fcn_data = curved_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
                                             (boundary_fcn, &global_ip_flux_params);

  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = quad->p.user_data;

        /* if (d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_LOBATTO")){ */
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
#if (P4EST_DIM)==3
        mesh_object.q[2] = ed->q[2];
#endif
        
        d4est_quadrature_apply_mass_matrix
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &f[ed->nodal_stride],
           ed->deg,
           ed->J_quad,
           ed->deg_quad,
           &prob_vecs->rhs[ed->nodal_stride]
          );
        
        /* } */
        /* else if(d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_QUADRATURE")){ */
          /* d4est_element_data_apply_fofufofvlj */
            /* ( */
             /* d4est_ops, */
             /* d4est_geom, */
             /* NULL, */
             /* NULL, */
             /* ed, */
             /* ed->deg_stiffness, */
             /* d4est_geom->geom_quad_type, */
             /* (P4EST_DIM), */
             /* &prob_vecs->rhs[ed->nodal_stride], */
             /* f_fcn_ext, */
             /* NULL, */
             /* NULL, */
             /* NULL */
            /* ); */
        /* } */
        /* else { */
          /* D4EST_ABORT("Should not happen\n"); */
        /* } */


        printf("elem_id, rhs sum = %d %.25f\n", ed->id, d4est_linalg_vec_sum(&prob_vecs->rhs[ed->nodal_stride], d4est_lgl_get_nodes((P4EST_DIM), ed->deg)));
        
        
      }
    }


  
  int local_nodes = prob_vecs->local_nodes;
  double* u_eq_0 = P4EST_ALLOC_ZERO(double, local_nodes);
  double* tmp = prob_vecs->u;

  prob_vecs->u = u_eq_0;
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad);
  d4est_linalg_vec_axpy(-1., prob_vecs->Au, prob_vecs->rhs, local_nodes);

  prob_vecs->u = tmp;
  P4EST_FREE(u_eq_0);

  prob_vecs->flux_fcn_data = curved_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
                                             (zero_fcn, &global_ip_flux_params);

#ifdef NASTY_DEBUG
  DEBUG_PRINT_ARR_DBL(prob_vecs->rhs, prob_vecs->local_nodes);
#endif
  
  P4EST_FREE(f);
}

/* static */
/* void problem_build_rhs_weakbc */
/* ( */
/*  p4est_t* p4est, */
/*  d4est_elliptic_data_t* prob_vecs, */
/*  d4est_elliptic_eqns_t* prob_fcns, */
/*  p4est_ghost_t* ghost, */
/*  d4est_element_data_t* ghost_data, */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_geometry_t* d4est_geom, */
/*  problem_input_t* input, */
/*  void* user */
/* ) */
/* { */
/*   double* f = P4EST_ALLOC(double, prob_vecs->local_nodes); */
/*   d4est_mesh_init_field */
/*     ( */
/*      p4est, */
/*      f, */
/*      f_fcn, */
/*      d4est_ops, */
/*      d4est_geom */
/*     ); */

/*   for (p4est_topidx_t tt = p4est->first_local_tree; */
/*        tt <= p4est->last_local_tree; */
/*        ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int Q = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int q = 0; q < Q; ++q) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q); */
/*         d4est_element_data_t* ed = quad->p.user_data; */
        
/*         if (d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_LOBATTO")){ */
/*           d4est_operators_apply_curved_mass_matrix(d4est_ops, */
/*                                           &f[ed->nodal_stride], */
/*                                           ed->deg, */
/*                                           ed->J_quad, */
/*                                           ed->deg_quad, */
/*                                           d4est_geom->geom_quad_type, */
/*                                           (P4EST_DIM), */
/*                                           &prob_vecs->rhs[ed->nodal_stride] */
/*                                          ); */
/*         } */
/*         else if(d4est_util_match(input->rhs_compute_method,"COMPUTE_RHS_ON_QUADRATURE")){ */
/*           d4est_element_data_apply_fofufofvlj */
/*             ( */
/*              d4est_ops, */
/*              d4est_geom, */
/*              NULL, */
/*              NULL, */
/*              ed, */
/*              ed->deg_stiffness, */
/*              d4est_geom->geom_quad_type, */
/*              (P4EST_DIM), */
/*              &prob_vecs->rhs[ed->nodal_stride], */
/*              f_fcn_ext, */
/*              NULL, */
/*              NULL, */
/*              NULL */
/*             ); */
/*         } */
/*         else { */
/*           D4EST_ABORT("Should not happen\n"); */
/*         } */


/*         printf("elem_id, rhs sum = %d %.25f\n", ed->id, d4est_linalg_vec_sum(&prob_vecs->rhs[ed->nodal_stride], d4est_lgl_get_nodes((P4EST_DIM), ed->deg))); */
        
        
/*       } */
/*     } */
  
/*   P4EST_FREE(f); */
/* } */



static
void
problem_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad
)
{
  problem_apply_lhs(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom, d4est_quad);
  d4est_linalg_vec_xpby(prob_vecs->rhs, -1., prob_vecs->Au, prob_vecs->local_nodes);
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
  
  d4est_elliptic_data_t prob_vecs;
  prob_vecs.rhs = rhs;
  prob_vecs.Au = Au;
  prob_vecs.u = u;
  prob_vecs.local_nodes = local_nodes;

  ip_flux_t* ip_flux = ip_flux_dirichlet_new(p4est, "[IP_FLUX]", input_file, zero_fcn);
  /* ip_flux_t* ip_flux = ip_flux_dirichlet_new(p4est, "[IP_FLUX]", input_file, boundary_fcn); */
  prob_vecs.flux_fcn_data = ip_flux->mortar_fcn_ptrs;
  global_ip_flux_params = *(ip_flux->ip_flux_params);
  prob_vecs.user = &input;

  d4est_elliptic_eqns_t prob_fcns;
  prob_fcns.build_residual = problem_build_residual;
  prob_fcns.apply_lhs = problem_apply_lhs;
  
  d4est_mesh_geometry_storage_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = d4est_quadrature_new(p4est, d4est_ops, d4est_geom, input_file, "quadrature", "[QUADRATURE]");


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

  }


  d4est_mesh_update
    (
     p4est,
     ghost,
     ghost_data,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     geometric_factors,
     INITIALIZE_QUADRATURE_DATA,
     INITIALIZE_GEOMETRY_DATA,
     INITIALIZE_GEOMETRY_ALIASES,
     problem_set_degrees,
     (void*)&input
    );
  

    
    local_nodes = d4est_mesh_get_local_nodes(p4est);

    Au = P4EST_REALLOC(Au, double, local_nodes);
    u = P4EST_REALLOC(u, double, local_nodes);
    rhs = P4EST_REALLOC(rhs, double, local_nodes);
    u_analytic = P4EST_REALLOC(u_analytic, double, local_nodes);
    error = P4EST_REALLOC(error, double, local_nodes);

    prob_vecs.Au = Au;
    prob_vecs.u = u;
    prob_vecs.rhs = rhs;
    prob_vecs.local_nodes = local_nodes;

    d4est_linalg_fill_vec(prob_vecs.u, 0., prob_vecs.local_nodes);


    /* d4est_mesh_init_field( */
    /*                                   p4est, */
    /*                                   u, */
    /*                                   u_poly, */
    /*                                   d4est_ops, */
    /*                                   d4est_geom */
    /* );   */
    
    d4est_mesh_init_field(
                                      p4est,
                                      u_analytic,
                                      analytic_solution_fcn,
                                      d4est_ops,
                                      d4est_geom
    );  
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    DEBUG_PRINT_ARR_DBL_SUM(error, local_nodes);
    DEBUG_PRINT_ARR_DBL_SUM(prob_vecs.u, local_nodes);
    DEBUG_PRINT_ARR_DBL_SUM(u_analytic, local_nodes);
    
    double initial_l2_norm_sqr_local = d4est_mesh_compute_l2_norm_sqr
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       error,
       local_nodes,
       DO_NOT_STORE_LOCALLY
      );
    printf("Initial local l2 norm = %.25f\n", sqrt(initial_l2_norm_sqr_local));


    
    /* problem_build_rhs_weakbc */
    /*   ( */
    /*    p4est, */
    /*    &prob_vecs, */
    /*    &prob_fcns, */
    /*    ghost, */
    /*    ghost_data, */
    /*    d4est_ops, */
    /*    d4est_geom, */
    /*    &input, */
    /*    NULL */
    /*   ); */

    problem_build_rhs
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       ghost,
       ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &input,
       NULL
      );

    
    clock_t begin = 0;
    clock_t end = -1;

    if (world_rank == 0){
      begin = clock();
    }

    /* double volume; */
    /* double surface_area; */
    /* d4est_mesh_compute_grid_volume_and_surface_area(p4est, &volume, &surface_area); */

    /* double R2 = global_cubed_sphere_attrs.R2; */
    /* double sqrl = 2*global_cubed_sphere_attrs.Clength; */
    /* double volume_theory = (4./3.)*M_PI*R2*R2*R2 - sqrl*sqrl*sqrl; */
    /* double surface_area_theory = (4)*M_PI*R2*R2 + sqrl*sqrl*6; */
    
    /* printf("[GRID_INFO]: volume numerical/theoretical = %.25f/%.25f\n", volume, volume_theory); */
    /* printf("[GRID_INFO]: surface_area numerical/theoretical = %.25f/%.25f\n", surface_area, surface_area_theory); */

    /* double jac_error_sum = 0.; */
    /* double drdx_error_sum = 0.; */
    /* double drdxjacdrdx_error_sum = 0.; */
    
 /*    for (p4est_topidx_t tt = p4est->first_local_tree; */
/*          tt <= p4est->last_local_tree; */
/*          ++tt) */
/*     { */
/*       p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
/*       sc_array_t* tquadrants = &tree->quadrants; */
/*       int QQ = (p4est_locidx_t) tquadrants->elem_count; */
/*       for (int qq = 0; qq < QQ; ++qq) { */
/*         p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq); */
/*         d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data); */

/*         int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM),ed->deg_quad); */
/*         double rst [(P4EST_DIM)]; */
/*         double drst_dxyz_i [(P4EST_DIM)][(P4EST_DIM)]; */
/*         double drdxjacdrdx_i_anal [(P4EST_DIM)][(P4EST_DIM)]; */
/*         double drdxjacdrdx_i_mult [(P4EST_DIM)][(P4EST_DIM)]; */

/*         d4est_rst_t rst_points */
/*           = d4est_operators_get_rst_points(d4est_ops, ed->deg_quad, (P4EST_DIM), GAUSS); */
  
/*         for (int i = 0; i < volume_nodes; i++){ */
/*           rst[0] = rst_points.r[i]; */
/*           rst[1] = rst_points.s[i]; */
/* #if (P4EST_DIM)==3 */
/*           rst[2] = rst_points.t[i]; */
/* #endif */

/*           if (d4est_geom->DRDX == NULL){ */
/*             D4EST_ABORT("d4est_geom->DRDX == NULL"); */
/*           } */
/*           d4est_geom->DRDX(d4est_geom, tt, ed->q, ed->dq, rst, drst_dxyz_i); */

/*           if (d4est_geom->JAC == NULL){ */
/*             D4EST_ABORT("d4est_geom->JAC == NULL"); */
/*           } */

/*           d4est_geom->JACDRDXDRDX(d4est_geom, tt, ed->q, ed->dq, rst, drdxjacdrdx_i_anal); */
/*           double jac; */
/*           d4est_geom->JAC(d4est_geom, tt, ed->q, ed->dq, rst, &jac); */
          
/*           for (int d1 = 0; d1 < (P4EST_DIM); d1++) { */
/*             for (int d2 = 0; d2 < (P4EST_DIM); d2++) { */
/*               drdxjacdrdx_i_mult[d1][d2] = 0.; */
/*               for (int k = 0; k < (P4EST_DIM); k++){ */
/*                 drdxjacdrdx_i_mult[d1][d2] += jac*ed->rst_xyz_quad[d1][k][i]*ed->rst_xyz_quad[d2][k][i]; */
/*               } */
/*               printf(" %d drdxjacdrdx_i_mult[%d][%d] drdxjacdrdx_i_anal[%d][%d] = %.15f %.15f\n", i, d1,d2,d1,d2,drdxjacdrdx_i_mult[d1][d2], drdxjacdrdx_i_anal[d1][d2]); */
/*               drdxjacdrdx_error_sum += fabs(drdxjacdrdx_i_mult[d1][d2] - drdxjacdrdx_i_anal[d1][d2]); */
/*               drdx_error_sum += fabs(ed->rst_xyz_quad[d1][d2][i] - drst_dxyz_i[d1][d2]); */
/*             }  */
/*           } */

/*           jac_error_sum += fabs(jac - ed->J_quad[i]); */


/*           printf(" %d jac_mult jac_anal = %.15f %.15f\n", i, ed->J_quad[i], jac); */
/*           printf(" %d drdx_mult drdx_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[0][0][i], drst_dxyz_i[0][0]); */
/*           printf(" %d drdy_mult drdy_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[0][1][i], drst_dxyz_i[0][1]); */
/*           printf(" %d drdz_mult drdz_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[0][2][i], drst_dxyz_i[0][2]); */
/*           printf(" %d dsdx_mult dsdx_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[1][0][i], drst_dxyz_i[1][0]); */
/*           printf(" %d dsdy_mult dsdy_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[1][1][i], drst_dxyz_i[1][1]); */
/*           printf(" %d dsdz_mult dsdz_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[1][2][i], drst_dxyz_i[1][2]); */
/*           printf(" %d dtdx_mult dtdx_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[2][0][i], drst_dxyz_i[2][0]); */
/*           printf(" %d dtdy_mult dtdy_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[2][1][i], drst_dxyz_i[2][1]); */
/*           printf(" %d dtdz_mult dtdz_anal = %.15f %.15f\n",i, ed->rst_xyz_quad[2][2][i], drst_dxyz_i[2][2]); */

/*           printf("drdx_error_sum jac_error_sum drdxjacdrdx_error_sum = %.15f %.15f %.15f\n", drdx_error_sum, jac_error_sum, drdxjacdrdx_error_sum); */
/*         } */
                
/*       } */
/*     } */



     /* if (input.solve_with_multigrid){ */
     
    /* int min_level, max_level; */

    /* multigrid_get_level_range(p4est, &min_level, &max_level); */
    /* printf("[min_level, max_level] = [%d,%d]\n", min_level, max_level); */

    /* need to do a reduce on min,max_level before supporting multiple proc */
    /* D4EST_ASSERT(proc_size == 1); */
    /* int num_of_levels = max_level + 1; */
     
      
    
    /* multigrid_smoother_t* smoother = multigrid_smoother_cheby_d4est_init */
    /*                                ( */
    /*                                 p4est, */
    /*                                 num_of_levels, */
    /*                                 input_file */
    /*                                ); */


    /* multigrid_smoother_t* smoother = multigrid_smoother_krylov_petsc_init(p4est, input_file); */
    
    /* multigrid_bottom_solver_t* bottom_solver = multigrid_bottom_solver_cg_d4est_init */
                                               /* ( */
                                                /* p4est, */
                                                /* input_file */
                                               /* ); */

    /* multigrid_bottom_solver_t* bottom_solver = multigrid_bottom_solver_krylov_petsc_init */
    /*                                            ( */
    /*                                             p4est, */
    /*                                             input_file */
    /*                                            ); */
    
    /* multigrid_bottom_solver_t* bottom_solver = multigrid_bottom_solver_cheby_d4est_init */
    /*                                            ( */
    /*                                             p4est, */
    /*                                             num_of_levels, */
    /*                                             input_file */
    /*                                            ); */
    
    /* multigrid_logger_t* logger = multigrid_logger_residual_init */
                                 /* ( */
                                 /* ); */

    /* multigrid_element_data_updater_t* updater */
      /* = multigrid_element_data_updater_curved_init */
      /* ( */
       /* num_of_levels, */
       /* &ghost, */
       /* &ghost_data, */
       /* geometric_factors, */
       /* d4est_geom, */
       /* set_deg_quad, */
       /* &input */
      /* ); */

    /* multigrid_user_callbacks_t* matrix_op_callbacks = multigrid_matrix_operator_init */
    /*                                                   ( */
    /*                                                    p4est, */
    /*                                                    num_of_levels, */
    /*                                                    d4est_ops, */
    /*                                                    d4est_element_data_get_local_matrix_nodes, */
    /*                                                    NULL */
    /*                                                   ); */
    
    /* multigrid_data_t* mg_data = multigrid_data_init(p4est, */
                                                    /* d4est_ops, */
                                                    /* num_of_levels, */
                                                    /* logger, */
                                                    /* NULL, //matrix_op_callbacks, */
                                                    /* matrix_op_callbacks, */
                                                    /* updater, */
                                                    /* input_file */
                                                   /* ); */




      /* prob_vecs.user = matrix_op_callbacks->user; */

      /* multigrid_matrix_curved_fofu_fofv_mass_operator_setup_deg_quad_eq_deg */
      /*   ( */
      /*    p4est, */
      /*    d4est_ops, */
      /*    d4est_geom, */
      /*    NULL, */
      /*    NULL, */
      /*    helmholtz_fcn, */
      /*    NULL, */
      /*    NULL, */
      /*    NULL, */
      /*    matrix_op_callbacks->user, */
      /*    set_deg_gauss, */
      /*    &input */
      /*   ); */
    
      /* multigrid_solve */
      /*   ( */
      /*    p4est, */
      /*    &prob_vecs, */
      /*    &prob_fcns, */
      /*    mg_data */
      /*   ); */

    /* if(input.use_mg_as_pc_for_ksp){ */
      /* krylov_pc_t* pc = krylov_pc_multigrid_create(mg_data, NULL); */
      /* krylov_petsc_params_t petsc_params; */
 
      /* krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &petsc_params); */
                
      /* krylov_info_t info = */
        /* krylov_petsc_solve */
        /* ( */
         /* p4est, */
         /* &prob_vecs, */
         /* (void*)&prob_fcns, */
         /* &ghost, */
         /* (void**)&ghost_data, */
         /* d4est_ops, */
         /* d4est_geom, */
         /* &petsc_params, */
         /* pc */
        /* ); */

      /* krylov_pc_multigrid_destroy(pc); */
    /* } */
    /* else { */
      /* multigrid_solve */
        /* ( */
         /* p4est, */
         /* &prob_vecs, */
         /* &prob_fcns, */
         /* mg_data */
        /* );      */
    /* } */
      /* multigrid_smoother_cheby_d4est_destroy(smoother); */
      /* multigrid_smoother_krylov_petsc_destroy(smoother); */

      /* multigrid_bottom_solver_cg_d4est_destroy(bottom_solver); */
      /* multigrid_bottom_solver_krylov_petsc_destroy(bottom_solver); */
      /* multigrid_bottom_solver_cheby_d4est_destroy(bottom_solver); */
      /* multigrid_logger_residual_destroy(logger); */
      /* multigrid_element_data_updater_curved_destroy(updater, num_of_levels); */
      /* multigrid_matrix_operator_destroy(matrix_op_callbacks); */
      /* multigrid_data_destroy(mg_data); */

     
     /* } */
     /* else { */
      /* krylov_petsc_params_t petsc_params; */
      /* krylov_petsc_input(p4est, input_file, "krylov_petsc", "[KRYLOV_PETSC]", &petsc_params);       */
       
      /* krylov_info_t info = */
      /*   krylov_petsc_solve */
      /*   ( */
      /*    p4est, */
      /*    &prob_vecs, */
      /*    (void*)&prob_fcns, */
      /*    &ghost, */
      /*    (void**)&ghost_data, */
      /*    d4est_ops, */
      /*    d4est_geom, */
      /*    d4est_quad, */
      /*    &petsc_params, */
      /*    NULL */
      /*   ); */


   d4est_cg_params_t cg_params;
    d4est_cg_input(p4est, input_file, "d4est_cg", "[D4EST_CG]", &cg_params);

    d4est_cg_solve
      (
       p4est,
       &prob_vecs,
       &prob_fcns,
       &ghost,
       (void**)&ghost_data,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       &cg_params
      );
    

      
     /* } */
    
    d4est_mesh_init_field(
                                      p4est,
                                      u_analytic,
                                      analytic_solution_fcn,
                                      d4est_ops,
                                      d4est_geom
                                     );
    d4est_linalg_vec_axpyeqz(-1., u, u_analytic, error, local_nodes);

    
    double local_l2_norm_sqr = d4est_mesh_compute_l2_norm_sqr
      (
       p4est,
       d4est_ops,
       d4est_geom,
       d4est_quad,
       error,
       local_nodes,
       DO_NOT_STORE_LOCALLY
      );



    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    /* DEBUGGING STARTS HERE */
    
    /* d4est_mesh_init_field( */
    /*                                   p4est, */
    /*                                   u, */
    /*                                   analytic_solution_fcn, */
    /*                                   d4est_ops, */
    /*                                   d4est_geom */
    /*                                  ); */

    /* prob_vecs.flux_fcn_data.bndry_fcn = boundary_fcn; */
    /* problem_apply_lhs(p4est, ghost, ghost_data, &prob_vecs, d4est_ops, d4est_geom); */


    /* double* Au_modal = P4EST_ALLOC(double, local_nodes); */
    /* int nodal_stride = 0; */
    /* for (p4est_topidx_t tt = p4est->first_local_tree; */
    /*      tt <= p4est->last_local_tree; */
    /*      ++tt) */
    /* { */
    /*   p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt); */
    /*   sc_array_t* tquadrants = &tree->quadrants; */
    /*   int QQ = (p4est_locidx_t) tquadrants->elem_count; */
    /*   for (int qq = 0; qq < QQ; ++qq) { */
    /*     p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, qq); */
    /*     d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data); */
    /*     int volume_nodes = d4est_lgl_get_nodes((P4EST_DIM),ed->deg); */
    /*     d4est_operators_convert_nodal_to_modal(d4est_ops, &Au[nodal_stride], */
    /*                                   (P4EST_DIM), ed->deg, &Au_modal[nodal_stride]); */
    /*     nodal_stride += volume_nodes; */
    /*   } */
    /* } */
    
    
    /* DEBUG_PRINT_2ARR_DBL(Au, Au_modal, local_nodes); */

    /* double* Au_error = P4EST_ALLOC(double, local_nodes); */
    /* double* jacobian_lgl = P4EST_ALLOC(double, local_nodes); */
    /* d4est_element_data_compute_jacobian_on_lgl_grid(p4est,d4est_geom, d4est_ops, jacobian_lgl); */
    /* d4est_linalg_copy_1st_to_2nd(Au, Au_error, local_nodes); */
    /* d4est_linalg_vec_fabs(Au_error, local_nodes); */

    /* d4est_linalg_vec_fabs(error, local_nodes); /\* not needed, except for vtk output *\/ */
    /* vtk_nodal_vecs_debug_t vtk_nodal_vecs_debug; */
    /* vtk_nodal_vecs_debug.u = u; */
    /* vtk_nodal_vecs_debug.Au = prob_vecs.Au; */
    /* vtk_nodal_vecs_debug.Au_error = Au_error; */
    /* vtk_nodal_vecs_debug.jacobian = jacobian_lgl; */

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
    /*    vtk_field_plotter_debug, */
    /*    (void*)&vtk_nodal_vecs_debug */
    /*   ); */

    /* P4EST_FREE(jacobian_lgl); */
    /* P4EST_FREE(deg_array); */
    /* P4EST_FREE(Au_error); */
    /* P4EST_FREE(Au_modal); */
    
    
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    /* DEBUGGING ENDS HERE */
    
    double* jacobian_lgl = P4EST_ALLOC(double, local_nodes);
    d4est_mesh_compute_jacobian_on_lgl_grid(p4est, d4est_ops, d4est_geom, jacobian_lgl);


    /* DEBUG_PRINT_ARR_DBL(jacobian_lgl, local_nodes); */
    
    d4est_linalg_vec_fabs(error, local_nodes); /* not needed, except for vtk output */
    vtk_nodal_vecs_t vtk_nodal_vecs;
    vtk_nodal_vecs.u = u;
    vtk_nodal_vecs.u_analytic = u_analytic;
    vtk_nodal_vecs.error = error;
    vtk_nodal_vecs.jacobian = jacobian_lgl;

    int* deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_mesh_get_array_of_degrees(p4est, deg_array);

    char save_as [500];
    sprintf(save_as, "%s_level_%d", "test_uniform_cubedsphere", input.num_unifrefs);
    
    /* vtk output */
    d4est_vtk_save_geometry_and_dg_fields
      (
       save_as,
       p4est,
       d4est_ops,
       deg_array,
       input_file,
       "d4est_vtk_geometry",
       vtk_field_plotter,
       (void*)&vtk_nodal_vecs
      );

    P4EST_FREE(jacobian_lgl);
    P4EST_FREE(deg_array);
    double local_nodes_dbl = (double)local_nodes;
    double local_reduce [2];
    double global_reduce [2];

    local_reduce[0] = local_l2_norm_sqr;
    local_reduce[1] = local_nodes_dbl;
    
    sc_reduce
      (
       &local_reduce[0],
       &global_reduce[0],
       2,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       0,
       sc_MPI_COMM_WORLD
      );

    double global_l2_norm_sqr = global_reduce[0];
    double global_nodes_dbl = global_reduce[1];

    if (world_rank == 0){
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf
      (
       "[HP_AMR]: %d, %d, %d, %.25f, %f\n",
       level,
       (int)p4est->global_num_quadrants,
       (int)global_nodes_dbl,
       sqrt(global_l2_norm_sqr),
       time_spent
      );
    }

  d4est_mesh_geometry_storage_destroy(geometric_factors);
  ip_flux_dirichlet_destroy(ip_flux);
  
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }

  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  
  P4EST_FREE(Au);
  P4EST_FREE(rhs);
  P4EST_FREE(u_analytic);
  P4EST_FREE(error);
  P4EST_FREE(u);
}
