#include <d4est_solver_cg.h>
#include <d4est_util.h>
#include <ini.h>
#include <d4est_linalg.h>
#include <sc_reduce.h>

static
int d4est_solver_cg_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_cg_params_t* pconfig = (d4est_solver_cg_params_t*)user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"atol")) {
    D4EST_ASSERT(pconfig->atol == -1);
    pconfig->atol = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"rtol")) {
    D4EST_ASSERT(pconfig->rtol == -1);
    pconfig->rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"monitor")) {
    D4EST_ASSERT(pconfig->monitor == -1);
    D4EST_ASSERT(atoi(value) == 0 || atoi(value) == 1);
    pconfig->monitor = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"iter")) {
    D4EST_ASSERT(pconfig->iter == -1);
    D4EST_ASSERT(atoi(value) >= 0);
    pconfig->iter = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_cg_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_solver_cg_params_t* input
)
{
  input->monitor = -1;
  input->iter = -1;
  input->rtol = -1;
  input->atol = -1;

  D4EST_ASSERT(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  
  if (ini_parse(input_file, d4est_solver_cg_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input->iter, -1);
  D4EST_CHECK_INPUT(input_section, input->monitor, -1);
  D4EST_CHECK_INPUT(input_section, input->atol, -1);
  D4EST_CHECK_INPUT(input_section, input->rtol, -1);
    
  if(p4est->mpirank == 0){
    printf("%s: iter = %d\n",printf_prefix, input->iter);
    printf("%s: atol = %.15f\n",printf_prefix, input->atol);
    printf("%s: rtol = %.15f\n",printf_prefix, input->rtol);
    printf("%s: monitor = %d\n",printf_prefix, input->monitor);
  }
}

void d4est_solver_cg_solve
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 p4est_ghost_t** ghost,
 d4est_element_data_t** ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 void* cg_params,
 krylov_pc_t* krylov_pc
)
{ 
  d4est_solver_cg_params_t* params = cg_params;
  int local_nodes;
  double delta_new, delta_0, delta_old, beta, alpha;

  double* Au;
  double* u;
  double* rhs;

  double* d;
  double* r;

  local_nodes = vecs->local_nodes;
  Au = vecs->Au;
  u = vecs->u;
  rhs = vecs->rhs;

  const int imax = params->iter;
  const double atol = params->atol;
  const double rtol = params->rtol;
  double d_dot_Au;

  d = P4EST_ALLOC(double, local_nodes);
  r = P4EST_ALLOC(double, local_nodes);

  /* first iteration data, store Au in r */
  d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     *ghost,
     *ghost_data,
     fcns,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );
  

  DEBUG_PRINT_ARR_DBL_SUM(vecs->u, vecs->local_nodes);
  DEBUG_PRINT_ARR_DBL_SUM(vecs->rhs, vecs->local_nodes);
  DEBUG_PRINT_ARR_DBL_SUM(vecs->Au, vecs->local_nodes);
  
  d4est_util_copy_1st_to_2nd(Au, r, local_nodes);

  /* r = f - Au ; Au is stored in r so r = rhs - r */
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);
  d4est_util_copy_1st_to_2nd(r, d, local_nodes);
  delta_new = d4est_linalg_vec_dot(r, r, local_nodes);

  double delta_new_global;
  double d_dot_Au_global;

  sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
               sc_MPI_COMM_WORLD);

  delta_new = delta_new_global;
  delta_0 = delta_new;

  /* start working on d */
  vecs->u = d;

  for (int i = 0; i < imax && (delta_new > atol*atol + delta_0 * rtol * rtol); i++) {
    
    d4est_elliptic_eqns_apply_lhs
      (
     p4est,
     *ghost,
     *ghost_data,
     fcns,
     vecs,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
      );
    
    d_dot_Au = d4est_linalg_vec_dot(d, Au, local_nodes);

    sc_allreduce(&d_dot_Au, &d_dot_Au_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);

    d_dot_Au = d_dot_Au_global;
    alpha = delta_new / d_dot_Au;

    d4est_linalg_vec_axpy(alpha, d, u, local_nodes);

    /* r = r - Au*alpha */
    d4est_linalg_vec_axpy(-alpha, Au, r, local_nodes);

    delta_old = delta_new;
    delta_new = d4est_linalg_vec_dot(r, r, local_nodes);

    sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                 sc_MPI_COMM_WORLD);
    delta_new = delta_new_global;

    beta = delta_new / delta_old;
    d4est_linalg_vec_xpby(r, beta, d, local_nodes);

    if (p4est->mpirank == 0 && params->monitor == 1)
      printf("[CG_SOLVER]: ITER %03d FNRMSQR %.30f\n", i, delta_new);

    /* DEBUG_PRINT_3ARR_DBL(u, r, Au, local_nodes); */
  }

  /* set back from d */
  vecs->u = u;

  P4EST_FREE(d);
  P4EST_FREE(r);
}
