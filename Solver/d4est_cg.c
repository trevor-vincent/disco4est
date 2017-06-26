#include <d4est_cg.h>
#include <util.h>
#include <ini.h>
#include <d4est_linalg.h>
#include <sc_reduce.h>

#define NASTY_DEBUG

static
int d4est_cg_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_cg_params_t* pconfig = (d4est_cg_params_t*)user;
  if (util_match_couple(section,pconfig->input_section,name,"atol")) {
    mpi_assert(pconfig->atol == -1);
    pconfig->atol = atof(value);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"rtol")) {
    mpi_assert(pconfig->rtol == -1);
    pconfig->rtol = atof(value);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"monitor")) {
    mpi_assert(pconfig->monitor == -1);
    mpi_assert(atoi(value) == 0 || atoi(value) == 1);
    pconfig->monitor = atoi(value);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"iter")) {
    mpi_assert(pconfig->iter == -1);
    mpi_assert(atoi(value) >= 0);
    pconfig->iter = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_cg_input
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_cg_params_t* input
)
{
  input->monitor = -1;
  input->iter = -1;
  input->rtol = -1;
  input->atol = -1;

  mpi_assert(sizeof(input->input_section) <= 50);
  snprintf (input->input_section, sizeof(input->input_section), "%s", input_section);
  
  if (ini_parse(input_file, d4est_cg_input_handler, input) < 0) {
    mpi_abort("Can't load input file");
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

void d4est_cg_solve
(
 p4est_t* p4est,
 problem_data_t* vecs,
 weakeqn_ptrs_t* fcns,
 p4est_ghost_t** ghost,
 void** ghost_data, 
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_cg_params_t* params
)
{ 
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
  fcns->apply_lhs(p4est, *ghost, *ghost_data, vecs, d4est_ops, d4est_geom, d4est_quad);

  /* DEBUG_PRINT_2ARR_DBL(vecs->u, vecs->Au, vecs->local_nodes); */
  
  d4est_linalg_copy_1st_to_2nd(Au, r, local_nodes);

  DEBUG_PRINT_2ARR_DBL(Au, rhs, local_nodes);

  /* r = f - Au ; Au is stored in r so r = rhs - r */
  d4est_linalg_vec_xpby(rhs, -1., r, local_nodes);
  d4est_linalg_copy_1st_to_2nd(r, d, local_nodes);
  delta_new = d4est_linalg_vec_dot(r, r, local_nodes);


  
  double delta_new_global;
  double d_dot_Au_global;

  sc_allreduce(&delta_new, &delta_new_global, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
               sc_MPI_COMM_WORLD);

  delta_new = delta_new_global;
  delta_0 = delta_new;

  /* start working on d */
  vecs->u = d;
  /* DEBUG_PRINT_2ARR_DBL(d, rhs, vecs->local_nodes); */
  /* printf("d sum = %.25f\n", d4est_linalg_vec_sum(d, vecs->local_nodes)); */
  
  int i;

#ifdef NASTY_DEBUG
  printf("imax = %d\n", imax);
  printf("monitor = %d\n", params->monitor);
  printf("delta_new > atol*atol + delta_0 * rtol * rtol = %d\n", delta_new > atol*atol + delta_0 * rtol * rtol);
  DEBUG_PRINT_2ARR_DBL(vecs->u, vecs->rhs, vecs->local_nodes);
#endif
  /* DEBUG_PRINT_2ARR_DBL(u, r, local_nodes); */
    
  for (i = 0; i < imax && (delta_new > atol*atol + delta_0 * rtol * rtol); i++) {

    
    /* Au = A*d; */
    fcns->apply_lhs(p4est, *ghost, *ghost_data, vecs, d4est_ops, d4est_geom, d4est_quad);
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
