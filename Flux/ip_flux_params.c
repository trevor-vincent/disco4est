#include <ip_flux_params.h>
#include <grid_functions.h>
#include <util.h>
#include <ini.h>

static double
ip_flux_params_penalty_meanp_sqr_over_meanh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double mean_p = .5*(deg_m + deg_p);
  double mean_p_sqr = mean_p*mean_p;
  double mean_h = .5*(h_m + h_p);
  return (penalty_prefactor*mean_p_sqr)/mean_h;
}

static double
ip_flux_params_penalty_maxp_sqr_over_minh
(
 int deg_m,
 double h_m,
 int deg_p,
 double h_p,
 double penalty_prefactor
)
{
  double max_deg = (deg_m > deg_p) ? deg_m : deg_p;
  double min_h = (h_m < h_p) ? h_m : h_p;
  return (penalty_prefactor*(max_deg)*(max_deg))/min_h;
}

static void
ip_flux_params_get_string_from_h_calc
(
 h_calc_method_t h_calc,
 char string [50]
)
{
  if (h_calc == H_EQ_J_DIV_SJ){
    string = "H_EQ_J_DIV_SJ";
  }
  else if (h_calc == H_EQ_J_DIV_SJ_MIN){
    string = "H_EQ_J_DIV_SJ_MIN";
  }
  else if (h_calc == H_EQ_VOLUME_DIV_AREA){
    string = "H_EQ_VOLUME_DIV_AREA";
  }
  else {
    string = "NOT_SET";
  }
}

static void
ip_flux_params_get_string_from_bc_eval
(
 bc_eval_t bc_eval,
 char string [50]
)
{
  if (bc_eval == BC_EVAL_ON_QUADRATURE_POINTS){
    string = "BC_EVAL_ON_QUADRATURE_POINTS";
  }
  else if (bc_eval == BC_EVAL_ON_LOBATTO_POINTS){
    string = "BC_EVAL_ON_LOBATTO_POINTS";
  }
  else {
    string = "NOT_SET";
  }
}

static void
ip_flux_params_get_string_from_penalty_fcn
(
 penalty_calc_t fcn,
 char string [50]
)
{
  if (fcn == ip_flux_params_penalty_maxp_sqr_over_minh){
    string = "maxp_sqr_over_minh";
  }
  else if (fcn == ip_flux_params_penalty_meanp_sqr_over_meanh){
    string = "meanp_sqr_over_meanh";
  }
  else {
    string = "not_set";
  }
}

static penalty_calc_t
ip_flux_params_get_penalty_fcn_from_string
(
 const char* string
)
{
  if (util_match(string,"maxp_sqr_over_minh")){
    return ip_flux_params_penalty_maxp_sqr_over_minh;
  }
  else if (util_match(string,"meanp_sqr_over_meanh")){
    return ip_flux_params_penalty_meanp_sqr_over_meanh;
  }
  else {
    mpi_abort("This ip flux penalty calculation fcn does not exist");
    return NULL;
  }
}

static
int ip_flux_params_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  ip_flux_params_t* pconfig = (ip_flux_params_t*)user;
  if (util_match_couple(section,"ip_flux_params",name,"name")) {
    mpi_assert(pconfig->name[0] == '*');
    snprintf (pconfig->name, sizeof(pconfig->name), "%s", value);
  }
  
  if (util_match_couple(section,"ip_flux_params",name,"ip_flux_penalty_prefactor")) {
    mpi_assert(pconfig->ip_flux_penalty_prefactor == -1);
    pconfig->ip_flux_penalty_prefactor = atof(value);
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_penalty_calculate_fcn")) {
    mpi_assert(pconfig->ip_flux_penalty_calculate_fcn == NULL);
    pconfig->ip_flux_penalty_calculate_fcn = ip_flux_params_get_penalty_fcn_from_string(value);
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_bc_eval")) {
    mpi_assert(pconfig->ip_flux_bc_eval == BC_EVAL_NOTSET);
    if(util_match(value, "BC_EVAL_ON_QUADRATURE_POINTS")){
      pconfig->ip_flux_bc_eval = BC_EVAL_ON_QUADRATURE_POINTS;
    }
    else if (util_match(value, "BC_EVAL_ON_LOBATTO_POINTS")){
      pconfig->ip_flux_bc_eval = BC_EVAL_ON_LOBATTO_POINTS;
    }
    else {
      printf("ip_flux_params_bc_eval = %s\n", value);
      mpi_abort("ip_flux_params_bc_eval is not set to BC_EVAL_ON_LOBATTO_POINTS or BC_EVAL_ON_QUADRATURE_POINTS\n");
    }
  }
  else if (util_match_couple(section,"ip_flux_params",name,"ip_flux_h_calc")) {
    mpi_assert(pconfig->ip_flux_h_calc == H_EQ_NOTSET);
    if(util_match(value, "H_EQ_J_DIV_SJ")){
      pconfig->ip_flux_h_calc = H_EQ_J_DIV_SJ;
    }
    else if(util_match(value, "H_EQ_J_DIV_SJ_MIN")){
      pconfig->ip_flux_h_calc = H_EQ_J_DIV_SJ_MIN;
    }
    else if (util_match(value, "H_EQ_VOLUME_DIV_AREA")){
      pconfig->ip_flux_h_calc = H_EQ_VOLUME_DIV_AREA;
    }
    else {
      printf("ip_flux_params_h_calc = %s\n", value);
      mpi_abort("ip_flux_params_h_calc is not set to H_EQ_J_DIV_SJ or H_EQ_VOLUME_DIV_AREA\n");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static void
ip_flux_params_input
(
 p4est_t* p4est,
 const char* printf_prefix,
 const char* input_file,
 ip_flux_params_t* input
)
{
  /* set defaults */
  input->name[0] = '*';
  input->ip_flux_bc_eval = BC_EVAL_NOTSET;
  input->ip_flux_h_calc = H_EQ_NOTSET;
  input->ip_flux_penalty_calculate_fcn = NULL;
  input->ip_flux_penalty_prefactor = -1.;

  if (ini_parse(input_file, ip_flux_params_input_handler, input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_penalty_prefactor, -1);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_penalty_calculate_fcn, NULL);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_h_calc, H_EQ_NOTSET);
  D4EST_CHECK_INPUT("ip_flux_params", input->ip_flux_bc_eval, BC_EVAL_NOTSET);
  D4EST_CHECK_INPUT("ip_flux_params", input->name[0], '*');
  
  char penalty_calculate_fcn [50];
  char h_eq [50];
  char bc_eval [50];  

  ip_flux_params_get_string_from_penalty_fcn (input->ip_flux_penalty_calculate_fcn,penalty_calculate_fcn);
  ip_flux_params_get_string_from_bc_eval (input->ip_flux_bc_eval,bc_eval);
  ip_flux_params_get_string_from_h_calc (input->ip_flux_h_calc,h_eq);
  
  if(p4est->mpirank == 0){
    printf("%s: ip_flux_params_penalty_prefactor = %f\n", printf_prefix, input->ip_flux_penalty_prefactor);
    printf("%s: ip_flux_params_penalty_calculate_fcn = %s\n", printf_prefix, penalty_calculate_fcn);
    printf("%s: ip_flux_params_h_calc = %s\n", printf_prefix, h_eq);
    printf("%s: ip_flux_params_bc_eval = %s\n", printf_prefix, bc_eval);
  }
}

ip_flux_params_t*
ip_flux_params_new
(
 p4est_t* p4est,
 const char* print_prefix,
 const char* input_file
)
{
  ip_flux_params_t* ip_flux_params = P4EST_ALLOC(ip_flux_params_t, 1); 
  ip_flux_params_input(p4est, print_prefix, input_file, ip_flux_params);
  return ip_flux_params;
}

void
ip_flux_params_destroy
(
 ip_flux_params_t* params
){
  P4EST_FREE(params);
}
