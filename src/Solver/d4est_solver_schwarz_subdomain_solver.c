#include <pXest.h>
#include <d4est_solver_schwarz_subdomain_solver_cg.h>
#include <d4est_solver_schwarz_subdomain_solver.h>
#include <ini.h>

static
int d4est_solver_schwarz_subdomain_solver_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_solver_schwarz_subdomain_solver"); // TODO: get from function argument
  d4est_solver_schwarz_subdomain_solver_t* pconfig = (d4est_solver_schwarz_subdomain_solver_t*)user;

  if (d4est_util_match_couple(section,pconfig->input_section, name, "subdomain_solver")) {
    if(d4est_util_match(value, "cg")){
      pconfig->solver_type = SUBDOMAIN_SOLVER_CG;

      
    }
    else {
      zlog_error(c_default, "%s is not a supported subdomain solver", value);
      D4EST_ABORT("");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_solver_schwarz_subdomain_solver_t*
d4est_solver_schwarz_subdomain_solver_init
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section
)
{
  zlog_category_t* c_default = zlog_get_category("d4est_solver_schwarz_subdomain_solver");
  d4est_solver_schwarz_subdomain_solver_t* schwarz_subdomain_solver 
    = P4EST_ALLOC(d4est_solver_schwarz_subdomain_solver_t, 1);

  schwarz_subdomain_solver->solver_fcn = NULL;
  schwarz_subdomain_solver->destroy_fcn = NULL;
  schwarz_subdomain_solver->input_section = input_section;
  
  if(
     ini_parse(input_file,
               d4est_solver_schwarz_subdomain_solver_input_handler,
               schwarz_subdomain_solver) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }
  
  D4EST_CHECK_INPUT(input_section, schwarz_subdomain_solver->solver_type, SUBDOMAIN_SOLVER_NOT_SET);

  if(schwarz_subdomain_solver->solver_type == SUBDOMAIN_SOLVER_CG){
    schwarz_subdomain_solver->solver_fcn = d4est_solver_schwarz_subdomain_solver_cg;
    schwarz_subdomain_solver->destroy_fcn = d4est_solver_schwarz_subdomain_solver_cg_destroy;
    schwarz_subdomain_solver->solver_ctx =
      d4est_solver_schwarz_subdomain_solver_cg_init
      (
       p4est,
       input_file,
       input_section
      );
  }
  else{
    D4EST_ABORT("Not a supported subdomain solver");
  }
    
  if(schwarz_subdomain_solver->solver_fcn == NULL){
    D4EST_ABORT("Subdomain solver fcn not set");
  }
  if(schwarz_subdomain_solver->destroy_fcn == NULL){
    D4EST_ABORT("Subdomain solver destroy fn not set");
  }
  
  return schwarz_subdomain_solver;
}

void
d4est_solver_schwarz_subdomain_solver_destroy
(
  d4est_solver_schwarz_subdomain_solver_t* schwarz_subdomain_solver 
){
  schwarz_subdomain_solver->destroy_fcn
    (
     schwarz_subdomain_solver->solver_ctx
    );
  P4EST_FREE(schwarz_subdomain_solver);
}
