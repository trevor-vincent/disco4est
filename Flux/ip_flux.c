#include <ip_flux_params.h>
#include <ip_flux.h>
#include <util.h>
#include <grid_functions.h>
#include <curved_Gauss_primal_sipg_kronbichler_flux_fcns.h>
#include <curved_Gauss_primal_sipg_kronbichler_uniform_flux_fcns.h>
#include <curved_Gauss_primal_sipg_hesthaven_flux_fcns.h>

ip_flux_t*
ip_flux_dirichlet_new
(
 p4est_t* p4est,
 const char* print_prefix,
 const char* input_file,
 grid_fcn_t bndry_fcn
)
{
  ip_flux_t* ip_flux = P4EST_ALLOC(ip_flux_t, 1); 
  ip_flux->ip_flux_params = ip_flux_params_new(p4est, print_prefix, input_file);

  /* From Hesthavens Nodal DG book */
  /* if (util_match(ip_flux->ip_flux_params->name,"hesthaven")) { */
  /*   ip_flux->curved_flux_fcn_ptrs =  */
  /*     curved_Gauss_primal_sipg_hesthaven_flux_dirichlet_fetch_fcns */
  /*     ( */
  /*      bndry_fcn, */
  /*      ip_flux->ip_flux_params */
  /*     ); */
  /* } */
  /*  */
  if (util_match(ip_flux->ip_flux_params->name,"kronbichler")) {
    ip_flux->curved_flux_fcn_ptrs =
      curved_primal_sipg_kronbichler_flux_dirichlet_fetch_fcns
      (
       bndry_fcn,
       ip_flux->ip_flux_params
      );
  }
  /* else if (util_match(ip_flux->ip_flux_params->name,"kronbichler_uniform")) { */
  /*   ip_flux->curved_flux_fcn_ptrs = */
  /*     curved_Gauss_primal_sipg_kronbichler_uniform_flux_dirichlet_fetch_fcns */
  /*     ( */
  /*      bndry_fcn, */
  /*      ip_flux->ip_flux_params */
  /*     ); */
  /* } */
  else {
    mpi_abort("ip_flux name must be hesthaven or kronbichler");
  }
  
  return ip_flux;
}

void
ip_flux_dirichlet_destroy
(
 ip_flux_t* ip_flux
)
{
  ip_flux_params_destroy(ip_flux->ip_flux_params);
  P4EST_FREE(ip_flux);
}
