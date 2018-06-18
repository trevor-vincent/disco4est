#include <d4est_power_method.h>
#include <d4est_linalg.h>
#include <sc_reduce.h>


double
d4est_power_method
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 d4est_ghost_t* ghost,
 d4est_ghost_data_t* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 double atol,
 double rtol,
 int imax,
 int imin
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_power_method");

  int nodes = vecs->local_nodes;
  
  double* u = P4EST_ALLOC(double, nodes);
  double* Au = P4EST_ALLOC(double, nodes);

  D4EST_ASSERT(imin <= imax);
  
  for (int i = 0; i < nodes; i++){
    u[i] = (double)i;
  }
  
  d4est_elliptic_data_t vecs_for_power_method;
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_for_power_method);
  vecs_for_power_method.u = u;
  vecs_for_power_method.Au = Au;

  double sqrt_global_norm = -1;
  double sqrt_global_norm_prev = -1;
  
  for(int k = 0;  k < imax; k++){

    d4est_elliptic_eqns_apply_lhs
    (
     p4est,
     ghost,
     ghost_data,
     fcns,
     &vecs_for_power_method,
     d4est_ops,
     d4est_geom,
     d4est_quad,
     d4est_factors
    );


    double local_norm = d4est_linalg_vec_dot(vecs_for_power_method.Au,
                                             vecs_for_power_method.Au, nodes);
    double global_norm;
    
    sc_allreduce
      (
       &local_norm,
       &global_norm,
       1,
       sc_MPI_DOUBLE,
       sc_MPI_SUM,
       sc_MPI_COMM_WORLD
      );

    sqrt_global_norm = sqrt(global_norm);
    if (p4est->mpirank == 0)
      zlog_info(c_default, "k,imax,eig = %d,%d,%.15f", k, imax, sqrt_global_norm);
    
    if (k != 0 && k > imin && fabs(sqrt_global_norm - sqrt_global_norm_prev) < atol + rtol*sqrt_global_norm_prev){
      break;
    }
    
    sqrt_global_norm_prev = sqrt_global_norm;

    for (int i = 0; i < nodes; i++){
      vecs_for_power_method.u[i] = vecs_for_power_method.Au[i]/sqrt_global_norm;
    }
  }

  return sqrt_global_norm;
}
