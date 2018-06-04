#include <d4est_power_method.h>
#include <d4est_linalg.h>
#include <sc_reduce.h>


double
d4est_power_method
(
 p4est_t* p4est,
 d4est_elliptic_data_t* vecs,
 d4est_elliptic_eqns_t* fcns,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 int imax
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_power_method");

  int nodes = vecs->local_nodes;
  
  double* u = P4EST_ALLOC(double, nodes);
  double* Au = P4EST_ALLOC(double, nodes);

  for (int i = 0; i < nodes; i++){
    u[i] = (double)i;
  }
  
  d4est_elliptic_data_t vecs_for_power_method;
  d4est_elliptic_data_copy_ptrs(vecs, &vecs_for_power_method);
  vecs_for_power_method.u = u;
  vecs_for_power_method.Au = Au;

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

    double sqrt_global_norm = sqrt(global_norm);

    for (int i = 0; i < nodes; i++){
      vecs_for_power_method.u[i] = vecs_for_power_method.Au[i]/sqrt_global_norm;
    }

    zlog_info(c_default, "%d/%d: eig_val = %.15f\n", k, imax, sqrt_global_norm);
  }


  double local_u_dot_u = d4est_linalg_vec_dot(vecs_for_power_method.u,vecs_for_power_method.u,nodes);
  double local_u_dot_Au = d4est_linalg_vec_dot(vecs_for_power_method.u,vecs_for_power_method.Au,nodes);
  double local_dots [] = {local_u_dot_u, local_u_dot_Au};
  double global_dots [2];
     
  sc_allreduce
    (
     &local_dots[0],
     &global_dots[0],
     2,
     sc_MPI_DOUBLE,
     sc_MPI_SUM,
     sc_MPI_COMM_WORLD
    );
  
  
  D4EST_FREE(u);
  D4EST_FREE(Au);
  return global_dots[1]/global_dots[0];
}
