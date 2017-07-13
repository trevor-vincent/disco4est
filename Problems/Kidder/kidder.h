#ifndef KIDDER_H
#define KIDDER_H 

typedef struct {

  d4est_poisson_flux_data_t* flux_data;
  kidder_params_t* kidder_params;
  
} problem_elliptic_eqn_params_t;

typedef struct {

  double kidder_a;
  double kidder_P;
  double kidder_E;
  
} kidder_params_t;

static
double
kidder_analytic_solution_fcn
(
 double x,
 double y,
 double z,
 void* user
)
{
  kidder_params_t* params = user;
  
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double r3 = r*r2;
  double r4 = r2*r2;

  double a = params->kidder_a;
  double E = params->kidder_E;
  
  double t1 = 1.;
  double t2 = 2.*E/r;
  double t3 = 6*a*a/(r2);
  double t4 = 2*a*a*E/(r3);
  double t5 = a*a*a*a/(r4);

  return pow((t1 + t2 + t3 + t4 + t5), .25);
}

static
double
kidder_boundary_fcn
(
 double x,
 double y,
 double z,
 vodi* user
)
{
  return kidder_analytic_solution_fcn(x,y,z,user);
}


static
double
kidder_residual_nonlinear_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  return kidder_helmholtz_fcn(x,y,z)*(1./pow(u,7));
}

static
double kidder_jacobian_nonlinear_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  return -7.*kidder_helmholtz_fcn(x,y,z)*(1./pow(u,8));
}


static
void kidder_apply_jac
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 d4est_element_data_t* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* user
)
{
  problem_elliptic_eqn_params_t* params = user;
  kidder_params_t* kidder_params = params->kidder_params;
  d4est_poisson_flux_data_t* flux_data = params->flux_data;
  flux_data->boundary_condition = zero_fcn;
  
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad);
  
  double* M_nonlinear_term = P4EST_ALLOC(double, prob_vecs->local_nodes);

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
        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
        mesh_object.q[2] = ed->q[2];        

        d4est_quadrature_apply_fofufofvlilj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           &prob_vecs->u0[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->xyz_quad,
           ed->J_quad,
           ed->deg_quad,
           &M_nonlinear_term[ed->nodal_stride],
           kidder_jacobian_nonlinear_term,
           kidder_params,
           NULL,
           NULL
          );        
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_nonlinear_term, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_nonlinear_term);
}

void
kidder_build_residual
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_elliptic_data_t* prob_vecs,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom
)
{
  problem_elliptic_eqn_params_t* params = user;
  kidder_params_t* kidder_params = params->kidder_params;
  d4est_poisson_flux_data_t* flux_data = params->flux_data;
  flux_data->boundary_condition = kidder_boundary_fcn;
  
  d4est_poisson_apply_aij(p4est, ghost, ghost_data, prob_vecs, d4est_ops, d4est_geom);

  double* M_neg_1o8_K2_psi_neg7_vec= P4EST_ALLOC(double, prob_vecs->local_nodes);
 
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



        d4est_quadrature_volume_t mesh_object;
        mesh_object.dq = ed->dq;
        mesh_object.tree = ed->tree;
        mesh_object.element_id = ed->id;
        mesh_object.q[0] = ed->q[0];
        mesh_object.q[1] = ed->q[1];
        mesh_object.q[2] = ed->q[2];

        d4est_quadrature_apply_fofufofvlj
          (
           d4est_ops,
           d4est_geom,
           d4est_quad,
           &mesh_object,
           QUAD_OBJECT_VOLUME,
           QUAD_INTEGRAND_UNKNOWN,
           &prob_vecs->u[ed->nodal_stride],
           NULL,
           ed->deg,
           ed->J_quad,
           ed->xyz_quad,
           ed->deg_quad,
           &M_neg_1o8_K2_psi_neg7_vec[ed->nodal_stride],
           kidder_neg_1o8_K2_psi_neg7,
           kidder_params,
           NULL,
           NULL
          );      
        
      }
    }

  d4est_linalg_vec_axpy(1.0,
                  M_neg_1o8_K2_psi_neg7_vec,
                  prob_vecs->Au,
                  prob_vecs->local_nodes);

  P4EST_FREE(M_neg_1o8_K2_psi_neg7_vec);
}



#endif
