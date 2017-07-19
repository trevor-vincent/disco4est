#ifndef OKENDON_FCNS_H
#define OKENDON_FCNS_H 

typedef struct {

  int deg;
  int deg_quad;
  double p;
  d4est_poisson_flux_data_t* flux_data_for_residual;
  d4est_poisson_flux_data_t* flux_data_for_jac;
  
} okendon_params_t;


static double
okendon_analytic_solution
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  okendon_params_t* input = user;
  double p = input->p;
#if (P4EST_DIM)==2
  double M = 1./(pow((2./(1. - p))*((2./(1. - p))), (1./(1.-p))));
  double r2 = x*x + y*y;
  return M*pow(r2,1/(1-p));
#elif (P4EST_DIM)==3
  double M = 1./(pow((2./(1. - p))*(1 + (2./(1. - p))), (1./(1.-p))));
  double r2 = x*x + y*y + z*z;
  return M*pow(r2,1/(1-p));
#else
  D4EST_ABORT("dim = 2 or 3");
#endif
}


static double
okendon_boundary_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  okendon_analytic_solution(x,y,z,user);
}


static
double okendon_residual_nonlinear_fofu_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  okendon_params_t* input = user;
  double p = input->p;
  return pow(u,p);
}


static
double okendon_jacobian_nonlinear_fofu0_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  okendon_params_t* input = user;
  double p = input->p;
  return p/pow(u,1-p);
}

static
void okendon_apply_jac
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
  okendon_params_t* params = user;
  d4est_poisson_flux_data_t* flux_data = params->flux_data_for_jac;
  
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad);
  
  double* M_fofu0_u_vec = P4EST_ALLOC(double, prob_vecs->local_nodes);

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
           &M_fofu0_u_vec[ed->nodal_stride],
           okendon_jacobian_nonlinear_fofu0_term,
           user,
           NULL,
           NULL
          );
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_fofu0_u_vec, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_fofu0_u_vec);
}

static void
okendon_build_residual
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
  okendon_params_t* params = user;

  d4est_poisson_flux_data_t* flux_data = params->flux_data_for_residual;
  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          flux_data,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad
                         );

  double* M_fof_vec= P4EST_ALLOC(double, prob_vecs->local_nodes);
 
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
           &M_fof_vec[ed->nodal_stride],
           okendon_residual_nonlinear_fofu_term,
           user,
           NULL,
           NULL
          );        
      }
    }

  d4est_linalg_vec_axpy(1.0,
                  M_fof_vec,
                  prob_vecs->Au,
                  prob_vecs->local_nodes);

  P4EST_FREE(M_fof_vec);
}


#endif
