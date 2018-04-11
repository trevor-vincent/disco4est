#ifndef BOYEN_YORK_MODEL_FCNS_H
#define BOYEN_YORK_MODEL_FCNS_H 

typedef struct {

  double a;
  double P;
  double E;
  d4est_poisson_flux_data_t* flux_data_for_residual;
  d4est_poisson_flux_data_t* flux_data_for_jac;
  
} boyen_york_model_params_t;

static
int boyen_york_model_params_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  boyen_york_model_params_t* pconfig = (boyen_york_model_params_t*)user;
  if (d4est_util_match_couple(section,"problem",name,"P")) {
    D4EST_ASSERT(pconfig->P == -1);
    pconfig->P = atof(value);
  }
  else if (d4est_util_match_couple(section,"problem",name,"a")) {
    D4EST_ASSERT(pconfig->a == -1);
    pconfig->a = atof(value);
  }
  else {
    return 0;
  }
  return 1;
}


static
boyen_york_model_params_t
boyen_york_model_params
(
 const char* input_file
)
{
  boyen_york_model_params_t input;
  input.a = -1;
  input.P = -1;
  
  if (ini_parse(input_file, boyen_york_model_params_handler, &input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT("problem", input.a, -1);
  D4EST_CHECK_INPUT("problem", input.P, -1);

  D4EST_ASSERT(input.a > 0);
  D4EST_ASSERT(input.P > 0);
  double P = input.P;
  double a = input.a;
  printf("[PROBLEM]: a = %f\n",input.a);
  printf("[PROBLEM]: P = %f\n",input.P);
  printf("[PROBLEM]: E = %f\n",sqrt(P*P+4*a*a));
  return input;
}

static
double
boyen_york_model_analytic_solution
(
 double x,
 double y,
 double z,
 void* user
)
{
  boyen_york_model_params_t* params = user;
  
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double r3 = r*r2;
  double r4 = r2*r2;

  double a = params->a;
  double P = params->P;
  double E = sqrt(P*P + 4*a*a);
  
  double t1 = 1.;
  double t2 = 2.*E/r;
  double t3 = 6*a*a/(r2);
  double t4 = 2*a*a*E/(r3);
  double t5 = a*a*a*a/(r4);

  return pow((t1 + t2 + t3 + t4 + t5), .25);
}

static
double
boyen_york_model_boundary_fcn
(
 double x,
 double y,
 double z,
 void* user
)
{
  return boyen_york_model_analytic_solution(x,y,z,user);
}

static
double
boyen_york_model_helmholtz_fcn
(
 double x,
 double y,
 double z,
 void* user
){
  boyen_york_model_params_t* params = user;
  double a = params->a;
  double P = params->P;
  double r2 = x*x + y*y + z*z;
  double r4 = r2*r2;
  double a2 = a*a;
  double P2 = P*P;
  
  return .75*(P2/r4)*(1-(a2/r2))*(1-(a2/r2));  
}


static
double
boyen_york_model_residual_nonlinear_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  /* we use u^2^3.5 since it is the same as u^7 if u is positive, which we assume it is */
  double term = boyen_york_model_helmholtz_fcn(x,y,z,user)*(1./pow(u*u,3.5));
  /* if (term > 1000.){ */
    /* printf("term, x, y, z, u = %.25f, %.25f, %.25f, %.25f, %.25f\n", term,x,y,z,u); */
  /* } */
  /* printf("no0nlinear res term = %.25f\n", term); */
  return term;
}

static
double boyen_york_model_jacobian_nonlinear_term
(
 double x,
 double y,
 double z,
 double u,
 void* user
)
{
  /* we use u^2^4 since it is the same as u^8 if u is positive, which we assume it is */
  return -7.*boyen_york_model_helmholtz_fcn(x,y,z,user)*(1./pow(u*u,4));
}


static
void boyen_york_model_apply_jac
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
  boyen_york_model_params_t* params = user;
  d4est_poisson_flux_data_t* flux_data = params->flux_data_for_jac;
  
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
           boyen_york_model_jacobian_nonlinear_term,
           params,
           NULL,
           NULL
          );        
      }
    }
  
  d4est_linalg_vec_axpy(1.0, M_nonlinear_term, prob_vecs->Au, prob_vecs->local_nodes);
  P4EST_FREE(M_nonlinear_term);
}

void
boyen_york_model_build_residual
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
  boyen_york_model_params_t* params = user;


  double* Abc= P4EST_ALLOC(double, prob_vecs->local_nodes);
  d4est_poisson_build_rhs_with_strong_bc(p4est, ghost, ghost_data, d4est_ops, d4est_geom, d4est_quad, prob_vecs, params->flux_data_for_residual, Abc, zero_fcn, NULL);



  d4est_poisson_apply_aij(p4est,
                          ghost,
                          ghost_data,
                          prob_vecs,
                          params->flux_data_for_jac,
                          d4est_ops,
                          d4est_geom,
                          d4est_quad
                         );


  for (int i = 0; i < prob_vecs->local_nodes; i++){
    prob_vecs->Au[i] -= Abc[i];
  }


  P4EST_FREE(Abc);
  
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
           boyen_york_model_residual_nonlinear_term,
           params,
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


double
boyen_york_model_initial_guess
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif
 void* user
)
{
  return 1000.;
}


#endif
