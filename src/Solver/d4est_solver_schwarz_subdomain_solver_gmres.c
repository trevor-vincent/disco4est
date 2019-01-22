#include <pXest.h>
#include <d4est_solver_schwarz_metadata.h>
#include <d4est_solver_schwarz_geometric_data.h>
#include <d4est_solver_schwarz_subdomain_solver_gmres.h>
#include <d4est_solver_schwarz_operators.h>
#include <d4est_solver_schwarz_apply_lhs.h>
#include <d4est_linalg.h>
#include <zlog.h>
#include <ini.h>
#include <time.h>

static
int d4est_solver_schwarz_subdomain_solver_gmres_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_solver_schwarz_subdomain_solver_gmres_t* pconfig = (d4est_solver_schwarz_subdomain_solver_gmres_t*)user;
  const char* input_section = pconfig->input_section;
  
  if (d4est_util_match_couple(section,input_section,name,"subdomain_atol")) {
    D4EST_ASSERT(pconfig->subdomain_atol == -1);
    pconfig->subdomain_atol = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_rtol")) {
    D4EST_ASSERT(pconfig->subdomain_rtol == -1);
    pconfig->subdomain_rtol = atof(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_iter")) {
    D4EST_ASSERT(pconfig->subdomain_iter == -1);
    pconfig->subdomain_iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_inner_iter")) {
    D4EST_ASSERT(pconfig->subdomain_inner_iter == -1);
    pconfig->subdomain_inner_iter = atoi(value);
  }
  else if (d4est_util_match_couple(section,input_section,name,"subdomain_verbose")) {
    D4EST_ASSERT(pconfig->verbose == -1);
    pconfig->verbose = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void
d4est_solver_schwarz_subdomain_solver_gmres_destroy
(
 void* gmres_params
){
  P4EST_FREE(gmres_params);
}


d4est_solver_schwarz_subdomain_solver_gmres_t*
d4est_solver_schwarz_subdomain_solver_gmres_init
(
 p4est_t* p4est,
 const char* input_file,
 const char* input_section
)
{

  d4est_solver_schwarz_subdomain_solver_gmres_t* solver_gmres =
    P4EST_ALLOC(d4est_solver_schwarz_subdomain_solver_gmres_t, 1);
  
  zlog_category_t* c_default = zlog_get_category("d4est_solver_schwarz_subdomain_solver_gmres");
  solver_gmres->subdomain_iter = -1;
  solver_gmres->subdomain_inner_iter = -1;
  solver_gmres->subdomain_rtol = -1;
  solver_gmres->subdomain_atol = -1;
  solver_gmres->verbose = -1;
  solver_gmres->input_section = input_section;

  if(
     ini_parse(input_file,
               d4est_solver_schwarz_subdomain_solver_gmres_input_handler,
               solver_gmres) < 0
  ){
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, solver_gmres->subdomain_iter, -1);
  D4EST_CHECK_INPUT(input_section, solver_gmres->subdomain_inner_iter, -1);
  D4EST_CHECK_INPUT(input_section, solver_gmres->subdomain_rtol, -1);
  D4EST_CHECK_INPUT(input_section, solver_gmres->subdomain_atol, -1);
  D4EST_CHECK_INPUT(input_section, solver_gmres->verbose, -1);

  if (solver_gmres->subdomain_iter <= 0 ||
      solver_gmres->subdomain_rtol <= 0 ||
      solver_gmres->subdomain_atol <= 0 
     ){
    D4EST_ABORT("Some subdomain solver options are <= 0");
  }

  return solver_gmres;
}


static
double **dmatrix ( int nrl, int nrh, int ncl, int nch )
{
  int i;
  double **m;
  int nrow = nrh - nrl + 1;
  int ncol = nch - ncl + 1;

  m= malloc ((sizeof(double*))*(nrow + 1));
  
  if ( !m ) 
  {
    printf ( "\n" );
    printf ( "DMATRIX - Fatal error!\n" );
    printf ( "  Failure allocating pointers to rows.\n");
    D4EST_ABORT("");
  }
  m = m + 1;
  m = m - nrl;
  m[nrl] = ( double * ) malloc ( (size_t) ( ( nrow * ncol + 1 ) * sizeof ( double ) ) );

  if ( !m[nrl] ) 
  {
    printf ( "\n" );
    printf ( "DMATRIX - Fatal error!\n" );
    printf ( "  Failure allocating rows.\n");
    D4EST_ABORT("");
  }
  m[nrl] = m[nrl] + 1;
  m[nrl] = m[nrl] - ncl;

  for ( i = nrl + 1; i <= nrh; i++ ) 
  { 
    m[i] = m[i-1] + ncol;
  }
  return m;
}

static
void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch )
{
  free ( ( char* ) ( m[nrl] + ncl - 1 ) );
  free ( ( char* ) ( m + nrl - 1 ) );
  return;
}

static
void mult_givens ( double c, double s, int k, double *g )
{
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k+1];
  g2 = s * g[k] + c * g[k+1];

  g[k]   = g1;
  g[k+1] = g2;

  return;
}

d4est_solver_schwarz_subdomain_solver_info_t
d4est_solver_schwarz_subdomain_solver_gmres
(
 p4est_t* p4est,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 d4est_mesh_data_t* d4est_factors,
 d4est_ghost_t* ghost,
 d4est_solver_schwarz_operators_t* schwarz_ops,
 d4est_solver_schwarz_metadata_t* schwarz_metadata,
 d4est_solver_schwarz_geometric_data_t* schwarz_geometric_data,
 d4est_solver_schwarz_apply_lhs_t* apply_lhs,
 double* du_restricted_field_over_subdomain,
 double* rhs_restricted_field_over_subdomain,
 int subdomain,
 void* params,
 int debug_amr,
 int debug_ksp,
 int debug_mg
)
{
  clock_t begin = clock();
  zlog_category_t* c_default = zlog_get_category("d4est_schwarz_subdomain");

  d4est_solver_schwarz_subdomain_solver_gmres_t* gmres_params
    = params;

  d4est_solver_schwarz_subdomain_metadata_t* sub_data =
    &schwarz_metadata->subdomain_metadata[subdomain];
  
  int itr_max = gmres_params->subdomain_iter;
  double tol_abs = gmres_params->subdomain_atol;
  double tol_rel = gmres_params->subdomain_rtol;
  int mr = gmres_params->subdomain_inner_iter;
  
  int nodes
    = schwarz_metadata->subdomain_metadata[subdomain].restricted_nodal_size;

  if (nodes < mr){
    D4EST_ABORT("nodes < mr in schwarz_subdomain_solver_gmres");
  }
  
  int n = nodes;
  double *x = du_restricted_field_over_subdomain;
  /* double* rhs = vecs->rhs; */
  
  double av;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double **h;
  double htmp;
  int i;
  int itr;
  int itr_used;
  int j;
  int k;
  int k_copy;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  double *s;
  double **v;
  double *y;

  itr_used = 0;

  c = ( double * ) malloc ( mr * sizeof ( double ) );
  g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
  h = dmatrix(0,mr,0,mr-1);
  r = ( double * ) malloc ( n * sizeof ( double ) );
  s = ( double * ) malloc ( mr * sizeof ( double ) );
  v = dmatrix(0,mr,0,n-1);
  y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );

  double* Au_restricted_field_over_subdomain = P4EST_ALLOC(double, nodes);
  
  for ( itr = 0; itr < itr_max; itr++ ) 
  {

    apply_lhs->apply_lhs_fcn
      (
       p4est,
       schwarz_ops->d4est_ops,
       d4est_geom,
       d4est_quad,
       d4est_factors,
       ghost,
       schwarz_ops,
       schwarz_metadata,
       schwarz_geometric_data,
       subdomain,
       du_restricted_field_over_subdomain,
       Au_restricted_field_over_subdomain,
       apply_lhs->apply_lhs_ctx
      );

    for ( i = 0; i < n; i++ )
    {
      r[i] = rhs_restricted_field_over_subdomain[i] - Au_restricted_field_over_subdomain[i];
    }

    rho = sqrt ( d4est_linalg_vec_dot ( r, r, n ) );
    
    if (gmres_params->verbose >= 2){
      zlog_info(c_default, "rank subdomain core_tree iters r2 %d %d %d %d %.15f", p4est->mpirank, subdomain, sub_data->core_tree, itr, rho);
    }
      
    if ( itr == 0 )
    {
      rho_tol = rho * tol_rel;
    }

    for ( i = 0; i < n; i++ )
    {
      v[0][i] = r[i] / rho;
    }

    g[0] = rho;
    for ( i = 1; i < mr + 1; i++ )
    {
      g[i] = 0.0;
    }

    for ( i = 0; i < mr + 1; i++ )
    {
      for ( j = 0; j < mr; j++ ) 
      {
        h[i][j] = 0.0;
      }
    }

    for ( k = 0; k < mr; k++ )
    {
      k_copy = k;

      apply_lhs->apply_lhs_fcn
        (
         p4est,
         schwarz_ops->d4est_ops,
         d4est_geom,
         d4est_quad,
         d4est_factors,
         ghost,
         schwarz_ops,
         schwarz_metadata,
         schwarz_geometric_data,
         subdomain,
         v[k],
         v[k+1],
         apply_lhs->apply_lhs_ctx
        );
      

      /* av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) ); */
      av = sqrt ( d4est_linalg_vec_dot (v[k+1], v[k+1], n ) );

      for ( j = 0; j < k+1; j++ )
      {
        /* h[j][k] = r8vec_dot ( n, v[k+1], v[j] ); */
        h[j][k] = d4est_linalg_vec_dot (v[k+1], v[j], n);
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
        }
      }

      /* h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) ); */
      h[k+1][k] = sqrt ( d4est_linalg_vec_dot ( v[k+1], v[k+1], n ) );

      if ( ( av + delta * h[k+1][k] ) == av )
      {
        for ( j = 0; j < k+1; j++ )
        {
          /* htmp = r8vec_dot ( n, v[k+1], v[j] ); */
          htmp = d4est_linalg_vec_dot( v[k+1], v[j], n );
          h[j][k] = h[j][k] + htmp;
          for ( i = 0; i < n; i++ ) 
          {
            v[k+1][i] = v[k+1][i] - htmp * v[j][i];
          }
        }
        /* h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) ); */
        h[k+1][k] = sqrt ( d4est_linalg_vec_dot ( v[k+1], v[k+1], n ) );
      }

      if ( h[k+1][k] != 0.0 )
      {
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] / h[k+1][k];
        }
      }

      if ( 0 < k )
      {
        for ( i = 0; i < k + 2; i++ )
        {
          y[i] = h[i][k];
        }
        for ( j = 0; j < k; j++ ) 
        {
          mult_givens ( c[j], s[j], j, y );
        }
        for ( i = 0; i < k + 2; i++ ) 
        {
          h[i][k] = y[i];
        }
      }

      mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );
      c[k] = h[k][k] / mu;
      s[k] = -h[k+1][k] / mu;
      h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
      h[k+1][k] = 0.0;
      mult_givens ( c[k], s[k], k, g );

      rho = fabs ( g[k+1] );

      itr_used = itr_used + 1;
      
      if (gmres_params->verbose == 3){
        zlog_info(c_default, "restart rank subdomain core_tree iters r2 %d %d %d %d %.15f", p4est->mpirank, subdomain, sub_data->core_tree, k, rho);
      }
      
      if ( rho <= rho_tol && rho <= tol_abs )
      {
        break;
      }
    }

    k = k_copy;

    y[k] = g[k] / h[k][k];
    for ( i = k - 1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = i+1; j < k + 1; j++ ) 
      {
        y[i] = y[i] - h[i][j] * y[j];
      }
      y[i] = y[i] / h[i][i];
    }

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k + 1; j++ )
      {
        x[i] = x[i] + v[j][i] * y[j];
      }
    }

    if ( rho <= rho_tol && rho <= tol_abs ) 
    {
      break;
    }
  }

  if (gmres_params->verbose >= 1){
    clock_t end = clock();
    double time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
    zlog_info(c_default, "rank subdomain core_tree iters r time %d %d %d %d %.15f %f", p4est->mpirank, subdomain, sub_data->core_tree, itr_used, rho, time_spent);
    /* printf( "rank subdomain core_tree iters r time %d %d %d %d %.15f %f", p4est->mpirank, subdomain, sub_data->core_tree, itr_used, rho, time_spent); */
  }
  
  free ( c );
  free ( g );
  free_dmatrix(h,0,mr,0,mr-1);
  free ( r );
  free ( s );
  free_dmatrix(v,0,mr,0,n-1);
  free ( y );
  P4EST_FREE(Au_restricted_field_over_subdomain);


  return (d4est_solver_schwarz_subdomain_solver_info_t){.final_iter = itr_used,
      .final_res = rho};
}
