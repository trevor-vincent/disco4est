#define _GNU_SOURCE
#include <d4est_quadrature.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_quadrature_lobatto.h>
#include <d4est_quadrature_compactified.h>
#include <util.h>
#include <d4est_linalg.h>
#include <ini.h>

typedef struct {

  const char* input_section;
  char* name;

} d4est_quadrature_input_t;


static
int d4est_quadrature_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_quadrature_input_t* pconfig = (d4est_quadrature_input_t*)user;

  if (util_match_couple(section,pconfig->input_section,name,"name")) {
    mpi_assert(pconfig->name == NULL);
    D4EST_ASPRINTF(pconfig->name,"%s",value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
d4est_quadrature_input_t
d4est_quadrature_input
(
 const char* input_file,
 const char* input_section,
 const char* printf_prefix
)
{
  
  d4est_quadrature_input_t input;
  input.input_section = input_section;
  input.name = NULL;
  
  if (ini_parse(input_file, d4est_quadrature_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input.name, NULL); 
  printf("%s: Loading %s quadrature\n", printf_prefix, input.name);

  return input;
}


d4est_quadrature_t*
d4est_quadrature_new
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix
)
{
  
  d4est_quadrature_input_t input = d4est_quadrature_input(input_file, input_section, printf_prefix);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);

  if (util_match(input.name,"legendre")) {
    d4est_quad->quad_type = QUAD_GAUSS_LEGENDRE;
    d4est_quadrature_legendre_new(d4est_quad, d4est_geom, input_file);
  }
  else if (util_match(input.name,"lobatto")) {
    d4est_quad->quad_type = QUAD_GAUSS_LEGENDRE_LOBATTO;
    d4est_quadrature_lobatto_new(d4est_quad, d4est_geom, input_file);
  }
  else if (util_match(input.name,"legendre_compactified")){
    d4est_quad->quad_type = QUAD_GAUSS_LEGENDRE_COMPACTIFIED;
    d4est_quadrature_compactified_new(p4est, d4est_ops, d4est_geom, d4est_quad, input_file);
  }
  else if (util_match(input.name,"none")){
  }
  else {
    printf("[D4EST_ERROR]: You tried to use %s quadrature\n", input.name);
    mpi_abort("[D4EST_ERROR]: this quadrature is currently not supported");
  }

  free(input.name);

  return d4est_quad;
}

void
d4est_quadrature_reinit(p4est_t* p4est,
                        d4est_operators_t* d4est_ops,
                        d4est_geometry_t* d4est_geom,
                        d4est_quadrature_t* d4est_quad)
{
  if (d4est_quad->user_reinit != NULL)
    d4est_quad->user_reinit(p4est, d4est_ops, d4est_geom, d4est_quad);
}
  
void
d4est_quadrature_destroy(p4est_t* p4est,
                        d4est_operators_t* d4est_ops,
                        d4est_geometry_t* d4est_geom,
                        d4est_quadrature_t* d4est_quad){
  if (d4est_quad->user_destroy != NULL)
    d4est_quad->user_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  P4EST_FREE(d4est_quad);
}


void d4est_quadrature_apply_galerkin_integral
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geometry,
 d4est_quadrature_t* d4est_quadrature,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* in_quad,
 int deg_lobatto,
 double* jac_quad,
 int deg_quad,
 double* out
){
  double* quad_weights [(P4EST_DIM)];
  double* interp_lobatto_to_quad_trans [(P4EST_DIM)];

  int dim = (object_type == D4EST_QUAD_FACE) ? ((P4EST_DIM)-1) : (P4EST_DIM);
  int nodes_lobatto = deg_lobatto + 1;
  int nodes_quad = deg_quad + 1;
  int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad);
  double* w_j_in_quad = P4EST_ALLOC(double, volume_nodes_quad);
  
  for (int dir = 0; dir < dim; dir++){
     quad_weights[dir] = d4est_quadrature->get_weights
                         (
                          d4est_ops,
                          d4est_geometry,
                          d4est_quadrature,
                          object,
                          object_type,
                          integrand_type,
                          deg_quad,
                          dir
                         );
     
     interp_lobatto_to_quad_trans[dir] = d4est_quadrature->get_interp_trans
                                       (
                                        d4est_ops,
                                        d4est_geometry,
                                        d4est_quadrature,
                                        object,
                                        object_type,
                                        integrand_type,
                                        deg_lobatto,
                                        deg_quad,
                                        dir
                                       );
  }

  if (dim == 3){
    d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_xy(quad_weights[2], quad_weights[1], quad_weights[0], jac_quad, in_quad, nodes_quad, nodes_quad, nodes_quad, w_j_in_quad);
    d4est_linalg_kron_A1A2A3x_nonsqr(out, interp_lobatto_to_quad_trans[2], interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], w_j_in_quad,
                               nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad);
  }
  else if (dim == 2) {
    d4est_linalg_kron_vec1_o_vec2_dot_xy(quad_weights[1], quad_weights[0], jac_quad, in_quad, nodes_quad, nodes_quad, w_j_in_quad);
    d4est_linalg_kron_A1A2x_nonsqr(out, interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], w_j_in_quad, nodes_lobatto, nodes_quad,
                             nodes_lobatto, nodes_quad);
  }
  else if (dim == 1){
    d4est_linalg_kron_vec_dot_xy(quad_weights[0], jac_quad, in_quad, nodes_quad, w_j_in_quad);
    d4est_linalg_matvec_plus_vec(1.0, interp_lobatto_to_quad_trans[0], w_j_in_quad, 0., out, nodes_lobatto, nodes_quad);
  }
  
  else {    
    P4EST_FREE(w_j_in_quad);
    P4EST_FREE(in_quad);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }  
  P4EST_FREE(w_j_in_quad);
}


/* void d4est_quadrature_compute_mass_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_quadrature_t* d4est_quadrature, */
/*  d4est_geometry_t* d4est_geometry, */
/*  p4est_qcoord_t q, */
/*  p4est_qcoord_t dq, */
/*  p4est_topidx_t tree, */
/*  int deg_lobatto, */
/*  int deg_quad, */
/*  int dim, */
/*  double* jac_quad, */
/*  double* M */
/* ) */
/* { */
/*   int volume_nodes_lobatto = d4est_operators_get_nodes(dim,deg_lobatto); */
  
/*   double* u = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */
/*   double* Mu = P4EST_ALLOC_ZERO(double, volume_nodes_lobatto); */

/*   for (int i = 0; i < volume_nodes_lobatto; i++){ */
/*     u[i] = 1.; */
/*     d4est_quad_apply_mass_matrix */
/*       ( */
/*        d4est_ops, */
/*        d4est_quadrature, */
/*        d4est_geometry, */
/*        q, */
/*        dq, */
/*        tree, */
/*        u, */
/*        deg_lobatto, */
/*        jac_quad, */
/*        deg_quad, */
/*        dim, */
/*        Mu */
/*       ); */
/*     d4est_linalg_set_column(M, Mu, i, volume_nodes_lobatto, volume_nodes_lobatto); */
/*     u[i] = 0.; */
/*   } */

/*   P4EST_FREE(Mu); */
/*   P4EST_FREE(u); */
/* } */



void d4est_quadrature_apply_stiffness_matrix
(
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quadrature,
 d4est_geometry_t* d4est_geometry,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* in,
 int deg_lobatto,
 double* jac_quad,
 double* rst_xyz [(P4EST_DIM)][(P4EST_DIM)], /* must be padded with NULL if you want lower dim */
 int deg_quad,
 double* out
){
  mpi_assert (object_type == D4EST_QUAD_VOLUME);
  /* for now assert this, can get rid of by p-prolonging then p-restricting */


  double* quad_weights [(P4EST_DIM)];
  double* interp_lobatto_to_quad_trans [(P4EST_DIM)];
  double* interp_lobatto_to_quad [(P4EST_DIM)]; 

  int dim = (P4EST_DIM);
  int nodes_lobatto = deg_lobatto + 1;
  int nodes_quad = deg_quad + 1;

    int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad);
    int volume_nodes_lobatto = d4est_operators_get_nodes(dim, deg_lobatto);
  
  for (int dir = 0; dir < dim; dir++){
    interp_lobatto_to_quad[dir] = d4est_quadrature->get_interp
                                (
                                 d4est_ops,
                                 d4est_geometry,
                                 d4est_quadrature,
                                 object,
                                 object_type,
                                 integrand_type,
                                 deg_lobatto,
                                 deg_quad,
                                 dir
                                );

   interp_lobatto_to_quad_trans[dir] = d4est_quadrature->get_interp_trans
                                       (
                                        d4est_ops,
                                        d4est_geometry,
                                        d4est_quadrature,
                                        object,
                                        object_type,
                                        integrand_type,
                                        deg_lobatto,
                                        deg_quad,
                                        dir
                                       );


   quad_weights[dir] = d4est_quadrature->get_weights
                                       (
                                        d4est_ops,
                                        d4est_geometry,
                                        d4est_quadrature,
                                        object,
                                        object_type,
                                        integrand_type,
                                        deg_quad,
                                        dir
                                       );
  }
  
  double* Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad);
  double* W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_quad);
  double* VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto);
  double* DTlp_VT_W_V_Dl_in = P4EST_ALLOC(double, volume_nodes_lobatto);
  
  d4est_linalg_fill_vec(out, 0., volume_nodes_lobatto);
  
  for (int k = 0; k < dim; k++){
    for (int lp = 0; lp < dim; lp++){
      for (int l = 0; l < dim; l++){
        /* Apply Dbar in lth direction */
        d4est_operators_apply_Dij(d4est_ops, in, dim, deg_lobatto, l, Dl_in);

        if ((P4EST_DIM) == 3) {

              d4est_linalg_kron_A1A2A3x_nonsqr(V_Dl_in, interp_lobatto_to_quad[2], interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], Dl_in,
                               nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto);
              
              d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_wxyz(quad_weights[2],quad_weights[1],quad_weights[0], jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, nodes_quad, nodes_quad, W_V_Dl_in);
          d4est_linalg_kron_A1A2A3x_nonsqr(VT_W_V_Dl_in, interp_lobatto_to_quad_trans[2], interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], W_V_Dl_in,
                                     nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad);
        }
        else if ((P4EST_DIM) == 2) {
          d4est_linalg_kron_A1A2x_nonsqr(V_Dl_in, interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], Dl_in, nodes_quad, nodes_lobatto,
                            nodes_quad, nodes_lobatto);
          
          d4est_linalg_kron_vec1_o_vec2_dot_wxyz(quad_weights[1], quad_weights[0], jac_quad, rst_xyz[l][k], rst_xyz[lp][k], V_Dl_in, nodes_quad, nodes_quad, W_V_Dl_in);
          d4est_linalg_kron_A1A2x_nonsqr(VT_W_V_Dl_in, interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], W_V_Dl_in, nodes_lobatto, nodes_quad,
                                   nodes_lobatto, nodes_quad);
        }
        else {
          P4EST_FREE(DTlp_VT_W_V_Dl_in);
          P4EST_FREE(VT_W_V_Dl_in);
          P4EST_FREE(W_V_Dl_in);
          P4EST_FREE(V_Dl_in);
          P4EST_FREE(Dl_in);
          mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
        }

        d4est_operators_apply_Dij_transpose(d4est_ops, VT_W_V_Dl_in, dim, deg_lobatto, lp, DTlp_VT_W_V_Dl_in);
        d4est_linalg_vec_axpy(1., DTlp_VT_W_V_Dl_in, out, volume_nodes_lobatto);
      }
    }
  }
  
  P4EST_FREE(DTlp_VT_W_V_Dl_in);
  P4EST_FREE(VT_W_V_Dl_in);
  P4EST_FREE(W_V_Dl_in);
  P4EST_FREE(V_Dl_in);
  P4EST_FREE(Dl_in);
}


void d4est_quadrature_apply_mass_matrix
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geometry,
 d4est_quadrature_t* d4est_quadrature,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* in,
 int deg_lobatto,
 double* jac_quad,
 int deg_quad,
 double* out
){

  mpi_assert(object_type == D4EST_QUAD_VOLUME);

  double* quad_weights [(P4EST_DIM)];
  double* interp_lobatto_to_quad_trans [(P4EST_DIM)];
  double* interp_lobatto_to_quad [(P4EST_DIM)]; 

  int dim = (P4EST_DIM);
  int nodes_lobatto = deg_lobatto+1;
  int nodes_quad = deg_quad+1;
  int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad);
  
  double* in_quad = P4EST_ALLOC(double, volume_nodes_quad);
  double* w_j_in_quad = P4EST_ALLOC(double, volume_nodes_quad);

  
  for (int dir = 0; dir < dim; dir++){
    interp_lobatto_to_quad[dir] = d4est_quadrature->get_interp
                                  (
                                   d4est_ops,
                                   d4est_geometry,
                                   d4est_quadrature,
                                   object,
                                   object_type,
                                   integrand_type,
                                   deg_lobatto,
                                   deg_quad,
                                   dir
                                  );

   interp_lobatto_to_quad_trans[dir] = d4est_quadrature->get_interp_trans
                                       (
                                        d4est_ops,
                                        d4est_geometry,
                                        d4est_quadrature,
                                        object,
                                        object_type,
                                        integrand_type,
                                        deg_lobatto,
                                        deg_quad,
                                        dir
                                       );


   quad_weights[dir] = d4est_quadrature->get_weights
                                       (
                                        d4est_ops,
                                        d4est_geometry,
                                        d4est_quadrature,
                                        object,
                                        object_type,
                                        integrand_type,
                                        deg_quad,
                                        dir
                                       );
  }
  
  if (dim == 2) {
    d4est_linalg_kron_A1A2x_nonsqr(in_quad, interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], in, nodes_quad, nodes_lobatto,
                             nodes_quad, nodes_lobatto);
    d4est_linalg_kron_vec1_o_vec2_dot_xy(quad_weights[1], quad_weights[0], jac_quad, in_quad, nodes_quad, nodes_quad, w_j_in_quad);
    d4est_linalg_kron_A1A2x_nonsqr(out, interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], w_j_in_quad, nodes_lobatto, nodes_quad,
                             nodes_lobatto, nodes_quad);
  }
  else if (dim == 3){
    d4est_linalg_kron_A1A2A3x_nonsqr(in_quad, interp_lobatto_to_quad[2], interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], in,
                               nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto);
    d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_xy(quad_weights[2], quad_weights[1], quad_weights[0], jac_quad, in_quad, nodes_quad, nodes_quad, nodes_quad, w_j_in_quad);
    d4est_linalg_kron_A1A2A3x_nonsqr(out, interp_lobatto_to_quad_trans[2], interp_lobatto_to_quad_trans[1], interp_lobatto_to_quad_trans[0], w_j_in_quad,
                               nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad);
  }
  else {
    P4EST_FREE(w_j_in_quad);
    P4EST_FREE(in_quad);
    mpi_abort("ERROR: Apply mass matrix ref space, wrong dimension.");
  }
  P4EST_FREE(w_j_in_quad);
  P4EST_FREE(in_quad);
}

/* void d4est_quadrature_form_fofufofvlilj_matrix */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_quadrature_t* d4est_quadrature, */
/*  d4est_geometry_t* d4est_geometry, */
/*  p4est_qcoord_t q, */
/*  p4est_qcoord_t dq, */
/*  p4est_topidx_t tree, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  double* jac_quad, */
/*  int deg_quad, */
/*  int dim, */
/*  double* out, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL) */
/*     fofu_fcn = identity_fcn; */
/*   if (fofv_fcn == NULL) */
/*     fofv_fcn = identity_fcn; */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   u, */
/*                   deg_lobatto, */
/*                   u_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   v, */
/*                   deg_lobatto, */
/*                   v_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
  
/*   double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv_jac[i] = jac_quad[i]; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */
/*   } */
  
/*   d4est_quadrature_compute_mass_matrix */
/*     ( */
/*      d4est_ops, */
/*      d4est_quadrature, */
/*      d4est_geometry, */
/*      q, */
/*      dq, */
/*      tree, */
/*      deg_lobatto, */
/*      deg_quad, */
/*      dim, */
/*      fofu_fofv_jac, */
/*      out */
/*     ); */
 
  
/*   if (u != NULL){ */
/*     P4EST_FREE(u_quad); */
/*   } */
/*   if (v != NULL){ */
/*     P4EST_FREE(v_quad); */
/*   } */
  
/* } */


/* void d4est_quadrature_apply_fofufofvlilj */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_quadrature_t* d4est_quadrature, */
/*  d4est_geometry_t* d4est_geometry, */
/*  d4est_geometry_object_t d4est_object, */
/*  double* vec, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  double* jac_quad, */
/*  int deg_quad, */
/*  int dim, */
/*  double* out, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL) */
/*     fofu_fcn = identity_fcn; */
/*   if (fofv_fcn == NULL) */
/*     fofv_fcn = identity_fcn; */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   u, */
/*                   deg_lobatto, */
/*                   u_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   v, */
/*                   deg_lobatto, */
/*                   v_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
  
/*   double* fofu_fofv_jac = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv_jac[i] = jac_quad[i]; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv_jac[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */
/*   } */
  
/*   d4est_quadrature_apply_mass_matrix */
/*     ( */
/*      d4est_ops, */
/*      d4est_quadrature, */
/*      d4est_geometry, */
/*      q, */
/*      dq, */
/*      tree, */
/*      vec, */
/*      deg_lobatto, */
/*      fofu_fofv_jac, */
/*      deg_quad, */
/*      out */
/*     ); */
  
/*   if (u != NULL){ */
/*     P4EST_FREE(u_quad); */
/*   } */
/*   if (v != NULL){ */
/*     P4EST_FREE(v_quad); */
/*   } */
  
/* } */

/* void d4est_quadrature_apply_fofufofvlj */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_quadrature_t* d4est_quadrature, */
/*  d4est_geometry_t* d4est_geometry, */
/*  d4est_geometry_object_t d4est_object, */
/*  double* u, */
/*  double* v, */
/*  int deg_lobatto, */
/*  double* jac_quad, */
/*  double* xyz_quad [(P4EST_DIM)], */
/*  int deg_quad, */
/*  int dim, */
/*  double* out, */
/*  grid_fcn_ext_t fofu_fcn, */
/*  void* fofu_ctx, */
/*  grid_fcn_ext_t fofv_fcn, */
/*  void* fofv_ctx */
/* ) */
/* { */
/*   if (fofu_fcn == NULL){ */
/*     fofu_fcn = identity_fcn; */
/*   } */
/*   if (fofv_fcn == NULL){ */
/*     fofv_fcn = identity_fcn; */
/*   } */
  
/*   double* u_quad = NULL; */
/*   double* v_quad = NULL; */

/*   int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad); */
  
/*   if (u != NULL){ */
/*     u_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   u, */
/*                   deg_lobatto, */
/*                   u_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
/*   if (v != NULL){ */
/*     v_quad = P4EST_ALLOC(double, volume_nodes_quad); */
/*     d4est_quadrature_interpolate( */
/*                   d4est_ops, */
/*                   d4est_quadrature, */
/*                   d4est_geometry, */
/*                   q, */
/*                   dq, */
/*                   tree, */
/*                   v, */
/*                   deg_lobatto, */
/*                   v_quad, */
/*                   deg_quad, */
/*                   dim */
/*                  ); */
/*   } */
  
/*   double* fofu_fofv = P4EST_ALLOC(double, volume_nodes_quad); */
/*   for (int i = 0; i < volume_nodes_quad; i++){ */
/*     fofu_fofv[i] = 1.0; */
/*     if (u != NULL || fofu_fcn != identity_fcn){ */
/*       fofu_fofv[i] *= fofu_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (u != NULL) ? u_quad[i] : 0, */
/*                                    fofu_ctx); */
/*       /\* printf("xyz_quad[0][i], xyz_quad[1][i], xyz_quad[2][i], fofu_fofv[i] = %f,%f,%f,%f\n",xyz_quad[0][i], xyz_quad[1][i], xyz_quad[2][i], fofu_fofv[i]); *\/ */
/*     } */
/*     if (v != NULL || fofv_fcn != identity_fcn){ */
/*       fofu_fofv[i] *= fofv_fcn(xyz_quad[0][i], */
/*                                    xyz_quad[1][i], */
/* #if (P4EST_DIM)==3 */
/*                                    xyz_quad[2][i], */
/* #endif */
/*                                    (v != NULL) ? v_quad[i] : 0, */
/*                                    fofv_ctx); */
/*     } */

/*   } */
 
/*   d4est_quadrature_apply_galerkin_quadral */
/*     ( */
/*      d4est_ops, */
/*      d4est_quadrature, */
/*      d4est_geometry, */
/*      q, */
/*      dq, */
/*      fofu_fofv, */
/*      deg_lobatto, */
/*      jac_quad, */
/*      deg_quad, */
/*      dim, */
/*      out */
/*     ); */

/*   P4EST_FREE(u_quad); */
/*   P4EST_FREE(v_quad); */
/*   P4EST_FREE(fofu_fofv); */
  
/* } */


d4est_rst_t
d4est_quadrature_get_rst_points
(
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quadrature,
 d4est_geometry_t* d4est_geometry,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 int degree
)
{
  int dim = (object_type == D4EST_QUAD_FACE) ? ((P4EST_DIM)-1) : (P4EST_DIM);
  d4est_rst_t rst_points;
  rst_points.r = d4est_quadrature->get_rst(d4est_ops, d4est_geometry, d4est_quadrature,object,object_type, integrand_type, degree, 0);
  rst_points.s = d4est_quadrature->get_rst(d4est_ops, d4est_geometry, d4est_quadrature,  object, object_type, integrand_type, degree, 1);
  if (dim == 3)
    rst_points.t = d4est_quadrature->get_rst(d4est_ops, d4est_geometry, d4est_quadrature,  object, object_type, integrand_type, degree, 2);
  else
    rst_points.t = NULL;
  return rst_points;
}


void
d4est_quadrature_interpolate
(
 d4est_operators_t* d4est_ops,
 d4est_quadrature_t* d4est_quadrature,
 d4est_geometry_t* d4est_geometry,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type,
 double* u_lobatto_in,
 int deg_lobatto,
 double* u_quad_out,
 int deg_quad
)
{
  /* mpi_assert(object_type == D4EST_QUAD_FACE || object_type == D4EST_QUAD_VOLUME); */
  double* interp_lobatto_to_quad [(P4EST_DIM)]; 
  int dim = (object_type == D4EST_QUAD_FACE) ? ((P4EST_DIM)-1) : (P4EST_DIM);
  int nodes_lobatto = deg_lobatto + 1;
  int nodes_quad = deg_quad + 1;

  for (int dir = 0; dir < dim; dir++){
    interp_lobatto_to_quad[dir] = d4est_quadrature->get_interp
                                (
                                 d4est_ops,
                                 d4est_geometry,
                                 d4est_quadrature,
                                 object,
                                 object_type,
                                 integrand_type, 
                                 deg_lobatto,
                                 deg_quad,
                                 dir
                                );
  }
  
  if (dim == 3){
    d4est_linalg_kron_A1A2A3x_nonsqr(u_quad_out, interp_lobatto_to_quad[2], interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], u_lobatto_in,
                               nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto, nodes_quad, nodes_lobatto);
  }
  else if (dim == 2){
    d4est_linalg_kron_A1A2x_nonsqr(u_quad_out, interp_lobatto_to_quad[1], interp_lobatto_to_quad[0], u_lobatto_in, nodes_quad, nodes_lobatto,
                             nodes_quad, nodes_lobatto);
  }
  else if (dim == 1){
    d4est_linalg_matvec_plus_vec(1.0,interp_lobatto_to_quad[0], u_lobatto_in, 0., u_quad_out, nodes_quad, nodes_lobatto);
  }
  else {
    mpi_abort("[D4EST_ERROR]: dim = 1, 2 or 3\n");
  }  
}

/* Computes \sum w_i u_i v_i */
double
d4est_quadrature_innerproduct
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 void* object,
 d4est_quadrature_object_type_t object_type,
 d4est_quadrature_integrand_type_t integrand_type, 
 double* u,
 double* v,
 double* jac_quad,
 int deg_quad
)
{

  double* quad_weights [(P4EST_DIM)];

  int dim = (object_type == D4EST_QUAD_FACE) ? ((P4EST_DIM)-1) : (P4EST_DIM);
  int nodes_quad = deg_quad+1;  
  
  for (int dir = 0; dir < dim; dir++)
    {
      quad_weights[dir] = d4est_quad->get_weights
                          (
                           d4est_ops,
                           d4est_geom,
                           d4est_quad,
                           object,
                           object_type,
                           integrand_type, 
                           deg_quad,
                           dir
                          );
    }
  
  mpi_assert(u != NULL);
  int volume_nodes_quad = d4est_operators_get_nodes(dim, deg_quad);
  double* wuv = P4EST_ALLOC(double, volume_nodes_quad);
  double wdotuv = 0.;

  if (dim == 3){
    if (v != NULL){
      d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_xy(quad_weights[2], quad_weights[1], quad_weights[0], u, v, nodes_quad, nodes_quad, nodes_quad, wuv);
    }
    else {
      d4est_linalg_kron_vec1_o_vec2_o_vec3_dot_x(quad_weights[2], quad_weights[1], quad_weights[0], u,nodes_quad,nodes_quad,nodes_quad,wuv);
    }
  }
  else if (dim == 2){
    if (v != NULL){
      d4est_linalg_kron_vec1_o_vec2_dot_xy(quad_weights[1], quad_weights[0],  u, v, nodes_quad, nodes_quad, wuv);
    }
    else {
      d4est_linalg_kron_vec1_o_vec2_dot_x(quad_weights[1], quad_weights[0], u, nodes_quad, nodes_quad, wuv);
    }
  }
  else if (dim == 1){
    if (v != NULL){
      d4est_linalg_kron_vec_dot_xy(quad_weights[0], u, v, deg_quad + 1, wuv);
    }
    else {
      d4est_linalg_kron_vec_dot_x(quad_weights[0], u, deg_quad + 1, wuv);
    }
  }
  else {
    mpi_abort("dim must be 1,2,3");
  }

  if (jac_quad == NULL){
    for (int i = 0; i < volume_nodes_quad; i++){
      wdotuv += wuv[i];
    }
  }
  else {
    for (int i = 0; i < volume_nodes_quad; i++){
      wdotuv += jac_quad[i]*wuv[i];
    }
  }
  P4EST_FREE(wuv);

  return wdotuv;
}
