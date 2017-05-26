#include <util.h>
#include <ini.h>
#include <d4est_geometry.h>
#include <d4est_operators.h>
#include <linalg.h>
#include <util.h>


#if (P4EST_DIM)==3
#include <d4est_geometry_cubed_sphere.h>
/* #include <d4est_geometry_shell.h> */
#endif

#if (P4EST_DIM)==2
#include <d4est_geometry_disk.h>
#endif

typedef struct {

  const char* input_section;
  char* name;
  geometric_quantity_compute_method_t DX_compute_method;
  geometric_quantity_compute_method_t JAC_compute_method;
  
} d4est_geometry_input_t;


static
int d4est_geometry_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_input_t* pconfig = (d4est_geometry_input_t*)user;

  if (util_match_couple(section,pconfig->input_section,name,"name")) {
    mpi_assert(pconfig->name == NULL);
    D4EST_ASPRINTF(pconfig->name,"%s",value);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"DX_compute_method")) {
    if(util_match(value, "numerical")){
      pconfig->DX_compute_method = GEOM_COMPUTE_NUMERICAL;
    }
    else if (util_match(value, "analytic")){
      pconfig->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
    }
    else {
      printf("[D4EST_ERROR]: You tried to use %s as a mapping type\n", value);
      mpi_abort("[D4EST_ERROR]: This mapping is not supported");
    }
  } 
  else if (util_match_couple(section,pconfig->input_section,name,"JAC_compute_method")) {
    if(util_match(value, "numerical")){
      pconfig->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;
    }
    else if (util_match(value, "analytic")){
      pconfig->JAC_compute_method = GEOM_COMPUTE_ANALYTIC;
    }
    else {
      printf("[D4EST_ERROR]: You tried to use %s as a mapping type\n", value);
      mpi_abort("[D4EST_ERROR]: This mapping is not supported");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static
d4est_geometry_input_t
d4est_geometry_input
(
 const char* input_file,
 const char* input_section,
 const char* printf_prefix
)
{
  
  d4est_geometry_input_t input;
  input.input_section = input_section;
  input.name = NULL;
  input.DX_compute_method = GEOM_COMPUTE_NOT_SET;
  input.JAC_compute_method = GEOM_COMPUTE_NOT_SET;
  
  if (ini_parse(input_file, d4est_geometry_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input.name, NULL);
  D4EST_CHECK_INPUT(input_section, input.DX_compute_method, GEOM_COMPUTE_NOT_SET);
  D4EST_CHECK_INPUT(input_section, input.JAC_compute_method, GEOM_COMPUTE_NOT_SET);
 
  printf("%s: Loading %s geometry\n", printf_prefix, input.name);
  if (input.DX_compute_method == GEOM_COMPUTE_NUMERICAL)
    printf("%s: Dx computation method = %s\n", printf_prefix, "numerical");
  if (input.DX_compute_method == GEOM_COMPUTE_ANALYTIC)
    printf("%s: Dx computation method = %s\n", printf_prefix, "analytic");
  if (input.JAC_compute_method == GEOM_COMPUTE_NUMERICAL)
    printf("%s: JAC computation method = %s\n", printf_prefix, "numerical");
  if (input.JAC_compute_method == GEOM_COMPUTE_ANALYTIC)
    printf("%s: JAC computation method = %s\n", printf_prefix, "analytic");


  free(input.name);
  
  return input;
}


d4est_geometry_t*
d4est_geometry_new(int mpirank,
                   const char* input_file,
                   const char* input_section,
                   const char* printf_prefix)
{

  d4est_geometry_input_t input = d4est_geometry_input(input_file, input_section, printf_prefix);
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->DX_compute_method = input.DX_compute_method;
  d4est_geom->JAC_compute_method = input.JAC_compute_method;
  
#if (P4EST_DIM)==3
  if (util_match(input.name,"cubed_sphere")) {
    d4est_geometry_cubed_sphere_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"cubed_sphere_7tree")) {
    d4est_geometry_cubed_sphere_7tree_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"cubed_sphere_with_cube_hole")) {
    d4est_geometry_cubed_sphere_with_cube_hole_new(mpirank, input_file, input_section, printf_prefix,  d4est_geom);
  }
  else if (util_match(input.name,"cubed_sphere_innerouter_shell")) {
    d4est_geometry_cubed_sphere_innerouter_shell_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"cubed_sphere_outer_shell_block")){
    d4est_geometry_cubed_sphere_outer_shell_block_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"cubed_sphere_inner_shell_block")){
    d4est_geometry_cubed_sphere_inner_shell_block_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"none")){}
  else {
    printf("[D4EST_ERROR]: You tried to use %s geometry\n", input.name);
    mpi_abort("[D4EST_ERROR]: this geometry is currently not supported");
  }

#endif
#if (P4EST_DIM)==2
  if (util_match(input.name,"disk")) {
    d4est_geometry_5treedisk_new(mpirank, input_file, input_section, printf_prefix,  d4est_geom);
  }
  else if (util_match(input.name,"disk_outer_wedge")){
    d4est_geometry_disk_outer_wedge_new(mpirank, input_file, input_section, printf_prefix, d4est_geom);
  }
  else if (util_match(input.name,"none")){}
  else {
    printf("[D4EST_ERROR]: You tried to use %s geometry\n", input.name);
    mpi_abort("[D4EST_ERROR]: this geometry is currently not supported");
  }
#endif
  free(input.name);

  if(d4est_geom->DX_compute_method == GEOM_COMPUTE_ANALYTIC && d4est_geom->DX == NULL){
    mpi_abort("[D4EST_ERROR]: This geometry does not support analytic derivatives");
  }

  if(d4est_geom->JAC_compute_method == GEOM_COMPUTE_ANALYTIC && d4est_geom->JAC == NULL){
    mpi_abort("[D4EST_ERROR]: This geometry does not support analytic jacobians");
  }
  
  return d4est_geom;
}

void
d4est_geometry_destroy(d4est_geometry_t* d4est_geom){
  p4est_connectivity_destroy(d4est_geom->p4est_conn);
  d4est_geom->destroy(d4est_geom);
  /* p4est_geometry_destroy (d4est_geom->p4est_geom); */
  /* P4EST_FREE(d4est_geom); */
}

void
d4est_geometry_octree_to_vertex (p8est_connectivity_t *connectivity,
                                 p4est_topidx_t which_tree,
                                 const double abc[3], double xyz[3])
{
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
            +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                  eta_x  * v[3 * vt[5] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                  eta_x  * v[3 * vt[7] + j]))
           );
    /* *INDENT-ON* */
  }
}


void
d4est_geometry_quadtree_to_vertex (p4est_connectivity_t *connectivity,
                                 p4est_topidx_t which_tree,
                                 const double abc[3], double xyz[3])
{
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
  eta_z = abc[2];

  /* bi/trilinear transformation */
  for (j = 0; j < 3; ++j) {
    /* *INDENT-OFF* */
    xyz[j] =
           ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                  eta_x  * v[3 * vt[1] + j]) +
                                  eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                  eta_x  * v[3 * vt[3] + j]))
           );
    /* *INDENT-ON* */
  }
}


void
d4est_geometry_compute_dxyz_drst_analytic
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_rst_t rst_points,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int deg,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg);
  double rst [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];

  
  for (int i = 0; i < volume_nodes; i++){
    rst[0] = rst_points.r[i];
    rst[1] = rst_points.s[i];
#if (P4EST_DIM)==3
    rst[2] = rst_points.t[i];
#endif

    d4est_geom->DX(d4est_geom, which_tree, q0, dq, rst, dxyz_drst_i);
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      for (int d2 = 0; d2 < (P4EST_DIM); d2++){
        dxyz_drst[d1][d2][i] = dxyz_drst_i[d1][d2];
      }
    }
    
  }
}

/* void */
/* d4est_geometry_compute_Jdrdxdrdx_analytic */
/* ( */
/*  p4est_topidx_t which_tree, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int deg, */
/*  quadrature_type_t quad_type, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_operators_t* d4est_ops, */
/*  double* Jdrdxdrdx [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* {  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg); */
/*   double rst [(P4EST_DIM)]; */
/*   double Jdrdxdrdx_i [(P4EST_DIM)][(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM), quad_type); */
  
/*   for (int i = 0; i < volume_nodes; i++){ */
/*     rst[0] = rst_points.r[i]; */
/*     rst[1] = rst_points.s[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[2] = rst_points.t[i]; */
/* #endif */

/*     d4est_geom->JACDRDXDRDX(d4est_geom, which_tree, q0, dq, rst, Jdrdxdrdx_i); */
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         Jdrdxdrdx[d1][d2][i] = Jdrdxdrdx_i[d1][d2]; */
/*       } */
/*     } */
    
/*   } */
/* } */


/* void */
/* d4est_geometry_compute_drst_dxyz_analytic */
/* ( */
/*  p4est_topidx_t which_tree, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int deg, */
/*  quadrature_type_t quad_type, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_operators_t* d4est_ops, */
/*  double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*   int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg); */
/*   double rst [(P4EST_DIM)]; */
/*   double drst_dxyz_i [(P4EST_DIM)][(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM), quad_type); */
  
/*   for (int i = 0; i < volume_nodes; i++){ */
/*     rst[0] = rst_points.r[i]; */
/*     rst[1] = rst_points.s[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[2] = rst_points.t[i]; */
/* #endif */

/*     d4est_geom->DRDX(d4est_geom, which_tree, q0, dq, rst, drst_dxyz_i); */
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         drst_dxyz[d1][d2][i] = drst_dxyz_i[d1][d2]; */
/*       } */
/*     } */
    
/*   } */
/* } */


/* void */
/* d4est_geometry_compute_jac_analytic */
/* ( */
/*  p4est_topidx_t which_tree, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int deg, */
/*  quadrature_type_t quad_type, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_operators_t* d4est_ops, */
/*  double* jac */
/* ) */
/* { */
/*   int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg); */
/*   double rst [(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM), quad_type); */
  
/*   for (int i = 0; i < volume_nodes; i++){ */
/*     rst[0] = rst_points.r[i]; */
/*     rst[1] = rst_points.s[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[2] = rst_points.t[i]; */
/* #endif */

/*      d4est_geom->JAC(d4est_geom, which_tree, q0, dq, rst, &jac[i]); */
    
/*   } */
/* } */

void
d4est_geometry_compute_dxyz_drst
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_rst_t rst_points,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int deg,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)] 
)
{
  if (d4est_geom->DX_compute_method == GEOM_COMPUTE_NUMERICAL){
    /* d4est_geometry_compute_dxyz_drst_isoparametric */
    /*   ( */
    /*    which_tree, */
    /*    q0, */
    /*    dq, */
    /*    deg, */
    /*    quad_type, */
    /*    d4est_geom, */
    /*    d4est_ops, */
    /*    dxyz_drst */
    /*   ); */
    mpi_abort("[D4EST_ERROR]: we currently do not support numerical/isoparametric mappings\n");
  }
  else if (d4est_geom->DX_compute_method == GEOM_COMPUTE_ANALYTIC) {
    d4est_geometry_compute_dxyz_drst_analytic
      (
       d4est_ops,
       d4est_geom,
       rst_points,
       which_tree,
       q0,
       dq,
       deg,
       dxyz_drst
      );
  }
  else {
    mpi_abort("derivative type must be ISOPARAMETRIC or ANALYTIC");
  }
}

/* void */
/* d4est_geometry_compute_dxyz_drst_isoparametric */
/* ( */
/*  p4est_topidx_t which_tree, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int deg, */
/*  quadrature_type_t quad_type, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_operators_t* d4est_ops, */
/*  double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*   double* xyz [(P4EST_DIM)]; */
/*   int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg); */
/*   for (int i = 0; i < (P4EST_DIM); i++){ */
/*     xyz[i] = P4EST_ALLOC(double, volume_nodes); */
/*   } */

/*   double rst [(P4EST_DIM)]; */
/*   double xyz_i [(P4EST_DIM)]; */
/*   double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)]; */
/*   double abc [(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM), QUAD_LOBATTO); */
  
/*   for (int i = 0; i < volume_nodes; i++){ */
/*     rst[0] = rst_points.r[i]; */
/*     rst[1] = rst_points.s[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[2] = rst_points.t[i]; */
/* #endif */

/*     d4est_geom->X(d4est_geom, which_tree, q0, dq, rst, COORDS_INTEG_RST, xyz_i); */
      
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d][i] = xyz_i[d]; */
/*     } */
/*   } */

/*   double* tmp = P4EST_ALLOC(double, volume_nodes); */
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       d4est_operators_apply_Dij(d4est_ops, &xyz[d][0], (P4EST_DIM), deg, d1, &dxyz_drst[d][d1][0]); */
/*     } */
/*   } */

/*   if (quad_type == QUAD_GAUSS){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*         d4est_operators_interp_Lobatto_to_Gauss(d4est_ops, &dxyz_drst[d][d1][0], deg, deg, tmp, (P4EST_DIM)); */
/*         linalg_copy_1st_to_2nd(tmp, &dxyz_drst[d][d1][0], volume_nodes); */
/*       } */
/*     } */
/*   } */
/*   P4EST_FREE(tmp); */
/*   for (int i = 0; i < (P4EST_DIM); i++){ */
/*     P4EST_FREE(xyz[i]); */
/*   } */
/* } */

void
d4est_geometry_compute_dxyz_drst_face_analytic
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_rst_t rst_points,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 int deg,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM),deg);
  int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg);
  double rst [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];

  d4est_operators_face_info_t face_info;
  d4est_operators_get_face_info(face, &face_info);
  
  for (int i = 0; i < face_nodes; i++){
    rst[face_info.a] = rst_points.r[i];
#if (P4EST_DIM)==3
    rst[face_info.b] = rst_points.s[i];
#endif
    rst[face_info.c] = face_info.sgn;

    d4est_geom->DX(d4est_geom, which_tree, q0, dq, rst, dxyz_drst_i);
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      for (int d2 = 0; d2 < (P4EST_DIM); d2++){
        dxyz_drst[d1][d2][i] = dxyz_drst_i[d1][d2];
      }
    }
  }
}

/* void */
/* d4est_geometry_compute_drst_dxyz_face_analytic */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg); */
/*   double rst [(P4EST_DIM)]; */
/*   double drst_dxyz_i [(P4EST_DIM)][(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM) - 1, quad_type); */

/*   d4est_operators_face_info_t face_info; */
/*   d4est_operators_get_face_info(face, &face_info); */

/*   if(d4est_geom->DRDX == NULL){ */
/*     mpi_abort("d4est_geom->DRDX == NULL"); */
/*   } */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     rst[face_info.a] = rst_points.r[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[face_info.b] = rst_points.s[i]; */
/* #endif */
/*     rst[face_info.c] = face_info.sgn; */

/*     d4est_geom->DRDX(d4est_geom, which_tree, q0, dq, rst, drst_dxyz_i); */
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         drst_dxyz[d1][d2][i] = drst_dxyz_i[d1][d2]; */
/*       } */
/*     } */
/*   } */
/* } */

/* void */
/* d4est_geometry_compute_drst_dxyz_times_jacobian_face_analytic */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* drst_dxyz_times_jac [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg); */
/*   double rst [(P4EST_DIM)]; */
/*   double drst_dxyz_times_jac_i [(P4EST_DIM)][(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM) - 1, quad_type); */

/*   d4est_operators_face_info_t face_info; */
/*   d4est_operators_get_face_info(face, &face_info); */

/*   if(d4est_geom->DRDX == NULL){ */
/*     mpi_abort("d4est_geom->DRDX == NULL"); */
/*   } */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     rst[face_info.a] = rst_points.r[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[face_info.b] = rst_points.s[i]; */
/* #endif */
/*     rst[face_info.c] = face_info.sgn; */

/*     d4est_geom->DRDX_JAC(d4est_geom, which_tree, q0, dq, rst, drst_dxyz_times_jac_i); */
/*     for (int d1 = 0; d1 < (P4EST_DIM); d1++){ */
/*       for (int d2 = 0; d2 < (P4EST_DIM); d2++){ */
/*         drst_dxyz_times_jac[d1][d2][i] = drst_dxyz_times_jac_i[d1][d2]; */
/*       } */
/*     } */
/*   } */
/* } */



/* void */
/* d4est_geometry_compute_jac_face_analytic */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* jac */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg); */
/*   double rst [(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM) - 1, quad_type); */

/*   d4est_operators_face_info_t face_info; */
/*   d4est_operators_get_face_info(face, &face_info); */


/*   if(d4est_geom->JAC == NULL){ */
/*     mpi_abort("d4est_geom->JAC == NULL"); */
/*   } */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     rst[face_info.a] = rst_points.r[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[face_info.b] = rst_points.s[i]; */
/* #endif */
/*     rst[face_info.c] = face_info.sgn; */

/*     d4est_geom->JAC(d4est_geom, which_tree, q0, dq, rst, &jac[i]); */
/*   } */
/* } */


/* void */
/* d4est_geometry_compute_jac_face_analytic */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_rst_t rst_points, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* jac */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg); */
/*   double rst [(P4EST_DIM)]; */

/*   d4est_operators_face_info_t face_info; */
/*   d4est_operators_get_face_info(face, &face_info); */


/*   if(d4est_geom->JAC == NULL){ */
/*     mpi_abort("d4est_geom->JAC == NULL"); */
/*   } */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     rst[face_info.a] = rst_points.r[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[face_info.b] = rst_points.s[i]; */
/* #endif */
/*     rst[face_info.c] = face_info.sgn; */

/*     d4est_geom->JAC(d4est_geom, which_tree, q0, dq, rst, &jac[i]); */
/*   } */
/* } */



/* void */
/* d4est_geometry_compute_xyz_face_analytic_old */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* xyz [(P4EST_DIM)] */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg); */
/*   double rst [(P4EST_DIM)]; */

/*   d4est_rst_t rst_points */
/*     = d4est_operators_get_rst_points(d4est_ops, deg, (P4EST_DIM) - 1, quad_type); */

/*   d4est_operators_face_info_t face_info; */
/*   d4est_operators_get_face_info(face, &face_info); */
  
/*   for (int i = 0; i < face_nodes; i++){ */
/*     rst[face_info.a] = rst_points.r[i]; */
/* #if (P4EST_DIM)==3 */
/*     rst[face_info.b] = rst_points.s[i]; */
/* #endif */
/*     rst[face_info.c] = face_info.sgn; */

/*     double xyz_i [(P4EST_DIM)]; */
/*     d4est_geom->X(d4est_geom, which_tree, q0, dq, rst, COORDS_INTEG_RST, xyz_i); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d][i] = xyz_i[d]; */
/*     } */
/*   } */
/* } */

void
d4est_geometry_compute_xyz_face_analytic
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_rst_t rst_points,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 int deg,
 double* xyz [(P4EST_DIM)]
)
{
  int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1,deg);
  double rst [(P4EST_DIM)];
  
  d4est_operators_face_info_t face_info;
  d4est_operators_get_face_info(face, &face_info);
  
  for (int i = 0; i < face_nodes; i++){
    rst[face_info.a] = rst_points.r[i];
#if (P4EST_DIM)==3
    rst[face_info.b] = rst_points.s[i];
#endif
    rst[face_info.c] = face_info.sgn;

    double xyz_i [(P4EST_DIM)];
    d4est_geom->X(d4est_geom, which_tree, q0, dq, rst, COORDS_INTEG_RST, xyz_i);
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d][i] = xyz_i[d];
    }
  }
}


/* void */
/* d4est_geometry_compute_dxyz_drst_face_isoparametric */
/* ( */
/*  d4est_operators_t* d4est_ops, */
/*  d4est_rst_t rst_points, */
/*  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*  p4est_qcoord_t dq, */
/*  int which_tree, */
/*  int face, */
/*  d4est_geometry_t* d4est_geom, */
/*  quadrature_type_t quad_type, */
/*  int deg, */
/*  double* dxyz_drst_on_face [(P4EST_DIM)][(P4EST_DIM)] */
/* ) */
/* { */
/*   int face_nodes = d4est_operators_get_nodes((P4EST_DIM)-1, deg); */
/*   int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), deg); */
  
/*   double* dxyz_drst_volume [(P4EST_DIM)][(P4EST_DIM)]; */
/*   D4EST_ALLOC_DBYD_MAT(dxyz_drst_volume, volume_nodes); */
  
/*   d4est_geometry_compute_dxyz_drst */
/*     ( */
/*      which_tree, */
/*      rst_points, */
/*      q0, */
/*      dq, */
/*      deg, */
/*      d4est_geom, */
/*      d4est_ops, */
/*      dxyz_drst_volume */
/*     ); */

/*   if (quad_type == QUAD_GAUSS){ */

/*     double* temp = P4EST_ALLOC(double, face_nodes); */
/*     for (int i = 0; i < (P4EST_DIM); i++){ */
/*       for (int j = 0; j < (P4EST_DIM); j++){ */
/*         d4est_operators_apply_slicer(d4est_ops, */
/*                             dxyz_drst_volume[i][j], */
/*                             (P4EST_DIM), */
/*                             face, */
/*                             deg, */
/*                             temp); */
      
/*         d4est_operators_interp_Lobatto_to_Gauss */
/*           ( */
/*            d4est_ops, */
/*            temp, */
/*            deg, */
/*            deg, */
/*            dxyz_drst_on_face[i][j], */
/*            (P4EST_DIM)-1 */
/*           ); */
/*       } */
/*     } */
/*     P4EST_FREE(temp); */
/*   } */
/*   else if (quad_type == QUAD_LOBATTO){ */
/*     for (int i = 0; i < (P4EST_DIM); i++){ */
/*       for (int j = 0; j < (P4EST_DIM); j++){ */
/*         d4est_operators_apply_slicer(d4est_ops, */
/*                             dxyz_drst_volume[i][j], */
/*                             (P4EST_DIM), */
/*                             face, */
/*                             deg, */
/*                             dxyz_drst_on_face[i][j]); */
/*       } */
/*     } */
/*   } */
/*   else { */
/*     mpi_abort("quad type not supported"); */
/*   } */

/*   D4EST_FREE_DBYD_MAT(dxyz_drst_volume); */
/* } */


void
d4est_geometry_compute_jacobian
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* jac,
 int volume_nodes
)
{
  for (int i = 0; i < volume_nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];

#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif
    
#if (P4EST_DIM) == 3
    jac[i] = xr*(ys*zt-zs*yt)
                - yr*(xs*zt-zs*xt)
                + zr*(xs*yt-ys*xt);
#elif (P4EST_DIM) == 2
    jac[i] = -xs*yr + xr*ys;
#else
    mpi_abort("DIM must be 2 or 3");
#endif
  }
}

void
d4est_geometry_compute_drst_dxyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* jac,
 double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)],
 int volume_nodes
)
{
  for (int i = 0; i < volume_nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];

#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif

    double J = jac[i];

    double* rx = &drst_dxyz[0][0][i];
    double* ry = &drst_dxyz[0][1][i];
    double* sx = &drst_dxyz[1][0][i];
    double* sy = &drst_dxyz[1][1][i];

#if (P4EST_DIM)==3
    double* rz = &drst_dxyz[0][2][i];
    double* sz = &drst_dxyz[1][2][i];
    
    double* tx = &drst_dxyz[2][0][i];
    double* ty = &drst_dxyz[2][1][i];
    double* tz = &drst_dxyz[2][2][i];
#endif

#if (P4EST_DIM)==3
    *rx =  (ys*zt - zs*yt)/(J);
    *ry = -(xs*zt - zs*xt)/(J);
    *rz =  (xs*yt - ys*xt)/(J);
    
    *sx = -(yr*zt - zr*yt)/(J);
    *sy =  (xr*zt - zr*xt)/(J);
    *sz = -(xr*yt - yr*xt)/(J);
    
    *tx =  (yr*zs - zr*ys)/(J);
    *ty = -(xr*zs - zr*xs)/(J);
    *tz =  (xr*ys - yr*xs)/(J);
#endif

#if (P4EST_DIM)==2
    *rx = ys/(J);
    *sx =-yr/(J);
    *ry =-xs/(J);
    *sy = xr/(J);
#endif
    
  }  
}

void
d4est_geometry_compute_drst_dxyz_times_jacobian
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* drst_dxyz_times_jac [(P4EST_DIM)][(P4EST_DIM)],
 int volume_nodes
)
{
  for (int i = 0; i < volume_nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];

#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif

    double* rx_jac = &drst_dxyz_times_jac[0][0][i];
    double* ry_jac = &drst_dxyz_times_jac[0][1][i];
    double* sx_jac = &drst_dxyz_times_jac[1][0][i];
    double* sy_jac = &drst_dxyz_times_jac[1][1][i];

#if (P4EST_DIM)==3
    double* rz_jac = &drst_dxyz_times_jac[0][2][i];
    double* sz_jac = &drst_dxyz_times_jac[1][2][i];
    
    double* tx_jac = &drst_dxyz_times_jac[2][0][i];
    double* ty_jac = &drst_dxyz_times_jac[2][1][i];
    double* tz_jac = &drst_dxyz_times_jac[2][2][i];
#endif

#if (P4EST_DIM)==3
    *rx_jac =  (ys*zt - zs*yt);
    *ry_jac = -(xs*zt - zs*xt);
    *sx_jac = -(yr*zt - zr*yt);
    *sy_jac =  (xr*zt - zr*xt);
    *rz_jac =  (xs*yt - ys*xt);
    *sz_jac = -(xr*yt - yr*xt);
    *tx_jac =  (yr*zs - zr*ys);
    *ty_jac = -(xr*zs - zr*xs);
    *tz_jac =  (xr*ys - yr*xs);
#endif
#if (P4EST_DIM)==2    
    *rx_jac = ys;
    *sx_jac =-yr;
    *ry_jac =-xs;
    *sy_jac = xr;
#endif
    
  }  
}

void
d4est_geometry_compute_qcoords_on_mortar
(
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 int num_faces_side, 
 int num_faces_mortar,
 int face,
 p4est_qcoord_t mortar_q0 [(P4EST_HALF)][(P4EST_DIM)],
 p4est_qcoord_t* mortar_dq
)
{
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = d4est_operators_is_child_left_or_right(c, d);
      mortar_q0[j][d] = e0_q[d] + cd*e0_dq;
    }
  }
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      dqa[dir][d] = (mortar_q0[(dir+1)][d] - mortar_q0[0][d]);
      if (num_faces_side != num_faces_mortar)
        dqa[dir][d] /= 2;
    }
  }    

  p4est_qcoord_t dq0mf0 [(P4EST_DIM)];
  for (int d = 0; d < (P4EST_DIM); d++){
    dq0mf0[d] = (mortar_q0[0][d] - e0_q[d])/2;
  }
    
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int c = 0; c < (P4EST_HALF); c++){
      mortar_q0[c][d] = e0_q[d];
      if (num_faces_side != num_faces_mortar)
        mortar_q0[c][d] += dq0mf0[d];
      for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
        int cd = d4est_operators_is_child_left_or_right(c, dir);
        mortar_q0[c][d] += cd*dqa[dir][d];
      }
    }
  }

  *mortar_dq = (num_faces_side == num_faces_mortar) ? e0_dq : e0_dq/2;
}



void
d4est_geometry_compute_xyz
(
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_rst_t rst_points,
 int which_tree,
 int deg,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 double* xyz [(P4EST_DIM)]
)
{    
  int volume_nodes = d4est_operators_get_nodes((P4EST_DIM), deg);
  
  double rst_i [(P4EST_DIM)]; 
  double xyz_i [(P4EST_DIM)];
  for (int i = 0; i < volume_nodes; i++){
    rst_i[0] = rst_points.r[i];
    rst_i[1] = rst_points.s[i];
#if (P4EST_DIM) == 3
    rst_i[2] = rst_points.t[i];
#endif
    
    d4est_geom->X(d4est_geom, which_tree, q, dq, rst_i, COORDS_INTEG_RST, xyz_i);
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d][i] = xyz_i[d];
    }
    /* printf("xyz[0], xyz[1], xyz[2], rst[0], rst[1], rst[2] = %f,%f,%f,%f,%f,%f\n", xyz_i[0], xyz_i[1], xyz_i[2], rst_i[0], rst_i[1], rst_i[2]); */
    
  }
}

/* void */
/* d4est_geometry_compute_geometric_data_on_mortar_TESTINGONLY */
/* ( */
/*  p4est_topidx_t e0_tree, */
/*  p4est_qcoord_t e0_q [(P4EST_DIM)], */
/*  p4est_qcoord_t e0_dq, */
/*  int num_faces_side, */
/*  int num_faces_mortar, */
/*  int* deg_mortar, */
/*  int face_side, */
/*  quadrature_type_t quad_type, */
/*  double* n [(P4EST_DIM)], */
/*  double* sj, */
/*  d4est_geometry_t* d4est_geom, */
/*  d4est_operators_t* d4est_ops, */
/*  double* xyz_storage [(P4EST_DIM)] */
/* ) */
/* { */

/*   double* xyz [(P4EST_DIM)]; */
/*   if (xyz_storage[0] != NULL){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d] = xyz_storage[d]; */
/*     } */
/*   } */
/*   else { */
/*     int total_face_mortar_nodes = 0; */
/*     for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */
/*       total_face_mortar_nodes += d4est_operators_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]); */
/*     } */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       xyz[d] = P4EST_ALLOC(double, total_face_mortar_nodes); */
/*     } */
/*   } */

/*   /\* Calculate the four "0" corners of  */
/*    * the mortar faces. In the case that */
/*    * there is only one mortar face, these */
/*    * will be the four corners of that face */
/*    *\/ */
  
/*   p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)]; */
  
/*   for (int j = 0; j < (P4EST_HALF); j++){ */
/*     int c = p4est_face_corners[face_side][j]; */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       int cd = d4est_operators_is_child_left_or_right(c, d); */
/*       q0[j][d] = e0_q[d] + cd*e0_dq; */
/*     } */
/*   } */

/*   /\* Calculate the vectors that span the face  */
/*    * there will be one in 2-D and two in 3-d *\/ */
  
/*   p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)]; */
  
/*   for (int d = 0; d < (P4EST_DIM); d++){ */
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]); */
/*       if (num_faces_side != num_faces_mortar) */
/*         dqa[dir][d] /= 2; */
/*     } */
/*   } */

/*   if (num_faces_side != num_faces_mortar){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       for (int c = 0; c < (P4EST_HALF); c++){ */
/*         q0[c][d] = q0[0][d]; */
/*         for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){ */
/*           int cd = d4est_operators_is_child_left_or_right(c, dir); */
/*           q0[c][d] += cd*dqa[dir][d]; */
/*         } */
/*       } */
/*     } */
/*   } */
  
/*   double* a [((P4EST_DIM)-1)]; */
/*   /\* double* xyz [(P4EST_DIM)]; *\/ */
/*   double* dxda [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double dxyz_drs_i [(P4EST_DIM)][((P4EST_DIM)-1)]; */
/*   double abc [] = {0.,0.,0.}; */
/*   double xyz_i [] = {0.,0.,0.}; */
/*   int face_mortar_nodal_stride = 0; */
/*   d4est_rst_t rst_points; */
  
/*   for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){ */


/*     int face_mortar_nodes = d4est_operators_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]); */
/*     /\* compute the LGL nodes in the directions of the face_mortar vectors *\/ */
/*     double* tmp = P4EST_ALLOC(double,face_mortar_nodes); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       /\* xyz[d] = P4EST_ALLOC(double, face_mortar_nodes); *\/ */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*         dxda[d][dir] = P4EST_ALLOC(double, face_mortar_nodes); */
/*     } */
     
/*     rst_points = d4est_operators_get_rst_points(d4est_ops, deg_mortar[face_mortar], (P4EST_DIM)-1, QUAD_LOBATTO); */
      
/*     for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*       if (dir == 0) */
/*         a[dir] = rst_points.r; */
/*       else */
/*         a[dir] = rst_points.s; */
/*     } */

    
/*     for (int i = 0; i < face_mortar_nodes; i++){ */
/*       if (xyz[0] != NULL){ */
/*         for (int d = 0; d < (P4EST_DIM); d++){ */
/*           /\* get "0" corner of this face_mortar *\/ */
/*           abc[d] = (double)q0[face_mortar][d]; */
       
/*           for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*             /\* add a fraction of the face_mortar vector in direction dir */
/*            * corresponding to the placement of the LGL node *\/ */
/*             double da = (a[dir][i] + 1.)/2.; */
/*             abc[d] += da*((double)dqa[dir][d]); */
/*             /\* rs[dir] = a[dir][i]; *\/ */
/*             /\* printf("abc[%d] = %f\n",d, abc[d]); *\/ */
/*           } */
        
/*           abc[d] /= (double)(P4EST_ROOT_LEN); */
/*         } */

/*         /\* convert vertex coords to physical coords *\/ */
/*         d4est_geom->X(d4est_geom, e0_tree, (p4est_qcoord_t [(P4EST_DIM)]){0},-1, abc, COORDS_TREE_UNITCUBE , xyz_i);      */
/*         for (int d = 0; d < (P4EST_DIM); d++){ */
/*           xyz[d][i] = xyz_i[d]; */
/*         } */
/*       } */

/*     } */

/*     /\* compute the tangent vectors in direction(s) "dir" *\/ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           d4est_operators_apply_Dij(d4est_ops, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], dir, dxda[d][dir]); */
/*         } */
/*       } */

/*     if (quad_type == QUAD_GAUSS){ */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         if (xyz[0] != NULL){ */
/*           d4est_operators_interp_Lobatto_to_Gauss(d4est_ops, xyz[d], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1); */
/*           linalg_copy_1st_to_2nd(tmp, xyz[d], face_mortar_nodes); */
/*         } */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){ */
/*           d4est_operators_interp_Lobatto_to_Gauss(d4est_ops, dxda[d][dir], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1); */
/*           linalg_copy_1st_to_2nd(tmp, dxda[d][dir], face_mortar_nodes); */
/*         } */
/*       } */
/*     } */
    
/*     /\* get the normal by taking the cross product of the tangent vectors */
/*      * in 2-d, we take the cross product of the tangent vector and zhat*\/ */
/*     for (int i = 0; i < face_mortar_nodes; i++){ */
/*       double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}}; */
/*       double n_i [] = {0.,0.,0.}; */
/*       for (int d = 0; d < (P4EST_DIM); d++) */
/*         for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*           vecs[dir][d] = dxda[d][dir][i]; */

/*       linalg_cross_prod */
/*         ( */
/*          vecs[0][0], */
/*          vecs[0][1], */
/*          vecs[0][2], */
/*          vecs[1][0], */
/*          vecs[1][1], */
/*          vecs[1][2], */
/*          &(n_i[0]), */
/*          &(n_i[1]), */
/*          &(n_i[2]) */
/*         ); */

/*       sj[i + face_mortar_nodal_stride] = 0.; */
/*       for (int d = 0; d < (P4EST_DIM); d++){ */
/*         /\* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) *\/ */
/*         if (face_side == 0 || face_side == 3 || face_side == 4){ */
/*           n_i[d] *= -1.; */
/*         } */
/*         sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d]; */
/*       } */
/*       sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]); */
/*       if (n[0] != NULL){ */
/*         for (int d = 0; d < (P4EST_DIM); d++){ */
/*           n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride]; */
/*         } */
/*       } */
/*     } */

/*     face_mortar_nodal_stride += d4est_operators_get_nodes((P4EST_DIM)-1, deg_mortar[face_mortar]); */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       /\* P4EST_FREE(xyz[d]); *\/ */
/*       for (int dir = 0; dir < ((P4EST_DIM)-1); dir++) */
/*         P4EST_FREE(dxda[d][dir]); */
/*     } */
/*     P4EST_FREE(tmp); */
/*   } */

/*   if (xyz_storage[0] == NULL){ */
/*     for (int d = 0; d < (P4EST_DIM); d++){ */
/*       P4EST_FREE(xyz[d]); */
/*     } */
/*   } */
  
/* } */

void
d4est_geometry_get_tree_coords_in_range_0_to_1
(
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double coords[3],
 coords_type_t coords_type,
 double tcoords[3]
)
{
  if (coords_type == COORDS_INTEG_RST){
    mpi_assert(coords[0] >= -1. && coords[0] <= 1.);
    mpi_assert(coords[1] >= -1. && coords[1] <= 1.);
    /* transform from coords in [-1,1]^3 to the subspace [q0, q0+dq]^3 in [0,1]^3*/
    tcoords[0] = d4est_operators_rtox(coords[0], (double)q0[0], (double)dq)/(double)P4EST_ROOT_LEN;
    tcoords[1] = d4est_operators_rtox(coords[1], (double)q0[1], (double)dq)/(double)P4EST_ROOT_LEN;
#if (P4EST_DIM)==3
    mpi_assert(coords[2] >= -1. && coords[2] <= 1.);
    tcoords[2] = d4est_operators_rtox(coords[2], (double)q0[2], (double)dq)/(double)P4EST_ROOT_LEN;
#endif
  }
  else if (coords_type == COORDS_P4EST_INT){
    tcoords[0] = coords[0]/(double)P4EST_ROOT_LEN;
    tcoords[1] = coords[1]/(double)P4EST_ROOT_LEN;
#if (P4EST_DIM)==3
    tcoords[2] = coords[2]/(double)P4EST_ROOT_LEN;
#endif
  }
  /* coords are already tree coords in range [0,1] */
  else if (coords_type == COORDS_TREE_UNITCUBE){
    tcoords[0] = coords[0];
    tcoords[1] = coords[1];
#if (P4EST_DIM)==3
    tcoords[2] = coords[2];
#endif
  }
  else {
    mpi_abort("[D4EST_ERROR]: You either did not set a coords type or used an unsupported coords type");
  }
}


int
d4est_geometry_does_element_touch_boundary
(
 p4est_t* p4est,
 p4est_quadrant_t* q,
 int which_tree
)
{
  int fbsum = 0;
  for (int face = 0; face < (P4EST_FACES); face++){
    fbsum += d4est_geometry_is_face_on_boundary(p4est, q, which_tree, face);
  }
  return (fbsum > 0);
}


int
d4est_geometry_is_face_on_boundary
(
 p4est_t* p4est,
 p4est_quadrant_t* q,
 int which_tree,
 int face
)
{
  p4est_qcoord_t      dh, xyz_temp;
  p4est_connectivity_t *conn = p4est->connectivity;
  int fbsum = 0;
  if (conn->tree_to_tree[P4EST_FACES * which_tree + face] != which_tree ||
      (int) conn->tree_to_face[P4EST_FACES * which_tree + face] != face) {
  }
  else {
    dh = P4EST_LAST_OFFSET (q->level);
    switch (face / 2) {
    case 0:
      xyz_temp = q->x;
      break;
    case 1:
      xyz_temp = q->y;
      break;
#ifdef P4_TO_P8
    case 2:
      xyz_temp = q->z;
      break;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    fbsum += (xyz_temp == ((face & 0x01) ? dh : 0));
  }

  return (fbsum > 0);

}
