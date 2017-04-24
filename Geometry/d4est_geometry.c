#include <util.h>
#include <ini.h>
#include <d4est_geometry.h>
#include <dgmath.h>
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
  mapping_type_t X_mapping_type;
  
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
  else if (util_match_couple(section,pconfig->input_section,name,"X_mapping_type")) {
    if(util_match(value, "isoparametric")){
      pconfig->X_mapping_type = MAP_ISOPARAMETRIC;
    }
    else if (util_match(value, "analytic")){
      pconfig->X_mapping_type = MAP_ANALYTIC;
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
  input.X_mapping_type = MAP_NONE;
  
  if (ini_parse(input_file, d4est_geometry_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, input.name, NULL);
  D4EST_CHECK_INPUT(input_section, input.X_mapping_type, MAP_NONE);
 
  printf("%s: Loading %s geometry\n", printf_prefix, input.name);
  if (input.X_mapping_type == MAP_ISOPARAMETRIC)
    printf("%s: Mapping type = %s\n", printf_prefix, "isoparametric");
  if (input.X_mapping_type == MAP_ANALYTIC)
    printf("%s: Mapping type = %s\n", printf_prefix, "analytic");
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
  d4est_geom->X_mapping_type = input.X_mapping_type;
  
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

  if(d4est_geom->X_mapping_type == MAP_ANALYTIC && d4est_geom->DX == NULL){
    mpi_abort("[D4EST_ERROR]:If X_mapping_type = analytic you must set DX function");
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
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int deg,
 quadrature_type_t quad_type,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes = dgmath_get_nodes((P4EST_DIM),deg);
  double rst [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];

  dgmath_rst_t rst_points
    = dgmath_get_rst_points(dgmath_jit_dbase, deg, (P4EST_DIM), quad_type);
  
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

void
d4est_geometry_compute_dxyz_drst
(

 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int deg,
 quadrature_type_t quad_type,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)] 
)
{
  if (d4est_geom->X_mapping_type == MAP_ISOPARAMETRIC){
    d4est_geometry_compute_dxyz_drst_isoparametric
      (
       which_tree,
       q0,
       dq,
       deg,
       quad_type,
       d4est_geom,
       dgmath_jit_dbase,
       dxyz_drst
      );
  }
  else if (d4est_geom->X_mapping_type == MAP_ANALYTIC) {
    d4est_geometry_compute_dxyz_drst_analytic
      (
       which_tree,
       q0,
       dq,
       deg,
       quad_type,
       d4est_geom,
       dgmath_jit_dbase,
       dxyz_drst
      );
  }
  else {
    mpi_abort("derivative type must be ISOPARAMETRIC or ANALYTIC");
  }
}

void
d4est_geometry_compute_dxyz_drst_isoparametric
(
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int deg,
 quadrature_type_t quad_type,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]
)
{
  double* xyz [(P4EST_DIM)];
  int volume_nodes = dgmath_get_nodes((P4EST_DIM),deg);
  for (int i = 0; i < (P4EST_DIM); i++){
    xyz[i] = P4EST_ALLOC(double, volume_nodes);
  }

  double rst [(P4EST_DIM)];
  double xyz_i [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];
  double abc [(P4EST_DIM)];

  dgmath_rst_t rst_points
    = dgmath_get_rst_points(dgmath_jit_dbase, deg, (P4EST_DIM), LOBATTO);
  
  for (int i = 0; i < volume_nodes; i++){
    rst[0] = rst_points.r[i];
    rst[1] = rst_points.s[i];
#if (P4EST_DIM)==3
    rst[2] = rst_points.t[i];
#endif

    d4est_geom->X(d4est_geom, which_tree, q0, dq, rst, COORDS_INTEG_RST, xyz_i);
      
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d][i] = xyz_i[d];
    }
  }

  double* tmp = P4EST_ALLOC(double, volume_nodes);
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      dgmath_apply_Dij(dgmath_jit_dbase, &xyz[d][0], (P4EST_DIM), deg, d1, &dxyz_drst[d][d1][0]);
    }
  }

  if (quad_type == GAUSS){
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int d1 = 0; d1 < (P4EST_DIM); d1++){
        dgmath_interp_GLL_to_GL(dgmath_jit_dbase, &dxyz_drst[d][d1][0], deg, deg, tmp, (P4EST_DIM));
        linalg_copy_1st_to_2nd(tmp, &dxyz_drst[d][d1][0], volume_nodes);
      }
    }
  }
  P4EST_FREE(tmp);
  for (int i = 0; i < (P4EST_DIM); i++){
    P4EST_FREE(xyz[i]);
  }
}

void
d4est_geometry_data_compute_dxyz_drst_face_analytic
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 d4est_geometry_t* d4est_geom,
 quadrature_type_t quad_type,
 int deg,
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int volume_nodes = dgmath_get_nodes((P4EST_DIM),deg);
  int face_nodes = dgmath_get_nodes((P4EST_DIM)-1,deg);
  double rst [(P4EST_DIM)];
  double dxyz_drst_i [(P4EST_DIM)][(P4EST_DIM)];

  dgmath_rst_t rst_points
    = dgmath_get_rst_points(dgmath_jit_dbase, deg, (P4EST_DIM) - 1, quad_type);

  dgmath_face_info_t face_info;
  dgmath_get_face_info(face, &face_info);
  
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


void
d4est_geometry_data_compute_xyz_face_analytic
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 d4est_geometry_t* d4est_geom,
 quadrature_type_t quad_type,
 int deg,
 double* xyz [(P4EST_DIM)]
)
{
  int face_nodes = dgmath_get_nodes((P4EST_DIM)-1,deg);
  double rst [(P4EST_DIM)];

  dgmath_rst_t rst_points
    = dgmath_get_rst_points(dgmath_jit_dbase, deg, (P4EST_DIM) - 1, quad_type);

  dgmath_face_info_t face_info;
  dgmath_get_face_info(face, &face_info);
  
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



void
d4est_geometry_data_compute_dxyz_drst_face_isoparametric
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int which_tree,
 int face,
 d4est_geometry_t* d4est_geom,
 quadrature_type_t quad_type,
 int deg,
 double* dxyz_drst_on_face [(P4EST_DIM)][(P4EST_DIM)]
)
{
  int face_nodes = dgmath_get_nodes((P4EST_DIM)-1, deg);
  int volume_nodes = dgmath_get_nodes((P4EST_DIM), deg);
  
  double* dxyz_drst_volume [(P4EST_DIM)][(P4EST_DIM)];
  D4EST_ALLOC_DBYD_MAT(dxyz_drst_volume, volume_nodes);
  
  d4est_geometry_compute_dxyz_drst
    (
     which_tree,
     q0,
     dq,
     deg,
     LOBATTO,
     d4est_geom,
     dgmath_jit_dbase,
     dxyz_drst_volume
    );

  if (quad_type == GAUSS){

    double* temp = P4EST_ALLOC(double, face_nodes);
    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
        dgmath_apply_slicer(dgmath_jit_dbase,
                            dxyz_drst_volume[i][j],
                            (P4EST_DIM),
                            face,
                            deg,
                            temp);
      
        dgmath_interp_GLL_to_GL
          (
           dgmath_jit_dbase,
           temp,
           deg,
           deg,
           dxyz_drst_on_face[i][j],
           (P4EST_DIM)-1
          );
      }
    }
    P4EST_FREE(temp);
  }
  else if (quad_type == LOBATTO){
    for (int i = 0; i < (P4EST_DIM); i++){
      for (int j = 0; j < (P4EST_DIM); j++){
        dgmath_apply_slicer(dgmath_jit_dbase,
                            dxyz_drst_volume[i][j],
                            (P4EST_DIM),
                            face,
                            deg,
                            dxyz_drst_on_face[i][j]);
      }
    }
  }
  else {
    mpi_abort("quad type not supported");
  }

  D4EST_FREE_DBYD_MAT(dxyz_drst_volume);
}

void
d4est_geometry_compute_jacobian_and_drst_dxyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* jac,
 double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)],
 double* drst_dxyz_times_jac [(P4EST_DIM)][(P4EST_DIM)],
 int volume_nodes
)
{
  for (int i = 0; i < volume_nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
#endif
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];
#if (P4EST_DIM)==3
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif

    double* rx;
    double* ry;
#if (P4EST_DIM)==3
    double* rz;
#endif
    double* sx;
    double* sy;
#if (P4EST_DIM)==3
    double* sz;
    double* tx;
    double* ty;
    double* tz;
#endif
    
    if (drst_dxyz != NULL){
     rx = &drst_dxyz[0][0][i];
     ry = &drst_dxyz[0][1][i];
#if (P4EST_DIM)==3
     rz = &drst_dxyz[0][2][i];
#endif
     sx = &drst_dxyz[1][0][i];
     sy = &drst_dxyz[1][1][i];
#if (P4EST_DIM)==3
     sz = &drst_dxyz[1][2][i];
    
     tx = &drst_dxyz[2][0][i];
     ty = &drst_dxyz[2][1][i];
     tz = &drst_dxyz[2][2][i];
#endif
    }

    double* rx_jac = NULL;
    double* ry_jac = NULL;
    double* rz_jac = NULL;


    double* sx_jac = NULL;
    double* sy_jac = NULL;
    double* sz_jac = NULL;


    double* tx_jac = NULL;
    double* ty_jac = NULL;
    double* tz_jac = NULL;
    
    if (drst_dxyz_times_jac != NULL){
     rx_jac = &drst_dxyz_times_jac[0][0][i];
     ry_jac = &drst_dxyz_times_jac[0][1][i];
#if (P4EST_DIM)==3
     rz_jac = &drst_dxyz_times_jac[0][2][i];
#endif
     sx_jac = &drst_dxyz_times_jac[1][0][i];
     sy_jac = &drst_dxyz_times_jac[1][1][i];
#if (P4EST_DIM)==3
     sz_jac = &drst_dxyz_times_jac[1][2][i];
    
     tx_jac = &drst_dxyz_times_jac[2][0][i];
     ty_jac = &drst_dxyz_times_jac[2][1][i];
     tz_jac = &drst_dxyz_times_jac[2][2][i];
#endif
    }
    
    double Jvar = 0;
    double* J = &Jvar;
    if (jac != NULL){
      J = &jac[i];
    }
    
#if (P4EST_DIM) == 3
    *J = xr*(ys*zt-zs*yt)
                - yr*(xs*zt-zs*xt)
                + zr*(xs*yt-ys*xt);

    if (drst_dxyz_times_jac != NULL){
    *rx_jac =  (ys*zt - zs*yt);
    *ry_jac = -(xs*zt - zs*xt);
    *rz_jac =  (xs*yt - ys*xt);
    *sx_jac = -(yr*zt - zr*yt);
    *sy_jac =  (xr*zt - zr*xt);
    *sz_jac = -(xr*yt - yr*xt);
    *tx_jac =  (yr*zs - zr*ys);
    *ty_jac = -(xr*zs - zr*xs);
    *tz_jac =  (xr*ys - yr*xs);
    }
    if (drst_dxyz != NULL && jac != NULL){
    *rx =  (ys*zt - zs*yt)/(*J);
    *ry = -(xs*zt - zs*xt)/(*J);
    *rz =  (xs*yt - ys*xt)/(*J);
    *sx = -(yr*zt - zr*yt)/(*J);
    *sy =  (xr*zt - zr*xt)/(*J);
    *sz = -(xr*yt - yr*xt)/(*J);
    *tx =  (yr*zs - zr*ys)/(*J);
    *ty = -(xr*zs - zr*xs)/(*J);
    *tz =  (xr*ys - yr*xs)/(*J);
    }
#elif (P4EST_DIM) == 2
    *J = -xs*yr + xr*ys;

    if (drst_dxyz_times_jac != NULL){
    *rx_jac = ys;
    *sx_jac =-yr;
    *ry_jac =-xs;
    *sy_jac = xr;
    }
    if (drst_dxyz != NULL && jac != NULL){
    *rx = ys/(*J);
    *sx =-yr/(*J);
    *ry =-xs/(*J);
    *sy = xr/(*J);
    }
#else
    mpi_abort("DIM must be 2 or 3");
#endif 
  }  
}




void
d4est_geometry_compute_geometric_data_on_mortar
(
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)], /* qcoord of first element of side */
 p4est_qcoord_t e0_dq, /* qcoord vector spanning first element of side */
 int num_faces_side, 
 int num_faces_mortar,
 int* deg_mortar_integ,
 int face,
 double* drst_dxyz_on_mortar_integ [(P4EST_DIM)][(P4EST_DIM)],
 double* sj_on_mortar_integ,
 double* n_on_mortar_integ [(P4EST_DIM)],
 double* j_div_sj_mortar_integ,
 quadrature_type_t quad_type,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 normal_compute_method_t n_compute_method
)
{
  
  /* double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)]; */
  double* dxyz_drst_on_face_integ [(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_on_face_integ[(P4EST_DIM)][(P4EST_DIM)];
  double* drst_dxyz_times_jac_on_face_integ[(P4EST_DIM)][(P4EST_DIM)];
  
  int max_deg = 0;
  for (int i = 0; i < num_faces_mortar; i++){
    max_deg = (deg_mortar_integ[i] > max_deg) ? deg_mortar_integ[i] : max_deg;
  }
  int volume_nodes_max = dgmath_get_nodes((P4EST_DIM), max_deg);
  int face_nodes_max = dgmath_get_nodes((P4EST_DIM)-1, max_deg);
  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      /* dxyz_drst[i][j] = P4EST_ALLOC(double, volume_nodes_max); */
      dxyz_drst_on_face_integ[i][j] = P4EST_ALLOC(double, face_nodes_max);
      drst_dxyz_times_jac_on_face_integ[i][j] = P4EST_ALLOC(double, face_nodes_max);
      
      if (drst_dxyz_on_mortar_integ == NULL){
        /* printf("drst_dxyz_on_mortar_integ is NULL\n"); */
        drst_dxyz_on_face_integ[i][j] = P4EST_ALLOC(double, face_nodes_max);
      }
    }

  double* temp = P4EST_ALLOC(double, volume_nodes_max);
  double* J_on_face_integ = P4EST_ALLOC(double, face_nodes_max);


  
  p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)];
  
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = dgmath_is_child_left_or_right(c, d);
      q0[j][d] = e0_q[d] + cd*e0_dq;
    }
  }
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
        dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]);
        if (num_faces_side != num_faces_mortar)
          dqa[dir][d] /= 2;
      }
    }    

    p4est_qcoord_t dq0mf0 [(P4EST_DIM)];
    for (int d = 0; d < (P4EST_DIM); d++){
      dq0mf0[d] = (q0[0][d] - e0_q[d])/2;
    }
    
    for (int d = 0; d < (P4EST_DIM); d++){
        for (int c = 0; c < (P4EST_HALF); c++){
          q0[c][d] = e0_q[d];
          if (num_faces_side != num_faces_mortar)
            q0[c][d] += dq0mf0[d];
          for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
            int cd = dgmath_is_child_left_or_right(c, dir);
            q0[c][d] += cd*dqa[dir][d];
          }
        }
    }

    p4est_qcoord_t q [(P4EST_DIM)];
    p4est_qcoord_t mortar_dq = (num_faces_side == num_faces_mortar) ? e0_dq : e0_dq/2;
 
  int face_mortar_integ_stride = 0;
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
  
    int face_mortar_integ_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar_integ[face_mortar]);
    /* double* xyz [(P4EST_DIM)]; */
    for (int d = 0; d < (P4EST_DIM); d++){
      q[d] = q0[face_mortar][d];
    }

    if (d4est_geom->X_mapping_type == MAP_ANALYTIC){
      d4est_geometry_data_compute_dxyz_drst_face_analytic(dgmath_jit_dbase,q, mortar_dq, e0_tree, face, d4est_geom, quad_type, deg_mortar_integ[face_mortar], dxyz_drst_on_face_integ);
    }
    else if (d4est_geom->X_mapping_type == MAP_ISOPARAMETRIC){
      /* printf("using MAP_ISOPARAMETRIC\n"); */
      d4est_geometry_data_compute_dxyz_drst_face_isoparametric(dgmath_jit_dbase,q, mortar_dq,e0_tree, face, d4est_geom, quad_type, deg_mortar_integ[face_mortar], dxyz_drst_on_face_integ);
    }
    else {
      mpi_abort("mapping type not supported");
    }


    if (drst_dxyz_on_mortar_integ != NULL){
      /* printf("drst_dxyz is not NULL\n"); */
      for (int i = 0; i < (P4EST_DIM); i++){
        for (int j = 0; j < (P4EST_DIM); j++){
          drst_dxyz_on_face_integ[i][j] = &drst_dxyz_on_mortar_integ[i][j][face_mortar_integ_stride];
        }
      }
    }
    
    d4est_geometry_compute_jacobian_and_drst_dxyz(dxyz_drst_on_face_integ, J_on_face_integ, drst_dxyz_on_face_integ, drst_dxyz_times_jac_on_face_integ, face_mortar_integ_nodes);

      int i0 = -1; 
      if (face == 0 || face == 1){
        i0 = 0;
      }
      else if (face == 2 || face == 3){
        i0 = 1;
      }
      else if (face == 4 || face == 5){
        i0 = 2;
      }
      else {
        mpi_abort("face must be < 6\n");
      }
      double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
      dgmath_face_info_t face_info;
      dgmath_get_face_info(face, &face_info);

      
      for (int i = 0; i < face_mortar_integ_nodes; i++){
        int is = face_mortar_integ_stride + i;
        double n_is [] = {0.,0.,0.};
        double sj_is = 0.;
        if(n_compute_method == COMPUTE_NORMAL_USING_JACOBIAN){
          double sgn = (face == 0 || face == 2 || face == 4) ? -1. : 1.;
          for (int d = 0; d < (P4EST_DIM); d++){
            /* n_is[d] = sgn*drst_dxyz_on_face_integ[i0][d][i]*J_on_face_integ[i]; */
            n_is[d] = sgn*drst_dxyz_times_jac_on_face_integ[i0][d][i];
            sj_is += n_is[d]*n_is[d];
          }
          sj_is = sqrt(sj_is);
          for (int d = 0; d < (P4EST_DIM); d++)
            n_is[d] /= sj_is;
        }
        
        else if (n_compute_method == COMPUTE_NORMAL_USING_CROSS_PRODUCT){
          double sgn = (face == 0 || face == 3 || face == 4) ? -1. : 1.;
          double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}};
          
          for (int d = 0; d < (P4EST_DIM); d++)
            for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
              vecs[dir][d] =  dxyz_drst_on_face_integ[d][dir == 0 ? face_info.a : face_info.b][i];
            }
          
          linalg_cross_prod
            (
             vecs[0][0],
             vecs[0][1],
             vecs[0][2],
             vecs[1][0],
             vecs[1][1],
             vecs[1][2],
             &(n_is[0]),
             &(n_is[1]),
             &(n_is[2])
            );


          sj_is = 0.;
          for (int d = 0; d < (P4EST_DIM); d++){
            /* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) */
            n_is[d] *= sgn;
            sj_is += n_is[d]*n_is[d];
          }
          sj_is = sqrt(sj_is);
          
          for (int d = 0; d < (P4EST_DIM); d++)
            n_is[d] /= sj_is;
          
        }
        else {
          mpi_abort("[D4EST_ERROR]: Only two ways to compute normal");
        }

        /* STORE COMPUTATIONS */
        if (n_on_mortar_integ != NULL){
          for (int d = 0; d < (P4EST_DIM); d++){
            n_on_mortar_integ[d][is] = n_is[d];
          }
        }
        if (sj_on_mortar_integ != NULL){
          sj_on_mortar_integ[is] = sj_is;
        }       
        if (j_div_sj_mortar_integ != NULL){
          j_div_sj_mortar_integ[is] = J_on_face_integ[i]/sj_is;
        }
      }  

    
    face_mortar_integ_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar_integ[face_mortar]);
  }

  P4EST_FREE(temp);
  P4EST_FREE(J_on_face_integ);

  for (int i = 0; i < (P4EST_DIM); i++)
    for (int j = 0; j < (P4EST_DIM); j++){
      /* P4EST_FREE(dxyz_drst[i][j]); */
      P4EST_FREE(dxyz_drst_on_face_integ[i][j]);
      if (drst_dxyz_on_mortar_integ == NULL){
        P4EST_FREE(drst_dxyz_on_face_integ[i][j]);
        P4EST_FREE(drst_dxyz_times_jac_on_face_integ[i][j]);
      }
    }
}




void
d4est_geometry_compute_xyz
(
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 d4est_geometry_t* d4est_geom,
 int which_tree,
 int deg,
 quadrature_type_t type,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 double* xyz [(P4EST_DIM)]
)
{  
  dgmath_rst_t rst_points = dgmath_get_rst_points(dgmath_jit_dbase,
                                                  deg,
                                                  (P4EST_DIM),
                                                  type);
  
  double* rst [(P4EST_DIM)] = {rst_points.r, rst_points.s, NULL};
#if (P4EST_DIM)==3
  rst[2] = rst_points.t;
#endif

  int volume_nodes = dgmath_get_nodes((P4EST_DIM), deg);
  
  double rst_i [(P4EST_DIM)]; 
  double xyz_i [(P4EST_DIM)];
  for (int i = 0; i < volume_nodes; i++){
    for (int d = 0; d < (P4EST_DIM); d++){
      rst_i[d] = rst[d][i];
    }
    
    d4est_geom->X(d4est_geom, which_tree, q, dq, rst_i, COORDS_INTEG_RST, xyz_i);
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d][i] = xyz_i[d];
    }
    /* printf("xyz[0], xyz[1], xyz[2], rst[0], rst[1], rst[2] = %f,%f,%f,%f,%f,%f\n", xyz_i[0], xyz_i[1], xyz_i[2], rst_i[0], rst_i[1], rst_i[2]); */
    
  }
}

void
d4est_geometry_compute_drst_dxyz
(
 double* dxyz_drst [(P4EST_DIM)][(P4EST_DIM)],
 double* drst_dxyz [(P4EST_DIM)][(P4EST_DIM)],
 int nodes
)
{
 for (int i = 0; i < nodes; i++){
    double xr = dxyz_drst[0][0][i];
    double xs = dxyz_drst[0][1][i];
#if (P4EST_DIM)==3
    double xt = dxyz_drst[0][2][i];
#endif
    
    double yr = dxyz_drst[1][0][i];
    double ys = dxyz_drst[1][1][i];
#if (P4EST_DIM)==3
    double yt = dxyz_drst[1][2][i];
    
    double zr = dxyz_drst[2][0][i];
    double zs = dxyz_drst[2][1][i];
    double zt = dxyz_drst[2][2][i];
#endif
    
    double* rx = &drst_dxyz[0][0][i];
    double* ry = &drst_dxyz[0][1][i];
#if (P4EST_DIM)==3
    double* rz = &drst_dxyz[0][2][i];
#endif
    double* sx = &drst_dxyz[1][0][i];
    double* sy = &drst_dxyz[1][1][i];
#if (P4EST_DIM)==3
    double* sz = &drst_dxyz[1][2][i];
    
    double* tx = &drst_dxyz[2][0][i];
    double* ty = &drst_dxyz[2][1][i];
    double* tz = &drst_dxyz[2][2][i];
#endif

#if (P4EST_DIM) == 3
    double J = xr*(ys*zt-zs*yt)
                 - yr*(xs*zt-zs*xt)
                 + zr*(xs*yt-ys*xt);
    *rx =  (ys*zt - zs*yt)/(J);
    *ry = -(xs*zt - zs*xt)/(J);
    *rz =  (xs*yt - ys*xt)/(J);
    *sx = -(yr*zt - zr*yt)/(J);
    *sy =  (xr*zt - zr*xt)/(J);
    *sz = -(xr*yt - yr*xt)/(J);
    *tx =  (yr*zs - zr*ys)/(J);
    *ty = -(xr*zs - zr*xs)/(J);
    *tz =  (xr*ys - yr*xs)/(J);
#elif (P4EST_DIM) == 2
    double J = -xs*yr + xr*ys;
    *rx = ys/(J);
    *sx =-yr/(J);
    *ry =-xs/(J);
    *sy = xr/(J);
#else
    mpi_abort("DIM must be 2 or 3");
#endif 
  }  
}



void
d4est_geometry_compute_geometric_data_on_mortar_TESTINGONLY
(
 p4est_topidx_t e0_tree,
 p4est_qcoord_t e0_q [(P4EST_DIM)],
 p4est_qcoord_t e0_dq,
 int num_faces_side,
 int num_faces_mortar,
 int* deg_mortar,
 int face_side,
 quadrature_type_t quad_type,
 double* n [(P4EST_DIM)],
 double* sj,
 d4est_geometry_t* d4est_geom,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 double* xyz_storage [(P4EST_DIM)]
)
{

  double* xyz [(P4EST_DIM)];
  if (xyz_storage[0] != NULL){
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d] = xyz_storage[d];
    }
  }
  else {
    int total_face_mortar_nodes = 0;
    for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){
      total_face_mortar_nodes += dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]);
    }
    for (int d = 0; d < (P4EST_DIM); d++){
      xyz[d] = P4EST_ALLOC(double, total_face_mortar_nodes);
    }
  }

  /* Calculate the four "0" corners of 
   * the mortar faces. In the case that
   * there is only one mortar face, these
   * will be the four corners of that face
   */
  
  p4est_qcoord_t q0 [(P4EST_HALF)][(P4EST_DIM)];
  
  for (int j = 0; j < (P4EST_HALF); j++){
    int c = p4est_face_corners[face_side][j];
    for (int d = 0; d < (P4EST_DIM); d++){
      int cd = dgmath_is_child_left_or_right(c, d);
      q0[j][d] = e0_q[d] + cd*e0_dq;
    }
  }

  /* Calculate the vectors that span the face 
   * there will be one in 2-D and two in 3-d */
  
  p4est_qcoord_t dqa [((P4EST_DIM)-1)][(P4EST_DIM)];
  
  for (int d = 0; d < (P4EST_DIM); d++){
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      dqa[dir][d] = (q0[(dir+1)][d] - q0[0][d]);
      if (num_faces_side != num_faces_mortar)
        dqa[dir][d] /= 2;
    }
  }

  if (num_faces_side != num_faces_mortar){
    for (int d = 0; d < (P4EST_DIM); d++){
      for (int c = 0; c < (P4EST_HALF); c++){
        q0[c][d] = q0[0][d];
        for (int dir = 0; dir < (P4EST_DIM) - 1; dir++){
          int cd = dgmath_is_child_left_or_right(c, dir);
          q0[c][d] += cd*dqa[dir][d];
        }
      }
    }
  }
  
  double* a [((P4EST_DIM)-1)];
  /* double* xyz [(P4EST_DIM)]; */
  double* dxda [(P4EST_DIM)][((P4EST_DIM)-1)];
  double dxyz_drs_i [(P4EST_DIM)][((P4EST_DIM)-1)];
  double abc [] = {0.,0.,0.};
  double xyz_i [] = {0.,0.,0.};
  int face_mortar_nodal_stride = 0;
  dgmath_rst_t rst_points;
  
  for (int face_mortar = 0; face_mortar < num_faces_mortar; face_mortar++){


    int face_mortar_nodes = dgmath_get_nodes((P4EST_DIM) - 1, deg_mortar[face_mortar]);
    /* compute the LGL nodes in the directions of the face_mortar vectors */
    double* tmp = P4EST_ALLOC(double,face_mortar_nodes);
    for (int d = 0; d < (P4EST_DIM); d++){
      /* xyz[d] = P4EST_ALLOC(double, face_mortar_nodes); */
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
        dxda[d][dir] = P4EST_ALLOC(double, face_mortar_nodes);
    }
     
    rst_points = dgmath_get_rst_points(dgmath_jit_dbase, deg_mortar[face_mortar], (P4EST_DIM)-1, LOBATTO);
      
    for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
      if (dir == 0)
        a[dir] = rst_points.r;
      else
        a[dir] = rst_points.s;
    }

    
    for (int i = 0; i < face_mortar_nodes; i++){
      if (xyz[0] != NULL){
        for (int d = 0; d < (P4EST_DIM); d++){
          /* get "0" corner of this face_mortar */
          abc[d] = (double)q0[face_mortar][d];
       
          for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
            /* add a fraction of the face_mortar vector in direction dir
           * corresponding to the placement of the LGL node */
            double da = (a[dir][i] + 1.)/2.;
            abc[d] += da*((double)dqa[dir][d]);
            /* rs[dir] = a[dir][i]; */
            /* printf("abc[%d] = %f\n",d, abc[d]); */
          }
        
          abc[d] /= (double)(P4EST_ROOT_LEN);
        }

        /* convert vertex coords to physical coords */
        d4est_geom->X(d4est_geom, e0_tree, (p4est_qcoord_t [(P4EST_DIM)]){0},-1, abc, COORDS_TREE_UNITCUBE , xyz_i);     
        for (int d = 0; d < (P4EST_DIM); d++){
          xyz[d][i] = xyz_i[d];
        }
      }

    }

    /* compute the tangent vectors in direction(s) "dir" */
      for (int d = 0; d < (P4EST_DIM); d++){
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
          dgmath_apply_Dij(dgmath_jit_dbase, xyz[d], ((P4EST_DIM))-1, deg_mortar[face_mortar], dir, dxda[d][dir]);
        }
      }

    if (quad_type == GAUSS){
      for (int d = 0; d < (P4EST_DIM); d++){
        if (xyz[0] != NULL){
          dgmath_interp_GLL_to_GL(dgmath_jit_dbase, xyz[d], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1);
          linalg_copy_1st_to_2nd(tmp, xyz[d], face_mortar_nodes);
        }
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++){
          dgmath_interp_GLL_to_GL(dgmath_jit_dbase, dxda[d][dir], deg_mortar[face_mortar], deg_mortar[face_mortar], tmp, (P4EST_DIM)-1);
          linalg_copy_1st_to_2nd(tmp, dxda[d][dir], face_mortar_nodes);
        }
      }
    }
    
    /* get the normal by taking the cross product of the tangent vectors
     * in 2-d, we take the cross product of the tangent vector and zhat*/
    for (int i = 0; i < face_mortar_nodes; i++){
      double vecs [2][3] = {{0.,0.,0.},{0.,0.,1.}};
      double n_i [] = {0.,0.,0.};
      for (int d = 0; d < (P4EST_DIM); d++)
        for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
          vecs[dir][d] = dxda[d][dir][i];

      linalg_cross_prod
        (
         vecs[0][0],
         vecs[0][1],
         vecs[0][2],
         vecs[1][0],
         vecs[1][1],
         vecs[1][2],
         &(n_i[0]),
         &(n_i[1]),
         &(n_i[2])
        );

      sj[i + face_mortar_nodal_stride] = 0.;
      for (int d = 0; d < (P4EST_DIM); d++){
        /* The normals are backwards for these 2(3) face_mortars in 2-d(3-d) */
        if (face_side == 0 || face_side == 3 || face_side == 4){
          n_i[d] *= -1.;
        }
        sj[i + face_mortar_nodal_stride] += n_i[d]*n_i[d];
      }
      sj[i + face_mortar_nodal_stride] = sqrt(sj[i + face_mortar_nodal_stride]);
      if (n[0] != NULL){
        for (int d = 0; d < (P4EST_DIM); d++){
          n[d][i + face_mortar_nodal_stride] = n_i[d]/sj[i + face_mortar_nodal_stride];
        }
      }
    }

    face_mortar_nodal_stride += dgmath_get_nodes((P4EST_DIM)-1, deg_mortar[face_mortar]);
    for (int d = 0; d < (P4EST_DIM); d++){
      /* P4EST_FREE(xyz[d]); */
      for (int dir = 0; dir < ((P4EST_DIM)-1); dir++)
        P4EST_FREE(dxda[d][dir]);
    }
    P4EST_FREE(tmp);
  }

  if (xyz_storage[0] == NULL){
    for (int d = 0; d < (P4EST_DIM); d++){
      P4EST_FREE(xyz[d]);
    }
  }
  
}

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
    tcoords[0] = dgmath_rtox(coords[0], (double)q0[0], (double)dq)/(double)P4EST_ROOT_LEN;
    tcoords[1] = dgmath_rtox(coords[1], (double)q0[1], (double)dq)/(double)P4EST_ROOT_LEN;
#if (P4EST_DIM)==3
    mpi_assert(coords[2] >= -1. && coords[2] <= 1.);
    tcoords[2] = dgmath_rtox(coords[2], (double)q0[2], (double)dq)/(double)P4EST_ROOT_LEN;
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
