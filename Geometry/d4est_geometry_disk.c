#include <ini.h>
#include <d4est_util.h>
#include <pXest.h>
#include <d4est_geometry.h>
#include <d4est_geometry_disk.h>
#include <p4est_connectivity.h>
#include <zlog.h>


static double
secant_fcn(double x){
  return 1./cos(x);
}


static
int d4est_geometry_disk_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_disk_attr_t* pconfig = (d4est_geometry_disk_attr_t*)user;
  
  if (d4est_util_match_couple(section,pconfig->input_section,name,"R0")) {
    D4EST_ASSERT(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R1")) {
    D4EST_ASSERT(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R2")) {
    D4EST_ASSERT(pconfig->R2 == -1);
    pconfig->R2 = atof(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"compactify_outer_wedge")) {
    D4EST_ASSERT(pconfig->compactify_outer_wedge == -1);
    pconfig->compactify_outer_wedge = atoi(value);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static d4est_geometry_disk_attr_t*
d4est_geometry_5treedisk_input
(
 const char* input_file,
  const char* input_section
)
{
  d4est_geometry_disk_attr_t* input = P4EST_ALLOC(d4est_geometry_disk_attr_t, 1);
  /* D4EST_ASSERT(sizeof(input->input_section) <= 50); */
  /* snprintf (input->input_section, sizeof(input->input_section), "%s", input_section); */
  input->input_section = input_section;
  input->R0 = -1;
  input->R1 = -1;
  input->R2 = -1;
  input->compactify_outer_wedge = -1;
  
  if (ini_parse(input_file, d4est_geometry_disk_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input->input_section, input->R0, -1);
  D4EST_CHECK_INPUT(input->input_section, input->R1, -1);
  
  return input;
}

static d4est_geometry_disk_attr_t*
d4est_geometry_disk_outer_wedge_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_disk_attr_t* input = P4EST_ALLOC(d4est_geometry_disk_attr_t, 1);
  /* D4EST_ASSERT(sizeof(input->input_section) <= 50); */
  /* snprintf (input->input_section, sizeof(input->input_section), "%s", input_section); */

  input->input_section = input_section;
  input->R0 = -1;
  input->R1 = -1;
  input->R2 = -1;
  input->compactify_outer_wedge = -1;


  if (ini_parse(input_file, d4est_geometry_disk_input_handler, input) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input->input_section, input->R1, -1);
  D4EST_CHECK_INPUT(input->input_section, input->R2, -1);
  D4EST_CHECK_INPUT(input->input_section, input->compactify_outer_wedge, -1);
  
  return input;
}



static
void d4est_geometry_5treedisk_linear_map
(
 double xref,
 double refmin,
 double refmax,
 double emin,
 double emax,
 double *x
){
  *x = emin + (emax-emin)*(xref - refmin)/(refmax - refmin);
}

static
void
d4est_geometry_5treedisk_map_cube_to_slab(
                    double xref,
                    double yref,
                    double cmin,
                    double cmax,
                    double emin,
                    double emax,
                    double* x,
                    double* y
){
  double xbar, ybar;
    /* map x from [0,1] to [emin, emax] */
  d4est_geometry_5treedisk_linear_map(xref, 0., 1., emin, emax, &xbar);
    /* map y from [0,1] to [-1,1] */
  d4est_geometry_5treedisk_linear_map(yref, 0., 1., -1., 1., &ybar);
  double xmin = (1. - cmin)*emin + emin*cmin/sqrt(1. + ybar*ybar);
  double xmax = (1. - cmax)*emax + emax*cmax/sqrt(1. + ybar*ybar);
  *x = xmin + (xmax - xmin)*(xbar - emin)/(emax - emin);
  *y = (*x)*ybar;
}



static void
d4est_geometry_disk_outer_wedge_DX(d4est_geometry_t* d4est_geom,
                                                 p4est_topidx_t which_tree,
                                                 p4est_qcoord_t q0 [2],
                                                 p4est_qcoord_t dq,
                                                 const double rst[2], /* [-1,1]^3 */
                                                 double dxyz_drst[2][2]
                                                )
{


  int compactify_outer_wedge = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->compactify_outer_wedge;
  double R1 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = bmin + 1.;
  bmax = bmax + 1.;

  if (!compactify_outer_wedge){
    dxyz_drst[0][0] = -((amax - amin)*M_PI*(R1*(-4 + bmax + bmin + bmax*s - bmin*s) - R2*(-2 + bmax + bmin + bmax*s - bmin*s)))/(16.*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[0][1] = -((bmax - bmin)*(R1 - R2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(2.*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][0] = ((amax - amin)*M_PI*(R1*(-4 + bmax + bmin + bmax*s - bmin*s) - R2*(-2 + bmax + bmin + bmax*s - bmin*s))*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(16.*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][1] = -((bmax - bmin)*(R1 - R2))/(2.*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  }
  else {
    dxyz_drst[0][0] = ((amax - amin)*M_PI*R1*R2)/(4.*(-(R2*(-4 + bmax + bmin + bmax*s - bmin*s)) + R1*(-2 + bmax + bmin + bmax*s - bmin*s))*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[0][1] = (-2*(bmax - bmin)*R1*(R1 - R2)*R2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),2)*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][0] = -((amax - amin)*M_PI*R1*R2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(4.*(-(R2*(-4 + bmax + bmin + bmax*s - bmin*s)) + R1*(-2 + bmax + bmin + bmax*s - bmin*s))*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][1] = (-2*(bmax - bmin)*R1*(R1 - R2)*R2)/(pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),2)*sqrt(pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  }
  
}
  
static void
d4est_geometry_disk_outer_wedge_X(d4est_geometry_t * geom,
                      p4est_topidx_t which_tree,
                      p4est_qcoord_t q0 [2],
                      p4est_qcoord_t dq,
                      const double coords[2],
                      coords_type_t coords_type,
                      double xyz[2]
                      )
{
  d4est_geometry_disk_attr_t* disk = geom->user;
  double tcoords [(P4EST_DIM)];

  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);

  /* transform topog coordinates in [0,1]^3 to cubed sphere vertex space [-1,1]^2 x [1,2]  */
  double abc [3];
  abc[0] = 2*tcoords[0] - 1;
  abc[1] = tcoords[1] + 1;

  double R = -1.;
  double x = tan (abc[0] * M_PI_4);
  if (disk->compactify_outer_wedge){
    double m = (2. - 1.)/((1./disk->R2) - (1./disk->R1));
    double t = (1.*disk->R1 - 2.*disk->R2)/(disk->R1 - disk->R2);
    R = m/(abc[1] - t);
  }
  else {
    R = disk->R1*(2. - abc[1]) + disk->R2*(abc[1] - 1.);
  }
  double q = R / sqrt (x * x + 1.);
  xyz[0] = +q * x;
  xyz[1] = q;
  xyz[2] = 0;
}

static void
d4est_geometry_5treedisk_X(d4est_geometry_t * geom,
                      p4est_topidx_t which_tree,
                      p4est_qcoord_t q0 [2],
                      p4est_qcoord_t dq,
                      const double coords[2],
                      coords_type_t coords_type,
                      double xyz[2]
                      )
{
  double* radii = (double*)geom->user;
  
  double R0 = radii[0];
  double R1 = radii[1];

  double tcoords[2];

  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  double xref = tcoords[0];
  double yref = tcoords[1];
  double x,y;
  
  if (which_tree == 0){
    /* bottom */
    d4est_geometry_5treedisk_map_cube_to_slab(yref, xref, 1., 0., -R1, -R0/sqrt(2), &y, &x);
    x *= -1.;
  }
  else if (which_tree == 1){
    /* left */
    d4est_geometry_5treedisk_map_cube_to_slab(xref, yref, 1., 0., -R1, -R0/sqrt(2), &x, &y);
    y *= -1;
  }
  else if (which_tree == 2){
    /* center */
    d4est_geometry_5treedisk_linear_map(xref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &x);
    d4est_geometry_5treedisk_linear_map(yref, 0., 1., -R0/sqrt(2), R0/sqrt(2), &y);
  }
  else if (which_tree == 3){
    /* right */
    d4est_geometry_5treedisk_map_cube_to_slab(xref, yref, 0, 1, R0/sqrt(2), R1, &x, &y);
  }
  else if (which_tree == 4){
    /* top */
    d4est_geometry_5treedisk_map_cube_to_slab(yref, xref, 0, 1, R0/sqrt(2), R1, &y, &x);
  }
  else{
    SC_ABORT_NOT_REACHED();
  }

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = 0.;
}

static void
d4est_geometry_disk_destroy
(
 d4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


void
d4est_geometry_5treedisk_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 zlog_category_t *c_default,
 d4est_geometry_t* d4est_geom
)
{
  D4EST_ASSERT((P4EST_DIM)==2);
  d4est_geometry_disk_attr_t* input = d4est_geometry_5treedisk_input(input_file, input_section);
  p4est_connectivity_t* conn = p4est_connectivity_new_disk();
  
  d4est_geom->user = input;
  d4est_geom->p4est_conn = conn;
  d4est_geom->destroy = d4est_geometry_disk_destroy;
  d4est_geom->X = d4est_geometry_5treedisk_X;
  d4est_geom->DX = NULL;
  d4est_geom->JAC = NULL;
  
  if (mpirank == 0){
    zlog_info(c_default, "NAME = 5treedisk");
    zlog_info(c_default, "R0 = %.25f", input->R0);
    zlog_info(c_default, "R1 = %.25f", input->R1);
  }
}

void
d4est_geometry_disk_outer_wedge_sj_analytic
(
 d4est_geometry_t * d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [2],
 p4est_qcoord_t dq,
 const double rst[2],
 double* sj
)
{
  int compactify_outer_wedge = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->compactify_outer_wedge;
  double R1 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = bmin + 1.;
  bmax = bmax + 1.;

  if (compactify_outer_wedge == 1){
    if (rst[0] == 1. || rst[0] == -1.){
      *sj = sqrt((4*pow(bmax - bmin,2)*pow(R1,2)*pow(R1 - R2,2)*pow(R2,2))/pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),4));
    }
    else if (rst[1] == 1. || rst[1] == -1.){
      *sj = sqrt((pow(amax - amin,2)*pow(M_PI,2)*pow(R1,2)*pow(R2,2))/(16.*pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),2)));
    }
    else {
      D4EST_ABORT("This is not a face");
    }
  }
  else {
    D4EST_ABORT("Code is not written yet");
  }
}

void
d4est_geometry_disk_outer_wedge_jac_analytic
(
 d4est_geometry_t * d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [2],
 p4est_qcoord_t dq,
 const double rst[2],
 double* j
)
{

  int compactify_outer_wedge = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->compactify_outer_wedge;
  double R1 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = bmin + 1.;
  bmax = bmax + 1.;
  
  if (compactify_outer_wedge == 1){
  *j = -((amax - amin)*(bmax - bmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2))/(2.*pow(-(R2*(-4 + bmax + bmin + bmax*s - bmin*s)) + R1*(-2 + bmax + bmin + bmax*s - bmin*s),3));
  }
  else {

  }
  
}

void
d4est_geometry_disk_outer_wedge_sj_div_jac_analytic
(
 d4est_geometry_t * d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [2],
 p4est_qcoord_t dq,
 const double rst[2],
 double* sj_div_jac
)
{

  int compactify_outer_wedge = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->compactify_outer_wedge;
  double R1 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_disk_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = bmin + 1.;
  bmax = bmax + 1.;
  
  if (compactify_outer_wedge == 1){
    if (rst[0] == 1. || rst[0] == -1.){
  *sj_div_jac = sqrt((16*pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),2))/(pow(amax - amin,2)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)));
    }
    else if (rst[1] == 1. || rst[1] == -1.){
  *sj_div_jac = sqrt(pow(R2*(-4 + bmax + bmin + bmax*s - bmin*s) - R1*(-2 + bmax + bmin + bmax*s - bmin*s),4)/(4.*pow(bmax - bmin,2)*pow(R1,2)*pow(R1 - R2,2)*pow(R2,2)));
    }
    else {
      D4EST_ABORT("This is not a face");
    }
  }
  else {
    D4EST_ABORT("Code is not written yet");
  }
}

void
d4est_geometry_disk_outer_wedge_new_aux
(
 d4est_geometry_t* d4est_geom,
 d4est_geometry_disk_attr_t* disk_attrs
)
{
  p4est_connectivity_t* conn = p4est_connectivity_new_unitsquare();
  d4est_geom->user = disk_attrs;
  d4est_geom->p4est_conn = conn;
  d4est_geom->destroy = d4est_geometry_disk_destroy;
  d4est_geom->X = d4est_geometry_disk_outer_wedge_X;
  d4est_geom->DX = d4est_geometry_disk_outer_wedge_DX;

}

void
d4est_geometry_disk_outer_wedge_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 zlog_category_t *c_default,
 d4est_geometry_t* d4est_geom
)
{
  D4EST_ASSERT((P4EST_DIM)==2);
  d4est_geometry_disk_attr_t* input = d4est_geometry_disk_outer_wedge_input(input_file, input_section);

  d4est_geometry_disk_outer_wedge_new_aux(d4est_geom, input);


  if (mpirank == 0){
    zlog_info(c_default, "NAME = disk_outer_wedge");
    zlog_info(c_default, "R1 = %.25f", input->R1);
    zlog_info(c_default, "R2 = %.25f", input->R2);
    zlog_info(c_default, "compactify outer wedge = %d", input->compactify_outer_wedge);
  }
}
