#include <pXest.h>
#include <d4est_geometry.h>
#include <ini.h>
#include <p8est_connectivity.h>
#include <d4est_util.h>
#include <petscsnes.h>
/* #include <arbquad.h> */

/* #define P4EST_DIM 3 */
#include <d4est_geometry_cubed_sphere.h>
#include <d4est_connectivity_cubed_sphere.h>

static inline double
secant_fcn(double x){
  return 1./cos(x);
}

static
int d4est_geometry_cubed_sphere_get_number_of_regions
(
 d4est_geometry_t* d4est_geom
){
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_7TREE){
    return 2;
  }
  else if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_13TREE){
    return 3;
  }
  else {
    D4EST_ABORT("Not supported yet");
    return -1;
  }
}

static
int d4est_geometry_cubed_sphere_get_region
(
 d4est_geometry_t* d4est_geom,
 p4est_qcoord_t q [(P4EST_DIM)],
 p4est_qcoord_t dq,
 int tree
){
  if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_7TREE){
    if (tree < 6) {   /* inner shell */
      return 0;
    }
    else {           /* center cube */
      return 1;
    }  
  }
  else if (d4est_geom->geom_type == GEOM_CUBED_SPHERE_13TREE){
    if (tree < 6) {         /* outer shell */
      return 0;
    }
    else if (tree < 12) {   /* inner shell */
      return 1;
    }
    else {                  /* center cube */
      return 2;
    }  
  }
  else {
    D4EST_ABORT("Not supported yet");
  }
}

static
int d4est_geometry_cubed_sphere_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_cubed_sphere_attr_t* pconfig = user;
  if (d4est_util_match_couple(section,pconfig->input_section,name,"R0")) {
    D4EST_ASSERT(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
    D4EST_ASSERT(pconfig->R0 > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R1")) {
    D4EST_ASSERT(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
    D4EST_ASSERT(pconfig->R1 > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"R2")) {
    D4EST_ASSERT(pconfig->R2 == -1);
    pconfig->R2 = atof(value);
    D4EST_ASSERT(pconfig->R2 > 0);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"compactify_outer_shell")) {
    D4EST_ASSERT(pconfig->compactify_outer_shell == -1);
    pconfig->compactify_outer_shell = atoi(value);
    D4EST_ASSERT(pconfig->compactify_outer_shell == 0 || pconfig->compactify_outer_shell == 1);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"compactify_inner_shell")) {
    D4EST_ASSERT(pconfig->compactify_inner_shell == -1);
    pconfig->compactify_inner_shell = atoi(value);
    D4EST_ASSERT(pconfig->compactify_inner_shell == 0 || pconfig->compactify_inner_shell == 1);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}





static
void d4est_geometry_cubed_sphere_inner_shell_block_X(
                                                d4est_geometry_t * geom,
                                                p4est_topidx_t which_tree,
                                                p4est_qcoord_t q0 [3],
                                                p4est_qcoord_t dq,
                                                const double coords[3],
                                                coords_type_t coords_type,
                                                double xyz[3]
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere = geom->user;
  double tcoords [3];
  double R;

  /* transform coords of element to [0,1]^3*/
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);

  /* transform topog coordinates in [0,1]^3 to cubed sphere vertex space [-1,1]^2 x [1,2]  */
  double abc [3];
  abc[0] = 2*tcoords[0] - 1;
  abc[1] = 2*tcoords[1] - 1;
  abc[2] = tcoords[2] + 1;

  if (sphere->compactify_inner_shell){
    double m = (2. - 1.)/((1./sphere->R1) - (1./sphere->R0));
    double t = (1.*sphere->R0 - 2.*sphere->R1)/(sphere->R0 - sphere->R1);
    R = m/(abc[2] - t);
  }
  else {
    R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
  }
  double p = 2. - abc[2];
  double tanx = tan (abc[0] * M_PI_4);
  double tany = tan (abc[1] * M_PI_4);
  double x = p * abc[0] + (1. - p) * tanx;
  double y = p * abc[1] + (1. - p) * tany;
  double q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  xyz[0] = +q * x;
  xyz[1] = +q * y;
  xyz[2] = +q;
}

/* static void  */
/* d4est_geometry_cubed_sphere_outer_shell_block_sj_sqr_div_J_sqr( */
/*                                                  d4est_geometry_t* d4est_geom, */
/*                                                  p4est_topidx_t which_tree, */
/*                                                  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*                                                  p4est_qcoord_t dq, */
/*                                                  const double rst[(P4EST_DIM)], /\* [-1,1]^3 *\/ */
/*                                                 double sj[(P4EST_DIM)] */
/* ) */
/* { */
/*   sj_sqr_div_J_sqr_vol[0] = (16*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(cos((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/(pow(amax - amin,2)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)); */

/*   sj_sqr_div_J_sqr_vol[1] = (16*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*(1 + pow(cos((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/(pow(bmax - bmin,2)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)); */

/*   sj_sqr_div_J_sqr_vol[2] = pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),4)/(4.*pow(cmax - cmin,2)*pow(R1,2)*pow(R1 - R2,2)*pow(R2,2)); */
/* } */
  


/* static void  */
/* d4est_geometry_cubed_sphere_outer_shell_block_sj( */
/*                                                  d4est_geometry_t* d4est_geom, */
/*                                                  p4est_topidx_t which_tree, */
/*                                                  p4est_qcoord_t q0 [(P4EST_DIM)], */
/*                                                  p4est_qcoord_t dq, */
/*                                                  const double rst[(P4EST_DIM)], /\* [-1,1]^3 *\/ */
/*                                                 double sj[(P4EST_DIM)] */
/* ) */
/* { */
/*   int compactify_inner_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_inner_shell; */
/*   double R0 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R0; */
/*   double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1; */
/*   /\* Element integration weight x-coordinates in [-1,1]^3 space of the element*\/ */
/*   double r = rst[0]; */
/*   double s = rst[1]; */
/*   double t = rst[2]; */

/*   /\* topological coordinates of element corners *\/ */
/*   double amin = q0[0]; */
/*   double amax = q0[0] + dq; */
/*   double bmin = q0[1]; */
/*   double bmax = q0[1] + dq; */
/*   double cmin = q0[2]; */
/*   double cmax = q0[2] + dq; */

/*   /\* transform element corners to [0,1]^3 topological space *\/ */
/*   amin /= (double)P4EST_ROOT_LEN; */
/*   amax /= (double)P4EST_ROOT_LEN; */
/*   bmin /= (double)P4EST_ROOT_LEN; */
/*   bmax /= (double)P4EST_ROOT_LEN; */
/*   cmin /= (double)P4EST_ROOT_LEN; */
/*   cmax /= (double)P4EST_ROOT_LEN; */

/*   /\* transform element corners to [-1,1]^2 x [1,2] topological space used in the cubed-sphere mapping*\/ */
/*   amin = 2.*amin - 1.; */
/*   amax = 2.*amax - 1.; */
/*   bmin = 2.*bmin - 1.; */
/*   bmax = 2.*bmax - 1.; */
/*   cmin = cmin + 1.; */
/*   cmax = cmax + 1.; */

/*   double sj_sqr_vol [3]; */
  
/*   sj_sqr_vol[0] = (pow(bmax - bmin,2)*pow(cmax - cmin,2)*pow(M_PI,2)*pow(R1,4)*pow(R1 - R2,2)*pow(R2,4)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),4))/(4.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),6)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2)); */

/*   sj_sqr_vol[1] = (pow(amax - amin,2)*pow(cmax - cmin,2)*pow(M_PI,2)*pow(R1,4)*pow(R1 - R2,2)*pow(R2,4)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),4)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(4.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),6)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2)); */

/*   sj_sqr_vol[2] = (pow(amax - amin,2)*pow(bmax - bmin,2)*pow(M_PI,4)*pow(R1,4)*pow(R2,4)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),4)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),4))/(256.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),4)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),3)); */

/*   sj[0] = sqrt(sj_sqr_vol[0]); */
/*   sj[1] = sqrt(sj_sqr_vol[1]); */
/*   sj[2] = sqrt(sj_sqr_vol[2]); */
  
/* } */

static void 
d4est_geometry_cubed_sphere_outer_shell_block_X_aux(
                                                    int compactify_outer_shell,
                                                    double R1,
                                                    double R2,
                                                    p4est_qcoord_t q0 [(P4EST_DIM)],
                                                    p4est_qcoord_t dq,
                                                    const double coords[(P4EST_DIM)],
                                                    coords_type_t coords_type,
                                                    double xyz[(P4EST_DIM)]
)
{
  double tcoords [(P4EST_DIM)];
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);

  /* transform topog coordinates in [0,1]^3 to cubed sphere vertex space [-1,1]^2 x [1,2]  */
  double abc [3];
  abc[0] = 2*tcoords[0] - 1;
  abc[1] = 2*tcoords[1] - 1;
  abc[2] = tcoords[2] + 1;

  double R = -1.;
  double x = tan (abc[0] * M_PI_4);
  double y = tan (abc[1] * M_PI_4);
  if (compactify_outer_shell){
    double m = (2. - 1.)/((1./R2) - (1./R1));
    double t = (1.*R1 - 2.*R2)/(R1 - R2);
    R = m/(abc[2] - t);
  }
  else {
    R = R1*(2. - abc[2]) + R2*(abc[2] - 1.);
  }
  double q = R / sqrt (x * x + y * y + 1.);  
  xyz[0] = +q * x;
  xyz[1] = +q * y;
  xyz[2] = +q;

}


static void 
d4est_geometry_cubed_sphere_outer_shell_block_X(
                                                d4est_geometry_t * d4est_geom,
                                                p4est_topidx_t which_tree,
                                                p4est_qcoord_t q0 [(P4EST_DIM)],
                                                p4est_qcoord_t dq,
                                                const double coords[(P4EST_DIM)],
                                                coords_type_t coords_type,
                                                double xyz[(P4EST_DIM)]
)
{
  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R0 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R0;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;
  
  d4est_geometry_cubed_sphere_outer_shell_block_X_aux
    (
     compactify_outer_shell,
     R1,
     R2,
     q0,
     dq,
     coords,
     coords_type,
     xyz
  );
}


static void
d4est_geometry_cubed_sphere_X(
                              d4est_geometry_t * geom,
                              p4est_topidx_t which_tree,
                              p4est_qcoord_t q0 [3],
                              p4est_qcoord_t dq,
                              const double coords[3],
                              coords_type_t coords_type,
                              double xyz[3]
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  double tcoords [3];
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom->p4est_conn, which_tree, tcoords, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    if (sphere->compactify_outer_shell){
      double m = (2. - 1.)/((1./sphere->R2) - (1./sphere->R1));
      double t = (1.*sphere->R1 - 2.*sphere->R2)/(sphere->R1 - sphere->R2);
      R = m/(abc[2] - t);
    }
    else {
      R = sphere->R1*(2. - abc[2]) + sphere->R2*(abc[2] - 1.);
    }
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;
    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4);
    tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);    
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }  
  switch (which_tree % 6) {
  case 0:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  case 1:                      /* top */
    xyz[0] = +q * x;
    xyz[1] = +q * y;
    xyz[2] = +q;
    break;
  case 2:                      /* back */
    xyz[0] = +q * x;
    xyz[1] = +q;
    xyz[2] = -q * y;
    break;
  case 3:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 5:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED();
  }
}

static void
d4est_geometry_cubed_sphere_with_sphere_hole_X(
                              d4est_geometry_t * geom,
                              p4est_topidx_t which_tree,
                              p4est_qcoord_t q0 [3],
                              p4est_qcoord_t dq,
                              const double coords[3],
                              coords_type_t coords_type,
                              double xyz[3]
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  double tcoords [3];
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom->p4est_conn, which_tree, tcoords, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    if (sphere->compactify_outer_shell){
      double m = (2. - 1.)/((1./sphere->R2) - (1./sphere->R1));
      double t = (1.*sphere->R1 - 2.*sphere->R2)/(sphere->R1 - sphere->R2);
      R = m/(abc[2] - t);
    }
    else {
      R = sphere->R1*(2. - abc[2]) + sphere->R2*(abc[2] - 1.);
    }
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    if (sphere->compactify_inner_shell){
      double m = (2. - 1.)/((1./sphere->R1) - (1./sphere->R0));
      double t = (1.*sphere->R0 - 2.*sphere->R1)/(sphere->R0 - sphere->R1);
      R = m/(abc[2] - t);
    }
    else {
      R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
    }
    q = R / sqrt (x * x + y * y + 1.);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }  
  switch (which_tree % 6) {
  case 0:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  case 1:                      /* top */
    xyz[0] = +q * x;
    xyz[1] = +q * y;
    xyz[2] = +q;
    break;
  case 2:                      /* back */
    xyz[0] = +q * x;
    xyz[1] = +q;
    xyz[2] = -q * y;
    break;
  case 3:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 5:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED();
  }
}
static void
d4est_geometry_cubed_sphere_7tree_X(
                              d4est_geometry_t * geom,
                              p4est_topidx_t which_tree,
                              p4est_qcoord_t q0 [3],
                              p4est_qcoord_t dq,
                              const double coords[3],
                              coords_type_t coords_type,
                              double xyz[3]
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere = geom->user;

  double              x, y, R, q;
  double              abc[3];

  double tcoords [3];
  d4est_geometry_get_tree_coords_in_range_0_to_1(q0, dq, coords, coords_type, tcoords);
  
  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom->p4est_conn, which_tree, tcoords, abc);

  /* assert that input points are in the expected range */
  if (which_tree < 6) {   /* inner shell */
    if (sphere->compactify_inner_shell){
      double m = (2. - 1.)/((1./sphere->R1) - (1./sphere->R0));
      double t = (1.*sphere->R0 - 2.*sphere->R1)/(sphere->R0 - sphere->R1);
      R = m/(abc[2] - t);
    }
    else {
      R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
    }
    double p = 2. - abc[2];
    double tanx = tan (abc[0] * M_PI_4);
    double tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }  
  switch (which_tree % 6) {
  case 0:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  case 1:                      /* top */
    xyz[0] = +q * x;
    xyz[1] = +q * y;
    xyz[2] = +q;
    break;
  case 2:                      /* back */
    xyz[0] = +q * x;
    xyz[1] = +q;
    xyz[2] = -q * y;
    break;
  case 3:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 5:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED();
  }
}

static void
d4est_geometry_cubed_sphere_inner_shell_block_DX(d4est_geometry_t* d4est_geom,
                                                 p4est_topidx_t which_tree,
                                                 p4est_qcoord_t q0 [3],
                                                 p4est_qcoord_t dq,
                                                 const double rst[3], /* [-1,1]^3 */
                                                 double dxyz_drst[3][3]
                                                )
{
  int compactify_inner_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_inner_shell;
  double R0 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R0;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space used in the cubed-sphere mapping*/
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;  

  /* FROM MATHEMATICA */
  
  if (compactify_inner_shell) {
    dxyz_drst[0][0] = -((amax - amin)*R0*R1*(M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)) + 8*(-2*(-4 + cmax + cmin + cmax*t - cmin*t) + (M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(32.*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][1] = -((bmax - bmin)*M_PI*R0*R1*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*(-((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(8.*(-(R1*(-4 + cmax + cmin + cmax*t - cmin*t)) + R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][2] = -((cmax - cmin)*R0*R1*((R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 4*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*(amax + amin + amax*r - amin*r - 2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.) + 16*(R0 - R1)*(-((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(8.*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][0] = -((amax - amin)*M_PI*R0*R1*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*(-((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.))/(8.*(-(R1*(-4 + cmax + cmin + cmax*t - cmin*t)) + R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][1] = -((bmax - bmin)*R0*R1*(M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)*((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)) + 8*(-2*(-4 + cmax + cmin + cmax*t - cmin*t) + (M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(32.*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][2] = -((cmax - cmin)*R0*R1*((R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 4*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t))*(bmax + bmin + bmax*s - bmin*s - 2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.) + 16*(R0 - R1)*(-((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(8.*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][0] = -((amax - amin)*M_PI*R0*R1*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(2.*sqrt(2)*(-(R1*(-4 + cmax + cmin + cmax*t - cmin*t)) + R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(-2*(-5 + cmax + cmin + cmax*t - cmin*t) + (-2 + cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + (-2 + cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2),1.5));
    dxyz_drst[2][1] = -((bmax - bmin)*M_PI*R0*R1*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*sqrt(2)*(-(R1*(-4 + cmax + cmin + cmax*t - cmin*t)) + R0*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(-2*(-5 + cmax + cmin + cmax*t - cmin*t) + (-2 + cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + (-2 + cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2),1.5));
    dxyz_drst[2][2] = -((cmax - cmin)*R0*R1*((-(R1*(-4 + cmax + cmin + cmax*t - cmin*t)) + R0*(-2 + cmax + cmin + cmax*t - cmin*t))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 4*(R0 - R1)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(2.*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R0*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
  }
  else {
    dxyz_drst[0][0] = ((amax - amin)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*(M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)) + 8*(-2*(-4 + cmax + cmin + cmax*t - cmin*t) + (M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(128.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][1] = -((bmax - bmin)*M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*(-((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][2] = ((cmax - cmin)*(-((R0*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t))*((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))) - 4*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*(amax + amin + amax*r - amin*r - 2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.) - 16*(R0 - R1)*(-((amax + amin + amax*r - amin*r)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][0] = -((amax - amin)*M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*(-((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][1] = ((bmax - bmin)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*(M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)*((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)) + 8*(-2*(-4 + cmax + cmin + cmax*t - cmin*t) + (M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(128.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][2] = ((cmax - cmin)*(-((R0*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t))*((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t) - 2*(-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))) - 4*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*(bmax + bmin + bmax*s - bmin*s - 2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.) - 16*(R0 - R1)*(-((bmax + bmin + bmax*s - bmin*s)*(-4 + cmax + cmin + cmax*t - cmin*t))/4. + ((-2 + cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][0] = -((amax - amin)*M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][1] = -((bmax - bmin)*M_PI*(-2 + cmax + cmin + cmax*t - cmin*t)*(-(R0*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(32.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][2] = ((cmax - cmin)*((R0*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t))*(-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 4*(R0 - R1)*(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/(8.*pow(5 + cmin*(-1 + t) - cmax*(1 + t) + ((-2 + cmax + cmin + cmax*t - cmin*t)*(pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
  }
}


static void
d4est_geometry_cubed_sphere_outer_shell_block_DRDX_JAC
(
 d4est_geometry_t* d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double rst[3], /* [-1,1]^3 */
 double drst_dxyz_times_jac[3][3]
)
{

  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;

  if(compactify_outer_shell){
 
drst_dxyz_times_jac[0][0] = -((bmax - bmin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

drst_dxyz_times_jac[0][1] = 0;

drst_dxyz_times_jac[0][2] = ((bmax - bmin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

drst_dxyz_times_jac[1][0] = 0;

drst_dxyz_times_jac[1][1] = -((amax - amin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

drst_dxyz_times_jac[1][2] = ((amax - amin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

drst_dxyz_times_jac[2][0] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(16.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));

drst_dxyz_times_jac[2][1] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));

drst_dxyz_times_jac[2][2] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(16.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));
 
  }
  else {

    drst_dxyz_times_jac[0][0] = ((bmax - bmin)*(cmax - cmin)*M_PI*(R1 - R2)*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(32.*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz_times_jac[0][1] = 0;

    drst_dxyz_times_jac[0][2] = -((bmax - bmin)*(cmax - cmin)*M_PI*(R1 - R2)*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(32.*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz_times_jac[1][0] = 0;

    drst_dxyz_times_jac[1][1] = ((amax - amin)*(cmax - cmin)*M_PI*(R1 - R2)*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/(32.*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz_times_jac[1][2] = -((amax - amin)*(cmax - cmin)*M_PI*(R1 - R2)*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(32.*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz_times_jac[2][0] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(256.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));

    drst_dxyz_times_jac[2][1] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(256.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));

    drst_dxyz_times_jac[2][2] = ((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(256.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2));

  }
  
}


static void
d4est_geometry_cubed_sphere_outer_shell_block_DRDX
(
 d4est_geometry_t* d4est_geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double rst[3], /* [-1,1]^3 */
 double drst_dxyz[3][3]
)
{

  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;

  if(compactify_outer_shell){

    drst_dxyz[0][0] = (4*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(cos((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((amax - amin)*M_PI*R1*R2);

    drst_dxyz[0][1] = 0;

    drst_dxyz[0][2] = (-2*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*sin((M_PI*(amax + amin + amax*r - amin*r))/4.)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((amax - amin)*M_PI*R1*R2);

    drst_dxyz[1][0] = 0;

    drst_dxyz[1][1] = (4*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(cos((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((bmax - bmin)*M_PI*R1*R2);

    drst_dxyz[1][2] = (-2*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*sin((M_PI*(bmax + bmin + bmax*s - bmin*s))/4.)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((bmax - bmin)*M_PI*R1*R2);

    drst_dxyz[2][0] = -(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(2.*(cmax - cmin)*R1*(R1 - R2)*R2*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz[2][1] = -(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*(cmax - cmin)*R1*(R1 - R2)*R2*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz[2][2] = -pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)/(2.*(cmax - cmin)*R1*(R1 - R2)*R2*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  }
  else {
    drst_dxyz[0][0] = (-16*pow(cos((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((amax - amin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t)));

    drst_dxyz[0][1] = 0;

    drst_dxyz[0][2] = (8*sin((M_PI*(amax + amin + amax*r - amin*r))/4.)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((amax - amin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t)));

    drst_dxyz[1][0] = 0;

    drst_dxyz[1][1] = (-16*pow(cos((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((bmax - bmin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t)));

    drst_dxyz[1][2] = (8*sin((M_PI*(bmax + bmin + bmax*s - bmin*s))/4.)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)))/((bmax - bmin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t)));

    drst_dxyz[2][0] = (-2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/((cmax - cmin)*(R1 - R2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz[2][1] = (-2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/((cmax - cmin)*(R1 - R2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drst_dxyz[2][2] = -2/((cmax - cmin)*(R1 - R2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  }
  
}


static void
d4est_geometry_cubed_sphere_outer_shell_block_drdxjacdrdx(
                                                          d4est_geometry_t* d4est_geom,
                                                          p4est_topidx_t which_tree,
                                                          p4est_qcoord_t q0 [3],
                                                          p4est_qcoord_t dq,
                                                          const double rst[3], /* [-1,1]^3 */
                                                          double drdxJdrdx [3][3]
                                                         )
{
  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;

  if(compactify_outer_shell){
    drdxJdrdx[0][0] = (-2*(bmax - bmin)*(cmax - cmin)*R1*(R1 - R2)*R2*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/((amax - amin)*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[0][1] = (-2*(cmax - cmin)*R1*(R1 - R2)*R2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[0][2] = 0;

    drdxJdrdx[1][0] = (-2*(cmax - cmin)*R1*(R1 - R2)*R2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[1][1] = (-2*(amax - amin)*(cmax - cmin)*R1*(R1 - R2)*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/((bmax - bmin)*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[1][2] = 0;

    drdxJdrdx[2][0] = 0;

    drdxJdrdx[2][1] = 0;

    drdxJdrdx[2][2] = -((amax - amin)*(bmax - bmin)*pow(M_PI,2)*R1*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(32.*(cmax - cmin)*(R1 - R2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    
  }
  else {

    drdxJdrdx[0][0] = -((bmax - bmin)*(cmax - cmin)*(R1 - R2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(2.*(amax - amin)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[0][1] = -((cmax - cmin)*(R1 - R2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[0][2] = 0;

    drdxJdrdx[1][0] = -((cmax - cmin)*(R1 - R2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[1][1] = -((amax - amin)*(cmax - cmin)*(R1 - R2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2))/(2.*(bmax - bmin)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));

    drdxJdrdx[1][2] = 0;

    drdxJdrdx[2][0] = 0;

    drdxJdrdx[2][1] = 0;

    drdxJdrdx[2][2] = -((amax - amin)*(bmax - bmin)*pow(M_PI,2)*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t),2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(128.*(cmax - cmin)*(R1 - R2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    
  }

}




static void
d4est_geometry_cubed_sphere_outer_shell_block_jac(d4est_geometry_t* d4est_geom,
                                                 p4est_topidx_t which_tree,
                                                 p4est_qcoord_t q0 [3],
                                                 p4est_qcoord_t dq,
                                                 const double rst[3], /* [-1,1]^3 */
                                                 double* jac
                                                )
{

  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;

  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;

  if(compactify_outer_shell){
    *jac = -((amax - amin)*(bmax - bmin)*(cmax - cmin)*pow(M_PI,2)*pow(R1,3)*(R1
               - R2)*pow(R2,3)*pow(secant_fcn((M_PI*(amax + amin + amax*r -
               amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s -
               bmin*s))/8.),2))/(8.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax
               + cmin + cmax*t - cmin*t),4)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s -
               bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
  }
  else {
    *jac = -((amax - amin)*(bmax - bmin)*(cmax - cmin)*pow(M_PI,2)*(R1
           - R2)*pow(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax
           + cmin + cmax*t - cmin*t),2)*pow(secant_fcn((M_PI*(amax + amin +
           amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s -
           bmin*s))/8.),2))/(512.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s -
           bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
  }
  
}


/* drst_dxyz_times_jac[0][0] = -((bmax - bmin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 */
/* - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + */
/* cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s */
/* - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))); */
/* drst_dxyz_times_jac[0][1] = 0; */
/* drst_dxyz_times_jac[0][2] = 0; */
/* drst_dxyz_times_jac[1][0] = -((amax - amin)*(cmax - cmin)*M_PI*pow(R1,2)*(R1 */
/* - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin + amax*r - */
/* amin*r))/8.),2))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + */
/* cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s */
/* - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))); */
/* drst_dxyz_times_jac[1][1] = ((bmax - bmin)*(cmax - */
/* cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(bmax */
/* + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - */
/* amin*r))/8.))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + */
/* cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s */
/* - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))); */
/* drst_dxyz_times_jac[1][2] = ((amax - amin)*(cmax - */
/* cmin)*M_PI*pow(R1,2)*(R1 - R2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax */
/* + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.))/(2.*pow(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + */
/* cmax + cmin + cmax*t - cmin*t),3)*(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s */
/* - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2))); */
/* drst_dxyz_times_jac[2][0] = ((amax - amin)*(bmax - */
/* bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin */
/* + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(16.*pow(R2*(-4 */
/* + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t */
/* - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2)); */
/* drst_dxyz_times_jac[2][1] = ((amax - amin)*(bmax - */
/* bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin */
/* + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(R2*(-4 */
/* + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t */
/* - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2)); */
/* drst_dxyz_times_jac[2][2] = ((amax - amin)*(bmax - */
/* bmin)*pow(M_PI,2)*pow(R1,2)*pow(R2,2)*pow(secant_fcn((M_PI*(amax + amin */
/* + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - */
/* bmin*s))/8.),2))/(16.*pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + */
/* cmax + cmin + cmax*t - cmin*t),2)*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s */
/* - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),2)); */


static void
d4est_geometry_cubed_sphere_outer_shell_block_DX_aux(int compactify_outer_shell,
                                                     double R1,
                                                     double R2,
                                                 p4est_qcoord_t q0 [3],
                                                 p4est_qcoord_t dq,
                                                 const double rst[3], /* [-1,1]^3 */
                                                 double dxyz_drst[3][3]
                                                )
{



  /* Element integration weight x-coordinates in [-1,1]^3 space of the element*/
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];

  /* topological coordinates of element corners */
  double amin = q0[0];
  double amax = q0[0] + dq;
  double bmin = q0[1];
  double bmax = q0[1] + dq;
  double cmin = q0[2];
  double cmax = q0[2] + dq;

  /* transform element corners to [0,1]^3 topological space */
  amin /= (double)P4EST_ROOT_LEN;
  amax /= (double)P4EST_ROOT_LEN;
  bmin /= (double)P4EST_ROOT_LEN;
  bmax /= (double)P4EST_ROOT_LEN;
  cmin /= (double)P4EST_ROOT_LEN;
  cmax /= (double)P4EST_ROOT_LEN;

  /* transform element corners to [-1,1]^2 x [1,2] topological space */
  amin = 2.*amin - 1.;
  amax = 2.*amax - 1.;
  bmin = 2.*bmin - 1.;
  bmax = 2.*bmax - 1.;
  cmin = cmin + 1.;
  cmax = cmax + 1.;  
  
  if(compactify_outer_shell){
    dxyz_drst[0][0] = ((amax - amin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[0][1] = -((bmax - bmin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[0][2] = (-2*(cmax - cmin)*R1*(R1 - R2)*R2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][0] = -((amax - amin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[1][1] = ((bmax - bmin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[1][2] = (-2*(cmax - cmin)*R1*(R1 - R2)*R2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[2][0] = -((amax - amin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[2][1] = -((bmax - bmin)*M_PI*R1*R2*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(4.*(-(R2*(-4 + cmax + cmin + cmax*t - cmin*t)) + R1*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[2][2] = (-2*(cmax - cmin)*R1*(R1 - R2)*R2)/(pow(R2*(-4 + cmax + cmin + cmax*t - cmin*t) - R1*(-2 + cmax + cmin + cmax*t - cmin*t),2)*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  }
  else {
    dxyz_drst[0][0] = -((amax - amin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[0][1] = ((bmax - bmin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[0][2] = -((cmax - cmin)*(R1 - R2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(2.*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[1][0] = ((amax - amin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[1][1] = -((bmax - bmin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[1][2] = -((cmax - cmin)*(R1 - R2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(2.*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
    dxyz_drst[2][0] = ((amax - amin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[2][1] = ((bmax - bmin)*M_PI*(R1*(-4 + cmax + cmin + cmax*t - cmin*t) - R2*(-2 + cmax + cmin + cmax*t - cmin*t))*pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/(16.*pow(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2),1.5));
    dxyz_drst[2][2] = -((cmax - cmin)*(R1 - R2))/(2.*sqrt(pow(secant_fcn((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2) + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2)));
  } 
}


static void
d4est_geometry_cubed_sphere_outer_shell_block_DX(d4est_geometry_t* d4est_geom,
                                                 p4est_topidx_t which_tree,
                                                 p4est_qcoord_t q0 [3],
                                                 p4est_qcoord_t dq,
                                                 const double rst[3], /* [-1,1]^3 */
                                                 double dxyz_drst[3][3]
                                                )
{
  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;

  d4est_geometry_cubed_sphere_outer_shell_block_DX_aux(compactify_outer_shell,
                                                       R1,
                                                       R2,
                                                       q0,
                                                       dq,
                                                       rst,
                                                       dxyz_drst);
  
}

static void
d4est_geometry_cubed_sphere_innerouter_shell_X
(
 d4est_geometry_t * geom,
 p4est_topidx_t which_tree,
 p4est_qcoord_t q0 [(P4EST_DIM)],
 p4est_qcoord_t dq,
 const double coords[(P4EST_DIM)],
 coords_type_t coords_type,
 double xyz[(P4EST_DIM)]
)
{
  if (which_tree == 0){
    d4est_geometry_cubed_sphere_outer_shell_block_X(geom,
                                                    which_tree,
                                                    q0,
                                                    dq,
                                                    coords,
                                                    coords_type,
                                                    xyz);
  }
  else if (which_tree == 1){
    d4est_geometry_cubed_sphere_inner_shell_block_X(geom,
                                                    which_tree,
                                                    q0,
                                                    dq,
                                                    coords,
                                                    coords_type,
                                                    xyz);
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: which_tree should only be 0 or 1 for innerouter_shell");
  }
}

static void
d4est_geometry_cubed_sphere_innerouter_shell_DX(d4est_geometry_t* d4est_geom,
                         p4est_topidx_t which_tree,
                         p4est_qcoord_t q0 [(P4EST_DIM)],
                         p4est_qcoord_t dq,
                         const double rst[(P4EST_DIM)], /* [-1,1]^3 */
                         double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
                                                )
{
  if (which_tree == 0){
    d4est_geometry_cubed_sphere_outer_shell_block_DX(d4est_geom,
                                                     which_tree,
                                                     q0,
                                                     dq,
                                                     rst,
                                                     dxyz_drst);
  }
  else if (which_tree == 1){
    d4est_geometry_cubed_sphere_inner_shell_block_DX(d4est_geom,
                                                     which_tree,
                                                     q0,
                                                     dq,
                                                     rst,
                                                     dxyz_drst);
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: which_tree should only be 0 or 1 for innerouter_shell");
  }
}

static void
d4est_geometry_cubed_sphere_DX_aux_rotate
(
 int which_tree,
 double dxyz_drst_top [(P4EST_DIM)][(P4EST_DIM)],
 double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
)
{
  switch (which_tree % 6) {
  case 0:                      /* front */
    /* xyz[0] = +q * x; */
    /* xyz[1] = -q; */
    /* xyz[2] = +q * y; */
    for (int d = 0; d < (P4EST_DIM); d++){
      dxyz_drst[0][d] = dxyz_drst_top[0][d];
      dxyz_drst[1][d] = -dxyz_drst_top[2][d];
      dxyz_drst[2][d] = dxyz_drst_top[1][d];
    }
    break;
 case 1:                      /* top */
   /* xyz[0] = +q * x; */
   /* xyz[1] = +q * y; */
   /* xyz[2] = +q; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_top[0][d];
     dxyz_drst[1][d] = dxyz_drst_top[1][d];
     dxyz_drst[2][d] = dxyz_drst_top[2][d];
    
   }
   break;
 case 2:                      /* back */
   /* xyz[0] = +q * x; */
   /* xyz[1] = +q; */
   /* xyz[2] = -q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_top[0][d];
     dxyz_drst[1][d] = dxyz_drst_top[2][d];
     dxyz_drst[2][d] = -dxyz_drst_top[1][d];
   }
   break;
 case 3:                      /* right */
   /* xyz[0] = +q; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = -q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_top[2][d];
     dxyz_drst[1][d] = -dxyz_drst_top[0][d];
     dxyz_drst[2][d] = -dxyz_drst_top[1][d];
   }
   break;
 case 4:                      /* bottom */
   /* xyz[0] = -q * y; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = -q; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = -dxyz_drst_top[1][d];
     dxyz_drst[1][d] = -dxyz_drst_top[0][d];
     dxyz_drst[2][d] = -dxyz_drst_top[2][d];      
   }
   break;
 case 5:                      /* left */
   /* xyz[0] = -q; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = +q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = -dxyz_drst_top[2][d];
     dxyz_drst[1][d] = -dxyz_drst_top[0][d];
     dxyz_drst[2][d] = dxyz_drst_top[1][d];
   }
   break;
  default:
    SC_ABORT_NOT_REACHED();
  }
}


static void
d4est_geometry_cubed_sphere_DX(d4est_geometry_t* d4est_geom,
                         p4est_topidx_t which_tree,
                         p4est_qcoord_t q0 [(P4EST_DIM)],
                         p4est_qcoord_t dq,
                         const double rst[(P4EST_DIM)], /* [-1,1]^3 */
                         double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
                                                )
{
  double dxyz_drst_temp [(P4EST_DIM)][(P4EST_DIM)];


  if (which_tree < 6){
    d4est_geometry_cubed_sphere_outer_shell_block_DX(d4est_geom,
                                                     which_tree,
                                                     q0,
                                                     dq,
                                                     rst,
                                                     dxyz_drst_temp);
  }
  else if (which_tree < 12){
    d4est_geometry_cubed_sphere_inner_shell_block_DX(d4est_geom,
                                                     which_tree,
                                                     q0,
                                                     dq,
                                                     rst,
                                                     dxyz_drst_temp);
  }
  else {
    for (int d1 = 0; d1 < (P4EST_DIM); d1++){
      for (int d2 = 0; d2 < (P4EST_DIM); d2++){
        dxyz_drst[d1][d2] = 0.;
      }
    }

    double Clength = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->Clength;
    
    /* topological coordinates of element corners */
    double amin = q0[0];
    double amax = q0[0] + dq;
    double bmin = q0[1];
    double bmax = q0[1] + dq;
    double cmin = q0[2];
    double cmax = q0[2] + dq;

    /* transform element corners to [0,1]^3 topological space */
    amin /= (double)P4EST_ROOT_LEN;
    amax /= (double)P4EST_ROOT_LEN;
    bmin /= (double)P4EST_ROOT_LEN;
    bmax /= (double)P4EST_ROOT_LEN;
    cmin /= (double)P4EST_ROOT_LEN;
    cmax /= (double)P4EST_ROOT_LEN;

    /* transform element corners to [-Clength,Clength]^3*/
    amin = 2.*Clength*amin - Clength;
    amax = 2.*Clength*amax - Clength;
    bmin = 2.*Clength*bmin - Clength;
    bmax = 2.*Clength*bmax - Clength;
    cmin = 2.*Clength*cmin - Clength;
    cmax = 2.*Clength*cmax - Clength;

    /* x = (amax-amin)*(r+1)/2 + amin */    
    /* y = (bmax-bmin)*(s+1)/2 + bmin */    
    /* z = (cmax-cmin)*(t+1)/2 + cmin */    
    dxyz_drst[0][0] = .5*(amax - amin);
    dxyz_drst[1][1] = .5*(bmax - bmin);
    dxyz_drst[2][2] = .5*(cmax - cmin);
    return;
  }

  d4est_geometry_cubed_sphere_DX_aux_rotate
    (
     which_tree,
     dxyz_drst_temp,
     dxyz_drst
    );
  
}

static void
d4est_geometry_cubed_sphere_with_sphere_hole_DX(d4est_geometry_t* d4est_geom,
                         p4est_topidx_t which_tree,
                         p4est_qcoord_t q0 [(P4EST_DIM)],
                         p4est_qcoord_t dq,
                         const double rst[(P4EST_DIM)], /* [-1,1]^3 */
                         double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
                                                )
{
  double dxyz_drst_temp [(P4EST_DIM)][(P4EST_DIM)];
  int compactify_outer_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_outer_shell;
  int compactify_inner_shell = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->compactify_inner_shell;
  double R0 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R0;
  double R1 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R1;
  double R2 = (( d4est_geometry_cubed_sphere_attr_t*)(d4est_geom->user))->R2;
  
  /* outer spherical shell wedges */
  if (which_tree < 6){
    d4est_geometry_cubed_sphere_outer_shell_block_DX_aux
      (
       compactify_outer_shell,
       R1,
       R2,
       q0,
       dq,
       rst,
       dxyz_drst_temp
      );
  }
  /* inner spherical shell wedges */
  else if (which_tree < 12){
    d4est_geometry_cubed_sphere_outer_shell_block_DX_aux
      (
       compactify_inner_shell,
       R0,
       R1,
       q0,
       dq,
       rst,
       dxyz_drst_temp
      );
  }

  d4est_geometry_cubed_sphere_DX_aux_rotate
    (
     which_tree,
     dxyz_drst_temp,
     dxyz_drst
    );
  
}



static void
d4est_geometry_cubed_sphere_7tree_DX(d4est_geometry_t* d4est_geom,
                         p4est_topidx_t which_tree,
                         p4est_qcoord_t q0 [(P4EST_DIM)],
                         p4est_qcoord_t dq,
                         const double rst[(P4EST_DIM)], /* [-1,1]^3 */
                         double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]
                                                )
{
  d4est_geometry_cubed_sphere_DX(d4est_geom,
                                 which_tree + 6,
                                 q0,
                                 dq,
                                 rst,
                                 dxyz_drst);
}

static
d4est_geometry_cubed_sphere_attr_t*
d4est_geometry_cubed_sphere_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = -1;
  sphere_attrs->R2 = -1;
  sphere_attrs->compactify_outer_shell = -1;
  sphere_attrs->compactify_inner_shell = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R2, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_outer_shell, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_inner_shell, -1);

  /* variables useful for the center cube */
  sphere_attrs->Clength = sphere_attrs->R0 / sqrt (3.);

  return sphere_attrs;
}


static
d4est_geometry_cubed_sphere_attr_t*
d4est_geometry_cubed_sphere_inner_shell_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = -1;
  sphere_attrs->R2 = -1;
  sphere_attrs->compactify_inner_shell = -1;
  sphere_attrs->compactify_outer_shell = -1;
    
  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_inner_shell, -1);
  /* D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_outer_shell, -1); */

  /* variables useful for the center cube */
  sphere_attrs->Clength = sphere_attrs->R0 / sqrt (3.);

  return sphere_attrs;
}


static
d4est_geometry_cubed_sphere_attr_t*
d4est_geometry_cubed_sphere_innerouter_shell_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = -1;
  sphere_attrs->R2 = -1;
  sphere_attrs->compactify_inner_shell = -1;
  sphere_attrs->compactify_outer_shell = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R2, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_inner_shell, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_outer_shell, -1);

  /* variables useful for the center cube */
  sphere_attrs->Clength = sphere_attrs->R0 / sqrt (3.);

  return sphere_attrs;
}

static
d4est_geometry_cubed_sphere_attr_t*
d4est_geometry_cubed_sphere_7tree_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = -1;
  /* sphere_attrs->R2 = -1; */
  sphere_attrs->compactify_inner_shell = -1;
  sphere_attrs->compactify_outer_shell = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1, -1);
  /* D4EST_CHECK_INPUT(input_section, sphere_attrs->R2, -1); */
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_inner_shell, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_outer_shell, -1);

  /* variables useful for the center cube */
  sphere_attrs->Clength = sphere_attrs->R0 / sqrt (3.);

  return sphere_attrs;
}


static
d4est_geometry_cubed_sphere_attr_t*
d4est_geometry_cubed_sphere_outer_shell_input
(
 const char* input_file,
 const char* input_section
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = P4EST_ALLOC(d4est_geometry_cubed_sphere_attr_t, 1);
  snprintf (sphere_attrs->input_section, sizeof(sphere_attrs->input_section), "%s", input_section);
  sphere_attrs->R0 = -1;
  sphere_attrs->R1 = -1;
  sphere_attrs->R2 = -1;
  sphere_attrs->compactify_outer_shell = -1;
  sphere_attrs->compactify_inner_shell = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  /* D4EST_CHECK_INPUT(input_section, sphere_attrs->R0, -1); */
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R1, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->R2, -1);
  D4EST_CHECK_INPUT(input_section, sphere_attrs->compactify_outer_shell, -1);


  return sphere_attrs;
}

static void
d4est_geometry_cubed_sphere_destroy
(
 d4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}

void
d4est_geometry_cubed_sphere_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_input(input_file, input_section);
  p4est_connectivity_t* conn = p8est_connectivity_new_sphere();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_X;
  d4est_geom->DX = d4est_geometry_cubed_sphere_DX;
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;
  
  if (mpirank == 0){
    printf("%s: NAME = cubed sphere\n", printf_prefix);
    printf("%s: R0 = %.25f\n", printf_prefix, sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix, sphere_attrs->R1);
    printf("%s: R2 = %.25f\n", printf_prefix, sphere_attrs->R2);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix, sphere_attrs->compactify_outer_shell);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix, sphere_attrs->compactify_inner_shell);
  }
}

void
d4est_geometry_cubed_sphere_7tree_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_7tree_input(input_file, input_section);
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_7tree();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_7tree_X;
  d4est_geom->DX = d4est_geometry_cubed_sphere_7tree_DX;
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;

  if (mpirank == 0){
    printf("%s: NAME = cubed sphere (7 tree)\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}


void
d4est_geometry_cubed_sphere_innerouter_shell_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_innerouter_shell_input(input_file, input_section);
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_innerouter_shell();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_innerouter_shell_X;
  d4est_geom->DX = d4est_geometry_cubed_sphere_innerouter_shell_DX;
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;
  
  if (mpirank == 0){
    printf("%s: NAME = cubed sphere innerouter shell\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: R2 = %.25f\n", printf_prefix , sphere_attrs->R2);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix , sphere_attrs->compactify_outer_shell);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}

void
d4est_geometry_cubed_sphere_inner_shell_block_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_inner_shell_input(input_file, input_section);
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_inner_shell_block_X;
  /* if (d4est_geom->X_mapping_type == MAP_ANALYTIC) */
  d4est_geom->DX = d4est_geometry_cubed_sphere_inner_shell_block_DX;
  /* else */
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;

  if (mpirank == 0){
    printf("%s: NAME = cubed sphere inner shell block\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}


void
d4est_geometry_cubed_sphere_outer_shell_block_new_aux
(
 d4est_geometry_t* d4est_geom,
 d4est_geometry_cubed_sphere_attr_t* sphere_attrs
)
{
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_outer_shell_block_X;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
  d4est_geom->DX = d4est_geometry_cubed_sphere_outer_shell_block_DX;
  d4est_geom->JAC = d4est_geometry_cubed_sphere_outer_shell_block_jac;

}

void
d4est_geometry_cubed_sphere_outer_shell_block_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_outer_shell_input(input_file, input_section);
  d4est_geometry_cubed_sphere_outer_shell_block_new_aux(d4est_geom, sphere_attrs);
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;
  
  if (mpirank == 0){
    printf("%s: NAME = cubed sphere outer shell block\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix , sphere_attrs->compactify_outer_shell);
  }
}

void
d4est_geometry_cubed_sphere_with_cube_hole_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_input(input_file, input_section);
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_with_hole();
  
  d4est_geom->p4est_conn = conn; 
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_X; /* same as cubed sphere, because
which_tree == 12 will never occur */
  d4est_geom->DX = d4est_geometry_cubed_sphere_DX; /* same as cubed sphere, because which_tree == 12 will never occur */
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy;
  
  if (mpirank == 0){
    printf("%s: NAME = cubed sphere with cube hole\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: R2 = %.25f\n", printf_prefix , sphere_attrs->R2);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix , sphere_attrs->compactify_outer_shell);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}


void
d4est_geometry_cubed_sphere_with_sphere_hole_new
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 const char* printf_prefix,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere_attrs = d4est_geometry_cubed_sphere_input(input_file, input_section);
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_with_hole();
  
  d4est_geom->p4est_conn = conn; 
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_with_sphere_hole_X;
  d4est_geom->DX = d4est_geometry_cubed_sphere_with_sphere_hole_DX; 
  d4est_geom->JAC = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy;
  d4est_geom->get_number_of_regions = d4est_geometry_cubed_sphere_get_number_of_regions;
  d4est_geom->get_region = d4est_geometry_cubed_sphere_get_region;
  
  if (mpirank == 0){
    printf("%s: NAME = cubed sphere with sphere hole\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: R2 = %.25f\n", printf_prefix , sphere_attrs->R2);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix , sphere_attrs->compactify_outer_shell);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}
