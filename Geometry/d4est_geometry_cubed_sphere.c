#include <pXest.h>
#include <d4est_geometry.h>
#include <ini.h>
#include <p8est_connectivity.h>
#include <util.h>
#include <petscsnes.h>

/* #define P4EST_DIM 3 */
#if (P4EST_DIM)==3
#include <d4est_geometry_cubed_sphere.h>

static inline double
secant_fcn(double x){
  return 1./cos(x);
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
  if (util_match_couple(section,pconfig->input_section,name,"R0")) {
    mpi_assert(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
    mpi_assert(pconfig->R0 > 0);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"R1")) {
    mpi_assert(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
    mpi_assert(pconfig->R1 > 0);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"R2")) {
    mpi_assert(pconfig->R2 == -1);
    pconfig->R2 = atof(value);
    mpi_assert(pconfig->R2 > 0);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"compactify_outer_shell")) {
    mpi_assert(pconfig->compactify_outer_shell == -1);
    pconfig->compactify_outer_shell = atoi(value);
    mpi_assert(pconfig->compactify_outer_shell == 0 || pconfig->compactify_outer_shell == 1);
  }
  else if (util_match_couple(section,pconfig->input_section,name,"compactify_inner_shell")) {
    mpi_assert(pconfig->compactify_inner_shell == -1);
    pconfig->compactify_inner_shell = atoi(value);
    mpi_assert(pconfig->compactify_inner_shell == 0 || pconfig->compactify_inner_shell == 1);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


static p4est_connectivity_t *
d4est_connectivity_new_sphere_with_cube_hole (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 12;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double        vertices[8 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
  };
  const p4est_topidx_t tree_to_vertex[12 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
  };
  const p4est_topidx_t tree_to_tree[12 * 6] = {
     5,  3,  4,  1,  6,  0,
     5,  3,  0,  2,  7,  1,
     5,  3,  1,  4,  8,  2,
     2,  0,  1,  4,  9,  3,
     2,  0,  3,  5, 10,  4,
     2,  0,  4,  1, 11,  5,
    11,  9, 10,  7, 6,  0,
    11,  9,  6,  8, 7,  1,
    11,  9,  7, 10, 8,  2,
     8,  6,  7, 10, 9,  3,
     8,  6,  9, 11, 10,  4,
     8,  6, 10,  7, 11,  5,
  };
  const int8_t        tree_to_face[12 * 6] = {
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  4,  4,
     9,  8,  3,  2,  4,  4,
     6,  0,  3,  6, 4,  4,
     1,  7,  7,  2, 4,  4,
     9,  8,  3,  2, 4,  4,
     6,  0,  3,  6,  4,  4,
  };

#if (P4EST_DIM)==3
  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
#else
  return NULL;
#endif
}


p4est_connectivity_t *
d4est_connectivity_new_sphere_7tree (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 16;
  const p4est_topidx_t num_trees = 7;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double vertices[16 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
    -1, -1, -1,
     1, -1, -1,
    -1,  1, -1,
     1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
  };
  const p4est_topidx_t tree_to_vertex[7 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    8,  9, 10, 11, 12, 13, 14, 15,
  };
  const p4est_topidx_t tree_to_tree[7 * 6] = {    
     5,  3, 4,  1, 6,  0, //
     5,  3,  0,  2, 6,  1,
     5,  3,  1, 4, 6,  2,
     2,  0,  1, 4, 6,  3,
     2,  0,  3, 5, 6,  4,
     2,  0, 4,  1, 6,  5,
     5,  3,  0,  2, 4,  1,
  };
  const int8_t tree_to_face[7 * 6] = {
     1,  7,  7,  2,  2,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6, 15,  5,
     1,  7,  7,  2, 19,  5,
     9,  8,  3,  2, 22,  5,
     6,  0,  3,  6,  6,  5,
    10, 22,  4, 16, 22,  5,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
}


static p4est_connectivity_t *
d4est_connectivity_new_sphere_innerouter_shell (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double        vertices[8 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
  };
  const p4est_topidx_t tree_to_vertex[2 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
  };
  const p4est_topidx_t tree_to_tree[2 * 6] = {
     0,  0,  0,  0,  1,  0, //0 CHANGE THESE to 0 and 0, i.e 0 -> 0 and 1 -> 0
     1,  1,  1,  1, 1,  0, //1
  };
  const int8_t        tree_to_face[12 * 6] = {
     0,  1,  2,  3,  5,  5, //1
     0,  1,  2,  3,  4,  4, //7
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
}



static void 
d4est_geometry_cubed_sphere_inner_shell_block_X(
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


static void 
d4est_geometry_cubed_sphere_outer_shell_block_X(
                                                d4est_geometry_t * geom,
                                                p4est_topidx_t which_tree,
                                                p4est_qcoord_t q0 [(P4EST_DIM)],
                                                p4est_qcoord_t dq,
                                                const double coords[(P4EST_DIM)],
                                                coords_type_t coords_type,
                                                double xyz[(P4EST_DIM)]
)
{
  d4est_geometry_cubed_sphere_attr_t* sphere = geom->user;
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
  if (sphere->compactify_outer_shell){
    double m = (2. - 1.)/((1./sphere->R2) - (1./sphere->R1));
    double t = (1.*sphere->R1 - 2.*sphere->R2)/(sphere->R1 - sphere->R2);
    R = m/(abc[2] - t);
  }
  else {
    R = sphere->R1*(2. - abc[2]) + sphere->R2*(abc[2] - 1.);
  }
  double q = R / sqrt (x * x + y * y + 1.);  
  xyz[0] = +q * x;
  xyz[1] = +q * y;
  xyz[2] = +q;
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
  d4est_geometry_cubed_sphere_X(geom,
                                which_tree + 6,
                                q0,
                                dq,
                                coords,
                                coords_type,
                                xyz
                               );
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
    mpi_abort("[D4EST_ERROR]: which_tree should only be 0 or 1 for innerouter_shell");
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
    mpi_abort("[D4EST_ERROR]: which_tree should only be 0 or 1 for innerouter_shell");
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
  
  switch (which_tree % 6) {
  case 0:                      /* front */
    /* xyz[0] = +q * x; */
    /* xyz[1] = -q; */
    /* xyz[2] = +q * y; */
    for (int d = 0; d < (P4EST_DIM); d++){
      dxyz_drst[0][d] = dxyz_drst_temp[0][d];
      dxyz_drst[1][d] = -dxyz_drst_temp[2][d];
      dxyz_drst[2][d] = dxyz_drst_temp[1][d];
    }
    break;
 case 1:                      /* top */
   /* xyz[0] = +q * x; */
   /* xyz[1] = +q * y; */
   /* xyz[2] = +q; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_temp[0][d];
     dxyz_drst[1][d] = dxyz_drst_temp[1][d];
     dxyz_drst[2][d] = dxyz_drst_temp[2][d];
    
   }
   break;
 case 2:                      /* back */
   /* xyz[0] = +q * x; */
   /* xyz[1] = +q; */
   /* xyz[2] = -q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_temp[0][d];
     dxyz_drst[1][d] = dxyz_drst_temp[2][d];
     dxyz_drst[2][d] = -dxyz_drst_temp[1][d];
   }
   break;
 case 3:                      /* right */
   /* xyz[0] = +q; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = -q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = dxyz_drst_temp[2][d];
     dxyz_drst[1][d] = -dxyz_drst_temp[0][d];
     dxyz_drst[2][d] = -dxyz_drst_temp[1][d];
   }
   break;
 case 4:                      /* bottom */
   /* xyz[0] = -q * y; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = -q; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = -dxyz_drst_temp[1][d];
     dxyz_drst[1][d] = -dxyz_drst_temp[0][d];
     dxyz_drst[2][d] = -dxyz_drst_temp[2][d];      
   }
   break;
 case 5:                      /* left */
   /* xyz[0] = -q; */
   /* xyz[1] = -q * x; */
   /* xyz[2] = +q * y; */
   for (int d = 0; d < (P4EST_DIM); d++){
     dxyz_drst[0][d] = -dxyz_drst_temp[2][d];
     dxyz_drst[1][d] = -dxyz_drst_temp[0][d];
     dxyz_drst[2][d] = dxyz_drst_temp[1][d];
   }
   break;
  default:
    SC_ABORT_NOT_REACHED();
  }
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
    mpi_abort("Can't load input file");
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
    mpi_abort("Can't load input file");
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
    mpi_abort("Can't load input file");
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
  sphere_attrs->R2 = -1;
  sphere_attrs->compactify_inner_shell = -1;
  sphere_attrs->compactify_outer_shell = -1;

  if (ini_parse(input_file, d4est_geometry_cubed_sphere_input_handler, sphere_attrs) < 0) {
    mpi_abort("Can't load input file");
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
    mpi_abort("Can't load input file");
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
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_DX;
  else
    d4est_geom->DX = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;

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
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_7tree_DX;
  else
    d4est_geom->DX = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;

  if (mpirank == 0){
    printf("%s: NAME = cubed sphere\n", printf_prefix );
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
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_innerouter_shell_DX;
  else
    d4est_geom->DX = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;

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
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_inner_shell_block_DX;
  else
    d4est_geom->DX = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;

  if (mpirank == 0){
    printf("%s: NAME = cubed sphere inner shell block\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
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
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();
  
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_outer_shell_block_X;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy; 
  d4est_geom->p4est_conn = conn;
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_outer_shell_block_DX;
  else
    d4est_geom->DX = NULL;

  
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
  p4est_connectivity_t* conn = d4est_connectivity_new_sphere_with_cube_hole();
  
  d4est_geom->p4est_conn = conn; 
  d4est_geom->user = sphere_attrs;
  d4est_geom->X = d4est_geometry_cubed_sphere_X; /* same as cubed sphere, because
which_tree == 12 will never occur */
  if (d4est_geom->X_mapping_type == MAP_ANALYTIC)
    d4est_geom->DX = d4est_geometry_cubed_sphere_DX; /* same as cubed sphere, because which_tree == 12 will never occur */
  else
    d4est_geom->DX = NULL;
  d4est_geom->destroy = d4est_geometry_cubed_sphere_destroy;
  
  if (mpirank == 0){
    printf("%s: NAME = compact_sphere_with_cube_hole\n", printf_prefix );
    printf("%s: R0 = %.25f\n", printf_prefix , sphere_attrs->R0);
    printf("%s: R1 = %.25f\n", printf_prefix , sphere_attrs->R1);
    printf("%s: R2 = %.25f\n", printf_prefix , sphere_attrs->R2);
    printf("%s: compactify_outer_shell = %d\n", printf_prefix , sphere_attrs->compactify_outer_shell);
    printf("%s: compactify_inner_shell = %d\n", printf_prefix , sphere_attrs->compactify_inner_shell);
  }
}

#endif
