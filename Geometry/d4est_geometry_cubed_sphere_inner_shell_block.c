#include <d4est_geometry.h>
#include <d4est_geometry_cubed_sphere_inner_shell_block.h>
#include <util.h>
#include <ini.h>

typedef struct {

  double R0;
  double R1;
  int compactify
  
} d4est_geometry_cubed_sphere_inner_shell_block_input_t;

static
int d4est_geometry_cubed_sphere_inner_shell_block_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  d4est_geometry_cubed_sphere_inner_shell_block_input_t* pconfig = (d4est_geometry_cubed_sphere_inner_shell_block_input_t*)user;
  if (util_match_couple(section,"geometry",name,"R0")) {
    mpi_assert(pconfig->R0 == -1);
    pconfig->R0 = atof(value);
  }
  else if (util_match_couple(section,"geometry",name,"R1")) {
    mpi_assert(pconfig->R1 == -1);
    pconfig->R1 = atof(value);
  }
  else if (util_match_couple(section,"geometry",name,"compactify")) {
    mpi_assert(pconfig->compactify == -1);
    pconfig->compactify = atoi(value);
    mpi_assert(pconfig->compactify == 0 || pconfig->compactify == 1);
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

d4est_geometry_cubed_sphere_inner_shell_block_input_t*
d4est_geometry_cubed_sphere_inner_shell_block_input
(
 const char* input_file
)
{
  
  d4est_geometry_cubed_sphere_inner_shell_block_input_t* input
    = P4EST_ALLOC(d4est_geometry_cubed_sphere_inner_shell_block_input_t, 1);
  
  input->R0 = -1;
  input->R1 = -1;
  
  if (ini_parse(input_file, d4est_geometry_cubed_sphere_inner_shell_block_input_handler, &input) < 0) {
    mpi_abort("Can't load input file");
  }

  D4EST_CHECK_INPUT("geometry", input->R0, -1);
  D4EST_CHECK_INPUT("geometry", input->R1, -1);
  
  return input;
}


void 
d4est_geometry_cubed_sphere_inner_shell_block_X(
                                                p4est_geometry_t * geom,
                                                p4est_topidx_t which_tree,
                                                const double rst[3], /* [0,1]^3  */
                                                double xyz[3])
{
  double* radii = (double*)geom->user;
  
  d4est_geometry_sphere_attr_t* sphere = geom->user;

  /* transform rst from [0,1]^3 to [0,1]^2 x [1,2]  */
  double abc [3];
  abc[0] = rst[0];
  abc[1] = rst[1];
  abc[2] = rst[2] + 1;

  if (sphere->compactify){
    double m = (2. - 1.)/((1./sphere->R1) - (1./sphere->R0));
    double t = (1.*sphere->R0 - 2.*sphere->R1)/(sphere->R0 - sphere->R1);
    R = m/(abc[2] - t);
  }
  else {
    R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);
  }
  x = p * abc[0] + (1. - p) * tanx;
  y = p * abc[1] + (1. - p) * tany;   
  q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  xyz[0] = +q * x;
  xyz[1] = +q * y;
  xyz[2] = +q;
}

void
d4est_geometry_cubed_sphere_inner_shell_block_DX(d4est_geometry_t* d4est_geom,
                                                 p4est_topidx_t which_tree,
                                                 p4est_qcoord_t q0 [3],
                                                 p4est_qcoord_t dq,
                                                 const double rst[3], /* [-1,1]^3 */
                                                 const double dxyz_drst[3][3]
                                                )
{
  int compactify = (( d4est_geometry_cubed_sphere_inner_shell_block_input_t*)(d4est_geom->p4est_geom->user))->compactify;
  double R0 = (( d4est_geometry_cubed_sphere_inner_shell_block_input_t*)(d4est_geom->p4est_geom->user))->R0;
  double R1 = (( d4est_geometry_cubed_sphere_inner_shell_block_input_t*)(d4est_geom->p4est_geom->user))->R1;
  double r = rst[0];
  double s = rst[1];
  double t = rst[2];
  double amin = q[0];
  double amax = q[0] + dq;
  double bmin = q[1];
  double bmax = q[1] + dq;
  double cmin = q[2];
  double cmax = q[2] + dq;

  /* FROM MATHEMATICA */
  
  if (compactify) {
    dxyz_drst[0][0] = ((amax - amin)*R0*R1*(-(M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*
          tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
            2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))) - 
       4*(-4*(-2 + cmax + cmin + cmax*t - cmin*t) + M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2))*
        (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
             (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (32.*(-2*R1 + cmin*(R0 - R1)*(-1 + t) - cmax*(R0 - R1)*(1 + t))*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
                                               (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));


    dxyz_drst[0][1] = -((bmax - bmin)*M_PI*R0*R1*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*
      (-((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
        ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/
   (8.*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5))

    dxyz_drst[0][2] = -((cmax - cmin)*R0*R1*(-((2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
           ((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
             2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*
           (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))) + 
        2*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
         (amax + amin + amax*r - amin*r - 2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*
         (-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
           (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 
        16*(R0 - R1)*(-((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
           ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*
         (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
              (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (8.*pow(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t),2)*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
                                               (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));

    dxyz_drst[1][0] = -((amax - amin)*M_PI*R0*R1*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*
      tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*(-((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
        ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.))/
   (8.*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
                                               (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    
    dxyz_drst[1][1] = ((bmax - bmin)*R0*R1*(-(M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*
          tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)*((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
            2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))) - 
       4*(-4*(-2 + cmax + cmin + cmax*t - cmin*t) + M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))*
        (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
             (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (32.*(-2*R1 + cmin*(R0 - R1)*(-1 + t) - cmax*(R0 - R1)*(1 + t))*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
                                               (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    
    dxyz_drst[1][2] = -((cmax - cmin)*R0*R1*(-((2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
           ((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
             2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*
           (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))) + 
        2*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
         (bmax + bmin + bmax*s - bmin*s - 2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*
         (-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
           (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 
        16*(R0 - R1)*(-((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
           ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.)*
         (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
              (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (8.*pow(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t),2)*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
                                               (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));

    dxyz_drst[2][0] = -((amax - amin)*M_PI*R0*R1*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*
      tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/
   (2.*Sqrt(2)*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
     pow(-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
       (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2),1.5));

    dxyz_drst[2][1] = -((bmax - bmin)*M_PI*R0*R1*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*
      tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/
   (2.*Sqrt(2)*(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
     pow(-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
       (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2),1.5));
    
    dxyz_drst[2][2] = -((cmax - cmin)*R0*R1*((2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t))*
         (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 
        4*(R0 - R1)*(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
              (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (2.*pow(2*R1 - cmin*(R0 - R1)*(-1 + t) + cmax*(R0 - R1)*(1 + t),2)*
     pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));  
  }
  else {
    dxyz_drst[0][0] = -((amax - amin)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      (-(M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*
           ((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
             2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))) - 
        4*(-4*(-2 + cmax + cmin + cmax*t - cmin*t) + M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2))*
         (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
              (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (128.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][1] = -((bmax - bmin)*M_PI*(cmax + cmin + cmax*t - cmin*t)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*
      (-((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
        ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[0][2] = ((cmax - cmin)*((-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
        ((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
          2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*
        (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 
       2*(R0*(-2 + cmax + cmin + cmax*t - cmin*t) - R1*(cmax + cmin + cmax*t - cmin*t))*
        (amax + amin + amax*r - amin*r - 2*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))*
        (-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
          (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 
       16*(R0 - R1)*(-((amax + amin + amax*r - amin*r)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
          ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/2.)*
        (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
             (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));

    dxyz_drst[1][0] = -((amax - amin)*M_PI*(cmax + cmin + cmax*t - cmin*t)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.)*
      (-((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
        ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][1] = -((bmax - bmin)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      (-(M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.)*
           ((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
             2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))) - 
        4*(-4*(-2 + cmax + cmin + cmax*t - cmin*t) + M_PI*(cmax + cmin + cmax*t - cmin*t)*pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2))*
         (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
              (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (128.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[1][2] = ((cmax - cmin)*((-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
        ((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t) - 
          2*(cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*
        (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) + 
       2*(R0*(-2 + cmax + cmin + cmax*t - cmin*t) - R1*(cmax + cmin + cmax*t - cmin*t))*
        (bmax + bmin + bmax*s - bmin*s - 2*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))*
        (-2*(-3 + cmax + cmin + cmax*t - cmin*t) + (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + 
          (cmax + cmin + cmax*t - cmin*t)*pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 
       16*(R0 - R1)*(-((bmax + bmin + bmax*s - bmin*s)*(-2 + cmax + cmin + cmax*t - cmin*t))/4. + 
          ((cmax + cmin + cmax*t - cmin*t)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/2.)*
        (3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
             (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));

    dxyz_drst[2][0] = -((amax - amin)*M_PI*(cmax + cmin + cmax*t - cmin*t)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      pow(sec((M_PI*(amax + amin + amax*r - amin*r))/8.),2)*tan((M_PI*(amax + amin + amax*r - amin*r))/8.))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][1] = -((bmax - bmin)*M_PI*(cmax + cmin + cmax*t - cmin*t)*(-(R0*(-2 + cmax + cmin + cmax*t - cmin*t)) + R1*(cmax + cmin + cmax*t - cmin*t))*
      pow(sec((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)*tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.))/
   (32.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    dxyz_drst[2][2] = ((cmax - cmin)*((R0*(-2 + cmax + cmin + cmax*t - cmin*t) - R1*(cmax + cmin + cmax*t - cmin*t))*
        (-2 + pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)) - 
       4*(R0 - R1)*(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
             (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.)))/
   (8.*pow(3 + cmin*(-1 + t) - cmax*(1 + t) + ((cmax + cmin + cmax*t - cmin*t)*
          (pow(tan((M_PI*(amax + amin + amax*r - amin*r))/8.),2) + pow(tan((M_PI*(bmax + bmin + bmax*s - bmin*s))/8.),2)))/2.,1.5));
    
  }
}


static void
d4est_geometry_cubed_sphere_inner_shell_block_destroy
(
 p4est_geometry_t* geom
)
{
  P4EST_FREE(geom->user);
  P4EST_FREE(geom);
}


void
d4est_geometry_cubed_sphere_inner_shell_block_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  d4est_geometry_cubed_sphere_inner_shell_block_input_t input = d4est_geometry_cubed_sphere_inner_shell_block_input(input_file);
  p4est_geometry_t *cubed_sphere_inner_shell_block_geometry = P4EST_ALLOC(p4est_geometry_t,1);
  p4est_connectivity_t* conn = p8est_connectivity_new_unitcube();

  double* radii = P4EST_ALLOC(double, 2);
  radii[0] = input.R0;
  radii[1] = input.R1;
  
  cubed_sphere_inner_shell_block_geometry->user = (void*)radii;
  cubed_sphere_inner_shell_block_geometry->destroy = d4est_geometry_cubed_sphere_inner_shell_block_destroy;
  cubed_sphere_inner_shell_block_geometry->X = d4est_geometry_cubed_sphere_inner_shell_block_X;

  d4est_geom->p4est_geom = cubed_sphere_inner_shell_block_geometry;
  d4est_geom->p4est_conn = conn;
  
  printf("[GEOMETRY_INFO]: NAME = cubed_sphere_inner_shell_block\n");
  printf("[GEOMETRY_INFO]: R0 = %.25f\n", radii[0]);
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", radii[1]);
}
