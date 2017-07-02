#include "../pXest/pXest.h"
#include "../GridFunctions/grid_functions.h"

double
zero_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,  
 double z
#endif
)
{
  return 0.;
}

double
sinpix_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,  
 double z
#endif
)
{
  double pi = 3.1415926535897932384626433832795;
  double sinpix = sin(pi*x);
  sinpix *= sin(pi*y);
#if (P4EST_DIM)==3
  sinpix *= sin(pi*z);
#endif
  return sinpix;
}

double identity_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif // (P4EST_DIM)==3
 double u,
 void* user
){
  return u;
}


