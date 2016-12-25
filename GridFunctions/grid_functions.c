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
