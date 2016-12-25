#ifndef GRID_FUNCTIONS_H
#define GRID_FUNCTIONS_H 

#include "../pXest/pXest.h"

/** 
 * Function pointer for 
 * functions of the form
 * F(x,y,z) or F(x,y)
 * where x,y,z is a nodal
 * grid point
 * 
 * @param grid_fcn_t 
 * 
 * @return 
 */
typedef double
(*grid_fcn_t)
(
 double,
 double
#if (P4EST_DIM)==3
 ,  
 double
#endif // (P4EST_DIM)==3
);

double zero_fcn
(
 double x,
 double y
#if (P4EST_DIM)==3
 ,  
 double z
#endif // (P4EST_DIM)==3
);

typedef double
(*grid_fcn_ext_t)
(
 double, //x
 double, //y
#if (P4EST_DIM)==3
 double, //z
#endif
 double,  //gradu or u
 void*
);


double identity_fcn
(
 double x,
 double y,
#if (P4EST_DIM)==3
 double z,
#endif // (P4EST_DIM)==3
 double u,
 void* user
);

#endif
