#ifndef D4EST_XYZ_FUNCTIONS_EXT_H
#define D4EST_XYZ_FUNCTIONS_EXT_H 

#include <d4est_geometry.h>
#include <d4est_element_data.h>

typedef double
(*d4est_xyz_fcn_ext_t)
(
 double,
 double,
#if (P4EST_DIM)==3
 double,
#endif
 void* user,
 d4est_geometry_t*,
 d4est_element_data_t*
);

#endif
