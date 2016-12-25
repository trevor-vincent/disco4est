#ifndef D4EST_GEOMETRY_H
#define D4EST_GEOMETRY_H 

#include <pXest.h>

typedef enum {INTERP_X_ON_LOBATTO, COMPUTE_DX_ON_LOBATTO, COMPUTE_DX_ON_GAUSS} dxdr_method_t;

typedef struct{

  p4est_geometry_t* p4est_geom;
  void(*dxdr)(p4est_geometry_t * geom,
              p4est_topidx_t which_tree,
              const double rst[(P4EST_DIM)], // \in [-1,1] such as GL or GLL points
              double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],
              p4est_qcoord_t q0 [(P4EST_DIM)],
              p4est_qcoord_t dq);




  void(*dxdr_face)(
                   p4est_geometry_t * geom,
                   p4est_topidx_t which_tree,
                   const double rs[(P4EST_DIM)-1], // \in [-1,1] such as GL or GLL points
                   double dxyz_drs[(P4EST_DIM)][(P4EST_DIM)-1],
                   p4est_qcoord_t q0 [(P4EST_DIM)],
                   p4est_qcoord_t dqa [(P4EST_DIM)-1][(P4EST_DIM)]
  );

  
  int interp_to_Gauss;
  dxdr_method_t dxdr_method;
  
} d4est_geometry_t;

#endif
