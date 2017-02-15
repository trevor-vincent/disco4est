#ifndef D4EST_GEOMETRY_H
#define D4EST_GEOMETRY_H 

#include <p8est_connectivity.h>
#include <pXest.h>

typedef enum {INTERP_X_ON_LOBATTO, COMPUTE_DX_ON_LOBATTO, COMPUTE_DX_ON_GAUSS} dxdr_method_t;

typedef struct{

  p4est_geometry_t* p4est_geom;
  p4est_connectivity_t* p4est_conn;

  /* if the analytic derivatives are known */
  void(*dxdr)(p4est_geometry_t * geom,
              p4est_topidx_t which_tree,
              const double rst[(P4EST_DIM)], // \in [-1,1] such as GL or GLL points
              double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],
              p4est_qcoord_t q0 [(P4EST_DIM)],
              p4est_qcoord_t dq);



  /* if the analytic derivatives are known */
  void(*dxdr_face)(
                   p4est_geometry_t * geom,
                   p4est_topidx_t which_tree,
                   const double rs[(P4EST_DIM)-1], // \in [-1,1] such as GL or GLL points
                   double dxyz_drs[(P4EST_DIM)][(P4EST_DIM)-1],
                   p4est_qcoord_t q0 [(P4EST_DIM)],
                   p4est_qcoord_t dqa [(P4EST_DIM)-1][(P4EST_DIM)]
  );

  
  /* int interp_to_Gauss; */
  dxdr_method_t dxdr_method;
  
} d4est_geometry_t;

void d4est_geometry_octree_to_vertex(p8est_connectivity_t *connectivity,p4est_topidx_t which_tree,const double abc[3],double xyz[3]);
void d4est_geometry_destroy(d4est_geometry_t *d4est_geom);
d4est_geometry_t *d4est_geometry_new(const char *input_file);

#endif
