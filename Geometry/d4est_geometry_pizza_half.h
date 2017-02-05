#ifndef D4EST_GEOMETRY_PIZZA_HALF_H
#define D4EST_GEOMETRY_PIZZA_HALF_H 

#include <pXest.h>
/* #include <p4est.h> */
/* #include <p4est_connectivity.h> */
/* #include <p4est_geometry.h> */

p4est_geometry_t *d4est_geometry_new_pizza_half(p4est_connectivity_t *conn, double R0,
                                          double R1);


void 
d4est_geometry_pizza_half_dxdr(p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double rst[3], // \in [-1,1] such as GL or GLL points
                               const p4est_qcoord_t dq,
                               double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);

void 
d4est_geometry_pizza_half_dxda(p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double rst[3], // \in [-1,1] such as GL or GLL points
                               double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);

void 
d4est_geometry_pizza_half_X(p4est_geometry_t * geom,
                            p4est_topidx_t which_tree,
                            const double rst[3], double xyz[3]);


void 
d4est_geometry_pizza_half_dxdr2(p4est_geometry_t * geom,
                                p4est_topidx_t which_tree,
                                const double rst[(P4EST_DIM)], // \in [-1,1] such as GL or GLL points
                                double dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],
                                p4est_qcoord_t q0 [(P4EST_DIM)],
                                p4est_qcoord_t dq);


void 
d4est_geometry_pizza_half_dxdr_face(
                                    p4est_geometry_t * geom,
                                    p4est_topidx_t which_tree,
                                    const double rs[(P4EST_DIM)-1], // \in [-1,1] such as GL or GLL points
                                    double dxyz_drs[(P4EST_DIM)][(P4EST_DIM)-1],
                                    p4est_qcoord_t q0 [(P4EST_DIM)],
                                    p4est_qcoord_t dqa [(P4EST_DIM)-1][(P4EST_DIM)]
);


#endif
