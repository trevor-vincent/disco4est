#ifndef D4EST_GEOMETRY_H
#define D4EST_GEOMETRY_H 

#include <p8est_connectivity.h>
#include <pXest.h>
#include <dgmath.h>


/** This object encapsulates a custom geometry transformation. */
typedef struct d4est_geometry d4est_geometry_t;

typedef enum {MAP_ISOPARAMETRIC, MAP_ANALYTIC, MAP_NONE} mapping_type_t;

/**
* There are three types of coords used to describe a point in an element.
* The first are the integration coords for gaussian quadrature. These are
* specified on the reference cube [-1,1]^DIM and represent a local coordinate
* system for the element. These are called COORDS_INTEG_RST. A point in an
* element can also have coordinates in the tree it is in, here d4est commonly
* used two different coords. The first is COORDS_P4EST_INT which is the built-in
* coords used by p4est, which represent points in the tree in the cube [0, P4EST_ROOT_LEN]^DIM
* where P4EST_ROOT_LEN is some big integer. These coords are stored as integers.
* The second is COORDS_TREE_UNITCUBE which are stored as dbls in the range [0,1]^DIM
* 
*/

typedef enum {COORDS_INTEG_RST, COORDS_P4EST_INT, COORDS_TREE_UNITCUBE} coords_type_t;

/** Forward transformation from the reference unit square to physical space.
 * Note that the two-dimensional connectivities have 3D vertex coordinates
 * that can be used in the transformation if so desired.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 */
typedef void        (*d4est_geometry_X_t) (d4est_geometry_t*,
                                           p4est_topidx_t ,
                                           p4est_qcoord_t [(P4EST_DIM)],
                                           p4est_qcoord_t,
                                           const double [(P4EST_DIM)],
                                           coords_type_t,
                                           double [(P4EST_DIM)]
                                          );

typedef void        (*d4est_geometry_DX_t) (d4est_geometry_t*,
                                            p4est_topidx_t ,
                                            p4est_qcoord_t [(P4EST_DIM)],
                                            p4est_qcoord_t,
                                            const double [(P4EST_DIM)],
                                            double [(P4EST_DIM)][(P4EST_DIM)]
                                            );

/** Destructor prototype for a user-allocated \a p4est_geometry_t.
 * It is invoked by p4est_geometry_destroy.  If the user chooses to
 * reserve the structure statically, simply don't call p4est_geometry_destroy.
 */
typedef void        (*d4est_geometry_destroy_t) (d4est_geometry_t * geom);


struct d4est_geometry {

  p4est_connectivity_t* p4est_conn;

  /* Mapping from [-1,1]^3 to grid coordinates*/
  d4est_geometry_X_t X;
  mapping_type_t X_mapping_type;
  
  /* Analytic derivatives of the mapping */
  d4est_geometry_DX_t DX;
  
  d4est_geometry_destroy_t destroy;
  void* user;
  
};

/* This file was automatically generated.  Do not edit! */
void
d4est_geometry_get_tree_coords_in_range_0_to_1
(
 p4est_qcoord_t q0 [3],
 p4est_qcoord_t dq,
 const double coords[3],
 coords_type_t coords_type,
 double tcoords[3]
);
void d4est_geometry_compute_geometric_data_on_mortar_TESTINGONLY(p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int *deg_mortar,int face_side,quadrature_type_t quad_type,double *n[(P4EST_DIM)],double *sj,d4est_geometry_t *d4est_geom,dgmath_jit_dbase_t *dgmath_jit_dbase,double *xyz_storage[(P4EST_DIM)]);
void d4est_geometry_compute_drst_dxyz(double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double *drst_dxyz[(P4EST_DIM)][(P4EST_DIM)],int nodes);
void d4est_geometry_compute_xyz(dgmath_jit_dbase_t *dgmath_jit_dbase,d4est_geometry_t *d4est_geom,int which_tree,int deg,quadrature_type_t type,p4est_qcoord_t q[(P4EST_DIM)],p4est_qcoord_t dq,double *xyz[(P4EST_DIM)]);
void d4est_geometry_compute_geometric_data_on_mortar(p4est_topidx_t e0_tree,p4est_qcoord_t e0_q[(P4EST_DIM)],p4est_qcoord_t e0_dq,int num_faces_side,int num_faces_mortar,int *deg_mortar_integ,int face,double *drst_dxyz_on_mortar_integ[(P4EST_DIM)][(P4EST_DIM)],double *sj_on_mortar_integ,double *n_on_mortar_integ[(P4EST_DIM)],double *j_div_sj_mortar_integ,quadrature_type_t quad_type,d4est_geometry_t *d4est_geom,dgmath_jit_dbase_t *dgmath_jit_dbase);
void d4est_geometry_compute_jacobian_and_drst_dxyz(double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)],double *jac,double *drst_dxyz[(P4EST_DIM)][(P4EST_DIM)],int volume_nodes);
void d4est_geometry_data_compute_dxyz_drst_face_isoparametric(dgmath_jit_dbase_t *dgmath_jit_dbase,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int which_tree,int face,d4est_geometry_t *d4est_geom,quadrature_type_t quad_type,int deg,double *dxyz_drst_on_face[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_data_compute_dxyz_drst_face_analytic(dgmath_jit_dbase_t *dgmath_jit_dbase,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int which_tree,int face,d4est_geometry_t *d4est_geom,quadrature_type_t quad_type,int deg,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst_isoparametric(p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,quadrature_type_t quad_type,d4est_geometry_t *d4est_geom,dgmath_jit_dbase_t *dgmath_jit_dbase,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst(p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,quadrature_type_t quad_type,d4est_geometry_t *d4est_geom,dgmath_jit_dbase_t *dgmath_jit_dbase,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_compute_dxyz_drst_analytic(p4est_topidx_t which_tree,p4est_qcoord_t q0[(P4EST_DIM)],p4est_qcoord_t dq,int deg,quadrature_type_t quad_type,d4est_geometry_t *d4est_geom,dgmath_jit_dbase_t *dgmath_jit_dbase,double *dxyz_drst[(P4EST_DIM)][(P4EST_DIM)]);
void d4est_geometry_quadtree_to_vertex(p4est_connectivity_t *connectivity,p4est_topidx_t which_tree,const double abc[3],double xyz[3]);
void d4est_geometry_octree_to_vertex(p8est_connectivity_t *connectivity,p4est_topidx_t which_tree,const double abc[3],double xyz[3]);
void d4est_geometry_destroy(d4est_geometry_t *d4est_geom);
d4est_geometry_t *d4est_geometry_new(int mpirank,const char *input_file);

#endif

