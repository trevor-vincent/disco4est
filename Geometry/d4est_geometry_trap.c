#include <d4est_geometry_trap.h>
#include <d4est_util.h>

static p4est_connectivity_t *
p4est_connectivity_new_trap (void)
{
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[4 * 3] = {
    0, 0, 0,
    1, -.5, 0,
    0, 1, 0,
    1, 1.5, 0,
  };
  const p4est_topidx_t tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[1 * 4] = {
    0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 4] = {
    0, 1, 2, 3,
  };

#if (P4EST_DIM)==2
  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
#else
  D4EST_ABORT("[ERROR]: trapezoid geometry only supports DIM=2");
  return NULL;
#endif
  
}

void
d4est_geometry_trap_new
(
 const char* input_file,
 d4est_geometry_t* d4est_geom
)
{
  p4est_connectivity_t* conn = p4est_connectivity_new_trap();
  p4est_geometry_t* geom = p4est_geometry_new_connectivity(conn);

  d4est_geom->p4est_conn = conn;
  d4est_geom->p4est_geom = geom;

  printf("[GEOMETRY_INFO]: NAME = trapezoid\n");
}
