#include <pXest.h>
#include <d4est_connectivity_cubed_sphere.h>


p4est_connectivity_t *
d4est_connectivity_new_sphere_7tree (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 16;
  const p4est_topidx_t num_trees = 7;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double vertices[16 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
    -1, -1, -1,
     1, -1, -1,
    -1,  1, -1,
     1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
  };
  const p4est_topidx_t tree_to_vertex[7 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    8,  9, 10, 11, 12, 13, 14, 15,
  };
  const p4est_topidx_t tree_to_tree[7 * 6] = {    
     5,  3, 4,  1, 6,  0, //
     5,  3,  0,  2, 6,  1,
     5,  3,  1, 4, 6,  2,
     2,  0,  1, 4, 6,  3,
     2,  0,  3, 5, 6,  4,
     2,  0, 4,  1, 6,  5,
     5,  3,  0,  2, 4,  1,
  };
  const int8_t tree_to_face[7 * 6] = {
     1,  7,  7,  2,  2,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6, 15,  5,
     1,  7,  7,  2, 19,  5,
     9,  8,  3,  2, 22,  5,
     6,  0,  3,  6,  6,  5,
    10, 22,  4, 16, 22,  4,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
}


 p4est_connectivity_t *
d4est_connectivity_new_sphere_innerouter_shell (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double        vertices[8 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
  };
  const p4est_topidx_t tree_to_vertex[2 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
  };
  const p4est_topidx_t tree_to_tree[2 * 6] = {
     0,  0,  0,  0,  1,  0, //0 CHANGE THESE to 0 and 0, i.e 0 -> 0 and 1 -> 0
     1,  1,  1,  1, 1,  0, //1
  };
  const int8_t        tree_to_face[12 * 6] = {
     0,  1,  2,  3,  5,  5, //1
     0,  1,  2,  3,  4,  4, //7
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
}

 p4est_connectivity_t *
d4est_connectivity_new_sphere_with_hole (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 12;
  const p4est_topidx_t ctt_offset = 0;
  const p4est_topidx_t ett_offset = 0;
  const double        vertices[8 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
  };
  const p4est_topidx_t tree_to_vertex[12 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
  };
  const p4est_topidx_t tree_to_tree[12 * 6] = {
     5,  3,  4,  1,  6,  0,
     5,  3,  0,  2,  7,  1,
     5,  3,  1,  4,  8,  2,
     2,  0,  1,  4,  9,  3,
     2,  0,  3,  5, 10,  4,
     2,  0,  4,  1, 11,  5,
    11,  9, 10,  7, 6,  0,
    11,  9,  6,  8, 7,  1,
    11,  9,  7, 10, 8,  2,
     8,  6,  7, 10, 9,  3,
     8,  6,  9, 11, 10,  4,
     8,  6, 10,  7, 11,  5,
  };
  const int8_t        tree_to_face[12 * 6] = {
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  4,  4,
     9,  8,  3,  2,  4,  4,
     6,  0,  3,  6, 4,  4,
     1,  7,  7,  2, 4,  4,
     9,  8,  3,  2, 4,  4,
     6,  0,  3,  6,  4,  4,
  };

#if (P4EST_DIM)==3
  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &ett_offset,
                                      NULL, NULL,
                                      NULL, &ctt_offset, NULL, NULL);
#else
  return NULL;
#endif
}
