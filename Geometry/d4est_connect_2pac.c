#include <d4est_connect_2pac.h>

p4est_connectivity_t *
p4est_connectivity_new_2pac_aligned (void)
{
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  /* const double        vertices[6 * 3] = { */
  /*   0, 0, 0, */
  /*   1, 0, 0, */
  /*   0, 1, 0, */
  /*   1, 1, 0, */
  /*   2, 0, 0, */
  /*   2, 1, 0, */
  /* }; */
  const double        vertices[6 * 3] = {
    -1., -1., 0,
    0., -.5, 0,
    -1., 1, 0,
    0., .5, 0,
    1, -1, 0,
    1, 1, 0,
  }; 
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 1, 2, 3,
    1, 4, 3, 5,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    0, 1, 0, 0,
    0, 1, 1, 1,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 0, 2, 3,
    1, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}


p4est_connectivity_t *
p4est_connectivity_new_2pac_aligned_CUBE (void)
{
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[6 * 3] = {
    -1., -1., 0,
    0.6, -1, 0,
    -1., 1, 0,
    0.6, 1, 0,
    1, -1, 0,
    1, 1, 0,
  }; 
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 1, 2, 3,
    1, 4, 3, 5,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    0, 1, 0, 0,
    0, 1, 1, 1,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 0, 2, 3,
    1, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}


p4est_connectivity_t *
p4est_connectivity_new_2pac_nonaligned (void)
{
 const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  /* const double        vertices[6 * 3] = { */
  /*   0, 0, 0, */
  /*   1, 0, 0, */
  /*   0, 1, 0, */
  /*   1, 1, 0, */
  /*   2, 0, 0, */
  /*   2, 1, 0, */
  /* }; */
  const double        vertices[6 * 3] = {
    -1., -1., 0,
    0., -.5, 0,
    -1., 1, 0,
    0., .5, 0,
    1, -1, 0,
    1, 1, 0,
  }; 
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 1, 2, 3,
    3, 1, 5, 4,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    0, 1, 0, 0,
    1, 1, 0, 1,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 6, 2, 3,
    0, 1, 5, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_2pac_nonaligned_CUBE (void)
{
 const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  /* const double        vertices[6 * 3] = { */
  /*   0, 0, 0, */
  /*   1, 0, 0, */
  /*   0, 1, 0, */
  /*   1, 1, 0, */
  /*   2, 0, 0, */
  /*   2, 1, 0, */
  /* }; */
  const double        vertices[6 * 3] = {
    -1., -1., 0,
    0., -1, 0,
    -1., 1, 0,
    0., 1, 0,
    1, -1, 0,
    1, 1, 0,
  }; 
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 1, 2, 3,
    3, 1, 5, 4,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    0, 1, 0, 0,
    1, 1, 0, 1,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 6, 2, 3,
    0, 1, 5, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

