#include <d4est_ghost.h>
#include <d4est_element_data.h>


static void
d4est_ghost_get_mirror_elements
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost
)
{
  size_t              zz;
  size_t              data_size;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      which_quad;
  p4est_quadrant_t   *mirror, *q;
  p4est_tree_t       *tree;

  data_size = p4est->data_size == 0 ? sizeof (void *) : p4est->data_size;

  for (zz = 0; zz < d4est_ghost->ghost->mirrors.elem_count; ++zz) {
    mirror = p4est_quadrant_array_index (&d4est_ghost->ghost->mirrors, zz);
    which_tree = mirror->p.piggy3.which_tree;
    D4EST_ASSERT (p4est->first_local_tree <= which_tree &&
                  which_tree <= p4est->last_local_tree);
    /* D4EST_ASSERT (prev_tree <= which_tree); */

    tree = p4est_tree_array_index (p4est->trees, which_tree);
    which_quad = mirror->p.piggy3.local_num - tree->quadrants_offset;
    D4EST_ASSERT (0 <= which_quad &&
                  which_quad < (p4est_locidx_t) tree->quadrants.elem_count);
    q = p4est_quadrant_array_index (&tree->quadrants, which_quad);
    d4est_ghost->mirror_elements[zz] = *(d4est_element_data_t*)q->p.user_data;
      /* p4est->data_size == 0 ? &q->p.user_data : q->p.user_data; */
      /* p4est->data_size == 0 ? &q->p.user_data : ; */
  } 
}

d4est_ghost_t*
d4est_ghost_init
(
 p4est_t* p4est
)
{
  d4est_ghost_t* d4est_ghost = D4EST_ALLOC(d4est_ghost_t, 1);
  
  d4est_ghost->ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  d4est_ghost->ghost_elements = D4EST_ALLOC (d4est_element_data_t, d4est_ghost->ghost->ghosts.elem_count);
  d4est_ghost->mirror_elements = D4EST_ALLOC(d4est_element_data_t, d4est_ghost->ghost->mirrors.elem_count);
  
  p4est_ghost_exchange_data(p4est, d4est_ghost->ghost, d4est_ghost->ghost_elements);
  d4est_ghost_get_mirror_elements(p4est, d4est_ghost);

  int stride = 0;
  for (int gid = 0; gid < d4est_ghost->ghost->ghosts.elem_count; gid++){
    d4est_element_data_t* ged = &d4est_ghost->ghost_elements[gid];
    int volume_nodes_lobatto = d4est_lgl_get_nodes((P4EST_DIM),ged->deg);
    ged->id = gid;
    ged->nodal_stride = stride;
    ged->quad_stride = -1;
    ged->sqr_nodal_stride = -1;
    /* for (int i = 0; i < (P4EST_FACES); i++){ */
    /*   ged->boundary_quad_vector_stride[i] = -1; */
    /*   ged->mortar_quad_vector_stride[i] = -1; */
    /*   ged->mortar_quad_scalar_stride[i] = -1; */
    /*   ged->mortar_quad_matrix_stride[i] = -1; */
    /* }     */
    stride += volume_nodes_lobatto;    
  }
  

  return d4est_ghost;
}

void
d4est_ghost_destroy
(
 d4est_ghost_t* d4est_ghost
)
{
  if (d4est_ghost != NULL){
    /* free the ghost_elements */
    D4EST_FREE(d4est_ghost->ghost_elements);
    D4EST_FREE(d4est_ghost->mirror_elements);

    /* free the ghost object */
    p4est_ghost_destroy(d4est_ghost->ghost);

    D4EST_FREE(d4est_ghost);
  }
}

void
d4est_ghost_update
(
 p4est_t* p4est,
 d4est_ghost_t* d4est_ghost
)
{
  p4est_ghost_destroy(d4est_ghost->ghost);
  D4EST_FREE(d4est_ghost->ghost_elements);
  D4EST_FREE(d4est_ghost->mirror_elements);
  
  d4est_ghost->ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  d4est_ghost->ghost_elements = D4EST_ALLOC(d4est_element_data_t, d4est_ghost->ghost->ghosts.elem_count);  
  d4est_ghost->mirror_elements = D4EST_ALLOC(d4est_element_data_t, d4est_ghost->ghost->ghosts.elem_count);  
  
  p4est_ghost_exchange_data(p4est, d4est_ghost->ghost, d4est_ghost->ghost_elements);
  d4est_ghost_get_mirror_elements(p4est, d4est_ghost);
}

