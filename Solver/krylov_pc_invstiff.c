#include <krylov_pc_invstiff.h>
#include <linalg.h>

krylov_pc_t*
krylov_pc_invstiff_create(p4est_t* p4est, dgmath_jit_dbase_t* dgmath_jit_dbase){

  krylov_pc_t* pc = P4EST_ALLOC(krylov_pc_t, 1);
  krylov_pc_invstiff_data_t* pc_data = P4EST_ALLOC(krylov_pc_invstiff_data_t, 1);
  
  pc_data->local_matrix_nodes = curved_element_data_get_local_matrix_nodes(p4est);
  pc_data->inv_stiff = P4EST_ALLOC(double, pc_data->local_matrix_nodes);
  pc_data->p4est = p4est;
  pc_data->dgmath_jit_dbase = dgmath_jit_dbase;
  pc->pc_data = (void*)pc_data;
  pc->pc_apply = krylov_pc_invstiff_apply;
  pc->pc_setup = krylov_pc_invstiff_setup;

  return pc;
}

void
krylov_pc_invstiff_setup
(
 krylov_pc_t* kpc
)
{
  krylov_pc_invstiff_data_t* invstiff_data = (krylov_pc_invstiff_data_t*)kpc->pc_data;
  int nodal_matrix_stride = 0;
  p4est_t* p4est = invstiff_data->p4est;
  dgmath_jit_dbase_t* dgmath_jit_dbase = invstiff_data->dgmath_jit_dbase;
  d4est_geometry_t* d4est_geom = kpc->pc_ctx->d4est_geom;
  
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);

        dgmath_compute_curved_inverse_stiffness_matrix
          (
           dgmath_jit_dbase,
           ed->deg,
           ed->J_integ,
           ed->rst_xyz_integ,
           ed->deg_integ,
           d4est_geom->geom_quad_type,
           (P4EST_DIM),
           &invstiff_data->inv_stiff[nodal_matrix_stride]
          );
        
        nodal_matrix_stride += volume_nodes*volume_nodes;
      }
    }
  
}

void
krylov_pc_invstiff_destroy(krylov_pc_t* pc){
  pc->pc_apply = NULL;
  pc->pc_setup = NULL;
  
  krylov_pc_invstiff_data_t* pc_data = (krylov_pc_invstiff_data_t*)pc->pc_data;
  P4EST_FREE(pc_data->inv_stiff);
  P4EST_FREE(pc_data); 
  P4EST_FREE(pc);
}


void
krylov_pc_invstiff_apply(krylov_pc_t* kpc, double* xp, double* yp)
{
  krylov_pc_invstiff_data_t* pc_data = kpc->pc_data;
  p4est_t* p4est = pc_data->p4est;
  
  int nodal_stride = 0;
  int nodal_matrix_stride = 0;
  for (p4est_topidx_t tt = p4est->first_local_tree;
       tt <= p4est->last_local_tree;
       ++tt)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        curved_element_data_t* ed = quad->p.user_data;
        int volume_nodes = dgmath_get_nodes((P4EST_DIM), ed->deg);

        linalg_matvec_plus_vec(1.0, &pc_data->inv_stiff[nodal_matrix_stride], &xp[nodal_stride],
                               0, &yp[nodal_stride], volume_nodes, volume_nodes);
        
        nodal_stride += volume_nodes;
        nodal_matrix_stride += volume_nodes*volume_nodes;
      }
    }
  
}





/* krylov_pc_invstiff.c ends here */
