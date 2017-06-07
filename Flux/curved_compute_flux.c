#include "../Flux/curved_compute_flux.h"
#include "../ElementData/d4est_element_data.h"
#include "../Utilities/util.h"

void curved_compute_flux_on_local_elements_aux(p4est_iter_face_info_t *info,
                                           void *user_data) {
  int i;
  int s_p, s_m;
  int mortar_side_id_m, mortar_side_id_p;
  int mortar_stride;

  curved_compute_flux_user_data_t *curved_compute_flux_user_data =
      (curved_compute_flux_user_data_t *)info->p4est->user_pointer;

  curved_flux_fcn_ptrs_t *flux_fcn_ptrs =
      (curved_flux_fcn_ptrs_t *)curved_compute_flux_user_data->flux_fcn_ptrs;
  
  d4est_operators_t *d4est_ops =
      (d4est_operators_t *)curved_compute_flux_user_data->d4est_ops;

  d4est_geometry_t* geom =
    (d4est_geometry_t *)curved_compute_flux_user_data->geom;

  d4est_quadrature_t* d4est_quad =
    (d4est_quadrature_t *)curved_compute_flux_user_data->d4est_quad;
  
  d4est_element_data_t *ghost_data = (d4est_element_data_t *)user_data;
  d4est_element_data_t *e_p[(P4EST_HALF)];
  d4est_element_data_t *e_m[(P4EST_HALF)];
  int e_m_is_ghost[(P4EST_HALF)];

  mortar_stride = curved_compute_flux_user_data->mortar_stride;
  sc_array_t *sides = &(info->sides);

  /* check if it's an interface boundary, otherwise it's a physical boundary */
  if (sides->elem_count == 2) {
    /* printf("tree_boundary = %d\n", info->tree_boundary); */
    /* printf("tree_orientation = %d\n", info->orientation); */
    p4est_iter_face_side_t *side[2];

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    /* iterate through both sides of the interface */
    for (s_m = 0; s_m < 2; s_m++) {
      /* initialize everything to non-ghost */
      for (i = 0; i < (P4EST_HALF); i++) e_m_is_ghost[i] = 0;

      s_p = (s_m == 0) ? 1 : 0;
      mortar_side_id_m = (s_m == 0) ? mortar_stride : mortar_stride + 1;
      mortar_side_id_p = (s_m == 0) ? mortar_stride + 1: mortar_stride;
        
      /* the minus side is hanging */
      if (side[s_m]->is_hanging) {
        for (i = 0; i < (P4EST_HALF); i++) {
          /* we skip ghosts, b.c those are handled by other processes */
          if (!side[s_m]->is.hanging.is_ghost[i]) {
            e_m[i] = (d4est_element_data_t *)side[s_m]
                         ->is.hanging.quad[i]
                         ->p.user_data;
            e_m_is_ghost[i] = 0;
          } else {
            /* e_m[i] = NULL; */
            e_m[i] = (d4est_element_data_t
                          *)&ghost_data[side[s_m]->is.hanging.quadid[i]];
            /* printf("ghost_id = %d, ghost e_m[i]->deg = %d\n",
             * side[s_m]->is.hanging.quadid[i],  e_m[i]->deg); */
            e_m_is_ghost[i] = 1;
          }
        }

        /* The + side must be full if the - side is hanging */
        /* We allow the + side to be a ghost, b.c we just copy data for this
         * case */
        if (side[s_p]->is.full.is_ghost)
          e_p[0] =
              (d4est_element_data_t *)&ghost_data[side[s_p]->is.full.quadid];
        else
          e_p[0] =
              (d4est_element_data_t *)side[s_p]->is.full.quad->p.user_data;

        int sum_ghost_array = util_sum_array_int(e_m_is_ghost, (P4EST_HALF));

        if (sum_ghost_array != 0 && sum_ghost_array != 4) {
          int dd;
          printf("sum_ghost_array != 0 || sum_ghost_array != 4, it = %d\n",
                 sum_ghost_array);
          for (dd = 0; dd < (P4EST_HALF); dd++) {
            printf(" ghost_id = %d\n", side[s_m]->is.hanging.quadid[dd]);
          }
        }

        /* unless every hanging face is a ghost, we calculate the flux */
        if (sum_ghost_array < (P4EST_HALF)) {
          flux_fcn_ptrs->flux_interface_fcn
            (
             e_m,
             (P4EST_HALF),
             side[s_m]->face,
             mortar_side_id_m,
             e_p,
             1,
             side[s_p]->face,
             mortar_side_id_p,
             e_m_is_ghost,
             info->orientation,
             d4est_ops,
             geom,
             d4est_quad,
             flux_fcn_ptrs->params
          );
        }

      }

      /* If the - side is full */
      else {
        /* We only do calculations if the side is not a ghost */
        if (!side[s_m]->is.full.is_ghost) {
          e_m[0] =
              (d4est_element_data_t *)side[s_m]->is.full.quad->p.user_data;
          // if + side is hanging
          if (side[s_p]->is_hanging) {
            for (i = 0; i < (P4EST_HALF); i++) {
              if (side[s_p]->is.hanging.is_ghost[i])
                e_p[i] = (d4est_element_data_t
                              *)&ghost_data[side[s_p]->is.hanging.quadid[i]];
              else
                e_p[i] = (d4est_element_data_t *)side[s_p]
                             ->is.hanging.quad[i]
                             ->p.user_data;
            }

            flux_fcn_ptrs->flux_interface_fcn
              (
               e_m,
               1,
               side[s_m]->face,
               mortar_side_id_m,
               e_p,
               (P4EST_HALF),
               side[s_p]->face,
               mortar_side_id_p,
               e_m_is_ghost,
               info->orientation,
               d4est_ops,
               geom,
               d4est_quad,
               flux_fcn_ptrs->params
              );

          }

          // if + side is full
          else {
            if (side[s_p]->is.full.is_ghost) {
              e_p[0] = (d4est_element_data_t
                            *)&ghost_data[side[s_p]->is.full.quadid];
            } else
              e_p[0] =
                  (d4est_element_data_t *)side[s_p]->is.full.quad->p.user_data;

            flux_fcn_ptrs->flux_interface_fcn(
                                              e_m, 1, side[s_m]->face, mortar_side_id_m, e_p, 1, side[s_p]->face, mortar_side_id_p, e_m_is_ghost,
                info->orientation, d4est_ops, geom, d4est_quad, flux_fcn_ptrs->params);
          }
        }
      }
    }
    curved_compute_flux_user_data->mortar_stride += 2;
  }

  /* Weak Enforcement of Boundary Conditions */
  else {
    p4est_iter_face_side_t *side;
    side = p4est_iter_fside_array_index_int(sides, 0);
    e_m[0] = (d4est_element_data_t *)side->is.full.quad->p.user_data;
    /* #ifndef NDEBUG */
    /* printf("e_m[0]->deg = %d\n", e_m[0]->deg); */
    /* #endif */
    mortar_side_id_m = mortar_stride;
    mortar_side_id_p = -1;
    
    flux_fcn_ptrs->flux_boundary_fcn(e_m[0], side->face, mortar_side_id_m,
                                     flux_fcn_ptrs->bndry_fcn, d4est_ops, geom, d4est_quad,
                                     flux_fcn_ptrs->params);

    curved_compute_flux_user_data->mortar_stride += 1;
  }
}


int
curved_compute_flux_on_local_elements
(
 p4est_t* p4est,
 p4est_ghost_t* ghost,
 void* ghost_data,
 d4est_operators_t* d4est_ops,
 d4est_geometry_t* d4est_geom,
 d4est_quadrature_t* d4est_quad,
 curved_flux_fcn_ptrs_t* fcn_ptrs,
 curved_compute_flux_exchange_data_option_t option
)
{
  if (option == EXCHANGE_GHOST_DATA)
    p4est_ghost_exchange_data(p4est,ghost,ghost_data);
  
  void* tmpptr = p4est->user_pointer;
  curved_compute_flux_user_data_t curved_compute_flux_user_data;
  curved_compute_flux_user_data.d4est_ops = d4est_ops;
  curved_compute_flux_user_data.geom = d4est_geom;
  curved_compute_flux_user_data.d4est_quad = d4est_quad;
  curved_compute_flux_user_data.mortar_stride = 0;
  
  p4est->user_pointer = &curved_compute_flux_user_data;
  curved_compute_flux_user_data.flux_fcn_ptrs = fcn_ptrs;

  p4est_iterate(p4est,
		ghost,
		(void*) ghost_data,
		NULL,
		curved_compute_flux_on_local_elements_aux,
#if (P4EST_DIM)==3
                NULL,
#endif
		NULL);

  p4est->user_pointer = tmpptr;
  return curved_compute_flux_user_data.mortar_stride;
}
