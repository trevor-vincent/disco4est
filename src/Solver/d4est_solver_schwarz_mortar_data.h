#ifndef D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H
#define D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H 

typedef struct {

  double* drst_dxyz_m_mortar_quad;
  double* drst_dxyz_p_mortar_quad_porder;
  double* sj_m_mortar_quad;
  double* n_m_mortar_quad;
  double* hm_mortar_quad;
  double* hp_mortar_quad;
  double* xyz_m_mortar_quad;
  double* xyz_m_mortar_lobatto;

  d4est_ghost_data_ext_t* mortar_side_ghost_data;

}d4est_solver_schwarz_mortar_data_t;


#endif
