#ifndef D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H
#define D4EST_SOLVER_SCHWARZ_MORTAR_DATA_H 

typedef struct {

  double* drst_dxyz_mortar_quad [(P4EST_DIM)][(P4EST_DIM)];
  double* sj_on_f_m_mortar_quad;
  double* n_on_f_m_mortar_quad [(P4EST_DIM)];

  d4est_ghost_data_ext_t* mortar_side_ghost_data;

}d4est_solver_schwarz_mortar_data_t;


#endif
