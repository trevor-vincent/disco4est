#ifndef D4EST_LAPLACIAN_AUX_H
#define D4EST_LAPLACIAN_AUX_H 

typedef enum {EVAL_BNDRY_FCN_NOT_SET, EVAL_BNDRY_FCN_ON_QUAD, EVAL_BNDRY_FCN_ON_LOBATTO} dirichlet_bndry_eval_method_t;

typedef enum {FLUX_SIPG, FLUX_NIPG, FLUX_IIPG, FLUX_NOT_SET} d4est_laplacian_flux_type_t;
typedef enum {BC_ROBIN, BC_DIRICHLET, BC_NOT_SET} d4est_laplacian_bc_t;

#endif
