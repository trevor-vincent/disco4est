#ifndef D4EST_LAPLACIAN_WITH_OPT_H
#define D4EST_LAPLACIAN_WITH_OPT_H 


/* This file was automatically generated.  Do not edit! */
void d4est_laplacian_flux_with_opt_destroy(d4est_laplacian_flux_with_opt_data_t *data);
d4est_laplacian_flux_with_opt_data_t *d4est_laplacian_flux_with_opt_new(p4est_t *p4est,const char *input_file,d4est_laplacian_bc_t bc_type,void *bc_data);
d4est_mortars_fcn_ptrs_t d4est_laplacian_flux_with_opt_fetch_fcns(d4est_laplacian_flux_with_opt_data_t *data);


#endif
