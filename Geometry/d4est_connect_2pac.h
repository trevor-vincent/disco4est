#ifndef D4EST_CONNECT_2PAC_H
#define D4EST_CONNECT_2PAC_H 

/* #include <pXest.h> */
#include <p4est_connectivity.h>

p4est_connectivity_t *
p4est_connectivity_new_2pac_aligned (void);
p4est_connectivity_t *
p4est_connectivity_new_2pac_nonaligned (void);
p4est_connectivity_t *
p4est_connectivity_new_2pac_aligned_CUBE (void);
p4est_connectivity_t *
p4est_connectivity_new_2pac_nonaligned_CUBE (void);
#endif
