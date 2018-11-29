/*  store laplacian mortar data e.g. the e_m and e_p shit in array
this will speed things up a lot, when we get to a ghost elmeent in a particular subdomain we will grab from d4est_solver_schwarz_mortar_data
 */

/* have is_ghost in d4est_laplacian_flux to check if ghost, if ghost then use */
