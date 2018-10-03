d4est_solver_schwarz_apply_stiffness_matrix
(
)
{
  for (int i = 0; i < sub_data->num_elements; i++){    
    d4est_mesh_data_on_element_t md_on_e = d4est_mesh_data_on_element
                                           (
                                            d4est_factors,
                                            ed
                                           );
    d4est_quadrature_apply_stiffness_matrix
      (
       d4est_ops,
       d4est_quad,
       d4est_geom,
       &mesh_vol,
       QUAD_OBJECT_VOLUME,
       QUAD_INTEGRAND_UNKNOWN,
       &u[which_field*local_nodes + ed->nodal_stride],           
       ed->deg,
       J_quad,
       md_on_e.rst_xyz_quad,
       ed->deg_quad,
       &Au[which_field*local_nodes + ed->nodal_stride]
      );     
  }
}

void
compute_Au_boundary_dirichlet
(
 double* u_field_over_subdomain,
 double* Au_field_over_subdomain
)
{

}


d4est_solver_schwarz_apply_mortar_matrices_dirichlet
(
)
{
 for (int i = 0; i < sub_data->num_elements; i++){
   for (int f = 0; f < (P4EST_FACES); f++){
     if (boundary){
       compute_Au_boundary_dirichlet(element, u_over_subdomains, Au_over_subdomains);
     }
     else {
       compute_Au_interface(element, u_over_subdomains, Au_over_subdomains);
     }
   }
 }

}


d4est_solver_schwarz_apply_laplacian
(
 
)
{

}
