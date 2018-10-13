#include <pXest.h>
#include <d4est_quadrature.h>
#include <d4est_element_data.h>
#include <d4est_quadrature_legendre.h>
#include <d4est_laplacian_flux_sipg.h>
#include <d4est_geometry.h>
#include <d4est_mesh.h>
#include <d4est_geometry_brick.h>
#include <petscsnes.h>
#include <d4est_linalg.h>
#include <d4est_mortars.h>
#include <d4est_amr.h>
#include <d4est_laplacian.h>
#include <d4est_laplacian_flux.h>
#include <d4est_geometry_brick.h>
#include <d4est_util.h>
#include <limits.h>


static void
problem_set_degrees_init
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg = 2;
  elem_data->deg_quad = 2;
}

static void
problem_set_degrees_amr
(
 d4est_element_data_t* elem_data,
 void* user_ctx
)
{
  elem_data->deg = 2;
  elem_data->deg_quad = 2;
}


static p4est_t*
problem_build_p4est
(
 sc_MPI_Comm mpicomm,
 p4est_connectivity_t* conn,
 p4est_locidx_t min_quadrants,
 int min_level, 
 int fill_uniform
)
{
  return p4est_new_ext
    (
     mpicomm,
     conn,
     min_quadrants,
     min_level,
     fill_uniform,
     sizeof(d4est_element_data_t),
     NULL,
     NULL
    );
}

int main(int argc, char *argv[])
{
  sc_MPI_Comm mpicomm;
  PetscInitialize(&argc,&argv,(char*)0,NULL);
  mpicomm = PETSC_COMM_WORLD;
  
  int proc_size;
  int proc_rank;
  MPI_Comm_size(mpicomm, &proc_size);
  MPI_Comm_rank(mpicomm, &proc_rank);
  p4est_init(NULL, SC_LP_ERROR);
  
  d4est_geometry_t* d4est_geom = P4EST_ALLOC(d4est_geometry_t, 1);
  d4est_geom->geom_type = GEOM_BRICK;
  d4est_geom->DX_compute_method = GEOM_COMPUTE_ANALYTIC;
  d4est_geom->JAC_compute_method = GEOM_COMPUTE_NUMERICAL;



  d4est_geometry_brick_attr_t attr;
  attr.X0 = 0.;
  attr.X1 = 1.;
  attr.Y0 = 0.;
  attr.Y1 = 1.;
  attr.Z0 = 0.;
  attr.Z1 = 1.;

  d4est_geometry_brick_new_aux(d4est_geom, &attr);
  
  p4est_t* p4est = problem_build_p4est
                   (
                    mpicomm,
                    d4est_geom->p4est_conn,
                    -1,
                    1,
                    1
                   );
  
  d4est_operators_t* d4est_ops = d4est_ops_init(20);
  d4est_mesh_data_t* geometric_factors = d4est_mesh_geometry_storage_init(p4est);
  d4est_quadrature_t* d4est_quad = P4EST_ALLOC(d4est_quadrature_t, 1);
  d4est_quad->quad_type = QUAD_TYPE_GAUSS_LEGENDRE;
  d4est_quadrature_legendre_new(d4est_quad, d4est_geom,"");
      
  p4est_ghost_t* ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  d4est_element_data_t* ghost_data = P4EST_ALLOC (d4est_element_data_t,
                                                   ghost->ghosts.elem_count);


  d4est_amr_t* d4est_amr_uniform = d4est_amr_init_uniform_h(p4est,7,5);

  int local_nodes = d4est_mesh_update
                    (
                     p4est,
                     ghost,
                     ghost_data,
                     d4est_ops,
                     d4est_geom,
                     d4est_quad,
                     geometric_factors,
                     INITIALIZE_QUADRATURE_DATA,
                     INITIALIZE_GEOMETRY_DATA,
                     INITIALIZE_GEOMETRY_ALIASES,
                     problem_set_degrees_init,
                     NULL
                    );
  
  double* poly_vec = P4EST_ALLOC(double, local_nodes);
  int same = 1;
  int same2 = 1;

  
  for (int level = 0; level < d4est_amr_uniform->num_of_amr_steps; ++level){

    printf("p4est->local_num_quadrants = %d\n", p4est->local_num_quadrants);
    
    local_nodes = d4est_mesh_update
                  (
                   p4est,
                   ghost,
                   ghost_data,
                   d4est_ops,
                   d4est_geom,
                   d4est_quad,
                   geometric_factors,
                   INITIALIZE_QUADRATURE_DATA,
                   INITIALIZE_GEOMETRY_DATA,
                   INITIALIZE_GEOMETRY_ALIASES,
                   problem_set_degrees_amr,
                   NULL
                  );

    int num_in_bin [(P4EST_CHILDREN)] = {0};

      for (p4est_topidx_t tt = p4est->first_local_tree;
           tt <= p4est->last_local_tree; tt++)
    {
      p4est_tree_t* tree = p4est_tree_array_index (p4est->trees, tt);
      sc_array_t* tquadrants = &tree->quadrants;
      int Q = (p4est_locidx_t) tquadrants->elem_count;
      for (int q = 0; q < Q; ++q) {
        p4est_quadrant_t* quad = p4est_quadrant_array_index (tquadrants, q);
        d4est_element_data_t* ed = (d4est_element_data_t*)(quad->p.user_data);

        double h = P4EST_QUADRANT_LEN(quad->level)/(double)P4EST_ROOT_LEN;
        double jacobian = pow(.5*h, (P4EST_DIM) ); 
        double surface_jacobian = pow(.5*h, (P4EST_DIM) - 1);
        double area = h*h;
        double volume = h*h*h;
        double diam_face = sqrt(h*h + h*h);
        double diam_volume = sqrt(h*h + h*h + h*h);
        double j_div_sj = jacobian/surface_jacobian;

        printf("element %d\n", ed->id);
        printf("volume = %f volume_theor = %f\n", ed->volume, volume);
        printf("diam = %f diam_theor = %f\n", ed->diam_volume, diam_volume);
        
        for (int f = 0; f < (P4EST_FACES); f++){
          printf("f %d, area = %f area_theor = %f\n", f, ed->area[f], area);
          printf("f %d, diam = %f diam_theor = %f\n", f, ed->diam_face[f], diam_face);
          printf("f %d, j_div_sj = %f j_div_sj_theor = %f\n", f, ed->j_div_sj_min[f], j_div_sj);
        }

        int bin = d4est_geometry_brick_which_child_of_root(ed->q, ed->dq);
        num_in_bin[bin] += 1;
      }
    }

      for (int i = 0; i < (P4EST_CHILDREN); i++){
        printf("num_in_bin[%d] = %d, theor = %d\n", i, num_in_bin[i], p4est->local_num_quadrants/(P4EST_CHILDREN));
      }
      
    
    d4est_amr_step
      (
       p4est,
       &ghost,
       &ghost_data,
       d4est_ops,
       d4est_amr_uniform,
       &poly_vec,
       NULL
      );
    
  }

    

    
  if (ghost) {
    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
  }
    
  d4est_mesh_geometry_storage_destroy(geometric_factors);
  d4est_quadrature_destroy(p4est, d4est_ops, d4est_geom, d4est_quad);
  d4est_amr_destroy(d4est_amr_uniform);
  d4est_ops_destroy(d4est_ops);
  p4est_destroy(p4est);
  d4est_geometry_destroy(d4est_geom);
  
  PetscFinalize();

}
