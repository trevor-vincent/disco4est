#include <d4est_reference.h>

 static const int d4est_reference_p8est_FToF_code[6][6] = {
  {0, 1, 1, 0, 0, 1}, {2, 0, 0, 1, 1, 0}, {2, 0, 0, 1, 1, 0},
  {0, 2, 2, 0, 0, 1}, {0, 2, 2, 0, 0, 1}, {2, 0, 0, 2, 2, 0}};

 static const int d4est_reference_p8est_code_to_perm[3][4] = {
  {1, 2, 5, 6}, {0, 3, 4, 7}, {0, 4, 3, 7}};

 static const int d4est_reference_p8est_perm_to_order[8][4] = {
  {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {1, 3, 0, 2},
  {2, 0, 3, 1}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};

int d4est_reference_is_child_left_or_right(int c, int dir) {
  if (dir == 0)
    if (c == 0 || c == 2 || c == 4 || c == 6)
      return 0;
    else
      return 1;
  else if (dir == 1)
    if (c == 0 || c == 1 || c == 4 || c == 5)
      return 0;
    else
      return 1;
  else if (dir == 2)
    if (c == 0 || c == 1 || c == 2 || c == 3)
      return 0;
    else
      return 1;
  else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_reference_is_child_left_or_right\n");
    return -1;
  }
}

void d4est_reference_child_to_parent_coords(int c, double* r) {
  *r *= .5;
  if (c == 0) {
    *r -= .5;
  }
  else if (c == 1) {
    *r += .5;
  }
  else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_reference_child_to_parent_coords\n");
  }
}

void d4est_reference_get_normal(int face, int dim, double* n) {
  D4EST_ASSERT((face < 2 * dim) && (dim == 2 || dim == 3));

  double* nx = &n[0];
  double* ny = &n[1];
  double* nz = NULL;

  *nx = 0.;
  *ny = 0.;

  if (dim == 3) {
    nz = &n[2];
    *nz = 0.;
  }
  if (face == 0) {
    *nx = -1.0;
  } else if (face == 1) {
    *nx = 1.0;
  } else if (face == 2) {
    *ny = -1.;
  } else if (face == 3) {
    *ny = 1.;
  } else if (face == 4) {
    *nz = -1.;
  } else if (face == 5) {
    *nz = 1.;
  } else {
    D4EST_ABORT("[D4EST_ERROR]: d4est_reference_get_normal");
  }
}

int d4est_reference_reorient_face_order
(
 int face_dim,
 int f_m,
 int f_p,
 int o,
 int i
)
{
  D4EST_ASSERT((face_dim == 1 && o < 2 && i < 2) || (face_dim == 2 && o < 4 && i < 4));

  if (face_dim == 1){
    if (o == 1){
      return (i == 0) ? 1 : 0;
    }
    else {
      return i;
    }
  }
  else if (face_dim == 2){
    int perm = d4est_reference_p8est_code_to_perm[d4est_reference_p8est_FToF_code[f_m][f_p]][o];
    return d4est_reference_p8est_perm_to_order[perm][i];
  }
  else {
    D4EST_ABORT("FACE DIM == 1, 2\n");
    return -1;
  }
}

int d4est_reference_corner_to_node
(
 int dim,
 int deg,
 int corner
)
{
  D4EST_ASSERT(
             (
              (dim == 3 && corner < 8)
              || (dim == 2 && dim < 4)
             )
            );
  
  int nodes_1d = deg+1;
  int nodes_2d = nodes_1d*nodes_1d;
  int nodes_3d = nodes_1d*nodes_1d*nodes_1d; 

#if (ORDERING == 1)
  if (dim == 2)  {
    if (corner == 0)
      return 0;
    else if (corner == 1)
      return nodes_2d - nodes_1d;
    else if (corner == 2)
      return nodes_1d - 1;
    else if (corner == 3)
      return nodes_2d - 1;
    else{
      D4EST_ABORT("[D4EST_ERROR]: corner < 4\n");
      return -1;
    }
  }

  else if (dim == 3){
    if (corner == 0)
      return 0;
    else if (corner == 1)
      return nodes_3d - nodes_2d;
    else if (corner == 2)
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_3d - nodes_1d;
    else if (corner == 4)
        return nodes_1d - 1;
    else if (corner == 5)
      return nodes_3d - nodes_2d + nodes_1d - 1;
    else if (corner == 6)
      return nodes_2d - 1;
    else if (corner == 7)
      return nodes_3d - 1;
    else{
      D4EST_ABORT("[D4EST_ERROR]: corner < 8\n");
      return -1;
    }
  }

  else {
    D4EST_ABORT("[D4EST_ERROR]: DIM == 3 or DIM == 2\n");
    return -1;
  }
#endif
#if (ORDERING == 3)
    if (dim == 2)  {
    if (corner == 0)
      return 0;
    else if (corner == 1) 
      return nodes_1d - 1;
    else if (corner == 2)
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_2d - 1;
    else{
      D4EST_ABORT("[D4EST_ERROR]: corner < 4\n");
      return -1;
    }
  }

  else if (dim == 3){
    if (corner == 0)
      return 0;
    else if (corner == 1)
        return nodes_1d - 1;
    else if (corner == 2) 
      return nodes_2d - nodes_1d;
    else if (corner == 3)
      return nodes_2d - 1;
    else if (corner == 4)
      return nodes_3d - nodes_2d;
    else if (corner == 5)
      return nodes_3d - nodes_2d + nodes_1d - 1;
    else if (corner == 6)
      return nodes_3d - nodes_1d;
    else if (corner == 7)
      return nodes_3d - 1;
    else{
      D4EST_ABORT("[D4EST_ERROR]: corner < 8\n");
      return -1;
    }
  }

  else {
    D4EST_ABORT("[D4EST_ERROR]: DIM == 3 or DIM == 2\n");
    return -1;
  }
#endif
}

double d4est_reference_rtox(double r, double xl, double h) {
  return h * (r + 1.) / 2. + xl;
}

void d4est_reference_rtox_array(double* r, double xl, double h, double* x, int nodes) {
  int i;
  for (i = 0; i < nodes; i++) {
    x[i] = d4est_reference_rtox(r[i], xl, h);
  }
}
