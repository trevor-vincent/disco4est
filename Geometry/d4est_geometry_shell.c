static void
p8est_geometry_shell_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double rst[3], double xyz[3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  d4est_geometry_octree_to_vertex (geom, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform abc[0] and y in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  R = sphere->R0*(2. - abc[2]) + sphere->R1*(abc[2] - 1.);  
  q = R / sqrt (x * x + y * y + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    xyz[0] = +q * y;
    xyz[1] = -q * x;
    xyz[2] = +q;
    break;
  case 2:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  case 1:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 0:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* back */
    xyz[0] = -q * x;
    xyz[1] = +q;
    xyz[2] = +q * y;
    break;
  case 5:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

p8est_geometry_t   *
p8est_geometry_new_shell (p8est_connectivity_t * conn, double R2, double R1)
{
  d4est_geometry_shell_input_t input = d4est_geometry_shell_input(input_file);
  p8est_geometry_t* shell_geom = P4EST_ALLOC(p8est_geometry_t, 1);
  p4est_connectivity_t* conn = p8est_connectivity_new_shell();
  
  d4est_geometry_shell_attr_t* shell_attrs = P4EST_ALLOC(d4est_geometry_shell_attr_t, 1);

  shell_attrs->conn = conn;
  shell_attrs->R2 = R2;
  shell_attrs->R1 = R1;
  shell_attrs->R2byR1 = R2 / R1;
  shell_attrs->R1sqrbyR2 = R1 * R1 / R2;
  shell_attrs->Rlog = log (R2 / R1);

  builtin->geom.name = "shell";
  builtin->geom.user = shell_attrs;
  builtin->geom.X = d4est_geometry_shell_X;

  printf("[GEOMETRY_INFO]: NAME = shell\n");
  printf("[GEOMETRY_INFO]: R1 = %.25f\n", shell_attrs->R1);
  printf("[GEOMETRY_INFO]: R2 = %.25f\n", shell_attrs->R2);
}


