#define _GNU_SOURCE
#include <pXest.h>
#include <d4est_vtk.h>
#include <d4est_geometry.h>
#include <d4est_util.h>
#include <d4est_mesh.h>
#include <d4est_util.h>
#include <ini.h>
#include <zlog.h>

#ifdef P4_TO_P8
#define D4EST_VTK_CELL_TYPE     11      /* VTK_VOXEL */
#else
#define D4EST_VTK_CELL_TYPE      8      /* VTK_PIXEL */
#endif /* !P4_TO_P8 */

#include <sc_io.h>

/* default parameters for the vtk context */
static const double d4est_vtk_scale = 0.95;
static const int    d4est_vtk_continuous = 0;

/* default parameters for d4est_vtk_write_file */
static const int    d4est_vtk_write_tree = 1;
static const int    d4est_vtk_write_level = 1;
static const int    d4est_vtk_write_rank = 1;
static const int    d4est_vtk_wrap_rank = 0;

/* #define D4EST_VTK_DEBUG */

#define D4EST_VTK_DOUBLES
#ifndef D4EST_VTK_DOUBLES
#define D4EST_VTK_FLOAT_NAME "Float32"
#define D4EST_VTK_FLOAT_TYPE float
#else
#define D4EST_VTK_FLOAT_NAME "Float64"
#define D4EST_VTK_FLOAT_TYPE double
#endif

/* #ifndef D4EST_VTK_BINARY */
/* #define D4EST_VTK_ASCII 1 */
/* #define D4EST_VTK_FORMAT_STRING "ascii" */
/* #else */
/* #define D4EST_VTK_FORMAT_STRING "binary */
/* #endif /\* D4EST_VTK_BINARY *\/ */

typedef enum {D4EST_VTK_INT, D4EST_VTK_FLOAT} d4est_vtk_storage_type_t;
typedef enum {D4EST_VTK_DG_GRID, D4EST_VTK_CORNER_GRID, D4EST_VTK_GRID_NOT_SET} d4est_vtk_grid_type_t;
typedef enum {D4EST_VTK_ASCII, D4EST_VTK_BINARY, D4EST_VTK_ZLIB_BINARY, D4EST_VTK_OUTPUT_NOT_SET} d4est_vtk_output_type_t;

/** Opaque context type for writing VTK output with multiple function calls.
 *
 * This structure holds all the information needed for the p4est vtk context.
 * It is used to relay necessary vtk information to the \b d4est_vtk_write_*
 * functions. This structure is initialized by \ref d4est_vtk_write_header and
 * destroyed by \b d4est_vtk_write_footer; it can also be destroyed manually
 * using the \b d4est_vtk_context_destroy function if necessary.
 *
 * The \a p4est member is a pointer to the local p4est.
 * The \a geom member is a pointer to the geometry used to create the p4est.
 * The \a num_points member holds the number of nodes present in the vtk output;
 * this is determined in \ref d4est_vtk_write_header using the \a scale parameter
 * and is used to assure the proper number of point variables are provided.
 * The \a filename member holds the vtk file basename: for error reporting.
 * The \a vtufilename, \a pvtufilename, and \a visitfilename members are the
 * vtk file names.
 * The \a vtufile, \a pvtufile, and \a visitfile members are the vtk file
 * pointers; opened by \ref d4est_vtk_write_header and closed by \b
 * d4est_vtk_write_footer.
 *
 */

struct d4est_vtk_context
{
  /* data passed initially */
  p4est_t            *p4est;       /**< The p4est structure must be alive. */
  char               *filename;    /**< Original filename provided is copied. */
  int* deg_array; /* for writing discontinuous Galerkin data */
  d4est_operators_t* d4est_ops;
  
  /* parameters that can optionally be set in a context */
  d4est_geometry_t   *geom;        /**< The geometry may be NULL. */
  double              scale;       /**< Parameter to shrink quadrants. */
  int                 continuous;  /**< Assume continuous point data? */

  /* internal context data */

  int                 writing;     /**< True after d4est_vtk_write_header. */
  p4est_locidx_t      num_corners; /**< Number of local element corners. */
  p4est_locidx_t      num_points;  /**< Number of VTK points written. */
  p4est_locidx_t      num_cells;  /**< Number of VTK points written. */
  p4est_locidx_t     *node_to_corner;     /**< Map a node to an element corner. */
  p4est_nodes_t      *nodes;       /**< NULL? depending on scale/continuous. */
  char                vtufilename[BUFSIZ];   /**< Each process writes one. */
  char                pvtufilename[BUFSIZ];  /**< Only root writes this one. */
  char                visitfilename[BUFSIZ]; /**< Only root writes this one. */
  FILE               *vtufile;     /**< File pointer for the VTU file. */
  FILE               *pvtufile;    /**< Paraview meta file. */
  FILE               *visitfile;   /**< Visit meta file. */

  d4est_vtk_output_type_t output_type;
  d4est_vtk_grid_type_t grid_type;
  char* input_section;
  char* geometry_section;
  char* folder;
  int write_tree;
  int write_level;
  int write_rank;
  int wrap_rank;
  int write_deg;
    
  
};


static int
d4est_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                        size_t byte_length, d4est_vtk_output_type_t output_type)
{
  D4EST_ABORT("Binary currently has a bug in it, so use ascii output_type");
  if (output_type == D4EST_VTK_BINARY){
    return sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
  }
#ifdef D4EST_USE_ZLIB
  else if (output_type == D4EST_VTK_ZLIB_BINARY){
    return sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
  }
#endif
  else {
    D4EST_ABORT("wrong vtk type");
    return 1;
  }
}

static const char*
d4est_vtk_get_format_string(d4est_vtk_output_type_t output_type){
  if (output_type == D4EST_VTK_ASCII){
    return "ascii";
  }
  else {
    return "binary";
  }
}


static
int d4est_vtk_input_handler
(
 void* user,
 const char* section,
 const char* name,
 const char* value
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");

  d4est_vtk_context_t * pconfig = (d4est_vtk_context_t *)user;

  if (d4est_util_match_couple(section,pconfig->input_section,name,"filename")) {
    D4EST_ASSERT(pconfig->filename == NULL);
    asprintf(&pconfig->filename,"%s",value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"write_tree")) {
    D4EST_ASSERT(pconfig->write_tree == -1);
    pconfig->write_tree = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"write_rank")) {
    D4EST_ASSERT(pconfig->write_rank == -1);
    pconfig->write_rank = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"wrap_rank")) {
    D4EST_ASSERT(pconfig->wrap_rank == -1);
    pconfig->wrap_rank = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"write_deg")) {
    D4EST_ASSERT(pconfig->write_deg == -1);
    pconfig->write_deg = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"write_level")) {
    D4EST_ASSERT(pconfig->write_level == -1);
    pconfig->write_level = atoi(value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"geometry_section")) {
    D4EST_ASSERT(pconfig->geometry_section == NULL);
    asprintf(&pconfig->geometry_section,"%s",value);
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"output_type")) {
    if(d4est_util_match(value, "ascii")){
      pconfig->output_type = D4EST_VTK_ASCII;
    }
    else if (d4est_util_match(value, "binary")){
      pconfig->output_type = D4EST_VTK_BINARY;
    }
    else if (d4est_util_match(value, "zlib_binary")){
      pconfig->output_type = D4EST_VTK_ZLIB_BINARY;
    }
    else {
      zlog_error(c_default, "You tried to use %s as an output type.", value);
      D4EST_ABORT("This output is not supported.");
    }
  }
  else if (d4est_util_match_couple(section,pconfig->input_section,name,"grid_type")) {
    if(d4est_util_match(value, "dg")){
      pconfig->grid_type = D4EST_VTK_DG_GRID;
    }
    else if (d4est_util_match(value, "corner")){
      pconfig->grid_type = D4EST_VTK_CORNER_GRID;
    }
    else {
      zlog_error(c_default, "You tried to use %s as a grid type.", value);
      D4EST_ABORT("This grid is not supported.");
    }
  }
  else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

static void
d4est_vtk_input
(
 int mpirank,
 const char* input_file,
 const char* input_section,
 d4est_vtk_context_t* cont
){
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  
  cont->filename = NULL;
  cont->geometry_section = NULL;
  cont->output_type = D4EST_VTK_OUTPUT_NOT_SET;
  cont->grid_type = D4EST_VTK_GRID_NOT_SET;
  cont->write_level = -1;
  cont->write_rank = -1;
  cont->write_deg = -1;
  cont->write_tree = -1;
  cont->wrap_rank = -1;
  asprintf(&cont->input_section,"%s",input_section);
  
  if (ini_parse(input_file, d4est_vtk_input_handler,cont) < 0) {
    D4EST_ABORT("Can't load input file");
  }

  D4EST_CHECK_INPUT(input_section, cont->filename, NULL);
  D4EST_CHECK_INPUT(input_section, cont->geometry_section, NULL);
  D4EST_CHECK_INPUT(input_section, cont->output_type, D4EST_VTK_OUTPUT_NOT_SET);
  D4EST_CHECK_INPUT(input_section, cont->grid_type, D4EST_VTK_GRID_NOT_SET);
  D4EST_CHECK_INPUT(input_section, cont->wrap_rank, -1);
  D4EST_CHECK_INPUT(input_section, cont->write_rank, -1);
  D4EST_CHECK_INPUT(input_section, cont->write_tree, -1);
  D4EST_CHECK_INPUT(input_section, cont->write_level, -1);
  D4EST_CHECK_INPUT(input_section, cont->write_deg, -1);

  if(mpirank == 0){
    zlog_debug(c_default, "Saving %s with geometry specified by %s", cont->filename, cont->geometry_section);
  if (cont->output_type == D4EST_VTK_ASCII)
    zlog_debug(c_default, "Saving %s in ASCII format", cont->filename);
  if (cont->output_type == D4EST_VTK_BINARY)
    zlog_debug(c_default, "Saving %s in BINARY format", cont->filename);
  if (cont->output_type == D4EST_VTK_ZLIB_BINARY)
    zlog_debug(c_default, "Saving %s in ZLIB_BINARY format", cont->filename);
  if (cont->grid_type == D4EST_VTK_DG_GRID)
    zlog_debug(c_default, "Saving %s in with DG grid", cont->filename);
  if (cont->grid_type == D4EST_VTK_CORNER_GRID)
    zlog_debug(c_default, "Saving %s in with CORNER grid", cont->filename);
  }
  
}


static d4est_vtk_context_t *
d4est_vtk_dg_context_new
(
  p4est_t *p4est,
  d4est_operators_t *d4est_ops,
  const char *input_file,
  const char *input_section
)
{
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");

#ifdef D4EST_VTK_DEBUG
  zlog_debug(c_default, "Starting d4est_vtk_context_new");
#endif
  d4est_vtk_context_t *cont;

  P4EST_ASSERT (p4est != NULL);
  /* Allocate, initialize the vtk context.  Important to zero all fields. */
  cont = P4EST_ALLOC_ZERO (d4est_vtk_context_t, 1);
  cont->d4est_ops = d4est_ops;
  cont->p4est = p4est;

  d4est_vtk_input(p4est->mpirank, input_file, input_section, cont);
  
  cont->deg_array = NULL;
  cont->scale = d4est_vtk_scale;
  cont->continuous = d4est_vtk_continuous;

  return cont;
}


static void
d4est_vtk_context_set_geom (d4est_vtk_context_t * cont,
                            d4est_geometry_t * geom)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_set_geom");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->geom = geom;
}

static void
d4est_vtk_context_set_deg_array (d4est_vtk_context_t * cont,
                                 int* deg_array)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_set_geom");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->deg_array = deg_array;
}

static void
d4est_vtk_context_set_scale (d4est_vtk_context_t * cont, double scale)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_set_scale");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);
  P4EST_ASSERT (0. < scale && scale <= 1.);

  cont->scale = scale;
}

static void
d4est_vtk_context_set_continuous (d4est_vtk_context_t * cont, int continuous)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_set_continuous");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->continuous = continuous;
}

static void
d4est_vtk_context_destroy (d4est_vtk_context_t * context)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_destroy");
#endif
  P4EST_ASSERT (context != NULL);
  P4EST_ASSERT (context->p4est != NULL);

  /* since this function is called inside write_header and write_footer, */
  /* we cannot assume a consistent state of all member variables */
  
  P4EST_ASSERT (context->filename != NULL);
  free(context->filename);
  free(context->folder);
  free(context->input_section);
  free(context->geometry_section);
  
  /* deallocate node storage */
  if (context->nodes != NULL) {
    p4est_nodes_destroy (context->nodes);
  }
  P4EST_FREE (context->node_to_corner);

  /* Close all file pointers. */
  if (context->vtufile != NULL) {
    if (fclose (context->vtufile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->vtufilename);
    }
    context->vtufile = NULL;
  }

  /* Close paraview master file */
  if (context->pvtufile != NULL) {
    /* Only the root process opens/closes these files. */
    P4EST_ASSERT (context->p4est->mpirank == 0);
    if (fclose (context->pvtufile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->pvtufilename);
    }
    context->pvtufile = NULL;
  }

  /* Close visit master file */
  if (context->visitfile != NULL) {
    /* Only the root process opens/closes these files. */
    P4EST_ASSERT (context->p4est->mpirank == 0);
    if (fclose (context->visitfile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->visitfilename);
    }
    context->visitfile = NULL;
  }

  /* Free context structure. */
  P4EST_FREE (context);
}

static d4est_vtk_context_t *
d4est_vtk_write_header (d4est_vtk_context_t * cont,
                           d4est_operators_t* d4est_ops,
                           d4est_vtk_output_type_t output_type)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_write_header");
#endif
  const double        intsize = 1.0 / P4EST_ROOT_LEN;
  int                 mpirank;
  int                 conti;
  double              scale;
  const char         *filename;
  const double       *v;
  const p4est_topidx_t *tree_to_vertex;
  p4est_topidx_t      first_local_tree, last_local_tree;
  p4est_locidx_t      Ncells, Ncorners;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  d4est_geometry_t   *geom;
  double              wx, wy, wz;
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
  int                 xi, yi, j, k;
#ifdef P4_TO_P8
  int                 zi;
#endif
  double              h2, eta_x, eta_y, eta_z = 0.;
  double              xyz[3], XYZ[3];   /* 3 not P4EST_DIM */
  size_t              num_quads, zz;
  p4est_topidx_t      jt;
  p4est_topidx_t      vt[P4EST_CHILDREN];
  p4est_locidx_t      quad_count, Npoints;
  p4est_locidx_t      sk, il, ntcid, *ntc;
  D4EST_VTK_FLOAT_TYPE *float_data;
  sc_array_t         *quadrants, *indeps;
  sc_array_t         *trees;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_nodes_t      *nodes;
  p4est_indep_t      *in;
  int* deg_array;
  /* check a whole bunch of assertions, here and below */
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  /* avoid uninitialized warning */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = -(k + 1);
  }
 
  /* from now on this context is officially in use for writing */
  cont->writing = 1;

  /* grab context variables */
  p4est = cont->p4est;
  filename = cont->filename;
  geom = cont->geom;
  scale = cont->scale;
  conti = cont->continuous;
  deg_array  = cont->deg_array;
  P4EST_ASSERT (filename != NULL);
  P4EST_ASSERT (deg_array != NULL);

  /* grab details from the forest */
  P4EST_ASSERT (p4est != NULL);
  mpirank = p4est->mpirank;
  connectivity = p4est->connectivity;
  P4EST_ASSERT (connectivity != NULL);
  v = connectivity->vertices;
  tree_to_vertex = connectivity->tree_to_vertex;
  if (geom == NULL) {
    SC_CHECK_ABORT (connectivity->num_vertices > 0,
                    "Must provide connectivity with vertex information");
    P4EST_ASSERT (v != NULL && tree_to_vertex != NULL);
  }
  trees = p4est->trees;
  first_local_tree = p4est->first_local_tree;
  last_local_tree = p4est->last_local_tree;

  Ncells = 0;
  for (int i = 0; i < p4est->local_num_quadrants; i++){
    Ncells += d4est_util_int_pow_int(deg_array[i], (P4EST_DIM));
  }
  /* Ncells = p4est->local_num_quadrants; */

  cont->num_cells = Ncells;
  
  cont->num_corners = Ncorners = P4EST_CHILDREN * Ncells;

  /* printf("cont->num_corners = %d\n",cont->num_corners); */
  /* printf("Ncells = %d\n\n",Ncells); */
  
  if (scale < 1. || !conti) {
    /* when we scale the quadrants we need each corner separately */
    cont->nodes = nodes = NULL;
    cont->num_points = Npoints = Ncorners;
    cont->node_to_corner = ntc = NULL;
    indeps = NULL;
  }
  else {
    D4EST_ABORT("We do not support scale = 1. yet");
  }

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s%s_%04d.vtu", cont->folder, filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", cont->vtufilename);
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
  if (output_type == D4EST_VTK_ZLIB_BINARY){
    fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
  }
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Npoints, (long long) Ncells);
  fprintf (cont->vtufile, "      <Points>\n");

  float_data = P4EST_ALLOC (D4EST_VTK_FLOAT_TYPE, 3 * Npoints);

  /* write point position data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, d4est_vtk_get_format_string(output_type));

  /* if (nodes == NULL) { */
  /* loop over the trees */

  int stride = 0;
  for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
    tree = p4est_tree_array_index (trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    /* retrieve corners of the tree */
    if (geom == NULL) {
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];
      }
    }
    else {
      /* provoke crash on logic bug */
      P4EST_ASSERT (vt[0] == -1);
      v = NULL;
    }
      
    /* loop over the elements in tree and calculate vertex coordinates */
    for (zz = 0; zz < num_quads; ++zz, ++quad_count) {


      double* vtk_rst = d4est_operators_fetch_vtk_rst
                        (
                         d4est_ops,
                         deg_array[quad_count],
                         (P4EST_DIM)
                        );


      quad = p4est_quadrant_array_index (quadrants, zz);
      h2 = .5 * intsize * P4EST_QUADRANT_LEN (quad->level);
      k = 0;

      int num_cells_in_element = d4est_util_int_pow_int(deg_array[quad_count], (P4EST_DIM));
      for (int ec = 0; ec < num_cells_in_element; ec++) {
        for (int corn = 0; corn < (P4EST_CHILDREN); corn++){

          double r_01 = d4est_reference_rtox(vtk_rst[0 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
          double s_01 = d4est_reference_rtox(vtk_rst[1 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
          double t_01 = d4est_reference_rtox(vtk_rst[2 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
          eta_x = intsize * quad->x + h2 * (1. + (r_01 * 2 - 1) * scale);
          eta_y = intsize * quad->y + h2 * (1. + (s_01 * 2 - 1) * scale);
          eta_z = 0.;
#if (P4EST_DIM==3)
          eta_z = intsize * quad->z + h2 * (1. + (t_01 * 2 - 1) * scale);
#endif
          if (geom != NULL) {
            xyz[0] = eta_x;
            xyz[1] = eta_y;
            xyz[2] = eta_z;
            geom->X (geom, jt, (p4est_qcoord_t [(P4EST_DIM)]){0}, -(P4EST_ROOT_LEN), xyz, COORDS_TREE_UNITCUBE, XYZ);
            /* printf("xyz, XYZ = %f, %f, %f, %f,%f,%f\n", xyz[0], xyz[1], xyz[2], XYZ[0], XYZ[1], XYZ[2]); */
            for (j = 0; j < 3; ++j) {
              float_data[stride + corn*3 + ec*3*(P4EST_CHILDREN) + j] =
                (D4EST_VTK_FLOAT_TYPE) XYZ[j];
            }
          }
          else {
            for (j = 0; j < 3; ++j) {
              /* *INDENT-OFF* */
              xyz[j] =
                ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                 eta_x  * v[3 * vt[1] + j]) +
                                 eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                           eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
                 +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                 eta_x  * v[3 * vt[5] + j]) +
                                 eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                           eta_x  * v[3 * vt[7] + j]))
#endif
                );
              /* *INDENT-ON* */
              float_data[stride + corn*3 + ec*3*(P4EST_CHILDREN) + j] =
                (D4EST_VTK_FLOAT_TYPE) XYZ[j];
            }
          }
        }
      }
      stride += num_cells_in_element*3*(P4EST_CHILDREN);
    }
  }

  if (output_type == D4EST_VTK_ASCII){
    
    for (il = 0; il < Npoints; ++il) {
      wx = float_data[3 * il + 0];
      wy = float_data[3 * il + 1];
      wz = float_data[3 * il + 2];

      fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
               "     %24.16e %24.16e %24.16e\n",
#else
               "          %16.8e %16.8e %16.8e\n",
#endif
               wx, wy, wz);
    }
  }
  else {
    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                     sizeof (*float_data) * 3 * Npoints,
                                     output_type);
    fprintf (cont->vtufile, "\n");
    if (retval) {
      d4est_vtk_context_destroy (cont);
      P4EST_FREE (float_data);
      return NULL;
    }
  }
  
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

  if (output_type == D4EST_VTK_ASCII){
    for (sk = 0, il = 0; il < Ncells; ++il) {
      fprintf (cont->vtufile, "         ");
      for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
        fprintf (cont->vtufile, " %lld", nodes == NULL ?
                 (long long) sk : (long long) nodes->local_nodes[sk]);
      }
      fprintf (cont->vtufile, "\n");
    }
  }
  else {
    fprintf (cont->vtufile, "          ");
    if (nodes == NULL) {
      locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncorners);
      for (il = 0; il < Ncorners; ++il) {
        locidx_data[il] = il;
      }
      retval =
        d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                sizeof (p4est_locidx_t) * Ncorners, output_type);
      P4EST_FREE (locidx_data);
    }
    else {
      retval =
        d4est_vtk_write_binary (cont->vtufile, (char *) nodes->local_nodes,
                                sizeof (p4est_locidx_t) * Ncorners, output_type);
    }
    fprintf (cont->vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

  }
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

  fprintf (cont->vtufile, "         ");
  if (output_type == D4EST_VTK_ASCII){
    for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
      if (!(sk % 8) && il != Ncells)
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
  }
  else {
    locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
    for (il = 1; il <= Ncells; ++il)
      locidx_data[il - 1] = P4EST_CHILDREN * il;  /* same type */

    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                     sizeof (p4est_locidx_t) * Ncells,output_type);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (locidx_data);

    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", d4est_vtk_get_format_string(output_type));
  fprintf (cont->vtufile, "         ");

  if (output_type == D4EST_VTK_ASCII){
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", D4EST_VTK_CELL_TYPE);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
  }
  else {
    uint8_data = P4EST_ALLOC (uint8_t, Ncells);
    for (il = 0; il < Ncells; ++il)
      uint8_data[il] = D4EST_VTK_CELL_TYPE;

    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                                     sizeof (*uint8_data) * Ncells, output_type);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (uint8_data);

    if (retval) {
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

  }
  
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s%s.pvtu", cont->folder,filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
    if (output_type == D4EST_VTK_ZLIB_BINARY){
      fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
    }
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             D4EST_VTK_FLOAT_NAME, d4est_vtk_get_format_string(output_type));
    fprintf (cont->pvtufile, "    </PPoints>\n");



    fprintf (cont->pvtufile, "    <PCells>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"connectivity\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "Int32", d4est_vtk_get_format_string(output_type));
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"offsets\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "Int32", d4est_vtk_get_format_string(output_type));
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"types\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "UInt8", d4est_vtk_get_format_string(output_type));
    fprintf (cont->pvtufile, "    </PCells>\n");
    
    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in d4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s%s.visit", cont->folder,filename);
    cont->visitfile = fopen (cont->visitfilename, "wb");
    if (!cont->visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->visitfilename);
      d4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  /* the nodes object is no longer needed */
  if (nodes != NULL) {
    p4est_nodes_destroy (cont->nodes);
    cont->nodes = NULL;
  }
  return cont;
}

static d4est_vtk_context_t *
d4est_vtk_write_cell_scalar
(
 d4est_vtk_context_t * cont,
 const char *scalar_name,
 void* values,
 d4est_vtk_storage_type_t storage_type,
 d4est_vtk_output_type_t output_type
)
{
  P4EST_ASSERT (cont != NULL && cont->writing);

  /* Write cell data. */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           (storage_type == D4EST_VTK_FLOAT) ? (D4EST_VTK_FLOAT_NAME) : "Int32",
           scalar_name,
           d4est_vtk_get_format_string(output_type));

  /* for (il = 0; il < Ncells; ++il) { */
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  const p4est_locidx_t Ncells = cont->num_cells;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
  p4est_topidx_t      jt;
  p4est_locidx_t      il;
  p4est_locidx_t sk;
  p4est_locidx_t tc;

  if (output_type == D4EST_VTK_ASCII){
  for (tc = 0, il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
    tree = p4est_tree_array_index (trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;
    for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
      int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
      for (int ec = 0; ec < num_cells_in_element; ec++, tc++){

        if(storage_type == D4EST_VTK_FLOAT){
        fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
                 "     %24.16e\n",
#else
                 /* "          %16.8e\n", */
#endif
                 ((double*)values)[il]);
        }
        else {
        fprintf (cont->vtufile,
                 "     %d\n",
                 ((int*)values)[il]);
        }
      }
    }
  }
  }
  else {
    D4EST_ABORT("Currently, binary is not supported");
    /* D4EST_VTK_FLOAT_TYPE* float_data = P4EST_ALLOC (D4EST_VTK_FLOAT_TYPE, Ncells); */
  /* for (tc = 0, il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) { */
    /* tree = p4est_tree_array_index (trees, jt); */
    /* quadrants = &tree->quadrants; */
    /* num_quads = quadrants->elem_count; */
    /* for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) { */
      /* int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM)); */
      /* for (int ec = 0; ec < num_cells_in_element; ec++, tc++){ */
      /* float_data[il] = */
        /* (D4EST_VTK_FLOAT_TYPE) (values)[il]; */
    /* } */
    /* } */
  /* } */
    /* fprintf (cont->vtufile, "          "); */
    /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
    /* int retval = d4est_vtk_write_binary (cont->vtufile, (char *) float_data, */
                                     /* sizeof (*float_data) * Ncells, output_type); */
    /* fprintf (cont->vtufile, "\n"); */

    /* P4EST_FREE (float_data); */

    /* if (retval) { */
      /* d4est_vtk_context_destroy (cont); */
      /* return NULL; */
    /* } */
  }
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell scalar file\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}


static int
d4est_vtk_write_footer (d4est_vtk_context_t * cont)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_write_footer");
#endif
  int                 p;
  int                 procRank = cont->p4est->mpirank;
  int                 numProcs = cont->p4est->mpisize;

  P4EST_ASSERT (cont != NULL && cont->writing);

  fprintf (cont->vtufile, "    </Piece>\n");
  fprintf (cont->vtufile, "  </UnstructuredGrid>\n");
  fprintf (cont->vtufile, "</VTKFile>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing footer\n");
    d4est_vtk_context_destroy (cont);
    return -1;
  }

  /* Only have the root write to the parallel vtk file */
  if (procRank == 0) {
    fprintf (cont->visitfile, "!NBLOCKS %d\n", numProcs);

    /* Write data about the parallel pieces into both files */
    for (p = 0; p < numProcs; ++p) {
      fprintf (cont->pvtufile,
               "    <Piece Source=\"%s_%04d.vtu\"/>\n", cont->filename, p);
      fprintf (cont->visitfile, "%s_%04d.vtu\n", cont->filename, p);
    }
    fprintf (cont->pvtufile, "  </PUnstructuredGrid>\n");
    fprintf (cont->pvtufile, "</VTKFile>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      d4est_vtk_context_destroy (cont);
      return -1;
    }

    if (ferror (cont->visitfile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      d4est_vtk_context_destroy (cont);
      return -1;
    }
  }

  /* Destroy context structure. */
  d4est_vtk_context_destroy (cont);

  return 0;
}


static double*
d4est_vtk_convert_nodal_to_vtk(
                               p4est_t* p4est,
                               d4est_vtk_context_t* cont,
                               d4est_operators_t* d4est_ops,
                               double* values,
                               d4est_vtk_grid_type_t grid_type
)
{

  double* vtk_array = P4EST_ALLOC(double, cont->num_points);
  sc_array_t         *quadrants, *indeps;
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  size_t              num_quads, zz;
  int vtk_stride = 0;
  int nodal_stride = 0;
  int il, sk, jt;
  for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
    tree = p4est_tree_array_index (trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    for (zz = 0; zz < num_quads; ++zz, ++il) {
      int num_nodes_in_element = d4est_util_int_pow_int(cont->deg_array[il] + 1, (P4EST_DIM));
      int num_points_in_element = (grid_type == D4EST_VTK_DG_GRID) ?
                                  d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM))*(P4EST_CHILDREN) :
                                  (P4EST_CHILDREN);
      if (grid_type == D4EST_VTK_DG_GRID){
        d4est_operators_convert_nodal_to_vtk
          (
           d4est_ops,
           &values[nodal_stride],
           (P4EST_DIM),
           cont->deg_array[il],
           &vtk_array[vtk_stride]
          );
      }
      else if (grid_type == D4EST_VTK_CORNER_GRID){
        for (int i = 0; i < (P4EST_CHILDREN); i++){
          vtk_array[vtk_stride + i]
            = values[d4est_reference_corner_to_node((P4EST_DIM), cont->deg_array[il],i)];
        }
      }
      vtk_stride += num_points_in_element;
      nodal_stride += num_nodes_in_element;
    }
  }

  return vtk_array;
}



static d4est_vtk_context_t *
d4est_vtk_write_point_scalar (d4est_vtk_context_t * cont,
                              const char *scalar_name,
                              double * values,
                              d4est_vtk_grid_type_t grid_type,
                              d4est_vtk_output_type_t output_type
                             )
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_write_point_scalar \n");
#endif
  p4est_locidx_t      il, ddl;
  int                 use_nodes;
// #ifdef P4EST_ENABLE_DEBUG
//   int                 Ncorners;
// #endif
  int                 Npoints;
  int                 retval;
  D4EST_VTK_FLOAT_TYPE *float_data;
  p4est_locidx_t     *ntc;

  P4EST_ASSERT (cont != NULL && cont->writing);
  Npoints = cont->num_points;
  ntc = cont->node_to_corner;
  use_nodes = 0;

  /* write point data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, scalar_name, d4est_vtk_get_format_string(output_type));

  double* vtk_array = d4est_vtk_convert_nodal_to_vtk
                      (
                       cont->p4est,
                       cont,
                       cont->d4est_ops,
                       values,
                       grid_type
                      );

  
  /* #ifdef D4EST_VTK_ASCII */
  if (output_type == D4EST_VTK_ASCII){
  for (il = 0; il < Npoints; ++il) {
    ddl = use_nodes ? ntc[il] : il;
    // P4EST_ASSERT (0 <= ddl && ddl < Ncorners);

    fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
             "     %24.16e\n",
#else
             "          %16.8e\n",
#endif
             vtk_array[il]);
  }
  }
  else {
  float_data = P4EST_ALLOC (D4EST_VTK_FLOAT_TYPE, Npoints);
  for (il = 0; il < Npoints; ++il) {
    ddl = use_nodes ? ntc[il] : il;
    // P4EST_ASSERT (0 <= ddl && ddl < Ncorners);
    float_data[il] =
      (D4EST_VTK_FLOAT_TYPE) vtk_array[il];
  }

  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = d4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * Npoints, output_type);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (float_data);

  if (retval) {
    d4est_vtk_context_destroy (cont);
    return NULL;
  }
  }

  
  P4EST_FREE(vtk_array);
  
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point scalar\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}


static d4est_vtk_context_t *
d4est_vtk_write_data
(
 d4est_vtk_context_t * cont,
 int num_point_scalars,
 const char** names,
 double** values,
 d4est_vtk_grid_type_t grid_type,
 d4est_vtk_output_type_t output_type
)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_context_destroy");
#endif
  const int           num_point_all = num_point_scalars;
  int                 mpirank;
  int                 retval;
  int                 i, all;
  int                 scalar_strlen;
  char                point_scalars[BUFSIZ];
  const char         *name;
  d4est_vtk_context_t *list_end;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (cont->p4est != NULL);

  /* This function needs to do nothing if there is no data. */
  if (!num_point_scalars) {
    return cont;
  }
  mpirank = cont->p4est->mpirank;


  /* Gather point data. */
  all = 0;
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    name = names[all];
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
  }


  fprintf (cont->vtufile, "      <PointData");
  fprintf (cont->vtufile, " Scalars=\"%s\"", point_scalars);
  fprintf (cont->vtufile, ">\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* now we are done counting and checking we write the fields */
  all = 0;
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    cont = d4est_vtk_write_point_scalar (cont, names[all], values[all], grid_type, output_type);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }

  fprintf (cont->vtufile, "      </PointData>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PPointData>\n");

    all = 0;
    for (i = 0; i < num_point_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], d4est_vtk_get_format_string(output_type));

    fprintf (cont->pvtufile, "    </PPointData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  return cont;
}


static d4est_vtk_context_t *
d4est_vtk_write_element_data
(
 d4est_vtk_context_t * cont,
 int write_tree,
 int write_level,
 int write_rank,
 int wrap_rank,
 int write_deg,
 int num_element_float_fields,
 const char ** float_names,
 double ** float_values,
 int num_element_int_fields,
 const char ** int_names,
 int ** int_values,
 d4est_vtk_output_type_t output_type
)
{
#ifdef D4EST_VTK_DEBUG
  zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  zlog_debug(c_default, "Starting d4est_vtk_write_cell_datav");
#endif
  /* This function needs to do nothing if there is no data. */
  if (!
      (write_tree || write_level || write_rank || wrap_rank || write_deg
       || num_element_int_fields || num_element_float_fields))
    return cont;

  const int           mpirank = cont->p4est->mpirank;
  int                 retval;
  int                 i, all = 0;
  int                 scalar_strlen, vector_strlen;
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  const p4est_locidx_t Ncells = cont->num_cells;
  char                cell_scalars[BUFSIZ], cell_vectors[BUFSIZ];
  const char         *name;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
  p4est_locidx_t      sk;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
  p4est_topidx_t      jt;
  p4est_locidx_t      il;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (wrap_rank >= 0);

  /* Gather cell data. */
  scalar_strlen = 0;
  cell_scalars[0] = '\0';
  for (i = 0; i < num_element_float_fields; ++all, ++i) {
    name = float_names[all];
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;
  }

  all = 0;
  for (i = 0; i < num_element_int_fields; ++all, ++i) {
    name = int_names[all];
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;
  }



  
  char                vtkCellDataString[BUFSIZ] = "";
  int                 printed = 0;

  if (write_tree)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed, "treeid");

  if (write_level)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",level" : "level");

  if (write_deg)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",degree" : "degree");
  

  if (write_rank)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",mpirank" : "mpirank");

  if (num_element_float_fields + num_element_int_fields)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_scalars);

  fprintf (cont->vtufile, "      <CellData Scalars=\"%s\">\n",
           vtkCellDataString);

  if (output_type != D4EST_VTK_ASCII){
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  }

  if (write_tree) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

    if (output_type == D4EST_VTK_ASCII){
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          fprintf (cont->vtufile, " %lld", (long long) jt);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    }
    else {
      
      /* for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) { */
    /*   tree = p4est_tree_array_index (trees, jt); */
    /*   num_quads = tree->quadrants.elem_count; */
    /*   for (zz = 0; zz < num_quads; ++zz, ++il) { */
    /*     int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM)); */
    /*     for (int ec = 0; ec < num_cells_in_element; ec++){ */
    /*     locidx_data[il] = (p4est_locidx_t) jt; */
    /*     } */
    /*   } */
    /* } */
    /* fprintf (cont->vtufile, "          "); */
    /* retval = d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data, */
    /*                                  sizeof (*locidx_data) * Ncells, output_type); */
    /* fprintf (cont->vtufile, "\n"); */
    /* if (retval) { */
    /*   P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n"); */
    /*   d4est_vtk_context_destroy (cont); */

    /*   /\* P4EST_FREE (values); *\/ */
    /*   /\* P4EST_FREE (names); *\/ */
    /*   P4EST_FREE (locidx_data); */
    /*   P4EST_FREE (uint8_data); */

    /*   return NULL; */
    /* } */
    }
    
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_level) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", d4est_vtk_get_format_string(output_type));

    if (output_type == D4EST_VTK_ASCII){
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          quad = p4est_quadrant_array_index (quadrants, zz);
          fprintf (cont->vtufile, " %d", (int) quad->level);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    }
    else {
    /* for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) { */
      /* tree = p4est_tree_array_index (trees, jt); */
    /*   quadrants = &tree->quadrants; */
    /*   num_quads = quadrants->elem_count; */
    /*   for (zz = 0; zz < num_quads; ++zz, ++il) { */
    /*     quad = p4est_quadrant_array_index (quadrants, zz); */

    /*             int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM)); */
    /*     for (int ec = 0; ec < num_cells_in_element; ec++){ */
    /*     uint8_data[il] = (uint8_t) quad->level; */
    /*     } */
    /*   } */
    /* } */

    /* fprintf (cont->vtufile, "          "); */
    /* retval = d4est_vtk_write_binary (cont->vtufile, (char *) uint8_data, */
    /*                                  sizeof (*uint8_data) * Ncells, output_type); */
    /* fprintf (cont->vtufile, "\n"); */
    /* if (retval) { */
    /*   P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n"); */
    /*   d4est_vtk_context_destroy (cont); */

    /*   P4EST_FREE (values); */
    /*   P4EST_FREE (names); */
    /*   P4EST_FREE (locidx_data); */

    /*   return NULL; */
    /* } */

    }
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_deg) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"degree\""
             " format=\"%s\">\n", "UInt8", d4est_vtk_get_format_string(output_type));
    if (output_type == D4EST_VTK_ASCII){
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          fprintf (cont->vtufile, " %d", (int) cont->deg_array[il]);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    }
    else {
    /* for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) { */
    /*   tree = p4est_tree_array_index (trees, jt); */
    /*   quadrants = &tree->quadrants; */
    /*   num_quads = quadrants->elem_count; */
    /*   for (zz = 0; zz < num_quads; ++zz, ++il) { */
    /*     int num_cells_in_element = d4est_util_int_pow_int(cont->deg_array[il], (P4EST_DIM)); */
    /*     for (int ec = 0; ec < num_cells_in_element; ec++){ */
    /*     uint8_data[il] = (uint8_t) cont->deg_array[il]; */
    /*     } */
    /*   } */
    /* } */

    /* fprintf (cont->vtufile, "          "); */
    /* retval = d4est_vtk_write_binary (cont->vtufile, (char *) uint8_data, */
    /*                                  sizeof (*uint8_data) * Ncells, output_type); */
    /* fprintf (cont->vtufile, "\n"); */

    /* P4EST_FREE (uint8_data); */

    /* if (retval) { */
    /*   P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n"); */
    /*   d4est_vtk_context_destroy (cont); */

    /*   P4EST_FREE (values); */
    /*   P4EST_FREE (names); */
    /*   P4EST_FREE (locidx_data); */

    /*   return NULL; */
    /* } */
    }
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

    if (output_type == D4EST_VTK_ASCII){
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
    }
    else {
    /* for (il = 0; il < Ncells; ++il) */
    /*   locidx_data[il] = (p4est_locidx_t) wrapped_rank; */

    /* fprintf (cont->vtufile, "          "); */
    /* retval = d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data, */
    /*                                  sizeof (*locidx_data) * Ncells, output_type); */
    /* fprintf (cont->vtufile, "\n"); */

    /* P4EST_FREE (locidx_data); */

    /* if (retval) { */
    /*   d4est_vtk_context_destroy (cont); */

    /*   P4EST_FREE (values); */
    /*   P4EST_FREE (names); */

    /*   return NULL; */
    /* } */
    }
    
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  all = 0;
  for (i = 0; i < num_element_float_fields; ++all, ++i) {
    cont = d4est_vtk_write_cell_scalar (cont, float_names[all], float_values[all], D4EST_VTK_FLOAT, output_type);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }

  all = 0;
  for (i = 0; i < num_element_int_fields; ++all, ++i) {
    cont = d4est_vtk_write_cell_scalar (cont, int_names[all], int_values[all], D4EST_VTK_INT, output_type);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }
  

  
  fprintf (cont->vtufile, "      </CellData>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PCellData Scalars=\"%s\">\n",
             vtkCellDataString);

    if (write_tree)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

    if (write_level)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", d4est_vtk_get_format_string(output_type));
    if (write_deg)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"deg\" format=\"%s\"/>\n",
               "UInt8", d4est_vtk_get_format_string(output_type));
    

    if (write_rank)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, d4est_vtk_get_format_string(output_type));

    all = 0;
    for (i = 0; i < num_element_float_fields; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, float_names[all], d4est_vtk_get_format_string(output_type));

    all = 0;
    for (i = 0; i < num_element_int_fields; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               "Int32", int_names[all], d4est_vtk_get_format_string(output_type));

    
    fprintf (cont->pvtufile, "    </PCellData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);

      return NULL;
    }
  }

  return cont;
}

void
d4est_vtk_save_aux
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 const char* input_file,
 const char* input_section,
 const char ** dg_field_names,
 double ** dg_fields,
 const char ** element_float_field_names,
 double ** element_float_fields,
  const char ** element_int_field_names,
 int ** element_int_fields,
 const char* folder,
 int sub_folder_number
){
 zlog_category_t *c_default = zlog_get_category("d4est_vtk");
  if (p4est->mpirank == 0)
    zlog_debug(c_default, "Saving mesh data to VTK file...");

  d4est_vtk_context_t* cont = d4est_vtk_dg_context_new(p4est, d4est_ops, input_file, input_section);

  cont->folder = d4est_util_add_cwd(folder);
  d4est_util_make_directory(cont->folder,0);
  if (sub_folder_number >= 0){
    asprintf(&cont->folder,"%s%d/", cont->folder, sub_folder_number);
  }
  d4est_util_make_directory(cont->folder,0);
  
  zlog_category_t *c_vtk_geom = zlog_get_category("d4est_vtk_geometry");
  d4est_geometry_t* geom_vtk = d4est_geometry_new
                               (
                                p4est->mpirank,
                                input_file,
                                cont->geometry_section,
                                c_vtk_geom
                               );
  
  d4est_vtk_context_set_geom(cont, geom_vtk);
  d4est_vtk_context_set_scale(cont, .99);

  int* deg_array = NULL;
  if(cont->grid_type == D4EST_VTK_DG_GRID){
    deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    d4est_mesh_get_array_of_degrees(p4est, (void*)deg_array, D4EST_INT);
  }
  else if (cont->grid_type == D4EST_VTK_CORNER_GRID){
    deg_array = P4EST_ALLOC(int, p4est->local_num_quadrants);
    for (int i = 0; i < p4est->local_num_quadrants; i++) deg_array[i] = 1;
  }
  else {
    D4EST_ABORT("Not an accepted grid type");
  }

  int num_dg_fields = 0;
  int num_element_float_fields = 0;
  int num_element_int_fields = 0;

  if (dg_field_names != NULL && dg_fields != NULL){
    for (int i = 0; dg_field_names[i] != NULL; i++){
      num_dg_fields++;
    }
  }

  if (element_float_field_names != NULL && element_float_fields != NULL){
    for (int i = 0; element_float_field_names[i] != NULL; i++){
      num_element_float_fields++;
    }
  }


  if (element_int_field_names != NULL && element_int_fields != NULL){
    for (int i = 0; element_int_field_names[i] != NULL; i++){
      num_element_int_fields++;
    }
  }
  
  d4est_vtk_context_set_deg_array(cont, deg_array);
  cont = d4est_vtk_write_header(cont, d4est_ops, cont->output_type);

  cont = d4est_vtk_write_data(cont, num_dg_fields, dg_field_names, dg_fields, cont->grid_type, cont->output_type);
  cont = d4est_vtk_write_element_data(cont, cont->write_tree, cont->write_level, cont->write_rank, cont->wrap_rank, cont->write_deg, num_element_float_fields, element_float_field_names, element_float_fields, num_element_int_fields, element_int_field_names, element_int_fields, cont->output_type);
    
  d4est_vtk_write_footer(cont);
  d4est_geometry_destroy(geom_vtk);
  P4EST_FREE(deg_array);
  
  if (p4est->mpirank == 0)
    zlog_debug(c_default, "VTK file saved.");

}

void
d4est_vtk_save
(
 p4est_t* p4est,
 d4est_operators_t* d4est_ops,
 const char* input_file,
 const char* input_section,
 const char ** dg_field_names,
 double ** dg_fields,
 const char ** element_float_field_names,
 double ** element_float_fields,
 const char ** element_int_field_names,
 int ** element_int_fields,
 int sub_folder_number
)
{
  d4est_vtk_save_aux
    (
     p4est,
     d4est_ops,
     input_file,
     input_section,
     dg_field_names,
     dg_fields,
     element_float_field_names,
     element_float_fields,
     element_int_field_names,
     element_int_fields,
     "VTK",
     sub_folder_number
    );
}
