/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/


#include <pXest.h>
#include "d4est_vtk.h"
#include "d4est_geometry.h"
#include <util.h>
#ifdef P4_TO_P8
/* #include <p8est_vtk.h> */
/* #include <p8est_nodes.h> */
#define D4EST_VTK_CELL_TYPE     11      /* VTK_VOXEL */
#else

/* #include <p4est_nodes.h> */
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

/** Write a cell scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref d4est_vtk_write_cell_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
d4est_vtk_context_t *d4est_vtk_write_cell_scalar (d4est_vtk_context_t *
                                                  cont,
                                                  const char *scalar_name,
                                                  sc_array_t * values);

/** Write a 3-vector cell field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref d4est_vtk_write_cell_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
d4est_vtk_context_t *d4est_vtk_write_cell_vector (d4est_vtk_context_t *
                                                  cont,
                                                  const char *vector_name,
                                                  sc_array_t * values);

/** Write a point scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref d4est_vtk_write_point_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
d4est_vtk_context_t *d4est_vtk_write_point_scalar (d4est_vtk_context_t *
                                                   cont,
                                                   const char *scalar_name,
                                                   sc_array_t * values);

/** Write a 3-vector point field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref d4est_vtk_write_point_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
d4est_vtk_context_t *d4est_vtk_write_point_vector (d4est_vtk_context_t *
                                                   cont,
                                                   const char *vector_name,
                                                   sc_array_t * values);

/* #ifdef P4_TO_P8 */
/* #define d4est_vtk_context               p8est_vtk_context */
/* #endif */

#ifndef D4EST_VTK_DOUBLES
#define D4EST_VTK_FLOAT_NAME "Float32"
#define D4EST_VTK_FLOAT_TYPE float
#else
#define D4EST_VTK_FLOAT_NAME "Float64"
#define D4EST_VTK_FLOAT_TYPE double
#endif

#ifndef D4EST_VTK_BINARY
#define D4EST_VTK_ASCII 1
#define D4EST_VTK_FORMAT_STRING "ascii"
#else
#define D4EST_VTK_FORMAT_STRING "binary"

static int
d4est_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                        size_t byte_length)
{
  printf("[D4EST_VTK]: Starting d4est_vtk_write_binary \n");
#ifndef D4EST_VTK_COMPRESSION
  return sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
#else
  return sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
#endif /* D4EST_VTK_COMPRESSION */
}

#endif /* D4EST_VTK_BINARY */

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
  dgmath_jit_dbase_t* dgbase;
  
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
};

d4est_vtk_context_t *
d4est_vtk_context_new (p4est_t * p4est, const char *filename)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_new \n");
#endif
  d4est_vtk_context_t *cont;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);

  /* Allocate, initialize the vtk context.  Important to zero all fields. */
  cont = P4EST_ALLOC_ZERO (d4est_vtk_context_t, 1);

  cont->p4est = p4est;
  cont->filename = P4EST_STRDUP (filename);
  cont->deg_array = NULL;
  
  cont->scale = d4est_vtk_scale;
  cont->continuous = d4est_vtk_continuous;

  return cont;
}

d4est_vtk_context_t *
d4est_vtk_dg_context_new (p4est_t * p4est, dgmath_jit_dbase_t* dgbase, const char *filename)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_new \n");
#endif
  d4est_vtk_context_t *cont;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);


  /* Allocate, initialize the vtk context.  Important to zero all fields. */
  cont = P4EST_ALLOC_ZERO (d4est_vtk_context_t, 1);
  cont->dgbase = dgbase;
  cont->p4est = p4est;
  cont->filename = P4EST_STRDUP (filename);
  cont->deg_array = NULL;
  
  cont->scale = d4est_vtk_scale;
  cont->continuous = d4est_vtk_continuous;

  return cont;
}


void
d4est_vtk_context_set_geom (d4est_vtk_context_t * cont,
                            d4est_geometry_t * geom)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_set_geom \n");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->geom = geom;
}

void
d4est_vtk_context_set_deg_array (d4est_vtk_context_t * cont,
                                 int* deg_array)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_set_geom \n");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->deg_array = deg_array;
}

void
d4est_vtk_context_set_scale (d4est_vtk_context_t * cont, double scale)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_set_scale \n");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);
  P4EST_ASSERT (0. < scale && scale <= 1.);

  cont->scale = scale;
}

void
d4est_vtk_context_set_continuous (d4est_vtk_context_t * cont, int continuous)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_set_continuous \n");
#endif
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->continuous = continuous;
}

void
d4est_vtk_context_destroy (d4est_vtk_context_t * context)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_destroy \n");
#endif
  P4EST_ASSERT (context != NULL);
  P4EST_ASSERT (context->p4est != NULL);

  /* since this function is called inside write_header and write_footer,
   * we cannot assume a consistent state of all member variables */

  P4EST_ASSERT (context->filename != NULL);
  P4EST_FREE (context->filename);

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

void
d4est_vtk_write_file (p4est_t * p4est, d4est_geometry_t * geom,
                      const char *filename)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_file \n");
#endif
  int                 retval;
  d4est_vtk_context_t *cont;

  /* allocate context and set parameters */
  cont = d4est_vtk_context_new (p4est, filename);
  d4est_vtk_context_set_geom (cont, geom);

  /* We do not write point data, so it is safe to set continuous to true.
   * This will not save any space though since the default scale is < 1. */
  d4est_vtk_context_set_continuous (cont, 1);

  /* write header, that is, vertex positions and quadrant-to-vertex map */
  cont = d4est_vtk_write_header (cont);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing header");

  /* write the tree/level/rank data */
  cont =
    d4est_vtk_write_cell_dataf (cont, d4est_vtk_write_tree,
                                d4est_vtk_write_level, d4est_vtk_write_rank,
                                d4est_vtk_wrap_rank, 0, 0, cont);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing cell data");

  /* properly write rest of the files' contents */
  retval = d4est_vtk_write_footer (cont);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

d4est_vtk_context_t *
d4est_vtk_write_header (d4est_vtk_context_t * cont)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_header \n");
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
#ifdef D4EST_VTK_ASCII
  double              wx, wy, wz;
#else
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
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
  P4EST_ASSERT (filename != NULL);

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
  Ncells = p4est->local_num_quadrants;

  cont->num_corners = Ncorners = P4EST_CHILDREN * Ncells;

  printf("cont->num_corners = %d\n",cont->num_corners);
  printf("Ncells = %d\n\n",Ncells);
  
  if (scale < 1. || !conti) {
    /* when we scale the quadrants we need each corner separately */
    cont->nodes = nodes = NULL;
    cont->num_points = Npoints = Ncorners;
    cont->node_to_corner = ntc = NULL;
    indeps = NULL;
  }
  else {
    /* if scale == 1. and the point data is continuous,
     * we can reuse shared quadrant corners */
    cont->nodes = nodes = p4est_nodes_new (p4est, NULL);
    indeps = &nodes->indep_nodes;
    cont->num_points = Npoints = nodes->num_owned_indeps;
    P4EST_ASSERT ((size_t) Npoints == indeps->elem_count);

    /* Establish a reverse lookup table from a node to its first reference.
     * It is slow to run twice through memory like this.  However, we also know
     * that writing data to disk is slower still, so we do not optimize.
     */
    cont->node_to_corner = ntc = P4EST_ALLOC (p4est_locidx_t, Npoints);
    memset (ntc, -1, Npoints * sizeof (p4est_locidx_t));
    for (sk = 0, il = 0; il < Ncells; ++il) {
      for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
        ntcid = nodes->local_nodes[sk];
        P4EST_ASSERT (0 <= ntcid && ntcid < Npoints);
        if (ntc[ntcid] < 0) {
          ntc[ntcid] = sk;
        }
      }
    }
#ifdef P4EST_ENABLE_DEBUG
    /* the particular version of nodes we call makes sure they are tight */
    for (ntcid = 0; ntcid < Npoints; ++ntcid) {
      P4EST_ASSERT (0 <= ntc[ntcid] && ntc[ntcid] < Ncorners);
    }
#endif
  }

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
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
#if defined D4EST_VTK_BINARY && defined D4EST_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
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
           D4EST_VTK_FLOAT_NAME, D4EST_VTK_FORMAT_STRING);

  if (nodes == NULL) {
    /* loop over the trees */
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
        quad = p4est_quadrant_array_index (quadrants, zz);
        h2 = .5 * intsize * P4EST_QUADRANT_LEN (quad->level);
        k = 0;
#ifdef P4_TO_P8
        for (zi = 0; zi < 2; ++zi) {
          eta_z = intsize * quad->z + h2 * (1. + (zi * 2 - 1) * scale);
#endif
          for (yi = 0; yi < 2; ++yi) {
            eta_y = intsize * quad->y + h2 * (1. + (yi * 2 - 1) * scale);
            for (xi = 0; xi < 2; ++xi) {
              P4EST_ASSERT (0 <= k && k < P4EST_CHILDREN);
              eta_x = intsize * quad->x + h2 * (1. + (xi * 2 - 1) * scale);
              if (geom != NULL) {
                xyz[0] = eta_x;
                xyz[1] = eta_y;
                xyz[2] = eta_z;
                geom->X (geom, jt, (p4est_qcoord_t [(P4EST_DIM)]){0}, -(P4EST_ROOT_LEN), xyz, COORDS_TREE_UNITCUBE, XYZ);
                for (j = 0; j < 3; ++j) {
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
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
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
                    (D4EST_VTK_FLOAT_TYPE) xyz[j];
                }
              }
              ++k;
            }
          }
#ifdef P4_TO_P8
        }
#endif
        P4EST_ASSERT (k == P4EST_CHILDREN);
      }
    }
    P4EST_ASSERT (P4EST_CHILDREN * quad_count == Npoints);
  }
 

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

  
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", nodes == NULL ?
               (long long) sk : (long long) nodes->local_nodes[sk]);
    }
    fprintf (cont->vtufile, "\n");
  }

  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");

  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", D4EST_VTK_FORMAT_STRING);
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", D4EST_VTK_CELL_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined D4EST_VTK_BINARY && defined D4EST_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             D4EST_VTK_FLOAT_NAME, D4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in d4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
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

/** Write VTK point data.
 *
 * This function exports custom point data to the vtk file; it is functionally
 * the same as \b d4est_vtk_write_point_dataf with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b d4est_vtk_write_point_dataf
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by \ref d4est_vtk_context_new.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static d4est_vtk_context_t *
d4est_vtk_write_point_datav (d4est_vtk_context_t * cont,
                             int num_point_scalars,
                             int num_point_vectors, va_list ap)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_destroy \n");
#endif
  const int           num_point_all = num_point_scalars + num_point_vectors;
  int                 mpirank;
  int                 retval;
  int                 i, all;
  int                 scalar_strlen, vector_strlen;
  char                point_scalars[BUFSIZ], point_vectors[BUFSIZ];
  const char         *name, **names;
  d4est_vtk_context_t *list_end;
  sc_array_t        **values;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (cont->p4est != NULL);

  /* This function needs to do nothing if there is no data. */
  if (!(num_point_scalars || num_point_vectors)) {
    return cont;
  }
  mpirank = cont->p4est->mpirank;

  /* Allocate storage to manage the data fields. */
  values = P4EST_ALLOC (sc_array_t *, num_point_all);
  names = P4EST_ALLOC (const char *, num_point_all);

  /* Gather point data. */
  all = 0;
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data type;"
                    " scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == (size_t) cont->num_corners,
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data count; see "
                    P4EST_STRING "_vtk.h for more details.");
  }

  /* keep variable all at current value */
  vector_strlen = 0;
  point_vectors[0] = '\0';
  for (i = 0; i < num_point_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data type;"
                    " vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == 3 * (size_t) cont->num_corners,
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data count; see "
                    P4EST_STRING "_vtk.h for more details.");
  }

  /* Check for pointer variable marking the end of variable data input. */
  list_end = va_arg (ap, d4est_vtk_context_t *);
  SC_CHECK_ABORT (list_end == cont,
                  P4EST_STRING "_vtk Error: the end of variable data must be"
                  " specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t pointer.  See " P4EST_STRING
                  "_vtk.h for more information.");

  fprintf (cont->vtufile, "      <PointData");
  fprintf (cont->vtufile, " Scalars=\"%s\"", point_scalars);
  fprintf (cont->vtufile, " Vectors=\"%s\"", point_vectors);
  fprintf (cont->vtufile, ">\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  /* now we are done counting and checking we write the fields */
  all = 0;
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    cont = d4est_vtk_write_point_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }

  for (i = 0; i < num_point_vectors; ++all, ++i) {
    cont = d4est_vtk_write_point_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point vectors");
  }

  fprintf (cont->vtufile, "      </PointData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PPointData>\n");

    all = 0;
    for (i = 0; i < num_point_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_point_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PPointData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}



d4est_vtk_context_t *
d4est_vtk_write_dg_header (d4est_vtk_context_t * cont, dgmath_jit_dbase_t* dgmath_jit_dbase)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_header \n");
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
#ifdef D4EST_VTK_ASCII
  double              wx, wy, wz;
#else
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
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
    Ncells += util_int_pow_int(deg_array[i], (P4EST_DIM));
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
    mpi_abort("[D4EST_ERROR]: We do not support scale = 1. yet");
  }

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
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
#if defined D4EST_VTK_BINARY && defined D4EST_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
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
           D4EST_VTK_FLOAT_NAME, D4EST_VTK_FORMAT_STRING);

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


        double* vtk_rst = dgmath_fetch_vtk_rst
                          (
                           dgmath_jit_dbase,
                           deg_array[quad_count],
                           (P4EST_DIM)
                          );


        quad = p4est_quadrant_array_index (quadrants, zz);
        h2 = .5 * intsize * P4EST_QUADRANT_LEN (quad->level);
        k = 0;

        int num_cells_in_element = util_int_pow_int(deg_array[quad_count], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++) {
          for (int corn = 0; corn < (P4EST_CHILDREN); corn++){

            double r_01 = dgmath_rtox(vtk_rst[0 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
            double s_01 = dgmath_rtox(vtk_rst[1 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
            double t_01 = dgmath_rtox(vtk_rst[2 + corn*3 + ec*3*(P4EST_CHILDREN)], 0., 1.);
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
    /* P4EST_ASSERT (k == P4EST_CHILDREN); */
    /* } */
  /*   } */
  /*   P4EST_ASSERT (P4EST_CHILDREN * quad_count == Npoints); */
  /* } */

  
  /* else { */
  /*   mpi_abort("[D4EST_ABORT]: We do not support nodes != NULL in vtk"); */
  /* } */
 

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

  
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", nodes == NULL ?
               (long long) sk : (long long) nodes->local_nodes[sk]);
    }
    fprintf (cont->vtufile, "\n");
  }

  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");

  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", D4EST_VTK_FORMAT_STRING);
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", D4EST_VTK_CELL_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined D4EST_VTK_BINARY && defined D4EST_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             D4EST_VTK_FLOAT_NAME, D4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");



    fprintf (cont->pvtufile, "    <PCells>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"connectivity\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "Int32", D4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"offsets\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "Int32", D4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"types\""
             " NumberOfComponents=\"1\" format=\"%s\"/>\n",
             "UInt8", D4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PCells>\n");
    
    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in d4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
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

d4est_vtk_context_t *
d4est_vtk_write_point_dataf (d4est_vtk_context_t * cont,
                             int num_point_scalars, int num_point_vectors,
                             ...)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_destroy \n");
#endif
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_point_scalars >= 0 && num_point_vectors >= 0);

  va_start (ap, num_point_vectors);
  cont = d4est_vtk_write_point_datav (cont, num_point_scalars,
                                      num_point_vectors, ap);
  va_end (ap);

  return cont;
}

/** Write VTK cell data.
 *
 * This function exports custom cell data to the vtk file; it is functionally
 * the same as \b d4est_vtk_write_cell_dataf with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b d4est_vtk_write_cell_dataf
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by \ref d4est_vtk_context_new.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static d4est_vtk_context_t *
d4est_vtk_write_cell_datav (d4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                            int num_cell_scalars,
                            int num_cell_vectors, va_list ap)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_cell_datav \n");
#endif
  /* This function needs to do nothing if there is no data. */
  if (!
      (write_tree || write_level || write_rank || wrap_rank
       || num_cell_vectors || num_cell_vectors))
    return cont;

  const int           mpirank = cont->p4est->mpirank;
  int                 retval;
  int                 i, all = 0;
  int                 scalar_strlen, vector_strlen;
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  char                cell_scalars[BUFSIZ], cell_vectors[BUFSIZ];
  const char         *name, **names;
  sc_array_t        **values;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
#ifdef D4EST_VTK_ASCII
  p4est_locidx_t      sk;
#else
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  p4est_topidx_t      jt;
  p4est_locidx_t      il;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (wrap_rank >= 0);

  values = P4EST_ALLOC (sc_array_t *, num_cell_scalars + num_cell_vectors);
  names = P4EST_ALLOC (const char *, num_cell_scalars + num_cell_vectors);

  /* Gather cell data. */
  scalar_strlen = 0;
  cell_scalars[0] = '\0';
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data type; scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    (size_t) cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data count; scalar data must contain exactly p4est->local_num_quadrants doubles.");
  }

  vector_strlen = 0;
  cell_vectors[0] = '\0';
  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data type; vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    3 * (size_t) cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data count; vector data must contain exactly 3*p4est->local_num_quadrants doubles.");
  }

  /* Check for pointer variable marking the end of variable data input. */
  d4est_vtk_context_t *end = va_arg (ap, d4est_vtk_context_t *);
  SC_CHECK_ABORT (end == cont, P4EST_STRING "_vtk Error: the end of variable "
                  "data must be specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t struct. See " P4EST_STRING
                  "_vtk.h for more information.");

  char                vtkCellDataString[BUFSIZ] = "";
  int                 printed = 0;

  if (write_tree)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed, "treeid");

  if (write_level)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",level" : "level");

  if (write_rank)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",mpirank" : "mpirank");

  if (num_cell_scalars)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_scalars);

  if (num_cell_vectors)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_vectors);

  fprintf (cont->vtufile, "      <CellData Scalars=\"%s\">\n",
           vtkCellDataString);

#ifndef D4EST_VTK_ASCII
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
#endif

  if (write_tree) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
#ifdef D4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        fprintf (cont->vtufile, " %lld", (long long) jt);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        locidx_data[il] = (p4est_locidx_t) jt;
      }
    }
    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (cont->vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);
      P4EST_FREE (uint8_data);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
    P4EST_ASSERT (il == Ncells);
  }

  if (write_level) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", D4EST_VTK_FORMAT_STRING);
#ifdef D4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        fprintf (cont->vtufile, " %d", (int) quad->level);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        uint8_data[il] = (uint8_t) quad->level;
      }
    }

    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                                     sizeof (*uint8_data) * Ncells);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (uint8_data);

    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
#ifdef D4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0; il < Ncells; ++il)
      locidx_data[il] = (p4est_locidx_t) wrapped_rank;

    fprintf (cont->vtufile, "          ");
    retval = d4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (locidx_data);

    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  all = 0;
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    cont = d4est_vtk_write_cell_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }

  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    cont = d4est_vtk_write_cell_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell vectors");
  }

  fprintf (cont->vtufile, "      </CellData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PCellData Scalars=\"%s\">\n",
             vtkCellDataString);

    if (write_tree)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

    if (write_level)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", D4EST_VTK_FORMAT_STRING);

    if (write_rank)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

    all = 0;
    for (i = 0; i < num_cell_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_cell_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PCellData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}

d4est_vtk_context_t *
d4est_vtk_write_cell_dataf (d4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                            int num_cell_scalars, int num_cell_vectors, ...)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_t *end = va_arg \n");
#endif
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_cell_scalars >= 0 && num_cell_vectors >= 0);

  va_start (ap, num_cell_vectors);
  cont = d4est_vtk_write_cell_datav (cont,
                                     write_tree, write_level,
                                     write_rank, wrap_rank,
                                     num_cell_scalars, num_cell_vectors, ap);
  va_end (ap);

  return cont;
}

d4est_vtk_context_t *
d4est_vtk_write_point_scalar (d4est_vtk_context_t * cont,
                              const char *scalar_name, sc_array_t * values)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_point_scalar \n");
#endif
  p4est_locidx_t      il, ddl;
  int                 use_nodes;
#ifdef P4EST_ENABLE_DEBUG
  int                 Ncorners;
#endif
  int                 Npoints;
#ifndef D4EST_VTK_ASCII
  int                 retval;
  D4EST_VTK_FLOAT_TYPE *float_data;
#endif
  p4est_locidx_t     *ntc;

  P4EST_ASSERT (cont != NULL && cont->writing);
#ifdef P4EST_ENABLE_DEBUG
  Ncorners = cont->num_corners;
#endif
  Npoints = cont->num_points;
  ntc = cont->node_to_corner;
  P4EST_ASSERT (values != NULL && values->elem_count == (size_t) Ncorners);
  if (ntc == NULL) {
    /* we are writing a discontinuous field, possibly due to vertex scaling */
    P4EST_ASSERT (cont->num_corners == cont->num_points);
    P4EST_ASSERT (cont->scale < 1. || !cont->continuous);
    use_nodes = 0;
  }
  else {
    /* we are definitely writing a continuous field, reusing corner values */
    P4EST_ASSERT (cont->scale == 1. && cont->continuous);
    use_nodes = 1;
  }

  /* write point data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, scalar_name, D4EST_VTK_FORMAT_STRING);

#ifdef D4EST_VTK_ASCII
  for (il = 0; il < Npoints; ++il) {
    ddl = use_nodes ? ntc[il] : il;
    P4EST_ASSERT (0 <= ddl && ddl < Ncorners);

    fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
             "     %24.16e\n",
#else
             "          %16.8e\n",
#endif
             *(double *) sc_array_index (values, ddl));
  }
#else
  float_data = P4EST_ALLOC (D4EST_VTK_FLOAT_TYPE, Npoints);
  for (il = 0; il < Npoints; ++il) {
    ddl = use_nodes ? ntc[il] : il;
    P4EST_ASSERT (0 <= ddl && ddl < Ncorners);
    float_data[il] =
      (D4EST_VTK_FLOAT_TYPE) * ((double *) sc_array_index (values, ddl));
  }

  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = d4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * Npoints);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (float_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point scalar\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}

d4est_vtk_context_t *
d4est_vtk_write_point_vector (d4est_vtk_context_t * cont,
                              const char *vector_name, sc_array_t * values)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_binary (cont->vtufile, \n");
#endif
  P4EST_ASSERT (cont != NULL && cont->writing);

  SC_ABORT (P4EST_STRING "_vtk_write_point_vector not implemented");
}

d4est_vtk_context_t *
d4est_vtk_write_cell_scalar (d4est_vtk_context_t * cont,
                             const char *scalar_name, sc_array_t * values)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_cell_scalar \n");
#endif
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  p4est_locidx_t      il;
#ifndef D4EST_VTK_ASCII
  int                 retval;
  D4EST_VTK_FLOAT_TYPE *float_data;
#endif

  P4EST_ASSERT (cont != NULL && cont->writing);

  /* Write cell data. */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, scalar_name, D4EST_VTK_FORMAT_STRING);

#ifdef D4EST_VTK_ASCII
  for (il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
             "     %24.16e\n",
#else
             "          %16.8e\n",
#endif
             *(double *) sc_array_index (values, il));
  }
#else
  float_data = P4EST_ALLOC (D4EST_VTK_FLOAT_TYPE, Ncells);
  for (il = 0; il < Ncells; ++il) {
    float_data[il] =
      (D4EST_VTK_FLOAT_TYPE) * ((double *) sc_array_index (values, il));
  }

  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = d4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * Ncells);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (float_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding scalar cell data\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell scalar file\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}


/* d4est_vtk_context_t * */
/* d4est_vtk_write_point_vector (d4est_vtk_context_t * cont, */
/*                               const char *vector_name, sc_array_t * values) */
/* { */
/* #ifdef D4EST_VTK_DEBUG */
/*   printf("[D4EST_VTK]: Starting d4est_vtk_write_binary (cont->vtufile, \n"); */
/* #endif */
/*   P4EST_ASSERT (cont != NULL && cont->writing); */

/*   SC_ABORT (P4EST_STRING "_vtk_write_point_vector not implemented"); */
/* } */

d4est_vtk_context_t *
d4est_vtk_write_dg_cell_scalar (d4est_vtk_context_t * cont,
                             const char *scalar_name, double* values)
{
  P4EST_ASSERT (cont != NULL && cont->writing);

  /* Write cell data. */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, scalar_name, D4EST_VTK_FORMAT_STRING);

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
  
  for (tc = 0, il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
    tree = p4est_tree_array_index (trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;
    for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
      int num_cells_in_element = util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
      for (int ec = 0; ec < num_cells_in_element; ec++, tc++){
        fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
                 "     %24.16e\n",
#else
                 "          %16.8e\n",
#endif
                 values[il]);
      }
    }
  }

  /* printf("tc = %d, Ncells = %d\n", tc, Ncells); */
  
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell scalar file\n");
    d4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}


d4est_vtk_context_t *
d4est_vtk_write_cell_vector (d4est_vtk_context_t * cont,
                             const char *vector_name, sc_array_t * values)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_binary (cont->vtufile, \n");
#endif
  P4EST_ASSERT (cont != NULL && cont->writing);

  SC_ABORT (P4EST_STRING "_vtk_write_cell_vector not implemented");
}

int
d4est_vtk_write_footer (d4est_vtk_context_t * cont)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_footer \n");
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


static d4est_vtk_context_t *
d4est_vtk_write_dg_cell_datav (d4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                               int write_deg,   
                            int num_cell_scalars,
                            int num_cell_vectors, va_list ap)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_cell_datav \n");
#endif
  /* This function needs to do nothing if there is no data. */
  if (!
      (write_tree || write_level || write_rank || wrap_rank || write_deg
       || num_cell_vectors || num_cell_vectors))
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
  const char         *name, **names;
  double        **values;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
#ifdef D4EST_VTK_ASCII
  p4est_locidx_t      sk;
#else
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  p4est_topidx_t      jt;
  p4est_locidx_t      il;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (wrap_rank >= 0);

  values = P4EST_ALLOC (double *, num_cell_scalars + num_cell_vectors);
  names = P4EST_ALLOC (const char *, num_cell_scalars + num_cell_vectors);

  /* Gather cell data. */
  scalar_strlen = 0;
  cell_scalars[0] = '\0';
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, double *);

    /* Validate input. */
    /* SC_CHECK_ABORT (values[all]->elem_size == sizeof (double), */
    /*                 P4EST_STRING */
    /*                 "_vtk: Error: incorrect cell scalar data type; scalar data must contain doubles."); */
    /* SC_CHECK_ABORT (values[all]->elem_count == */
    /*                 (size_t) cont->p4est->local_num_quadrants, */
    /*                 P4EST_STRING */
    /*                 "_vtk: Error: incorrect cell scalar data count; scalar data must contain exactly p4est->local_num_quadrants doubles."); */
  }

  vector_strlen = 0;
  cell_vectors[0] = '\0';
  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, double *);

  }

  /* Check for pointer variable marking the end of variable data input. */
  d4est_vtk_context_t *end = va_arg (ap, d4est_vtk_context_t *);
  SC_CHECK_ABORT (end == cont, P4EST_STRING "_vtk Error: the end of variable "
                  "data must be specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t struct. See " P4EST_STRING
                  "_vtk.h for more information.");

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

  if (num_cell_scalars)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_scalars);

  if (num_cell_vectors)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_vectors);

  fprintf (cont->vtufile, "      <CellData Scalars=\"%s\">\n",
           vtkCellDataString);

#ifndef D4EST_VTK_ASCII
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
#endif

  if (write_tree) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          fprintf (cont->vtufile, " %lld", (long long) jt);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_level) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", D4EST_VTK_FORMAT_STRING);
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          quad = p4est_quadrant_array_index (quadrants, zz);
          fprintf (cont->vtufile, " %d", (int) quad->level);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_deg) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"degree\""
             " format=\"%s\">\n", "UInt8", D4EST_VTK_FORMAT_STRING);
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        int num_cells_in_element = util_int_pow_int(cont->deg_array[il], (P4EST_DIM));
        for (int ec = 0; ec < num_cells_in_element; ec++){
          quad = p4est_quadrant_array_index (quadrants, zz);
          fprintf (cont->vtufile, " %d", (int) cont->deg_array[il]);
          if (!(sk % 20) && il != (Ncells - 1))
            fprintf (cont->vtufile, "\n         ");
        }
      }
    }
    fprintf (cont->vtufile, "\n");
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  all = 0;
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    cont = d4est_vtk_write_dg_cell_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }

  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    /* cont = d4est_vtk_write_cell_vector (cont, names[all], values[all]); */
    /* SC_CHECK_ABORT (cont != NULL, */
                    /* P4EST_STRING "_vtk: Error writing cell vectors"); */
  }

  fprintf (cont->vtufile, "      </CellData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PCellData Scalars=\"%s\">\n",
             vtkCellDataString);

    if (write_tree)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

    if (write_level)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", D4EST_VTK_FORMAT_STRING);
    if (write_deg)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"deg\" format=\"%s\"/>\n",
               "UInt8", D4EST_VTK_FORMAT_STRING);
    

    if (write_rank)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, D4EST_VTK_FORMAT_STRING);

    all = 0;
    for (i = 0; i < num_cell_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_cell_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PCellData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}


d4est_vtk_context_t *
d4est_vtk_write_dg_cell_dataf (d4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                               int write_rank, int wrap_rank, int write_deg,
                            int num_cell_scalars, int num_cell_vectors, ...)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_t *end = va_arg \n");
#endif
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_cell_scalars >= 0 && num_cell_vectors >= 0);

  va_start (ap, num_cell_vectors);
  cont = d4est_vtk_write_dg_cell_datav(cont,
                                       write_tree,
                                       write_level,
                                       write_rank,
                                       wrap_rank,
                                       write_deg,
                                       num_cell_scalars,
                                       num_cell_vectors,
                                       ap);
  va_end (ap);

  return cont;
}

double*
d4est_vtk_convert_nodal_to_vtk(p4est_t* p4est,
                               d4est_vtk_context_t* cont,
                               dgmath_jit_dbase_t* dgbase,
                               double* values)
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
      int num_points_in_element = util_int_pow_int(cont->deg_array[il], (P4EST_DIM))*(P4EST_CHILDREN);
      int num_nodes_in_element = util_int_pow_int(cont->deg_array[il] + 1, (P4EST_DIM));
      dgmath_convert_nodal_to_vtk(
                                  dgbase,
                                  &values[nodal_stride],
                                  (P4EST_DIM),
                                  cont->deg_array[il],
                                  &vtk_array[vtk_stride]
                                 );
      vtk_stride += num_points_in_element;
      nodal_stride += num_nodes_in_element;
    }
  }  

  return vtk_array;
}



d4est_vtk_context_t *
d4est_vtk_write_dg_point_scalar (d4est_vtk_context_t * cont,
                              const char *scalar_name, double * values)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_write_point_scalar \n");
#endif
  p4est_locidx_t      il, ddl;
  int                 use_nodes;
#ifdef P4EST_ENABLE_DEBUG
  int                 Ncorners;
#endif
  int                 Npoints;
#ifndef D4EST_VTK_ASCII
  int                 retval;
  D4EST_VTK_FLOAT_TYPE *float_data;
#endif
  p4est_locidx_t     *ntc;

  P4EST_ASSERT (cont != NULL && cont->writing);
#ifdef P4EST_ENABLE_DEBUG
  Ncorners = cont->num_corners;
#endif
  Npoints = cont->num_points;
  ntc = cont->node_to_corner;
  /* P4EST_ASSERT (values != NULL && values->elem_count == (size_t) Ncorners); */
  /* if (ntc == NULL) { */
    /* we are writing a discontinuous field, possibly due to vertex scaling */
  P4EST_ASSERT (cont->num_corners == cont->num_points);
  P4EST_ASSERT (cont->scale < 1. || !cont->continuous);
  use_nodes = 0;

    /* } */
  /* else { */
    /* we are definitely writing a continuous field, reusing corner values */
    /* P4EST_ASSERT (cont->scale == 1. && cont->continuous); */
    /* use_nodes = 1; */
  /* } */

  /* write point data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           D4EST_VTK_FLOAT_NAME, scalar_name, D4EST_VTK_FORMAT_STRING);

  double* vtk_array = d4est_vtk_convert_nodal_to_vtk
                      (
                       cont->p4est,
                       cont,
                       cont->dgbase,
                       values
                      );

  
/* #ifdef D4EST_VTK_ASCII */
  for (il = 0; il < Npoints; ++il) {
    ddl = use_nodes ? ntc[il] : il;
    P4EST_ASSERT (0 <= ddl && ddl < Ncorners);

    fprintf (cont->vtufile,
#ifdef D4EST_VTK_DOUBLES
             "     %24.16e\n",
#else
             "          %16.8e\n",
#endif
             vtk_array[il]);
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
d4est_vtk_write_dg_point_datav (d4est_vtk_context_t * cont,
                             int num_point_scalars,
                             int num_point_vectors, va_list ap)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_destroy \n");
#endif
  const int           num_point_all = num_point_scalars + num_point_vectors;
  int                 mpirank;
  int                 retval;
  int                 i, all;
  int                 scalar_strlen, vector_strlen;
  char                point_scalars[BUFSIZ], point_vectors[BUFSIZ];
  const char         *name, **names;
  d4est_vtk_context_t *list_end;
  double        **values;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (cont->p4est != NULL);

  /* This function needs to do nothing if there is no data. */
  if (!(num_point_scalars || num_point_vectors)) {
    return cont;
  }
  mpirank = cont->p4est->mpirank;

  /* Allocate storage to manage the data fields. */
  values = P4EST_ALLOC (double *, num_point_all);
  names = P4EST_ALLOC (const char *, num_point_all);

  /* Gather point data. */
  all = 0;
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, double *);

    /* /\* Validate input. *\/ */
    /* SC_CHECK_ABORT (values[all]->elem_size == sizeof (double), */
    /*                 P4EST_STRING */
    /*                 "_vtk: Error: incorrect point scalar data type;" */
    /*                 " scalar data must contain doubles."); */
    /* SC_CHECK_ABORT (values[all]->elem_count == (size_t) cont->num_corners, */
    /*                 P4EST_STRING */
    /*                 "_vtk: Error: incorrect point scalar data count; see " */
    /*                 P4EST_STRING "_vtk.h for more details."); */
  }

  /* keep variable all at current value */
  vector_strlen = 0;
  point_vectors[0] = '\0';
  for (i = 0; i < num_point_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, double *);

    /* Validate input. */
    /* SC_CHECK_ABORT (values[all]->elem_size == sizeof (double), */
                    /* P4EST_STRING */
                    /* "_vtk: Error: incorrect point vector data type;" */
                    /* " vector data must contain doubles."); */
    /* SC_CHECK_ABORT (values[all]->elem_count == 3 * (size_t) cont->num_corners, */
                    /* P4EST_STRING */
                    /* "_vtk: Error: incorrect point vector data count; see " */
                    /* P4EST_STRING "_vtk.h for more details."); */
  }

  /* Check for pointer variable marking the end of variable data input. */
  list_end = va_arg (ap, d4est_vtk_context_t *);
  SC_CHECK_ABORT (list_end == cont,
                  P4EST_STRING "_vtk Error: the end of variable data must be"
                  " specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t pointer.  See " P4EST_STRING
                  "_vtk.h for more information.");

  fprintf (cont->vtufile, "      <PointData");
  fprintf (cont->vtufile, " Scalars=\"%s\"", point_scalars);
  fprintf (cont->vtufile, " Vectors=\"%s\"", point_vectors);
  fprintf (cont->vtufile, ">\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  /* now we are done counting and checking we write the fields */
  all = 0;
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    cont = d4est_vtk_write_dg_point_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }

  for (i = 0; i < num_point_vectors; ++all, ++i) {
    mpi_abort("We do not support dg vectors yet");
    /* cont = d4est_vtk_write_dg_point_vector (cont, names[all], values[all]); */
    /* SC_CHECK_ABORT (cont != NULL, */
                    /* P4EST_STRING "_vtk: Error writing point vectors"); */
  }

  fprintf (cont->vtufile, "      </PointData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    d4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }



  

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PPointData>\n");

    all = 0;
    for (i = 0; i < num_point_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_point_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               D4EST_VTK_FLOAT_NAME, names[all], D4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PPointData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      d4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}


d4est_vtk_context_t *
d4est_vtk_write_dg_point_dataf (d4est_vtk_context_t * cont,
                                int num_point_scalars, int num_point_vectors,
                                ...)
{
#ifdef D4EST_VTK_DEBUG
  printf("[D4EST_VTK]: Starting d4est_vtk_context_destroy \n");
#endif
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_point_scalars >= 0 && num_point_vectors >= 0);

  va_start (ap, num_point_vectors);
  cont = d4est_vtk_write_dg_point_datav (cont, num_point_scalars,
                                      num_point_vectors, ap);
  va_end (ap);

  return cont;
}


void
d4est_vtk_save_geometry_and_dg_fields
(
 const char* save_as_filename,
 p4est_t* p4est,
 dgmath_jit_dbase_t* dgmath_jit_dbase,
 int* deg_array,
 const char* input_file,
 const char* input_section,
 d4est_vtk_user_fcn_t d4est_vtk_user_fcn,
 void* vtk_user
){


    d4est_geometry_t* geom_vtk = d4est_geometry_new
                                 (
                                  p4est->mpirank,
                                  input_file,
                                  input_section,
                                  "[D4EST_VTK_GEOMETRY]"
                                 );

    d4est_vtk_context_t* vtk_ctx = d4est_vtk_dg_context_new(p4est, dgmath_jit_dbase, save_as_filename);
    d4est_vtk_context_set_geom(vtk_ctx, geom_vtk);
    d4est_vtk_context_set_scale(vtk_ctx, .99);
    d4est_vtk_context_set_deg_array(vtk_ctx, deg_array);
    vtk_ctx = d4est_vtk_write_dg_header(vtk_ctx, dgmath_jit_dbase);    

    d4est_vtk_user_fcn(
                       vtk_ctx,
                       vtk_user
                      );

    d4est_vtk_write_footer(vtk_ctx);

    d4est_geometry_destroy(geom_vtk);
}
