/*
  This file is part of d4est.
  d4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  d4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  d4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with d4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file d4est_vtk.h
 *
 * Routines for printing a forest and associated fields to VTK format.
 *
 * \ingroup d4est
 */

#ifndef D4EST_VTK_H
#define D4EST_VTK_H

#include <dgmath.h>
#include <pXest.h>
#include <d4est_geometry.h>
/* SC_EXTERN_C_BEGIN; */

/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct d4est_vtk_context d4est_vtk_context_t;

/** Write the d4est in VTK format.
 *
 * This is a convenience function for the special case of writing out
 * the tree id, quadrant level, and MPI rank only.
 * One file is written per MPI rank, and one meta file on rank 0.
 * The quadrants are scaled to length .95; see \ref d4est_vtk_write_header.
 * This function will abort if there is a file error.
 *
 * \param [in] d4est    The d4est to be written.
 * \param [in] geom     A p4est_geometry_t structure or NULL for vertex space
 *                      as defined by d4est->connectivity.
 * \param [in] filename The first part of the file name which will have the
 *                      MPI rank appended to it: The output file will be
 *                      filename_rank.vtu, and the meta file filename.pvtu.
 */
void                d4est_vtk_write_file (p4est_t * d4est,
                                          d4est_geometry_t * geom,
                                          const char *filename);

/** The first call to write a VTK file using individual functions.
 *
 * Writing a VTK file is split into multiple functions that keep a context.
 * This is the first function that allocates the opaque context structure.
 * After allocation, further parameters can be set for the context.
 * Then, the header, possible data fields, and the footer must be written.
 * The process can be aborted any time by destroying the context.  In this
 * case, open files are closed cleanly with only partially written content.
 *
 * \param d4est     The d4est to be written.
 *                  If no geometry is specified in
 *                  \ref d4est_vtk_context_set_geom, we require
 *                  \b d4est->connectivity to have valid vertex arrays.
 * \param filename  The first part of the name which will have the processor
 *                  number appended to it (i.e., the output file will be
 *                  filename_rank.vtu).  The parallel meta-files for Paraview
 *                  and Visit use this basename too.
 *                  We copy this filename to internal storage, so it is not
 *                  needed to remain alive after calling this function.
 * \return          A VTK context fur further use.
 */
d4est_vtk_context_t *
d4est_vtk_context_new (p4est_t * p4est,
                       const char *filename);

/** Modify the geometry transformation registered in the context.
 * After \ref d4est_vtk_context_new, it is at the default NULL.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref d4est_vtk_write_header.
 * \param geom      A \ref p4est_geometry_t structure, or NULL for vertex space.
 *                  If NULL, \b d4est->connectivity->vertices and
 *                  \b tree_to_vertex must be non-NULL.
 */
void                d4est_vtk_context_set_geom (d4est_vtk_context_t * cont,
                                                d4est_geometry_t * geom);

/** Modify the context parameter for scaling the quadrants.
 * After \ref d4est_vtk_context_new, it is at the default 0.95.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref d4est_vtk_write_header.
 * \param [in] scale            Scale parameter must be in (0, 1].
 */
void                d4est_vtk_context_set_scale (d4est_vtk_context_t * cont,
                                                 double scale);

/** Modify the context parameter for expecting continuous point data.
 * If set to true, the point data is understood as a continuous field.
 * In this case, we can significantly reduce the file size when scale == 1.
 * For discontinuous point data, it should be set to false.
 * After \ref d4est_vtk_context_new, it is at the default false.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref d4est_vtk_write_header.
 * \param [in] continuous       Boolean parameter.
 */
void                d4est_vtk_context_set_continuous (d4est_vtk_context_t *
                                                      cont, int continuous);

/** Cleanly destroy a \ref d4est_vtk_context_t structure.
 *
 * This function closes all the file pointers and frees the context.
 * It can be called even if the VTK output
 * has only been partially written, the files' content will be incomplete.
 *
 * \param[in] context     The VTK file context to be destroyed.
 */
void                d4est_vtk_context_destroy (d4est_vtk_context_t * context);

/** Write the VTK header.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 *     vtk_context = d4est_vtk_context_new (d4est, "output");
 *     d4est_vtk_context_set_* (vtk_context, parameter);
 *     vtk_context = d4est_vtk_write_header (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     vtk_context = d4est_vtk_write_cell_data (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     vtk_context = d4est_vtk_write_point_data (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     retval = d4est_vtk_write_footer (vtk_context);
 *     if (retval) { error; }
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 *                         None of the vtk_write functions must have been called.
 *                         This context is the return value if no error occurs.
 *
 * \return          On success, an opaque context (d4est_vtk_context_t) pointer
 *                  that must be passed to subsequent d4est_vtk calls.  It is
 *                  required to call \ref d4est_vtk_write_footer eventually with
 *                  this value.  Returns NULL on error.
 */
d4est_vtk_context_t *d4est_vtk_write_header (d4est_vtk_context_t * cont);

/** Write VTK cell data.
 *
 * There are options to have this function write
 * the tree id, quadrant level, or MPI rank without explicit input data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues), followed
 * by a final argument of the VTK context cont (same as the first argument).
 * The cell scalar pairs come first, followed by the cell vector pairs, then cont.
 * Each 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables.  The number of
 * doubles in each sc_array must be exactly \a d4est->local_num_quadrants for
 * scalar data and \a 3*d4est->local_num_quadrants for vector data.
 *
 * \note The current d4est_vtk_context_t structure, \a cont, must be the first
 * and the last argument
 * of any call to this function; this argument is used to validate that the
 * correct number of variable arguments have been provided.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
d4est_vtk_context_t *d4est_vtk_write_cell_dataf (d4est_vtk_context_t * cont,
                                                 int write_tree,
                                                 int write_level,
                                                 int write_rank,
                                                 int wrap_rank,
                                                 int num_cell_scalars,
                                                 int num_cell_vectors, ...);

/** This is an alternate version of the varargs function.
 * Works exactly the same otherwise.
 * TODO: implement, also for vectors and point data.
 */
d4est_vtk_context_t *d4est_vtk_write_cell_data (d4est_vtk_context_t * cont,
                                                int write_tree,
                                                int write_level,
                                                int write_rank,
                                                int wrap_rank,
                                                int num_cell_scalars,
                                                int num_cell_vectors,
                                                const char *filenames[],
                                                sc_array_t * values[]);

/** Write VTK point data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont    A VTK context created by \ref d4est_vtk_context_new.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues) where
 * the point scalar pairs come first, followed by the point vector pairs.  Each
 * 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables. The number of
 * doubles in each sc_array must be exactly the number of components (1 for
 * scalar and 3 for vector) times 4 times number of elements.
 *
 * \note The current
 * d4est_vtk_context_t structure, cont, must be the last argument of any call
 * to this function; this argument is used to validate that the correct number
 * of variable arguments have been provided.
 *
 * \note The number of point scalar data in each
 * sc_array must be exactly \a D4EST_CHILDREN*local_num_quadrants, and the
 * number of point vector data must be exactly \a
 * 3*D4EST_CHILDREN*local_num_quadrants. I.e. there must be data for every
 * corner of every quadrant in the \a d4est, even if the corner is shared by
 * multiple quadrants.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
d4est_vtk_context_t *d4est_vtk_write_point_dataf (d4est_vtk_context_t * cont,
                                                  int num_point_scalars,
                                                  int num_point_vectors, ...);

/** Write the VTU footer and clean up.
 *
 * Writing a VTK file is split into a few routines.
 * This function writes the footer information to the VTK file and cleanly
 * destroys the VTK context.
 *
 * \param [in] cont Context is deallocated before the function returns.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 d4est_vtk_write_footer (d4est_vtk_context_t * cont);

d4est_vtk_context_t *
d4est_vtk_write_dg_header (d4est_vtk_context_t * cont, dgmath_jit_dbase_t* dgmath_jit_dbase);
d4est_vtk_context_t *
d4est_vtk_write_dg_cell_dataf (d4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                               int write_rank, int wrap_rank, int write_deg,
                               int num_cell_scalars, int num_cell_vectors, ...);

void
d4est_vtk_context_set_deg_array (d4est_vtk_context_t * cont,
                                 int* deg_array);

d4est_vtk_context_t *
d4est_vtk_write_dg_point_dataf (d4est_vtk_context_t * cont,
                                int num_point_scalars, int num_point_vectors,
                                ...);
d4est_vtk_context_t *
d4est_vtk_dg_context_new (p4est_t * p4est, dgmath_jit_dbase_t* dgbase, const char *filename);

#endif /* !D4EST_VTK_H */
