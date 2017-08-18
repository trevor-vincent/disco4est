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

#include <d4est_operators.h>
#include <pXest.h>
#include <d4est_geometry.h>
/* SC_EXTERN_C_BEGIN; */



/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct d4est_vtk_context d4est_vtk_context_t;


typedef void
(*d4est_vtk_user_fcn_t)
(
 d4est_vtk_context_t*,
 void*
);


/* This file was automatically generated.  Do not edit! */
void d4est_vtk_save_geometry_and_cell_fields(const char *save_as_filename,p4est_t *p4est,d4est_operators_t *d4est_ops,const char *input_file,const char *input_section,d4est_vtk_user_fcn_t d4est_vtk_user_fcn,void *vtk_user);
void d4est_vtk_save_geometry_and_dg_fields(const char *save_as_filename,p4est_t *p4est,d4est_operators_t *d4est_ops,int *deg_array,const char *input_file,const char *input_section,d4est_vtk_user_fcn_t d4est_vtk_user_fcn,void *vtk_user);
d4est_vtk_context_t *d4est_vtk_write_dg_point_dataf(d4est_vtk_context_t *cont,int num_point_scalars,int num_point_vectors,...);
d4est_vtk_context_t *d4est_vtk_write_dg_point_scalar(d4est_vtk_context_t *cont,const char *scalar_name,double *values);
double *d4est_vtk_convert_nodal_to_vtk(p4est_t *p4est,d4est_vtk_context_t *cont,d4est_operators_t *d4est_ops,double *values);
d4est_vtk_context_t *d4est_vtk_write_dg_cell_dataf(d4est_vtk_context_t *cont,int write_tree,int write_level,int write_rank,int wrap_rank,int write_deg,int num_cell_scalars,int num_cell_vectors,...);
d4est_vtk_context_t *d4est_vtk_write_dg_cell_scalar(d4est_vtk_context_t *cont,const char *scalar_name,double *values);
d4est_vtk_context_t *d4est_vtk_write_point_dataf(d4est_vtk_context_t *cont,int num_point_scalars,int num_point_vectors,...);
d4est_vtk_context_t *d4est_vtk_write_dg_header(d4est_vtk_context_t *cont,d4est_operators_t *d4est_ops);
int d4est_vtk_write_footer(d4est_vtk_context_t *cont);
d4est_vtk_context_t *d4est_vtk_write_cell_dataf(d4est_vtk_context_t *cont,int write_tree,int write_level,int write_rank,int wrap_rank,int num_cell_scalars,int num_cell_vectors,...);
d4est_vtk_context_t *d4est_vtk_write_header(d4est_vtk_context_t *cont);
void d4est_vtk_write_file(p4est_t *p4est,d4est_geometry_t *geom,const char *filename);
void d4est_vtk_context_destroy(d4est_vtk_context_t *context);
void d4est_vtk_context_set_continuous(d4est_vtk_context_t *cont,int continuous);
void d4est_vtk_context_set_scale(d4est_vtk_context_t *cont,double scale);
void d4est_vtk_context_set_deg_array(d4est_vtk_context_t *cont,int *deg_array);
void d4est_vtk_context_set_geom(d4est_vtk_context_t *cont,d4est_geometry_t *geom);
d4est_vtk_context_t *d4est_vtk_dg_context_new(p4est_t *p4est,d4est_operators_t *d4est_ops,const char *filename);
d4est_vtk_context_t *d4est_vtk_context_new(p4est_t *p4est,const char *filename);
d4est_vtk_context_t *d4est_vtk_write_point_vector(d4est_vtk_context_t *cont,const char *vector_name,sc_array_t *values);
d4est_vtk_context_t *d4est_vtk_write_point_scalar(d4est_vtk_context_t *cont,const char *scalar_name,sc_array_t *values);
d4est_vtk_context_t *d4est_vtk_write_cell_vector(d4est_vtk_context_t *cont,const char *vector_name,sc_array_t *values);
d4est_vtk_context_t *d4est_vtk_write_cell_scalar(d4est_vtk_context_t *cont,const char *scalar_name,sc_array_t *values);


#endif /* !D4EST_VTK_H */
