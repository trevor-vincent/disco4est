#ifndef PXEST_H
#define PXEST_H 

#include <global_macros.h>

#ifndef DIM
#error "DIM (dimension) not defined"
#elif (DIM>3||DIM<2)
#error "DIM > 3 or DIM < 2: Dimension not supported"
#elif DIM==3
#include <p4est_to_p8est.h>
#endif

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_nodes.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_nodes.h>
#endif

#endif
