/* PPF_Printf.c */
/* John May -- 23 Oct 96 */
/* Fortran wrapper for parallel print function */ 
/* Copyright (c) 1996-1998 The Regents of the University of California
 *                      and Lawrence Livermore National Laboratory
 * All rights reserved.
 */

/* Modification history */
/*
 *	 6 Nov 97	johnmay Added portability features for different MPIs
 *				and Fortran naming conventions
 *	13 Jan 98	johnmay	Change names to follow Ptools conventions
 */

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <ptools_ppf.h>

/* This is were all the Fortran to C poratbility problems live.
 * Three issues:  what do Fortran external names look like?
 * How are strings passed?  How do we convert Fortran communicators
 * to C communicators?
 * First issue was addressed by MPICH, and we've adopted their
 * solution: figure it out with a configure script and use macros.
 * Second issue seems to work OK with linkage below on all
 * systems tested so far, though I've heard Crays do things
 * differently.  I don't have one to try out.
 * Third issue is toughest.  MPI-2 will add MPI_Comm_f2c, which
 * does the conversion we want.  If that's not available, we
 * have two choices (so far): assume no conversion is needed,
 * or convert using MPICH's MPIR_ToPointer function.  This is
 * only needed on MPICH-based systems that have 64-bit pointers.
 * (DEC's MPI is an example).  So configure script tests for
 * these characteristics and sets flags as appropriate.
 */

/* Tested and works on:
	IBM MPI (PSSP 2.2 and 2.3)
	MPICH 1.1 on IBM 
	DEC native and MPICH 1.1 (Digital Unix 4.0B)
		-DEC native requires MPIR_ToPointer and FORTRANUNDERSCORE
		-DEC MPICH 1.1 sets FORTRANUNDERSCORE automatically in mpicc
	Meiko MPICH 1.0.13 uses FORTRANUNDERSCORE, pointers are 32-bit
 */

#ifndef PPF_USE_OLD_NAMES
# ifdef FORTRANCAPS
#  define ppf_print_ PPF_PRINT
# elif defined(FORTRANDOUBLEUNDERSCORE)
#  define ppf_print_ ppf_print__
# elif !defined(FORTRANUNDERSCORE)
#  define ppf_print_ ppf_print
# endif /* FORTRANCAPS */
#else /* old names... */
# ifdef FORTRANCAPS
#  define ppf_print_ PTC_PRINT
# elif defined(FORTRANDOUBLEUNDERSCORE)
#  define ppf_print_ ptc_print__
# elif !defined(FORTRANUNDERSCORE)
#  define ppf_print_ ptc_print
# else
#  define ppf_print_ ptc_print_
# endif/* FORTRANCAPS */
/* and this for the function call to PPF_Print... */
# define PPF_Print PTC_Print
#endif /*PPF_USE_OLD_NAMES */


/* Complicated tricks for converting communicators:  if we have
 * the MPI2 function, use it.  If we don't, then we can do a
 * direct copy if either 1) we're using MPICH 1.1 or later or
 * 2) we're using an older MPICH and pointers are 32 bits.  If
 * this is old MPICH and pointers are big, we call MPIR_ToPointer.
 * If we don't have that, we just guess that a direct copy is OK
 * because we have no other choice.
 */
#ifdef HAVE_MPI_COMM_F2C
#define FPTR_TO_COMM( pcomm ) (MPI_Comm_f2c(*pcomm))
typedef MPI_Comm * FCOMM;
#else /* No MPI_Comm_f2c */
#if !(defined(MPI_VERSION) && MPI_VERSION == 1 && MPI_SUBVERSION == 1) && \
	defined(POINTER_64_BITS) && defined(HAVE_MPIR_TOPOINTER)
/* we'll try a pointer conversion */
#define FPTR_TO_COMM( pcomm ) ( (MPI_Comm)MPIR_ToPointer(*(int*)(pcomm)) )
typedef MPI_Comm FCOMM;
#else /* try direct copy */
#define FPTR_TO_COMM( pcomm ) ( *pcomm )
typedef MPI_Comm * FCOMM;
#endif /* do pointer conversion */
#endif /* HAVE_MPI_COMM_F2C */


void ppf_print_( comm, string, fstatus, strlength )
	FCOMM comm;
	char * string;
	int * fstatus;
	int strlength;
{

	int size;
	int end = (strlength > 254) ? 254 : strlength;
	char strcopy[256];
	MPI_Comm c_comm;

	c_comm = FPTR_TO_COMM( comm );
#ifdef DEBUG
	if( MPI_Comm_size( c_comm, &size ) != MPI_SUCCESS ) {
		printf("Failed to convert communicator %d from Fotran to C!\n",
			comm);
	}
#endif

	/* Treat a one-space string as a NULL */
	if( strlength == 1 && string[0] == ' ' ) {
		*fstatus = PPF_Print( c_comm, NULL );
	} else {
		strncpy( strcopy, string, end );
		strcopy[end] = '\n';
		strcopy[end + 1] = 0;
		*fstatus = PPF_Print( c_comm, strcopy );
	}
	
}
