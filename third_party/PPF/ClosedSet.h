/*==============================================================================
 *
 * Name:        ClosedSet.h
 *
 * Description:
 *	Declarations for functions that maintain a closed set (a set
 * 	with a fixed maximum number of members).
 *
 * Traceability:
 *      Version         Author          Date            Description
 *      -------         ------          ----            -----------
 *      0.1             johnmay         05/03/95        Initial version
 *	0.2		johnmay		05/04/95	Added IsFull & IsEmpty
 *	1.0		johnmay		06/06/95	Initial stable version
 *	2.0		linda		11/10/95	initial phase 2
 *	2.1		linda		04/04/96	fixed FreeSet
 *	2.1x		johnmay		10/01/96	add macro to compute
 *							size to allow externally
 *							allocated memory;
 *							replace MakeSet with
 *							InitSet
 *
 *  Notes:
 *	 Copyright (c) 1995, 1996 The Regents of the Univerity of California
 *			and Lawrence Livermore National Laboratory
 *       All rights reserved.
 *
 *----------------------------------------------------------------------------*/

#ifndef _CLOSED_SET_H
#define _CLOSED_SET_H

typedef char *	ClosedSet;

void		InitSet( int length, ClosedSet mem );
void		AddSet( ClosedSet theSet, int position );
int		CheckSet( ClosedSet theSet, int position );
void		DropSet( ClosedSet theSet, int position );
void		ClearSet( ClosedSet theSet );
void		UnionSet( ClosedSet source, ClosedSet dest );
int		IsFull( ClosedSet theSet );
int		IsEmpty( ClosedSet theSet );

#define		LOGBITS	3
#define		CLOSED_SET_BYTES(entries) \
			(((entries - 1) >> LOGBITS) + 1 + sizeof(int))
#endif
