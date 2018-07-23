/* PPF_Print.c */
/* John May -- 23 Oct 96 */
/* Prototype parallel print function */
/* 
 *	 Copyright (c) 1996-1998 The Regents of the Univerity of California
 *			and Lawrence Livermore National Laboratory
 *       All rights reserved.
 */

/* Modification history */
/*
 *	31 Oct 96	johnmay	Accept NULL as input string and don't produce
 *				any output (even node number) for that node.
 *	13 Nov 96	johnmay Changed values of NODESUB and PERCENTSUB so
 *				they work more portably.
 *	 6 Nov 97	johnmay Fixed reduce to order merged strings correctly
 *	13 Jan 98	johnmay	Change names to follow Ptools conventions
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "ClosedSet.h"

#ifndef PPF_USE_OLD_NAMES
#include "ptools_ppf.h"
#else
#include "PTCPrint.h"
#define PPF_Print PTC_Print
#endif

typedef char LineData;

#define TEXT 0
#define SET PPF_MAXLINE		/* PPF_MAXLINE is defined in ptools_ppf.h */

#define MAXNODELIST 0x10000	/* good for at least 10000 nodes */
#define NODECHAR 'N'
#define NODESUB '\01'
#define PERCENTSUB '\02'

void PrintResult( LineData * dest_lines, int comm_size, size_t linesize );

void PrintReduceOp(	void *		invec,
			void *		inoutvec,
			int *		len,
			MPI_Datatype *	datatype )
{
	LineData *	src_data =
				(LineData *)((char*)invec + 2 * sizeof(int));
	LineData *	dest_data =
				(LineData *)((char*)inoutvec + 2 * sizeof(int));
	LineData *	psrcline, * pdestline;
	int		i, j;
	int		count = ((int *)invec)[0];
	int		linesize = ((int *)invec)[1];

	datatype = datatype;
	len = len;

	/* Run through each (valid) entry in dest_data and compare it with
	 * each valid entry in src_data.
	 */

	for( pdestline = dest_data, i = 0; i < count;
			pdestline += linesize, i++ ) {

		/* An empty set indicates that the string is blank; skip it */
		if( IsEmpty( pdestline+SET ) ) continue;

		/* We have a valid string; compare it with each of the
		 * strings in the incoming vector.
		 */
		for( psrcline = src_data, j = 0; j < count;
			psrcline += linesize, j++ ) {

			/* Skip blank entries */
			if( IsEmpty( psrcline+SET ) ) continue;

			/* If they match, remove this string (because it
			 * won't match any other entries), and mark the
			 * nodes where it appeared.  Once we've found a
			 * match to this dest string, there won't be any
			 * more in the source array because it should
			 * contain only one copy of each string, so we
			 * can quit this inner loop after the first match.
			 */
			if( strcmp( pdestline+TEXT, psrcline+TEXT ) == 0 ) {
				UnionSet( pdestline+SET, psrcline+SET );
				ClearSet( psrcline+SET );
				/* We want to put the result in the
				 * lowest-numbered position in the destination
				 * array so items will come out in the order
				 * of the lowest-numbered node in the set when
				 * we print the sets out.
				 */
				if( j < i ) {
					/* Move to lower-numbered position */
					LineData * lower
						 = dest_data + (j * linesize);
					UnionSet( lower+SET, pdestline+SET );
					ClearSet( pdestline+SET );
					strcpy( lower+TEXT, pdestline+TEXT );
				}
				break;
			} 
		}
	}

	/* All the matching strings should now be copied into the source.
	 * Those remaining are unique, so we can copy them all into the
	 * dest array.  Since each node sends only one string, there's no
	 * danger of overwriting an entry in dest when we copy an entry
	 * from source.
	 */

	for( pdestline = dest_data, psrcline = src_data, i = 0; i < count;
			pdestline += linesize,
			psrcline += linesize, i++ ) {

		/* If we find a string whose entry hasn't been cleared,
		 * copy it and its list of nodes.
		 */
		if( !IsEmpty( psrcline+SET ) ) {
			UnionSet( pdestline+SET, psrcline+SET );
			strcpy( pdestline+TEXT, psrcline+TEXT );
		}
	}

}

int PPF_Print( MPI_Comm comm, const char * string, ... )
{
	va_list		args;

	static int	inited;
	static MPI_Op	PTC_Print_reduce;
	static size_t	linesize, blocksize;
	static MPI_Datatype	type;
	static char *	src_block;
	static char *	dest_block;
	static int	rank, size;

	int 		i;
	int		localrank;
	int 		status;
	LineData *	src_lines;
	LineData *	dest_lines;
	LineData *	pline;
	char *		ptoken, * plast;
	char		temp[PPF_MAXLINE];

	if( !inited ) {
		status = MPI_Op_create( PrintReduceOp, 1, &PTC_Print_reduce );
		if( status != MPI_SUCCESS ) {
			fprintf( stderr, "Couldn't create Print Reduce op %d\n",
				status );
			return status;
		}

		/* Use comm world so that output is identified in terms of
		 * global rank.  (Should this be optional?)
		 */
		MPI_Comm_rank( MPI_COMM_WORLD, &rank );
		MPI_Comm_size( MPI_COMM_WORLD, &size );

		/* Get enough memory to hold size line entries + 
		 * count field + and linesize field
		 */
		linesize = CLOSED_SET_BYTES( size ) + PPF_MAXLINE;
		blocksize = linesize * size + 2 * sizeof(int);

		src_block = malloc( blocksize );
		dest_block = malloc( blocksize );

		/* Store the number of entries and linesize (since set size
		 * varies from run to run)
		 */
		((int*)src_block)[0] = size;
		((int*)src_block)[1] = linesize;

		MPI_Type_contiguous( blocksize, MPI_BYTE, &type );
		MPI_Type_commit( &type );

		inited = 1;
	}

	/* Initialize the line entries */
	src_lines = (LineData *)(src_block + 2 * sizeof(int));
	dest_lines = (LineData *)(dest_block + 2 * sizeof(int));
	
#if 0
	/* Truncate the string if necessary */
	if( strlen( string ) >= PPF_MAXLINE ) {
		string[PPF_MAXLINE - 1] = '\0';
	}
#endif

	/* Search for the pattern that means "put the node numbers here".
	 * Since it will begin with at %, we have to convert it to
	 * something else so sprintf doesn't try to interpret it.  When
	 * we finally do the printout, we'll replace the placeholder
	 * with the node list.  Likewise, since we're processing
	 * this string twice through printf-like functions, any %%
	 * patterns the user put in will get eaten up.  To prevent
	 * that, we change %% to something else before the first
	 * pass through sprintf and change it back before the second.
	 */

	if( string != NULL ) {	/* Don't try formatting an empty string */
		/* First copy into a buffer that we can modify */
		strncpy( temp, string, PPF_MAXLINE );
		temp[PPF_MAXLINE-1] = '\0';

		for( ptoken = strchr( temp, '%' );
				ptoken != NULL;
				ptoken = strchr( ptoken + 2, '%' ) ) {
			if( ptoken[1] == NODECHAR ) {
				*ptoken = PERCENTSUB;
				ptoken[1] = NODESUB;
			} else if( ptoken[1] == '%' ) {
				*ptoken = ptoken[1] = PERCENTSUB;
			}
		}
	}
		
	/* Initialize all the strings and sets in our list to be
	 * empty, except the one for this node.  That gets
	 * initialized with the formatted string that the user
	 * requested and the nodeset consisting of this node.
	 */
	va_start( args, string );
	for( pline = src_lines, i = 0; i < size;
		pline += linesize, i++ ) {
		InitSet( size, pline+SET );
		if( i == rank ) {
			if( string != NULL ) {
				vsprintf( pline+TEXT, temp, args );
				AddSet( pline+SET, rank );
			}
			/* If string is NULL, this node won't appear in
			 * the set for any string entry, so it won't
			 * produce any output.
			 */
		}
	}
	va_end( args );


	MPI_Reduce( src_block, dest_block, 1, type, PTC_Print_reduce, 0, comm );

	/* For this test, we use to input comm rank so that we're sure
	 * "node 0" is in this group.
	 */
	MPI_Comm_rank( comm, &localrank );
	if( localrank == 0 ) {
		PrintResult( dest_lines, size, linesize );
	}

	return MPI_SUCCESS;
}

#define INVALID -1
#define STARTRANGE( start, end, node ) if(start==INVALID) start=j; end=j;
#define ENDRANGE( start, end, firsttime ) \
	if(!firsttime) {sprintf(pnext,",%n",&numc); pnext+=numc;} \
	else firsttime=0; \
	if(start==end) sprintf(pnext,"%d%n",start,&numc); \
	else sprintf(pnext,"%d-%d%n",start,end,&numc); \
	pnext+=numc; \
	start = INVALID;

void PrintResult( LineData * dest_lines, int comm_size, size_t linesize )
{
	LineData *	pline;
	char	nodelist[MAXNODELIST];
	char *	pnext;
	char *	ptoken;
	int	numc;
	int	i;

	/* Look for lines that survived the merge and print them*/
	for( pline = dest_lines, i = 0; i < comm_size;
			pline += linesize, i++ ) {
		if( !IsEmpty( pline+SET ) ) {
			/* Create a compact list of nodes for this string */
			int j, start = INVALID, end, firsttime = 1;
			int hasnodepattern = 0;
			pnext = nodelist;

			for( j = 0; j < comm_size; j++ ) {
				if( CheckSet( pline+SET, j ) ) {
					STARTRANGE(start,end,j);
				} else if ( start != INVALID ) {
					ENDRANGE(start,end,firsttime);
				}
			}
			if( start != INVALID ) {
				ENDRANGE(start,end,firsttime);
			}

			/* Reconvert the node list substitution and the
			 * literal percent patterns.
			 */
			for( ptoken = strchr( pline+TEXT, PERCENTSUB );
				ptoken != NULL;
				ptoken = strchr( ptoken + 2, PERCENTSUB ) ) {

				if( ptoken[1] == NODESUB ) {
					/* Only allow one copy of the pattern;
					 * otherwise, make it a percent.
					 * (because we can't pass variable
					 * number of copies of the nodelist).
					 */
					*ptoken = '%';
					if( !hasnodepattern ) ptoken[1] = 's';
					else ptoken[1] = '%';
					hasnodepattern = 1;
				} else if( ptoken[1] == PERCENTSUB ) {
					*ptoken = ptoken[1] = '%';
				}
			}
	
			if( !hasnodepattern ) {
				/* Use separate printfs so the output line
				 * will be processed by print once in either
				 * case.
				 */
				printf( "%s ", nodelist );
				printf( pline+TEXT );
			} else {
				printf( pline+TEXT, nodelist );
			}
		}
	}
}
