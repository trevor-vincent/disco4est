/*==============================================================================
 *
 * Name:        ClosedSet.c
 *
 * Functions:
 *	InitSet
 *	AddSet
 *	CheckSet
 *	DropSet
 *	ClearSet
 *	UnionSet
 *	IsFull
 *	IsEmpty
 *
 * Description:
 *	Maintains a bit set of predetermined size; user can add
 *	and removed specific elements from the set, check for the
 *	presence of a given element, and check to see if a set is
 *	full or empty.
 *
 * Traceability:
 *      Version         Author          Date            Description
 *      -------         ------          ----            -----------
 *      0.1             johnmay         05/03/95        Initial code
 *	0.2		johnmay		05/04/95	Added full & empty check
 *	1.0		johnmay		06/06/95	Initial stable version
 *	1.1		linda		07/20/95	free => SIOF_FREE and
 *							malloc.h => stdlib.h
 *	1.2		johnmay		08/04/95	fix func name in comment
 *      2.0             linda           11/10/95        initial phase 2
 *	2.1		linda		04/04/96	fixed FreeSet
 *	2.1x		johnmay		10/01/96	removed SIOF dependence;
 *							modified for use with
 *							preallocated field
 *							(replaced MakeSet with
 *							InitSet).  Replaced
 *							FreeSet with ClearSet
 *							and added UnionSet
 *
 *  Notes:
 *	 Copyright (c) 1995, 1996 The Regents of the Univerity of California
 *			and Lawrence Livermore National Laboratory
 *       All rights reserved.
 *
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "ClosedSet.h"

/* Predefined masks for bits 0-7 of a byte */
static char mask[8] = {'\001','\002','\004','\010',
			'\020','\040','\100','\200'};

static int TestSet( ClosedSet theSet, int state );

#define	BITSPERBYTE 8
#define LOGBITS 3
#define	ENDMASK 7
#define ALLBITS 0xff

static int templen;
#define SETSIZE( set ) (memcpy(&templen,set,sizeof(int)),templen)

/*==============================================================================
 *
 * Function:	InitSet	- Initialize a set of a fixed maximum length
 *
 * Synopsis:
 *
 *	ClosedSet InitSet(
 *		int length )	** IN - maximum number of elements in set
 *
 * Description:
 *	This function initializes a set whose maximum size is given by the
 *	length parameter.  Potential elements are numbered 0 to length-1.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	Return Value:
 *		Object representing the set.
 *
 * Interfaces:
 *
 * Resources Used:
 *	None
 * 
 * Limitations:
 *
 * Assumptions:
 *	Assumes that memory passed in is propery allocated and large enough
 *	to hold the set (size should have been calculated with
 *	CLOSED_SET_BYTES macro).
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/

void InitSet( int length, ClosedSet mem )
{
	int 	bytes;

	if( length <= 0 ) return;

	bytes = CLOSED_SET_BYTES( length );

	/* Initialize set to be empty */
	memset( mem, 0, bytes );
	
	/* Store size at beginning of array */
	memcpy( mem, &length, sizeof(int) );

}

/*==============================================================================
 *
 * Function:	AddSet	- Insert an item in a set initialized by InitSet
 *
 * Synopsis:
 *
 *	void AddSet(
 *		ClosedSet	theSet,		** IN/OUT - Set to be changed
 *		int		position )	** IN - Element to be added
 *
 * Description:
 *	Inserts an element in a set initialized by InitSet.  No effect if
 *	the element was already in the set.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	theSet	Modified to contain the new element.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet.  No bounds checking is
 *	done, so specifiying a position < 0 or >= set size is will have
 *	undefined results.
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
void AddSet( ClosedSet theSet, int position )
{
	/* Get the bit position, skipping over length value stored at start */
	int 	byte = (position >> LOGBITS) + sizeof(int);	/* which byte */
	int	bit = position & ENDMASK;	/* which bit of that byte */

	theSet[byte] |= mask[bit];
}


/*==============================================================================
 *
 * Function:	CheckSet- Determine if a given element is in a set.
 *
 * Synopsis:
 *
 *	int CheckSet(
 *		ClosedSet	theSet,		** IN - set to check
 *		int		position )	** IN - element to check
 *
 * Description:
 *	Determine whether a given element has previously been inserted
 *	in the specified set.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	Return Value:
 *		Non-zero value if element is in set; zero otherwise
 *		Returns zero for out-of-bounds values.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet.  
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
int CheckSet( ClosedSet theSet, int position )
{
	int	length;
	int 	byte;
	int	bit;

	length = SETSIZE( theSet );
	if( position < 0 || position >= length )
		return 0;

	byte = (position >> LOGBITS) + sizeof(int);
	bit = position & ENDMASK;

	return theSet[byte] & mask[bit];
}


/*==============================================================================
 *
 * Function:	DropSet- Remove an element from a set, if it was there
 *
 * Synopsis:
 *
 *	void DropSet(
 *		ClosedSet	theSet,		** IN/OUT - set to change
 *		int		position )	** IN - element to remove
 *
 * Description:
 *	Remove an element from a set, if it was there.  If it was not
 *	in the specified set, no change is made.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	The set passed in is modified.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet.  No bounds checking is
 *	done, so specifiying a position < 0 or >= set size is will have
 *	undefined results.
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
void DropSet( ClosedSet theSet, int position )
{
	int 	byte = (position >> LOGBITS) + sizeof(int);
	int	bit = position & ENDMASK;
	
	theSet[byte] &= ~(mask[bit]);
}

/*==============================================================================
 *
 * Function:	ClearSet - make the set empty
 *
 * Synopsis:
 *
 *	void ClearSet(
 *		ClosedSet	*theSet)	** IN/OUT - set to clear
 *
 * Description:
 *	Makes the set empty.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	None.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
void ClearSet( ClosedSet theSet )
{
	int length;

	length = SETSIZE( theSet );

	memset( theSet + sizeof(int), 0,
		CLOSED_SET_BYTES(length) - sizeof(int) );
}

/*==============================================================================
 *
 * Function:	UnionSet - computes the union of two sets
 *
 * Synopsis:
 *
 *	void ClearSet(
 *		ClosedSet	*dest,		** IN/OUT - set to be modified
 *		ClosedSet	*source)	** IN 	  - set to be added
 *
 * Description:
 *	Computes the union of source and dest and stores the result in
 *	source.  If the maximum length of source is smaller that the
 *	maximum length of dest, the union works correctly (elements
 *	in dest beyond source's maximum are untouched).  However if
 *	the maximum length of source is greater than the maximum
 *	length of dest, dest is *not* expanded to include any elements
 *	beyond its maximum length.  Only elements up dest's maximum
 *	are included.  In that sense, this function doesn't compute
 *	a true union of the sets.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	None.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *	See description; true union operation is not performed if
 *	source contains elements beyond dest's maximum length.
 *
 * Assumptions:
 *	Both sets should have been initialized by InitSet
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
void UnionSet( ClosedSet dest, ClosedSet source )
{
	int bytes, lsrc, ldest, i;
	char * psrc, * pdest;

	lsrc = CLOSED_SET_BYTES(SETSIZE( source )) - sizeof(int);
	ldest = CLOSED_SET_BYTES(SETSIZE( dest )) - sizeof(int);
	bytes = ( lsrc < ldest ) ? lsrc : ldest;

	for( psrc = source + sizeof(int), pdest = dest + sizeof(int), i = 0;
		i < bytes; i++, psrc++, pdest++ ) {
		*pdest |= *psrc;
	}
}
		
/*==============================================================================
 *
 * Function:	IsEmpty - determines whether a set contains any elements
 *
 * Synopsis:
 *
 *	int IsEmpty(
 *		ClosedSet	theSet)	** IN - set to examine
 *
 * Description:
 *	Determines whether a set contains any elements.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	Returns non-zero if the set is empty; zero otherwise.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
int IsEmpty( ClosedSet theSet )
{
	return TestSet( theSet, 0 );
}

/*==============================================================================
 *
 * Function:	IsFull - determines whether a set is full
 *
 * Synopsis:
 *
 *	int IsFull(
 *		ClosedSet	theSet)	** IN - set to examine
 *
 * Description:
 *	Determines if all potential elements of a set have been
 *	inserted into it.
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	Returns non-zero if all elements are present; zero otherwise.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/
int IsFull( ClosedSet theSet )
{
	return TestSet( theSet, 1 );
}

/*==============================================================================
 *
 * Function:	TestSet - examines membership of a set
 *
 * Synopsis:
 *
 *	int TestSet(
 *		ClosedSet	theSet,	** IN - set to examine
 *		int		state )	** IN - state to test for
 *
 * Description:
 *	Utility function that IsEmpty and IsFull use to determine
 *	whether a set is empty or full.  Value of  state  parameter
 *	selects the test: 0 for emptiness, non-zero for fullness.
 *	Works by comparing all the bytes representing groups of
 *	set elements (except the last one) with a mask that has
 *	all its bits set to 1 (for fullness test) or 0 (for
 *	emptiness test).  Depending on the number of elements in
 *	the set, the last byte may not be fully used, so a
 *	special mask must be constructed to test it for fullness.
 *	(We assume the unused bits will be zero.)
 *
 * Other Inputs:
 *	None.
 *
 * Outputs:
 *	Returns non-zero if selected test is true; zero otherwise.
 *
 * Interfaces:
 *
 * Resources Used:
 * 
 * Limitations:
 *
 * Assumptions:
 *	theSet should have been initialized by InitSet
 *
 * Notes:
 *
 *----------------------------------------------------------------------------*/

static int TestSet( ClosedSet theSet, int state )
{
	int	length;
	int	bytes;
	int	bits;
	int	firstmask = (state ? ALLBITS : 0 );
	int	lastbyte;
	int	lastmask = 0;
	int	i;

	length = SETSIZE( theSet );

	bytes = ((length - 1) >> LOGBITS) + 1;
	bits = length & ENDMASK;
	if( !bits ) bits = BITSPERBYTE;	/* last byte uses all its bits */
	lastbyte = bytes + sizeof(int) - 1;

	/* See if all bits in first length-1 bytes match our mask
	 * (after skipping length field)
	 */
	for( i = sizeof(int); i < lastbyte; i++ ) {
		if( (char)(theSet[i]) != (char)(firstmask) ) return 0;
	}
	
	/* Cosntruct mask to check bits in final byte, since
	 * not all bits may be used.
	 */
	if( state ) {
		for( i = 0; i < bits; i++ ) {
			lastmask |= mask[i];
		}
	}
	/* Otherwise, for emptiness test, leave lastmask at 0 */

	/* Do the final bits match? */
	return ( (char)(theSet[lastbyte]) == (char)(lastmask));
}
