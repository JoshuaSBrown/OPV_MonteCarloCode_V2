#ifndef _MIDPOINT_H_
#define _MIDPOINT_H_

#include <stdarg.h>
#include <stdlib.h>

/* MidPoint structure are points between 
	 sites in a lattice. They contain the information
	 regarding the order of magnitude of the hopping 
	 rates between the two sites. 
	 Information is only stored in a mid point when both
	 hop rates from site1 to site2 and from site2 to site1
	 are on the same order of magnitude.
*/
typedef struct _MidPoint * MidPoint;
typedef struct _MidPoint const * const_MidPoint;
//Tools for accessing Mid Points
/* newMidPoint creates a mid point between two neigboring
	 sites. It stores the order of magnitude of the hoprate
	 between the two sites as well.
	 Function will return NULL if it is unable to allocate
	 memory for the point or if the other parameters are 
	 inappropriate
*/
MidPoint newMidPoint(int order, int Mid_ID, int nei1, int nei2);

/* Delete MidPoint simply deletes the Mid point passed
	 Will return -1 if mp is NULL
	 Will return 0 if successful
*/
int deleteMidPoint(MidPoint mp);

/* For printing the contents of a mid point
	 Will return 0 if successful and -1 if
	 mp is NULL
 */
int printMP(const_MidPoint mp);

/* Grabs the order of magnitude of the midpoint
*/
int getMP_order(const_MidPoint mp);

/* Grab the Mid Point Id return as an integer
	 or returns -1 if mp is NULL
*/
int getMP_id(const_MidPoint mp);

/* Grab the id of the nei1 of the mp
 */
int getMP_nei1(const_MidPoint mp);

/* Grab the id of the nei2 of the mp
 */
int getMP_nei2(const_MidPoint mp);

/* Get the next midpoint
 */
MidPoint getMP_next(const_MidPoint mp);

/* Set the Mid Point id will return 0 if success
	 and -1 if failed
*/
int setMP_id(MidPoint mp, int ID);

/* Links the midpoint mp with nextmp
 */
int setMP_nextMP(MidPoint mp, MidPoint nextmp);

/* This function compares the neighbors of two mid points
	 if the mid points share neighboring sites than the 
	 function returns a 1. If they do not share sites then
	 a 0 is returned. 
*/
int CompareNeiMidPoint(const_MidPoint mp1, MidPoint mp2);

/* Adds a midpoint to the end of the midpoing link list
 * will not add the midpoint if has the same id as one 
 * already in the list
 * returns:
 *  0 - success
 * -1 - if there is an error malformed input
 * -2 - if another mp with same id is foudn in the list
 */
int addMPToEnd(MidPoint mp1, MidPoint mp2);

/* Grab the next MP
 */
MidPoint getNextMP(const_MidPoint mp);
#endif
