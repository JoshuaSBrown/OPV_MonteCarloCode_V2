#ifndef _ARBARRAY_H_
#define _ARBARRAY_H_

#include <stdarg.h>
#include <stdlib.h>

#include "../MIDPOINT/midpoint.h"
#include "../OMLL/omll.h"

/* This type is an array where each element points
	 to either a MidPoint an OrderMagLL or a ClusterLL
	 The type of data held in the array is recorded as
	 a type attribute in the array. i.e.
	 Arb->type=0 - Stores OrderMagLL
	 Arb->type=1 - Stores ClusterLL
	 Arb->type=2 - Stores MidPoints
*/
typedef struct _ArbArray * ArbArray;
typedef struct _ArbArray const * const_ArbArray;

/////////////////////////////////////////////////////////
//Tools for accessing Arbirtrary array

/* Creates a new Arbitrary Array the size of the array is 
	 input as an integer and stored as a reserved attribute
	 the number of elements used in the array is defined by
	 the used attribute.
	 The elements range between 0 and < the total reserved
	 The array is defined as a specefic type. Where the 
	 elements are a certain data structure type i.e.
	 type = 0 - Order of Magnitude Link List Array
	 type = 1 - Cluster Link List Array
	 type = 2 - Mid Point Array
*/
ArbArray newArbArray(int len, int type);

/* The purpose of this function is to ensure that stack
	 overflow did not occur will return a -1 if the type 
	 of ArbArray has been changed to a number that is not
	 0, or 1 (Only checks to ensure ArbArray is an Array
	 of Link list or Order of magnitude link list not 
	 a mid point array which would be of type 2). 
*/
int ArbArrayCheck(const_ArbArray Arb);

int deleteAllMidPointArray(ArbArray * Arb);

/* This function deletes the arbitray array according to
	 it's type if it is of type 0 it only deletes the Order
	 of magnitude data structure and does not delete the 
	 Mid Points. 
	 If it is of type 1 it deletes all the nodes and the
	 link list.
	 If it is of type 2 it does not delete the mid points
	 only the Arbitrary array itself.
	 Returns a 0 if succesful and -1 if failed
*/
int deleteArbArray(ArbArray * Arb);

/* Prints the contents of the arbitrary array according
	 to it's type, a second parameter should be submitted
	 if the array is of type 1. The lowest order of 
	 magnitude should be submitted as a second argument:
	 printArbArray(ArbArray Arb, int orderLow);
*/
int printArbArray(const_ArbArray Arb, ...);

/*Sets the element of Arb to NULL
*/
int NullArbElement(ArbArray Arb, int element);

/* Append Arb element to the end, this is not the fastest
 * way to do this but we simply create a new Arb Array 
 * structure and copy over the elements from the last one
 * then we delete the old Arb Array and set it's pointer
 * equal to the new one. 
 */
int appendArbElement(ArbArray * Arb, void * ptr);

/* Sets an element of the ArbArray, increments the number
	 of elements used if it is NULL
*/
int setArbElement(ArbArray Arb, int element, void * ptr);

/* Grabs the elements of the ArbArray does not remove it
	 from the array though.
*/
void * getArbElement(const_ArbArray Arb, int element);

/* Returns the number of elements used
*/
int getElementsUsed(const_ArbArray Arb);

/* Returns the number of elements resereved
*/
int getElementsReserved(const_ArbArray Arb);

/* Removes the elements from the array and makes it NULL
*/
int removeArbElement(ArbArray Arb, int element);

/* Grabs the Midpoint at element 
 */
MidPoint getMP(const_ArbArray Arb,int element);

/* Grab neighbor one of midpoint Mid_ID specified in the
	 array mpA
	 will return the nei1 id or -1 if failed
*/
int getMPnei1(const_ArbArray Arb, int Mid_ID);

/* Grab neighbor two of midpoint Mid_ID specified in the 
	 array mpA
	 will return the nei1 id or -1 if failed
*/
int getMPnei2(const_ArbArray Arb, int Mid_ID);

/* Get the order of magnitude of the hopping rate for mid
	 Mid_ID from midpoint array mpA.
*/
int getMPOrder(const_ArbArray Arb, int Mid_ID);

/* Add a Mid Point to the Order of Magnitude Link list 
	 array specefied as Arb->Num[element]. This link list 
	 contains elements that point to the order of mangitude
	 link lists. 
	 If sucessful return 0
	 If Mid point was previously added returns -1
*/
int addToOrLL(ArbArray Arb, int element, MidPoint * mp);

/* This function gets the "OrderMagLL" and returns the 
	 pointer.
*/
OrderMagLL getOrderLL(const_ArbArray Arb, int element);

/* This function sets the default values for the order of 
	 magnitude LL in Arb->Num[element] it will only work
	 if the Arbitrary array used is of type 1
	 If a new LL is added and the elements was not NULL than
	 the size will be incremented.
*/
int setDefaultArbElem(ArbArray Arb, int element, int orderMag);

#endif
