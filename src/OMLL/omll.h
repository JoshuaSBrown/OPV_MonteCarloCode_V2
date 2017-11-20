#ifndef _OMLL_H_
#define _OMLL_H_

#include <stdarg.h>
#include <stdlib.h>

#include "../MIDPOINT/midpoint.h"

/* Multiple Link list of this type are used to sort the
	 mid pts such that all the mid points of a given order
	 exists in a given link list. 
	 Mid_ID ranges between 0 and < total reserved
*/
typedef struct _OrderMagLL * OrderMagLL;
typedef struct _OrderMagLL const * const_OrderMagLL;

////////////////////////////////////////////////////////
//Tools for accessing Order Magnitude Link List
/* This function creates an order of magnitude link list
	 starting at the mid point mp. The order of magnitude
	 is also defined and stored in the orderMag attribute
*/
OrderMagLL newOrLL(int orderMag);

/* Checks to see if all the MidPoints linked with mp
 * are of the same order of magnitude as OMLL
 */
int checkNewOrLL(OrderMagLL OMLL, MidPoint mp);

/* This function only deletes the OrderMagnitude data
	 type. It does not delete the MidPoints referenced
	 by the link List. 
*/
int deleteOrLL(OrderMagLL * OMLL);

/* Only deletes the OMLL the midpoints if there
 * are any are not touched
 */
int deleteOrLLwithOutDeletingMP(OrderMagLL * OMLL);

/* Delete the Midpoint that is located at the start
 * of the OMLL
 */
int deleteOMLL_startMP(OrderMagLL OMLL);

/* This function deletes the Midpoints held in the link
	 list as well as the actual OMLL.
*/
int deleteAllOrLL(OrderMagLL OMLL);

/* Print mid points from link list
*/
int printOrLL( const_OrderMagLL OMLL);

/* Grab size of OMLL
*/
int getOMLL_size(const_OrderMagLL OMLL);

/* Grans the value of the order
*/
int getOMLL_order(const_OrderMagLL OMLL);

/* Function grabs the starting Midpoint of the OMLL
*/
MidPoint getOMLLstartMP(const_OrderMagLL OMLL);

/* Setting the start of the OMLL with the first midpoing
 */
int setOMLL_startMP(OrderMagLL OMLL, MidPoint mp);

/* Adds MidPoint to an OMLL provided that there is not
 * already a midpoint with the same id in the OMLL
 * returns a value:
 * -2 - if another midpoint with the same id already exists
 *  0 - sucess
 * -1 - malformed input
 */
int addMPToOMLL(OrderMagLL OMLL, MidPoint mp);
#endif
