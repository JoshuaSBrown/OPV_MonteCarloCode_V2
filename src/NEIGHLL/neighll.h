#ifndef _NEIGHLL_H_
#define _NEIGHLL_H_

#include <stdarg.h>
#include <stdlib.h>

#include "../NEIGHNODE/neighnode.h"

typedef struct _NeighLL * NeighLL;
typedef struct _NeighLL const * const_NeighLL;
///////////////////////////////////////////////////////////
/* Creates a NeighLL datastructure contains the size, id and
	 the neighnodes attached to a cluster
*/
NeighLL newNeighLL(void);

/* Deletes the neighLL
*/
int deleteNeighLL(NeighLL neighLL);

/* Deletes the neighLL and all the NeighNodes in the neighLL
*/
int deleteNeighLLAll(NeighLL neighLL);

/* Prints contents of neighLL
*/
int printNeighLL(const_NeighLL neighLL);

/* Gets the number of neighbors
 */
int getNeighLL_numNeigh(NeighLL Nei);

int setNeighLL_start(NeighLL Nei, NeighNode NeighNod);

NeighNode getNeighLL_start(NeighLL neighll);

/* Adds a Neighbor Node to the end of the NeighLL
   note that you cannot have two neighboring nodes
   with the same id that would defeat the purpose. 
 */
int setNeighLL_addNeighNode(NeighLL NeiLL, NeighNode Nei);
#endif
