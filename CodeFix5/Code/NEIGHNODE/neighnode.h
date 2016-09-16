#ifndef _NEIGHNODE_H_
#define _NEIGHNODE_H_

#include <stdarg.h>
#include <stdlib.h>

typedef struct _NeighNode * NeighNode;
typedef struct _NeighNode const * const_NeighNode;

/* Each NeighborNode can have multiple hops that
	 end on it. The data structure Hop is used
	 to keep track of how many hops are allowed
	 as well as the time it takes when a hop is 
	 made and the probability that the specific 
	 hop is made.
*/
typedef struct _Hop * Hop;
typedef struct _Hop const * const_Hop;

/////////////////////////////////////////////////////////
//Tools for accessing Neighbor Nodes
NeighNode newNeighNode(int N_ID);

/* Deletes a node
*/
int deleteNeighNode(NeighNode NeighNod);

/* Prints contents of node
*/
int printNeighNode(const_NeighNode NeighNod);

/* Get the NeighNode ID and set as an integer
*/
int getNeighNode_id(const_NeighNode NeighNod);

int setNextNeighNode(NeighNode * NeighNod, NeighNode * NeighNod2);

/* Set the Next NeighNode to a NULL value
 */
int setNextNeighNodeToNULL(NeighNode * NeighNod);

/* Gets the next Neighbor in line
*/
NeighNode getNextNeigh(const_NeighNode NeighNod);

/* Sets the probability of moving to this neighboring node
	 for a new hop
*/
int setNeighNodeNew_p(NeighNode NeighNod, double pval);

/* Set the probability of moving to this neighboing node
	 taking the hop route specified by Elem
*/
int setNeighNode_p(NeighNode NeighNod, double pval, int Elem);

/* Sets the time of moving to this neighboring node
	 given that it takes the hop route specified by 
	 Elem.
*/
int setNeighNode_t(NeighNode NeighNod, double time, int Elem);

/* Gets the value of probability of hopping to NeighNod from
	 the hop route speciefied by Elem.
*/
double getNeighNode_p(const_NeighNode NeighNod, int Elem);

/* Gets the length of the number of hops to a given neighboring
	 site off the cluster
*/
int getNeighNode_hoplength(NeighNode NeighNod);

/* Gets the value of the time of hopping to NeighNod from
	 the hop route specified by Elem.
*/
double getNeighNode_t(NeighNode NeighNod, int Elem);

/* Set the NeighNode Id
*/
int setNeighNode_id(NeighNode NeighNod, int ID);

/* Sets the starting hop from the NeighNod
*/
int setNeighNode_hopstart(NeighNode NeighNod, Hop h);

/////////////////////////////////////////////////////////
//Tools for accessing Hop data structure (no getters puroposefully)
Hop newHop(void);

//Deletes all Hops linked to the hop h
int deleteAllHop(Hop h);

int printHop(Hop h);

int setHop_t(Hop h, double t);

int setHop_p(Hop h, double p);

int setHop_next(Hop h, Hop h2);
#endif
