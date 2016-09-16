#ifndef _NODE_H_
#define _NODE_H_

#include <stdarg.h>
#include <stdlib.h>

typedef struct _Node * Node;
typedef struct _Node const * const_Node;

/////////////////////////////////////////////////////////
//Tools for accessing Nodes
/* Function creates a Node. The Node acts like a
	 Link list it contains the id of a node and points to
	 another Node that is also in the same cluster
*/
Node newNode(int N_ID);

/* Deletes a node
*/
int deleteNode(Node Nod);

/* Deletes all nodes that are connected to Nod
 * as well as itself
 */
int deleteNodeAll(Node * Nod);

/* Prints contents of node attributes
*/
int printNode(const_Node Nod);

/* get the node id and return as an integer
*/
int getNode_id(const_Node Nod);

/* set the next node
 */
int setNextNode(Node * Nod, Node Nod2);

/* get the next node
*/
Node getNextNode(const_Node Nod);

/* Set the Node id
*/
int setNode_id(Node Nod, int ID);

/* Set Node's probability value
*/
int setNode_p(Node Nod, double pval);

/* Get Node's probability value
*/
double getNode_p(const_Node Nod);

/* Gets the last Node attached to Nod
 */
Node getNode_lastNode(Node Nod);

/* Set flags
	 If a node is within the 
	 cluster it is set to 1
	 if it is outside it is
	 left as a 0
*/
int setFlagFro(Node Nod);
int setFlagBeh(Node Nod);
int setFlagLef(Node Nod);
int setFlagRig(Node Nod);
int setFlagAbo(Node Nod);
int setFlagBel(Node Nod);

int getFlagFro(const_Node Nod);
int getFlagBeh(const_Node Nod);
int getFlagLef(const_Node Nod);
int getFlagRig(const_Node Nod);
int getFlagAbo(const_Node Nod);
int getFlagBel(const_Node Nod);

int getFlag(const_Node Nod, int index);
int setFlag(Node Nod, int index);

#endif
