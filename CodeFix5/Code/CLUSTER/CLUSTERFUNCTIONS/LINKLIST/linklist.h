
#ifndef _LINKLIST_H_
#define _LINKLIST_H_

/* Link list contains the length of the list
	 as well as the starting LLNode
*/
typedef struct _linklist * linklist;
typedef struct _linklist const * const_linklist;

/* LLNode contains the id of the node as well
	 as the pointer to the next node in line
*/
typedef struct _LLNode * LLNode;
typedef struct _LLNode const * const_LLNode;

/* Creates a new LLNode with id as the id
	 and returns it. The next LLNode that it
	 points to is initialized to NULL
*/
LLNode newLLNode(int id);

/* Creates a new link list starting with the
	 LLNode of id and length 1
*/
linklist newLinkList(int id);

/* Creates a blank new link list
*/
linklist newBlankLinkList();

/* Deletes the LL and all the nodes
	 within it
*/
int deleteLL( linklist * LL);

/* Deletes only the LLNode node
*/
int deleteLLNode(LLNode node);

/* Gets the nextLLNode attached
	 to node
*/
LLNode nextLLNode(LLNode node);

/* Adds a node with an id to the linklist
	 LL if a node does not already exist with
	 that identity. Also increments the length
	 of the LL
*/
int addLLNode(linklist LL, int id);

/* Removes the LLNode with identity id
	 then splices the two LLNodes that are
	 adjacent if there are any. Reinitializes
	 the starting LLnode of LL if the node removed
	 is at the start of the LL.
*/
int removeLLNode(linklist LL, int id);

/* Gets the length of the linklist LL
*/
int getLLlength(linklist LL);

/* Grabs the id of the starting node in the linklist
*/
int getLLstartID(linklist LL);

/* Prints the length of the LL
	 and the ids of the nodes contained
	 within in theLL
*/
int printLL(linklist LL);

#endif
