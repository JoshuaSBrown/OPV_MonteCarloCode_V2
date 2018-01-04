
#ifndef _LINKLIST2_H_
#define _LINKLIST2_H_

/* Unlike LINKLIST, LINKLIST2 allows you to
 * give the node any id interger not just
 * positive values 
 */

/* Link list contains the length of the list
	 as well as the starting LLNode
*/
typedef struct _linklist2 * linklist2;
typedef struct _linklist2 const * const_linklist2;

/* LLNode contains the id of the node as well
	 as the pointer to the next node in line
*/
typedef struct _LLNode2 * LLNode2;
typedef struct _LLNode2 const * const_LLNode2;

/* Creates a new LLNode with id as the id
	 and returns it. The next LLNode that it
	 points to is initialized to NULL
*/
LLNode2 newLLNode2(int id);

/* Creates a new link list starting with the
	 LLNode of id and length 1
*/
linklist2 newLinkList2(int id);

/* Creates a blank new link list
*/
linklist2 newBlankLinkList2();

/* Deletes the LL and all the nodes
	 within it
*/
int deleteLL2( linklist2 * LL);

/* Deletes only the LLNode node
*/
int deleteLLNode2(LLNode2 node);

/* Gets the nextLLNode attached
	 to node
*/
LLNode2 nextLLNode(LLNode2 node);

/* Adds a node with an id to the linklist
	 LL if a node does not already exist with
	 that identity. Also increments the length
	 of the LL
*/
int addLLNode2(linklist2 LL, int id);

/* Removes the LLNode with the sequence seq
	 then splices the two LLNodes that are
	 adjacent if there are any. Reinitializes
	 the starting LLnode of LL if the node removed
	 is at the start of the LL.
*/
int removeLLNode2(linklist2 LL, int seq);

/* Removes the LLNodes between and including seq1
 * and seq2 than joins the nodes on either side
 */
int removeLLNodes2(linklist2 LL, int seq1, int seq2);

/* Gets the length of the linklist LL
*/
int getLLlength2(linklist2 LL);

/* Grabs the id of the starting node in the linklist
*/
int getLLstartID2(linklist2 LL);

/* Prints the length of the LL
	 and the ids of the nodes contained
	 within in theLL
*/
int printLL2(linklist2 LL);

/* This function finds looks through the link list and
 * find the position in the link list of the last match
 * and returns it. If there are no matches returns -1
 */
int getLLLastMatch2(linklist2 LL, int match);

/* This function finds the number of nodes in the linklist
 * that have the same id as that which is submitted. If 
 * there are no matches returns 0 
 */
int getLLNumberMatch2(linklist2 LL, int match);

/* Function grabs the id of the Node listed at position seq
 * in the LL
 */
int getLLNodeID2(linklist2 LL, int seq);
#endif
