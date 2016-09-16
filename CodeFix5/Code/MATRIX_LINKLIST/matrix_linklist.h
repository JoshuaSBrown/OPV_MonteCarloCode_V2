
#ifndef _MATRIX_LINKLIST_H_
#define _MATRIX_LINKLIST_H_

#include "../MATRIX/matrix.h"

/* Unlike LINKLIST and LINKLIST2, MATRIX_LINKLIST
 * allows you to give a node an array of doubles
 * of any size. The node stores the double as an
 * attribute of type matrix integer
 */

/* matrix_linklist contains the length of the list
	 as well as the starting LL_MNode
*/
typedef struct _matrix_linklist * matrix_linklist;
typedef struct _matrix_linklist const * const_matrix_linklist;

/* LL_MNode contains a matrix as well as a node
 * that is connected to it. 
*/
typedef struct _LL_MNode * LL_MNode;
typedef struct _LL_MNode const * const_LL_MNode;

/* Creates a new LL_MNode when an array of 
 * doubles is passed to it. The next LLNode 
 * that it points to is initialized to NULL
*/
LL_MNode newLL_MNode(double data[],int rows);

/* Creates a new link list starting with the
	 LL_MNode with the matrix the initial node
   initialized with the double array.
*/
matrix_linklist newMatrix_LinkList(double data[], int rows);

/* Creates a blank new matrix_linklist
*/
matrix_linklist newBlankMatrixLinkList();

/* Deletes the LL and all the nodes
	 within it
*/
int deleteMatrixLL( matrix_linklist * LL);

/* Deletes only the LLNode node
*/
int deleteLL_MNode(LL_MNode * node);

/* Gets the nextLLNode attached
	 to node
*/
LL_MNode nextLL_MNode(LL_MNode node);

/* Adds a node with the double array data to 
 * the matrix_linklist LL at the end.It also \
 * increments the length of the LL.
 * data is an array of doubles
 * rows describes how long data is
*/
int addLL_MNode(matrix_linklist * LL, double data[],int rows);

/* This function is similar to the above but
 * instead of adding the node to the end it adds
 * it to the begginning of the matrix linklist
 * data is an array of doubles
 * rows describes how long data is
 */
int addLL_MNodeBegin(matrix_linklist * LL, double data[],int rows);

/* Removes the LLNode with the sequence seq
	 then splices the two LLNodes that are
	 adjacent if there are any. Reinitializes
	 the starting LLnode of LL if the node removed
	 is at the start of the LL.
*/
int removeLL_MNode(matrix_linklist * LL, int seq);

/* Removes the last node on the LL
 */
int removeLL_MNodeEnd(matrix_linklist * LL);

/* Removes the LLNodes between and including seq1
 * and seq2 than joins the nodes on either side
 */
int removeLL_MNodes(matrix_linklist * LL, int seq1, int seq2);

/* Gets the length of the linklist LL
*/
int getMatrixLLlength(matrix_linklist LL);

/* Grabs the row R of the of the starting node in 
 * the linklist
*/
double getMatrixLLstartAtRow(matrix_linklist LL, int R);

/* Sets the value at seq, and R 
 */
int setMatrixLLElem(matrix_linklist LL, int seq, int R, double val);

/* This function will set every elem of each node with value
 * val at row R. 
 */
int setMatrixLLValueAtRow(matrix_linklist LL, int R, double val);

/* Prints the length of the LL
	 and the ids of the nodes contained
	 within in the LL
*/
int printMatrixLL(matrix_linklist LL);

int printMNode(LL_MNode node);

/* This function finds looks through the link list and
 * find the position in the link list of the last match
 * at Row R and returns it. If there are no matches it
 * returns -1
 */
int getMatrixLLLastMatchAtRow(matrix_linklist LL, double match, int R);

/* This function finds the number of nodes in the linklist
 * that have the same number as that which is submitted. If 
 * there are no matches returns 0 
 */
int getMatrixLLNumberMatchAtRow(matrix_linklist LL, double match, int R);

/* Function returns the Elem at row 1 that are a match with row
 * R
 */
matrix getMatrixLLMatchAtRow(matrix_linklist LL, double match, int R,  int * numMatches);

/* Determines the number of elements in the matrix_linklist
 * at row R that have a value greater than match. Checks each
 * node. The return values is the total count. A value of -1 
 * is returned if there is an error.
 */
int getMatrixLLNumberElemGreaterThanMatchAtRow(matrix_linklist LL, double match, int R);

/* Function grabs the elem of the Node listed at position seq
 * in the LL and row R in the data matrix.
 */
double getMatrixLLNodeElem(matrix_linklist LL, int seq, int R);

/* Function moves the node defined at position seq in the 
 * matrix linklist to the start of the linklist. 
 */
int moveMatrixLLNodeToStart(matrix_linklist LL, int seq);

/* Function moves the node defined at position seq in the 
 * matrix linklist to the end of the linklist
 */
int moveMatrixLLNodeToEnd(matrix_linklist LL, int seq);

/* This funciton rearranges the Nodes so that the node at
 * position seq1 is exchanged with the node at position 
 * seq2.
 */
int exchangeMatrixLLNodes(matrix_linklist LL, int seq1, int seq2);

/* This function will increment the chosen element provided
 * it is within the scope of the matrix linklist
 */
int incMatrixLLElem(matrix_linklist LL, int seq, int row);

#endif
