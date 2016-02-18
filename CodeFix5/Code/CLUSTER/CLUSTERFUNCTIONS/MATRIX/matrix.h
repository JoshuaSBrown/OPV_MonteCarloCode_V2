#ifndef _MATRIX_H_
#define _MATRIX_H_

/*the type declartion of the ADT. */
typedef struct _matrix * matrix;
typedef struct _matrix const * const_matrix;
/* Creates a rows by columns matrix 
	 Returns NULL if rows <=0 or cols <=0 
	 and otherwise a pointer to the new 
	 matrix is returned. All the elements
	 are set equal to 0.0
 */
matrix newMatrix(int rows, int cols);

/* Does the same thing as newMatrix but
	 sets all the elements equal to set 
	 instead of 0.0
 */
matrix newMatrixSet(int rows, int cols, double set);

/* Duplicates an entirely new matrix
*/
matrix duplicateMatrix(const_matrix mtx);

/* Deletes a matrix Returns 0 if successful
	 and -1 if mtx is NULL.*/
int deleteMatrix(matrix * mtx);

int printMatrix(const_matrix mtx);

/* If matrix has three cols of data will
	 print the values in three seperate files.
	 Returns -1 if the matrix doesn't have 3 columns
	 or is NULL. Returns 0 if successful. 
 */
int printRandomEnergyFilesijk(const_matrix mtx);

/* Returns the number of rows of matrix mtx
 * if successful. Returns -1 if mtx is NULL.
 */
int getRows(const_matrix mtx);

/* Returns the number of columns of matrix mtx
	 if successful. Returns -1 if mtx is NULL.
 */
int getCols(const_matrix mtx);

/* resizes the matrix to length Row
	 but keeps the column widith of the
	 original.
	 If the new matrix has fewer rows it
	 will not copy all the data from
	 mtx.
	 Will return NULL if failed
 */
int resizeRow(matrix * mtx, int Row);

/* setAll is used to set all the elements
	 in a matrix mtx to a value of val
	 returns 0 if succesful and -1 if failed
 */
int setAll(matrix mtx, double val);

/* Sets all the rows in column col to 
	 the value specefied by val.
*/
int setAllRowsInCol(matrix mtxProbNeighDwell,int col, double val);

/*Sets the (row, col) element of mtx to 
	a val. Returns 0 if successful, -1 if
	mtx is NULL, and -2 if row or col are
	outside of the dimensions of mtx*/
int setE(matrix mtx, int row, int col, double val);

/*Gets the reference val of the (row, col) 
 * element of mtx.Returns 0 if successful, 
 * -1 if either mtx is NULL, and -2 
 * if row or col are outside of the dimenisons of 
 the mtx. */
double getE(const_matrix mtx, int row, int col);

/* Checks to see if there is a number within the
	 matrix mtx with value match. If there is ret-
	 urns a 1 if there isn't returns a 0
*/
int matchExist(const_matrix mtx, double match);

/* Checks to see if there is a number within the 
	 matrix in the column col that matches the value
	 of match. If there is returns the row number
	 else it returns a 0.
*/
int matchExistCol(const_matrix mtx, int col, double match);

/* Searches through the matrix mtx and find any
	 elements with the same value as match
	 these elements are replaced with the value
	 replace
*/
int matchReplace(matrix mtx, double match, double replace);

/* Searches through the column col in the matrix
	 mtx finds the first value that matches 
	 and returns the row. 
*/
int FindRowOfMatchInCol(const_matrix mtx, double match, int col);

/* Function adds up all the columns of data in a 
	 given row and returns the value
*/
double SumOfRow(matrix mtx, int row);

/* Adds all the rows in a given column col and re-
	 turns the value 
*/
double SumOfCol(matrix mtx, int col);

/* Divides each element in the matrix by the value
	 val
*/
int DivideEachElement(matrix * mtx, double val);

/* Divide each element in the column col by the 
	 value val
*/
int DivideEachElementCol(matrix * mtx, int col, double val);

#endif
