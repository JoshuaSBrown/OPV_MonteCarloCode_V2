#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//#include "../../../MEM/mem.h"
#include "matrix.h"

int main() {//still need DivideEachElement in Matrix and FindRowOfmatchincol
	
	//mem_init();
	printf("Testing: functions with mtx as null \n");
	matrix mtx = NULL; 
	int rv;
	int row=4;
	int col=3;
	int val = -999; 

	printf("Testing: printMatrix with mtx as null \n");
	rv = printMatrix(mtx);
	assert(rv==-1);
	printf("Testing: deleteMatrix with mtx as null \n");
	rv = deleteMatrix(&mtx);
	assert(rv==-1);
	printf("Testing: setAll with mtx as null \n");
	rv = setAll(mtx, val);
	assert(rv==-1);
	printf("Testing: printRandomEnergyFilesijk with mtx as null \n");
	rv= printRandomEnergyFilesijk(mtx); 
	assert(rv==-1);
	printf("Testing: getE with mtx as null \n");
	rv = getE(mtx,row, col);
	assert(rv==-1);
	printf("Testing: getRows with mtx as null \n");	
	rv = getRows(mtx);
	assert(rv==-1);
	printf("Testing: getCols with mtx as null \n");
	rv = getCols(mtx);
	assert(rv==-1);
	printf("Testing: setE with mtx as null \n");	
	rv = setE(mtx,row, col, val);
	assert(rv==-1);
	printf("Testing: deleteMatrix with mtx as null \n");
	rv=deleteMatrix(&mtx);
	assert(rv==-1);

	printf("Testing: functions with functional mtx\n");
	matrix mtxrv; 
	mtx = newMatrix(row, col);

	printf("Testing: newMatrix\n");
	mtxrv=newMatrix(0,1);
	assert(mtxrv==NULL);
	mtxrv=newMatrix(1,0);
	assert(mtxrv==NULL);
	printf("Expecting matrix of 0's\n");
	mtxrv= newMatrix(5,5); 
	assert (mtxrv!=NULL);

	printf("Testing: printMatrix\n");
	rv=printMatrix(NULL);
	assert(rv==-1);
	rv= printMatrix(mtxrv);
	assert(rv==0);

	printf("Testing: setE with functional mtx\n");
	rv=setE(mtx, 2, 5, 12); 
	assert(rv==-2);
	rv= setE(mtx, 2, 3, 343.1);
	assert(rv==0);
	rv= setE(mtx,3,2,902.8);
	assert(rv==0);
	printMatrix(mtx);
	printf("Testing: setE returns null properly\n");
	rv = setE(mtx,-1, col, val);
	assert(rv==-2);
	rv=setE(NULL,-1, col, val);
	rv = setE(mtx,15, col, val);
	assert(rv==-2);
	rv = setE(mtx,row, -1, val);
	assert(rv==-2);
	rv = setE(mtx,row, 15, val);
	assert(rv==-2);

	printf("Testing: getE with functional mtx\n");
	double drv;
	drv= getE(mtx, 2, 3);
	printMatrix(mtx);
	printf("drv %g\n",drv);
	assert(drv==343.1);
	drv= getE(mtx,3,2);
	//printf("%g\n",drv);
	assert(drv==902.8);
	printf("Testing: getE returns null properly\n");
	drv = getE(mtx,-1, col);
	assert(drv==-2);
	drv = getE(mtx,15, col);
	assert(drv==-2);
	drv = getE(mtx, row, -1);
	assert(drv==-2);
	drv = getE(mtx,row, 15);
	assert(drv==-2);
	drv=getE(NULL,row, col);
	assert(drv==-1); 
		
	printf("Testing: resizeRow\n");
	matrix mtx4 = NULL;
	rv = resizeRow(&mtx4,3); //why the ampersand?
	assert(rv==-1);
	printf("Test 1 past\n");
	mtx4 = newMatrix(5,3);
	printMatrix(mtx4);
	rv = resizeRow(&mtx4,4);
	assert(rv==0);
	printf("Test 2 past\n");
	rv = getRows(mtx4);
	printf("rows %d\n",rv);
	assert(rv==4);
	printf("Test 3 past\n");
	rv = resizeRow(&mtx4,0);
	assert(rv==-1);
	printf("Test 4 past\n");
	setE(mtx4,1,1,34);
	setE(mtx4,2,1,31);
	printMatrix(mtx4);
	rv = resizeRow(&mtx4,6);
	assert(rv==0);
	printf("Test 5 past\n");
	/*setE(mtx4,6,2,-12.3);
	printMatrix(mtx4);
*/

	matrix mtx3 = newMatrix(row,15);
	printf("Testing: printRandomEnergyFilesijk\n");
	rv= printRandomEnergyFilesijk(mtx3); 
	assert(rv==-1);
	printf("Testing: getRows\n"); 
	rv= getRows(mtx3); 
	assert(rv==4);
	printf("Testing: getCols\n"); 
	rv=getCols(mtx3); 
	assert(rv==15);
	resizeRow(&mtx3,3);
	rv=printRandomEnergyFilesijk(mtx3);
	assert(rv==-1); 

	printf("Testing: deleteMatrix\n");
	rv=deleteMatrix(NULL);
	assert(rv=-1);
	rv=deleteMatrix(&mtx);
	assert(rv==0);

	printf("Testing: duplicateMatrix");
	matrix mtx5=duplicateMatrix(NULL);
	assert(mtx5==NULL);
	mtx5=duplicateMatrix(mtx4);
	assert(mtx5!=NULL);
	printMatrix(mtx4);
	printMatrix(mtx5); 

	printf("Testing: newMatrixSet\n");
	mtxrv=newMatrixSet(0,1,5);
	assert(mtxrv==NULL);
	mtxrv=newMatrixSet(1,0,4);
	assert(mtxrv==NULL);
	mtxrv=newMatrixSet(row,col,14);
	assert(mtxrv!=NULL);
	printf("Expect matrix of 14's: ");
	printMatrix(mtxrv);

	printf("Testing:setAllRowsinCol");
	rv=setAllRowsInCol(mtxrv,0,val);
	assert(rv==-1);
	rv=setAllRowsInCol(mtxrv,10,val);
	assert(rv==-1);
	rv=setAllRowsInCol(NULL,2,val);
	assert(rv==-1);
	printf("\nRows in matrix %d cols %d\nValue of row %d and col %d\n",getRows(mtxrv),getCols(mtxrv),row,col);
	rv=setAllRowsInCol(mtxrv,2,2);
	assert(rv==0);
	printMatrix(mtxrv);

	printf("Testing:DivideEachElementInCol"); //may be problems due to *matrix
	rv=DivideEachElementCol(NULL,2,val);
	assert(rv==-1);
	rv=DivideEachElementCol(&mtxrv,0,val);
	assert(rv==-1);
	rv=DivideEachElementCol(&mtxrv,15,val);
	assert(rv==-1);
	rv=DivideEachElementCol(&mtxrv,2,2);
	assert(rv==0);
	printMatrix(mtxrv); 

	printf("Testing: setAll");
	rv=setAll(NULL,val);
	assert(rv==-1);
	rv=setAll(mtxrv,100);
	assert(rv==0);
	
	printf("Testing: DivideEachElement"); //may be problems due to *matrix
	rv=DivideEachElement(&mtxrv,25);
	assert(rv==0);
	rv=DivideEachElement(NULL,25);
	assert(rv==-1);
	//rv=DivideEachElement(mtxrv,25); //how to test !*mtx without NULL?
	assert(rv==-1);

	printf("Testing: SumOfCol");
	drv=SumOfCol(NULL,val);
	assert(drv==-1);
	drv=SumOfCol(mtxrv,0);
	assert(drv==-1);
	drv=SumOfCol(mtxrv,15);
	assert(drv==-1);
	drv=SumOfCol(mtxrv,3);
	assert(drv!=-1);
	
	printf("Testing: SumOfRow");
	drv=SumOfRow(NULL,val);
	assert(drv==-1);
	drv=SumOfRow(mtxrv,0);
	assert(drv==-1);
	drv=SumOfRow(mtxrv,15);
	assert(drv==-1);
	drv=SumOfRow(mtxrv,3);
	assert(drv!=-1);
	
	setE(mtxrv,2,2,13);
	printf("Testing: FindRowOfMatchInCol"); //may be problem w/ cont_matrix mtx
	rv = FindRowOfMatchInCol(NULL,13,2);  
	assert(rv==-1);
	rv = FindRowOfMatchInCol(mtxrv,13,15);
	assert(rv==-1);
	rv = FindRowOfMatchInCol(mtxrv,0,15);
	assert(rv==-1);
	rv = FindRowOfMatchInCol(mtxrv,13,1);
	assert(rv==-1);
	rv = FindRowOfMatchInCol(mtxrv,13,2);
	assert(rv==2);


	setE(mtxrv,1,2,15);
	printf("Testing: matchReplace");
	rv=matchReplace(NULL,15,4);
	assert(rv==-1);
	rv=matchReplace(mtxrv,15,4);
	assert(rv==0);

	printf("Testing: matchExistCol");
	rv=matchExistCol(NULL,col,4);
	assert(rv==-1);
	rv=matchExistCol(mtxrv,0,4);
	assert(rv==-1);
	rv=matchExistCol(mtxrv,5,4);
	assert(rv==-1);
	rv=matchExistCol(mtxrv,col,val);
	assert(rv==0);
	rv=matchExistCol(mtxrv,3,4);
	printMatrix(mtxrv);
	assert(rv==1); //or at least !0

	printf("Testing: matchExist");
	rv=matchExist(NULL,val);
	assert(rv==-1);
	rv=matchExist(mtxrv,4);
	assert(rv==1);

	//atexit(mem_term);
	return 0;
}
