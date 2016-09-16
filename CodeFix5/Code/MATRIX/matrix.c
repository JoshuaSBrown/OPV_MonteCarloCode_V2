#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "../ERROR/error.h"

#define ELEM(mtx, row, col) \
  mtx->data[(col-1) * mtx->rows + (row-1)]

struct _matrix {
	int rows;
	int cols;
	matrix Extra;
	double data[0];
};


matrix  newMatrix(int rows, int cols) {
  #ifdef _ERROR_CHECKING_ON_
  if(rows <= 0 ){
     #ifdef _ERROR_ 
     fprintf(stderr,"ERROR rows is less than 1 in newMatrix\n");
     #endif
     #ifdef _FORCE_HARD_CRASH_
     exit(1);
     #endif
     return NULL;
  }
  if(cols <= 0 ){
     #ifdef _ERROR_ 
     fprintf(stderr,"ERROR cols is less than 1 in newMatrix\n");
     #endif
     #ifdef _FORCE_HARD_CRASH_
     exit(1);
     #endif
     return NULL;
  }
  #endif
	int i;
	int j;
 	//printf("new Matrix\n"); 
	//allocate a matrix structure
	int rows1 = rows;
	int rowsNotAllocated = rows;
	int totalAllocated = 0;
	int Flag = 0;

	matrix m;
	matrix mtx;
	matrix * mtx2;

	while( rowsNotAllocated != 0){

		m = (matrix) malloc(sizeof(struct _matrix)+sizeof(double)*rows1*cols);

		if( m!=NULL ){
			
			totalAllocated+=rows1;
			//printf("m not null value of rows1 %d total Allocated %d\n",rows1,totalAllocated);
			if(Flag == 0){
				Flag = 1;
				mtx = m;
				mtx2 = &mtx;
			}else{
				(*mtx2)->Extra = m;
				(*mtx2) = (*mtx2)->Extra;
			}
			m->rows = rows1;
			m->cols = cols;
			m->Extra = NULL;

			for (i = 1; i<=rows1; i++){
				for (j=1;j<=cols;j++) {
					ELEM(m,i,j)=0.000;
				}
			}
			rowsNotAllocated = rowsNotAllocated-rows1;
			rows1 = rowsNotAllocated;

			m = m->Extra;
		}else{

			rows1 = rows1/2;

		}
	}

  #ifdef _ERROR_CHECKING_ON_
	if(mtx==NULL){
    #ifdef _ERROR_
		fprintf(stderr,"Rows*Cols %d rows %d cols %d\n",rows*cols,rows,cols);
		fprintf(stderr,"ERROR out of memory cannot create matrix in newMatrix!\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
	}
  #endif
  
	//set dimensions
  return mtx;
}

matrix  newMatrixSet(int rows, int cols, double set) {
  #ifdef _ERROR_CHECKING_ON_
  if(rows <= 0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR rows are less than 1 in newMatrixSet\n");
    #endif
    #ifdef _FORCE_HARD_FAIL_
    exit(1);
    #endif
    return NULL;
  }
  if( cols <= 0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR cols are less than 1 in newMatrixSet\n");
    #endif
    #ifdef _FORCE_HARD_FAIL_
    exit(1);
    #endif
    return NULL;
  }
  #endif
	//printf("newMatrixSet\n");
	int i;
	int j;
  
	//allocate a matrix structure
	int rows1 = rows;
	int rowsNotAllocated = rows;
	int Flag = 0;

	matrix m;
	matrix mtx;
	matrix * mtx2;

	while( rowsNotAllocated != 0){

		m = (matrix) malloc(sizeof(struct _matrix)+sizeof(double)*rows1*cols);


		if( m!=NULL ){
			
			//printf("m not null value of rows1 %d\n",rows1);
			if(Flag == 0){
				Flag = 1;
				mtx = m;
				mtx2 = &mtx;
			}else{
				(*mtx2)->Extra = m;
				(*mtx2) = (*mtx2)->Extra;
			}

			m->rows = rows1;
			m->cols = cols;
			m->Extra = NULL;

			for (i = 1; i<=rows1; i++){
				for (j=1;j<=cols;j++) {
					ELEM(m,i,j)=set;
				}
			}

			rowsNotAllocated = rowsNotAllocated-rows1;
			rows1 = rowsNotAllocated;

			m = m->Extra;
		}else{

			rows1 = rows1/2;

		}
	}

  #ifdef _ERROR_CHECKING_ON_
	if(mtx==NULL){
    #ifdef _ERROR_
		fprintf(stderr,"Rows*Cols %d rows %d cols %d\n",rows*cols,rows,cols);
		fprintf(stderr,"ERROR out of memory cannot create matrix!\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
	}
  #endif
	//set dimensions
	return mtx;
}

matrix duplicateMatrix(const_matrix mtx){

  #ifdef _ERROR_CHECKING_ON_	
  if(mtx==NULL){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR matrix NULL in duplicateMatrix!\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	matrix mtxNew = newMatrix(getRows(mtx),getCols(mtx));
//	matrix * mtxNew2;
//	matrix mtxtemp;
//	matrix mtx2;

	if(mtxNew==NULL){
		printf("ERROR no memory to duplicate matrix\n");
		return NULL;
	}
	//set Dimenions

	int i;
	int j;
	double val;
	for( i= 1;i<=getRows(mtx);i++){
		for(j=1;j<=getCols(mtx);j++){
			val = getE(mtx,i,j);
			//printf("i %d j %d val %g\n",i,j,val);
			setE(mtxNew,i,j,val);
		}
	}
/*
	mtxNew->Extra = NULL;
	mtxNew2 = &mtxNew;
	mtx2 = mtx;

	while(mtx2->Extra!=NULL){
		(*mtxNew2) = (*mtxNew2)->Extra;
		mtx2 = mtx2->Extra;

		mtxtemp = (matrix) malloc(sizeof(struct _matrix)+sizeof(double)*(mtx2->rows)*(mtx2->cols));
		if(mtxtemp==NULL){
			printf("ERROR no memory to duplicate matrix\n");
			return NULL;
		}

		(*mtxNew2)=mtxtemp;
		//set Dimenions
		(*mtxNew2)->rows = mtx2->rows;
		(*mtxNew2)->cols = mtx2->cols;

		int i;
		int j;
		double val;
		for( i= 1;i<=(*mtxNew2)->rows;i++){
			for(j=1;j<=(*mtxNew2)->cols;j++){
				val = getE(mtx2,i,j);
				setE((*mtxNew2),i,j,val);
			}
		}

	(*mtxNew2)->Extra = NULL;
	
	}
*/
	return mtxNew;

}

int deleteMatrix( matrix * mtx) {
 
  #ifdef _ERROR_CHECKING_ON_ 
	if (!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in deleteMatrix\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if (!(*mtx)){
     #ifdef _ERROR_
    fprintf(stderr,"ERROR *mtx is NULL in deleteMatrix\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	
  matrix mtx2;
  while( (*mtx)->Extra!=NULL){
			mtx2 = (*mtx);
			(*mtx) = (*mtx)->Extra;
			free(mtx2);
	}
	
	free(*mtx);

	*mtx = NULL;
  return 0;
}

int printMatrix(matrix mtx) {

  #ifdef _ERROR_CHECKING_ON_
	if (!mtx){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR Matrix does not appear to exist in printMatrix!\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif
	int i;
	int j;

	int totalRows = getRows(mtx);
	int totalCols = mtx->cols;

	matrix mtxtemp;

	printf("Rows %d Columns %d\n",totalRows, totalCols);

	mtxtemp = mtx;
	while(mtxtemp!=NULL){
		for (i=1;i<=mtxtemp->rows;i++) {
			for(j=1;j<=mtxtemp->cols;j++){
				printf("%.6g \t",getE(mtxtemp,i,j));			
			}
			printf("\n");
		}
		mtxtemp=mtxtemp->Extra;
	}
	return 0;
}

int printRandomEnergyFilesijk(matrix mtx) {
	#ifdef _ERROR_CHECKING_ON_
  if (!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in printRandomEnergyFilesijk\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if (mtx->cols!=3){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx does not have 3 columns in printRandomEnergyFilesijk\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	printf("\nCols: %d\n",mtx->cols);
	int i;
	int num;
	char *RandomFile="RandomEnergyi.txt";
	FILE * RandOut;

	int totalRows = getRows(mtx);
	matrix mtxtemp;

	if(( RandOut=fopen(RandomFile,"w"))==NULL) {
		printf("Error! unable to open RandomEnergyi.txt\n");
	}else{
		fprintf(RandOut,"%d\n\n",totalRows);
		num=1;
		mtxtemp=mtx;
		while(mtxtemp!=NULL){
			for(i=1;i<=mtxtemp->rows;i++){
				fprintf(RandOut,"%d \t %g\n", num, ELEM(mtxtemp,num,1));
				num++;
			}
			mtxtemp=mtxtemp->Extra;
		}
	}
	fclose(RandOut);

	char *RandomFile2="RandomEnergyj.txt";
	FILE * RandOut2;
	
	if(( RandOut2=fopen(RandomFile2,"w"))==NULL) {
		printf("Error! unable to open RandomEnergyj.txt\n");
	}else{
		fprintf(RandOut2,"%d\n\n",totalRows);
		num=1;
		mtxtemp=mtx;
		while(mtxtemp!=NULL){
			for(i=1;i<=mtxtemp->rows;i++){
				fprintf(RandOut2,"%d \t %g\n", num, ELEM(mtxtemp,num,2));
				num++;
			}
			mtxtemp=mtxtemp->Extra;
		}
	}
	fclose(RandOut2);

	char *RandomFile3="RandomEnergyk.txt";
	FILE * RandOut3;
	
	if(( RandOut3=fopen(RandomFile3,"w"))==NULL) {
		printf("Error! unable to open RandomEnergyj.txt\n");
	}else{
		fprintf(RandOut3,"%d\n\n",totalRows);
		num=1;
		mtxtemp=mtx;
		while(mtxtemp!=NULL){
			for(i=1;i<=mtxtemp->rows;i++){
				fprintf(RandOut3,"%d \t %g\n", num, ELEM(mtxtemp,num,3));
				num++;
			}
		mtxtemp=mtxtemp->Extra;
		}
	}
	fclose(RandOut3); 
	return 0;
}

int resizeRow(matrix * mtx, int Row){

  #ifdef _ERROR_CHECKING_ON_
	if(!(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in resizeRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if(!(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *mtx is NULL in resizeRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if(Row<=0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Row is less than 1 in resizeRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	//Find if we are adding or removing from the matrix
	//printf("ResizeRow\n");	
	//If flag equals 1 means decreasing the size of the
	//matrix if it is 0 it means we are increasing the size
	int i;
	int j;
	int cols = (*mtx)->cols;
	int flag = 0;
	int currentRow = 0;
	int rowResize = Row;
	matrix * mtxtemp;
	matrix * mtxprev;
	mtxprev = NULL;
	mtxtemp = mtx;

	if(Row>((*mtx)->rows+currentRow)){
		while((*mtxtemp)->Extra!=NULL && flag==0){
			currentRow+=(*mtxtemp)->rows;
			rowResize=Row-currentRow;
			(*mtxprev) = (*mtxtemp);
			(*mtxtemp)=(*mtxtemp)->Extra;
			if(Row<=((*mtxtemp)->rows+currentRow)){
				flag=1;
			}
		}
	}else{
		flag=1;
	}

	if(flag==0){
		//Increasing
		//If we are increasing mtxtemp->Extra should be NULL
		if((*mtxtemp)->Extra!=NULL){
			printf("ERROR should not be increasing size of matrix\n");
			return -1;
		}

		//Determine if it is the first matrix in the list
		if((*mtx)->Extra==NULL){
			//First matrix in the list
			matrix temp = duplicateMatrix((*mtx));
			matrix mtxnew = newMatrix(rowResize,cols);
			
			for(i=1;i<=rowResize;i++){
				for(j=1;j<=cols;j++){
					if(i<=temp->rows){
						setE(mtxnew,i,j,getE(temp,i,j));
					}else{
						setE(mtxnew,i,j,0);
					}
				}
			}

			deleteMatrix(mtx);
			deleteMatrix(&temp);
			(*mtx) = mtxnew;

		}else{
			//Not the first matrix in the list
			matrix temp = duplicateMatrix((*mtxtemp));
			matrix mtxnew = newMatrix(rowResize,12);

			for(i=1;i<=rowResize;i++){
				for(j=1;j<=cols;j++){
					if(i<=temp->rows){
						setE(mtxnew,i,j,getE(temp,i,j));
					}else{
						setE(mtxnew,i,j,0);
					}
				}
			}

			deleteMatrix(mtxtemp);
			deleteMatrix(&temp);
			(*mtxprev)->Extra = mtxnew;

		}

	}else{
    //Decreasing
		//mtxtemp is the matrix that is being reduced in size
		//All matrices linked after mtxtemp are to be deleted
		
    if((*mtxtemp)->Extra!=NULL){
			matrix mtxRemove = (*mtxtemp)->Extra;
			deleteMatrix(&mtxRemove);
			(*mtxtemp)->Extra=NULL;
      printf("Part 1\n");
		}
    

		//Now we just need to resize mtxtemp
	  printf("Creating new matrix full of zeros\n");
    matrix mtxnew = newMatrix(rowResize,cols);
		if(!mtxnew) {
			printf("ERROR unable to create new matrix returned NULL\n");
			return -1;
		}
		mtxnew->rows=rowResize;
		mtxnew->cols=cols;
		mtxnew->Extra=NULL;
      
		for (i=1;i<=rowResize;i++){
			for (j=1;j<=cols;j++){
				if(i<=mtxnew->rows){
					setE(mtxnew,i,j,getE(*mtx,i,j));
				}else{
					printf("ERROR should not have a value in the new matrix that is greater than the old one\n");
					return -1;
				}
			}
		}
    printf("Finished setting values equal\n");
		deleteMatrix(mtx);
		(*mtx) = mtxnew;


	}

	return 0;
}

/*
int resizeRow(matrix * mtx, int Row){

	if(!(*mtx)) return -1;
	if(Row<=0 ) return -1;

	int i;
	int j;
	int cols = (*mtx)->cols;
	int rows = (*mtx)->rows;
	int rowNew;
	int flag = 0;
	int totalRows = getRows(*mtx);
	int currentRow = 0;
	matrix mtxtemp = (*mtx);
	matrix mtxprev = NULL;
	matrix mtxRemove;

	if((*mtx)->Extra!=NULL){
		while( mtxtemp!=NULL && flag==0 ){

			if(Row>(mtxtemp->rows+currentRow)){
				currentRow = mtxtemp->rows;
				mtxprev = mtxtemp;
				mtxtemp=mtxtemp->Extra;
			}else{
				//This means the row we are resizing too
				//is smaller than the total size of the matrices
				flag=1;
			}
		}

		if(flag==1){
			//Go ahead and delete the excess matrices
			mtxRemove = mtxtemp->Extra;
			if(mtxRemove!=NULL){
				deleteMatrix(&mtxRemove);
				mtxtemp->Extra=NULL;
			}
		}
	}
	rowNew = Row-currentRow;

	if(mtxprev==NULL){
		//This means linked matrices were not used we are only dealing with one matrix
		//can therefore solve the problem as if it is only one matrix
		matrix temp = duplicateMatrix(*mtx); 
		matrix mtxnew = (matrix) realloc(*mtx, sizeof(struct _matrix)+sizeof(double)*Row*cols);
		if(!mtxnew) {
			printf("ERROR unable to realloc matrix returned NULL\n");
			return -1;
		}
		mtxnew->rows=Row;
		mtxnew->cols=cols;

		for (i=1;i<=Row;i++){
			for (j=1;j<=cols;j++){
				if(i<=rows){
					setE(mtxnew,i,j,getE(temp,i,j));
				}else{
					setE(mtxnew,i,j,0);
				}
			}
		}

		deleteMatrix(&temp);
		*mtx = mtxnew;

	}else{
	
		//If flag is 1 we are simply decreasing the size of matrix mtxtemp
		if(flag==1){
			matrix temp = duplicateMatrix(mtxtemp);
			matrix mtxnew = (matrix) realloc(mtxtemp,sizeof(struct _matrix)+sizeof(double)*rowNew*cols);
			if(!mtxnew){
				printf("ERROR unable to reallocate matrix returned NULL\n");
				return -1;
			}
			mtxnew->rows=rowNew;
			mtxnew->cols=cols;
			
			for (i=1;i<=rowNew;i++){
				for (j=1;j<=cols;j++){
					if(i<=temp->rows){
						setE(mtxnew,i,j,getE(temp,i,j));
					}else{
						//This if statement should not be triggered
						printf("ERROR trying to assign excess values to reallocated matrix\n");
						return -1;
					}
				}
			}	
			
			deleteMatrix(&temp);
			mtxnew->Extra = NULL;
			mtxprev->Extra = mtxnew;
			

		}else{
	
			//This is if we are reallocated the matrix to a larger size
			matrix temp = duplicateMatrix(mtxtemp);
			matrix mtxnew = newMatrix(rowNew,cols);
			
			if(!mtxnew){
				printf("ERROR unable to reallocate matrix returned NULL\n");
				return -1;
			}

			for (i=1;i<=rowNew;i++){
				for(j=1;j<=cols;j++){
					if(i<=temp->rows){
						setE(mtxnew,i,j,getE(temp,i,j));
					}else{
						setE(mtxnew,i,j,0);
					}
				}
			}

			deleteMatrix(&temp);
			mtxprev->Extra = mtxnew;
		}
	}

	return 0;
}
*/

int setAll(matrix mtx, double val){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in setAll\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
 	//printf("setAll\n"); 
	//set all data to val
  int i;
	int j;
	while(mtx!=NULL){
		for (i = 1; i<=mtx->rows; i++){
			for (j=1;j<=mtx->cols;j++) {
				ELEM(mtx,i,j)=val;
			}
		}
		mtx=mtx->Extra;
	}
	return 0;
}

int setAllRowsInCol(matrix mtx, int col, double val){

  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in setAllRowsInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if(col<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in setAllRowsInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(col>mtx->cols){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than mtx->cols in setAllRowsInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
	matrix mtxtemp;
	int i;
	mtxtemp=mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows;i++){
			ELEM(mtxtemp,i,col)=val;
		}
		mtxtemp=mtxtemp->Extra;
	}
	return 0;
}

int setE(matrix mtx, int row, int col, double val) {
  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in setE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }	
  if(row>getRows(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is greater than getRows(mtx) in setE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(row<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is less than 1 in setE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(col>getCols(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in setE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(col<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in setE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif

  matrix mtxtemp;
	matrix mtxprev;
	int rowtemp;
	int currentRow = 0;
	mtxtemp = mtx;
	mtxprev = mtx;
	while(row>(currentRow+mtxtemp->rows)){
		currentRow+=mtxtemp->rows;
		mtxprev = mtxtemp;
		mtxtemp = mtxtemp->Extra;
		if(mtxtemp==NULL){
			printf("ERROR exceeded size of matrix without finding the correct row!\n");
			return -1;
		}
	}
	
	rowtemp = row-currentRow;
  ELEM(mtxprev, rowtemp, col) = val;
  return 0;
}

int getRows(const_matrix mtx){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in getRows\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif 
    return -1;
  }

  #endif
	int totalRows = 0;
	const_matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		totalRows+=mtxtemp->rows;
		mtxtemp=mtxtemp->Extra;
	}
	return totalRows;
}

int getCols(const_matrix mtx){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in getCols\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	return mtx->cols;
}

double getE(const_matrix mtx, int row, int col) {
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in getE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(row>getRows(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is greater than getRows(mtx) in getE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(row<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is less than 1 in getE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(col>getCols(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in getE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(col<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in getE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	#endif
  int rowtemp;
	int totalRows = getRows(mtx);
	int totalCols = mtx->cols;
	if(row <= 0 || row > totalRows ||
			col <= 0 || col > totalCols){
		return -2;
	}
	
	int currentRow = 0;
	const_matrix mtxtemp = mtx;
	const_matrix mtxprev = mtx;
	while(row>(currentRow+mtxtemp->rows)){
		currentRow+=mtxtemp->rows;
		mtxprev = mtxtemp;
		mtxtemp = mtxtemp->Extra;

		if(mtxtemp==NULL){
			printf("ERROR exceeded size of matrix without finding the correct row!\n");
			return -1;
		}
	}
	
	rowtemp = row-currentRow;

	return (double) ELEM(mtxprev,rowtemp,col);
}

int matchExist(const_matrix mtx, double match){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in matchExist\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	#endif
  int i;
	int j;
	const_matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows; i++){
			for(j=1;j<=mtxtemp->cols;j++){

				if(ELEM(mtxtemp,i,j)==match){
					return 1;
				}
			}
		}
		mtxtemp=mtxtemp->Extra;
	}
	return 0;
}

int matchExistCol(const_matrix mtx, int col, double match){

  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in matchExistCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(col>getCols(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in matchExistCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(col<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in matchExistCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	#endif

	int i;
	const_matrix mtxtemp=mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows; i++){
			if(ELEM(mtxtemp,i,col)==match){
				return i;
			}
		}
	mtxtemp=mtxtemp->Extra;
	}
	return 0;
}

int matchReplace(matrix mtx, double match, double replace){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in matchReplace\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	
  int i;
	int j;
	
	matrix * mtxtemp = &mtx;
	while((*mtxtemp)!=NULL){
		for(i=1;i<=(*mtxtemp)->rows; i++){
			for(j=1;j<=(*mtxtemp)->cols;j++){

				if(ELEM((*mtxtemp),i,j)==match){
					setE((*mtxtemp),i,j,replace);
				}
			}
		}
		(*mtxtemp)=(*mtxtemp)->Extra;
	}
	return 0;
}

int FindRowOfMatchInCol(const_matrix mtx, double match, int col){
	
  #ifdef _ERROR_CHECKING_ON_
  if(!mtx){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in FindRowOfMatchInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  if(col<1){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in FindRowOfMatchInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  if(col>getCols(mtx)){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in FindRowOfMatchInCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif
	int i;
	int matchRow = 1;
	const_matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows; i++){
			if(ELEM(mtxtemp,i,col)==match){
				return matchRow;
			}
			matchRow++;
		}
		mtxtemp=mtxtemp->Extra;
	}

	return -1; 
}

double SumOfRow(matrix mtx, int row){
	
  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in SumOfRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(row<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is less than 1 in SumOfRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(row>getRows(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is greater than getRows(mtx) in SumOfRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
	
	int currentRow = 0;
	int rowtemp;
	int j;
	double val=0;
	matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		if(row<(currentRow+mtxtemp->rows)){
			rowtemp = row-currentRow;
			for(j=1;j<=mtxtemp->cols;j++){
				val+=ELEM(mtxtemp,rowtemp,j);
			}
		}else{
			currentRow+=mtxtemp->rows;
		}
		mtxtemp=mtxtemp->Extra;
	}
	return val;
}

double SumOfCol(matrix mtx, int col){
  
  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in SumOfCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(col<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in SumOfCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(col>getCols(mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in SumOfCols\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
	
	int i;
	double val=0;
	matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows; i++){
			val+=ELEM(mtxtemp,i,col);
		}
		mtxtemp=mtxtemp->Extra;
	}
	return val;
}

int DivideEachElement(matrix * mtx, double val){
	
  #ifdef _ERROR_CHECKING_ON_
	if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in DivideEachElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(!(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *mtx is NULL in DivideEachElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	int i;
	int j;

	for(i=1;i<=(*mtx)->rows;i++){
		for(j=1;j<=(*mtx)->cols;j++){
			setE((*mtx),i,j,ELEM((*mtx),i,j)/val);
		}
	}

	matrix mtxtemp = (*mtx)->Extra;
	while((mtxtemp)!=NULL){
		for(i=1;i<=(mtxtemp)->rows; i++){
			for(j=1;j<=(mtxtemp)->cols;j++){
				setE((mtxtemp),i,j,ELEM((mtxtemp),i,j)/val);					
			}
		}
		(mtxtemp)=(mtxtemp)->Extra;
	}

	if((*mtx)==NULL){
		printf("Error mtx is NULL\n");
		exit(1);
	}

	return 0;
}

int DivideEachElementCol(matrix * mtx, int col, double val){
	
  #ifdef _ERROR_CHECKING_ON_
  if(mtx==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in DivideEachElementCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
	if(!(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *mtx is NULL in DivideEachElementCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(col<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in DivideEachElementCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(col>getCols(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in DivideEachElementCol\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif

	int i;
	matrix mtxtemp = (*mtx);
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows;i++){
			setE(mtxtemp,i,col,ELEM(mtxtemp,i,col)/val);	
		}
		mtxtemp=mtxtemp->Extra;
	}

	return 0;
}

int incrementE(matrix * mtx, int row, int col){
	#ifdef _ERROR_CHECKING_ON_
  if(!mtx){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mtx is NULL in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(!(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *mtx is NULL in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}

  if(row>getRows(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is greater than getRows(mtx) in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(row<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR row is less than 1 in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(col>getCols(*mtx)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is greater than getCols(mtx) in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(col<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR col is less than 1 in incrementE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	#endif

  matrix mtxtemp;
	matrix mtxprev;
	int rowtemp;

	int currentRow = 0;
	mtxtemp = *mtx;
	mtxprev = *mtx;
	while(row>(currentRow+mtxtemp->rows)){
		currentRow+=mtxtemp->rows;
		mtxprev = mtxtemp;
		mtxtemp = mtxtemp->Extra;
		if(mtxtemp==NULL){
			printf("ERROR exceeded size of matrix without finding the correct row!\n");
			return -1;
		}
	}
	
	rowtemp = row-currentRow;
  ELEM(mtxprev, rowtemp, col) = ELEM(mtxprev, rowtemp, col)+1;
  return 0;

}
