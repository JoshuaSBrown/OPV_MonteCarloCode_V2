#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
//#include "../../../MEM/mem.h"

#define ELEM(mtx, row, col) \
  mtx->data[(col-1) * mtx->rows + (row-1)]

struct _matrix {
	int rows;
	int cols;
	matrix Extra;
	double data[0];
};


matrix  newMatrix(int rows, int cols) {
  if(rows <= 0 || cols <= 0 ) return NULL;

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

	if(mtx==NULL){
		printf("Rows*Cols %d rows %d cols %d\n",rows*cols,rows,cols);
		printf("ERROR out of memory cannot create matrix!\n");
		return NULL;
	}
  
	//set dimensions
  return mtx;
}

matrix  newMatrixSet(int rows, int cols, double set) {
 if(rows <= 0 || cols <= 0 ) return NULL;
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

	if(mtx==NULL){
		printf("Rows*Cols %d rows %d cols %d\n",rows*cols,rows,cols);
		printf("ERROR out of memory cannot create matrix!\n");
		return NULL;
	}
  
	//set dimensions

	return mtx;
}

matrix duplicateMatrix(const_matrix mtx){
	//printf("duplicateMatrix\n");
	if(mtx==NULL){
		printf("ERROR matrix NULL!\n");
		return NULL;
	}

	matrix mtxNew = newMatrix(getRows(mtx),getCols(mtx));
//	matrix * mtxNew2;
//	matrix mtxtemp;
//	matrix mtx2;

	if(mtxNew==NULL){
		//printf("rows %d cols %d\n",getRows(mtx),getCols(mtx));
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
  
	if (!mtx) return -1;
	if (!(*mtx)) return -1;
	//printf("deleteMatrix\n");
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

int printMatrix(const_matrix mtx) {
	if (!mtx){
		printf("Matrix does not appear to exist!\n");
		return -1;
	}
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

int printRandomEnergyFilesijk(const_matrix mtx) {
	if (!mtx || mtx->cols!=3) return -1; //null test seg faults
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

	if(!(*mtx)) return -1;
	if(Row<=0 ) return -1;
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
		}

		//Now we just need to resize mtxtemp
		matrix temp = duplicateMatrix((*mtxtemp)); 
		matrix mtxnew = (matrix) realloc((*mtxtemp), sizeof(struct _matrix)+sizeof(double)*rowResize*cols);
		if(!mtxnew) {
			printf("ERROR unable to realloc matrix returned NULL\n");
			return -1;
		}
		mtxnew->rows=rowResize;
		mtxnew->cols=cols;
		mtxnew->Extra=NULL;

		for (i=1;i<=rowResize;i++){
			for (j=1;j<=cols;j++){
				if(i<=temp->rows){
					setE(mtxnew,i,j,getE(temp,i,j));
				}else{
					printf("ERROR should not have a value in the new matrix that is greater than the old one\n");
					return -1;
				}
			}
		}

		deleteMatrix(&temp);
		(*mtxtemp) = mtxnew;


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
	if(!mtx) return -1;
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

	if(!mtx) return -1;

	matrix mtxtemp;

	if(col<=0 || col>mtx->cols)
		return -1;

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
  if(!mtx) {
		return -1;
	}
	//printf("SetE\n");
	matrix mtxtemp;
	matrix mtxprev;
	int rowtemp;
	int totalCols = mtx->cols;
	int totalRows = getRows(mtx);

  if (row <= 0 || row > totalRows ||
      col <= 0 || col > totalCols ){
    
		return -2;
	}

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
	if(!mtx) return -1;

	int totalRows = 0;
	matrix mtxtemp = mtx;
	while(mtxtemp!=NULL){
		totalRows+=mtxtemp->rows;
		mtxtemp=mtxtemp->Extra;
	}
	return totalRows;
}

int getCols(const_matrix mtx){
	if(!mtx) return -1;
	return mtx->cols;
}

double getE(const_matrix mtx, int row, int col) {
	if(!mtx){
		printf("ERROR matrix NULL!\n");
		return -1;
	}
	//printf("getE\n");
	int rowtemp;
	int totalRows = getRows(mtx);
	int totalCols = mtx->cols;
	if(row <= 0 || row > totalRows ||
			col <= 0 || col > totalCols){
		return -2;
	}
	
	int currentRow = 0;
	matrix mtxtemp = mtx;
	matrix mtxprev = mtx;
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
	if(!mtx){
		return -1;
	}
	//printf("matchExist\n");
	int i;
	int j;
	matrix mtxtemp = mtx;
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

	if(!mtx) return -1;
	
	if(col<1 || col>mtx->cols){
		return -1;
	}

	//printf("matchExistCol\n");

	int i;
	matrix mtxtemp=mtx;
	while(mtxtemp!=NULL){
		for(i=1;i<=mtxtemp->rows; i++){
			if(ELEM(mtxtemp,i,col)==match){
				return i;
			}
		}
	mtxtemp=mtxtemp->Extra;
	}
	return 0; //why not a neg, so we know it doesn't exist?
}

int matchReplace(matrix mtx, double match, double replace){
	if(!mtx){
		return -1;
	}
	//printf("matchReplace\n");
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
	if(!mtx || col<1 || col>mtx->cols){
		return -1;
	}
	//printf("FindRowOfMatchInCol\n");

	int i;
	int matchRow = 1;
	matrix mtxtemp = mtx;
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
	
	
	if(!mtx){
		return -1;
	}
	//printf("SumOfRow\n");
	int totalRows = getRows(mtx);
	if( row<1 || row>totalRows){
		return -1;
	}
	
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
	if(!mtx || col<1 || col>mtx->cols)
		return -1;
	
	//printf("SumOfCol\n");
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
	if(mtx==NULL)
		return -1;
	
	if(!(*mtx))
		return -1;
	//printf("DivideEachElement\n");	

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
	if(mtx==NULL){
		return -1;
	}
	//printf("DivideEachElementCol\n");	
	if(!(*mtx)){
		return -1;
	}

	if( col<1 || col>(*mtx)->cols){
		return -1;
	}


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
