#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "matrix.h"

int main(){


	//matrix mtx = newMatrix(10000000,12);
	//assert(mtx!=NULL);
	//sleep(3);
	/*matrix mtx2 = newMatrix(5000000,12);
	assert(mtx2!=NULL);
	sleep(3);
	matrix mtx3 = newMatrix(5000000,12);
	assert(mtx3!=NULL);
	sleep(3);
	*/
	//deleteMatrix(&mtx);
	//deleteMatrix(&mtx2);
	//deleteMatrix(&mtx3);
	int rv;
	double val;

	matrix mtx = newMatrix(10000000,12);
	assert(mtx!=NULL);
	
	rv = setE(mtx,1000,1,5);
	assert(rv==0);
	rv = setE(mtx,6000000,1,5);
	assert(rv==0);
	
	val = getE(mtx,1000,1);
	assert(val==5);
	val = getE(mtx,6000000,1);
	assert(val==5);
	

	deleteMatrix(&mtx);


	return 0;

}
