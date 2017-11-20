#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "arbarray.h"
#include "../MIDPOINT/midpoint.h"

int main() {
	int rv;
  void * Temp;
  void * Temp2;
//////////////////////////////////////////////////////////////////////
	printf("Testing:newArbArray\n");
	ArbArray Arb = newArbArray(-1, 0);
	assert(Arb==NULL);
	ArbArray Arb2 = newArbArray(3, -1);
	assert(Arb2==NULL);
	ArbArray Arb3 = newArbArray(3, 3);
	assert(Arb3==NULL);
	ArbArray Arb4 = newArbArray(4,1);
	assert(Arb4!=NULL);

	printf("Testing:ArbArrayCheck\n");
	rv = ArbArrayCheck(NULL);
	assert(rv==-1);
	rv = ArbArrayCheck(Arb4);
	assert(rv==0);

	printf("Testing:deleteArbArray\n");
	rv = deleteArbArray(NULL);
	assert(rv==-1);
	rv = deleteArbArray(&Arb4);
	assert(rv==0);
	ArbArray Arb5 = newArbArray(5,0);
	rv = deleteArbArray(&Arb5);
	assert(rv==0);
	ArbArray Arb6 = newArbArray(5,2);
	rv = deleteArbArray(&Arb6);
	assert(rv==0);

	Arb6 = newArbArray(5,2);
	printf("Testing:getArbElement\n");
	Temp = getArbElement(NULL,2);
	assert(Temp==NULL);
	Temp = getArbElement(Arb6,-1);
	assert(Temp==NULL);
	Temp = getArbElement(Arb6,5);
	assert(Temp==NULL);
	Temp = getArbElement(Arb6,4);
	assert(Temp==NULL);

	printf("Testing:setArbElement\n");
	Temp2 = (void *) newMidPoint(-5, 1, 99, 9);
	rv = setArbElement(NULL,1, Temp2);
	assert(rv==-1);
	rv = setArbElement(Arb6,5, Temp2);
	assert(rv==-1);
	rv = setArbElement(Arb6,-1, Temp2);
	assert(rv==-1);
	rv = setArbElement(Arb6,1, Temp);
	assert(rv==-1);
	rv = setArbElement(Arb6,1, Temp2);
	assert(rv==0);
	rv = setArbElement(Arb6,1, Temp2);
	assert(rv==0);
	void * Temp3 = getArbElement(Arb6, 1);
	MidPoint mp34 = (MidPoint) Temp3;
	assert(getMP_id(mp34)==1);
	
	deleteMidPoint((MidPoint) Temp2);
	deleteArbArray(&Arb);
	deleteArbArray(&Arb2);
	deleteArbArray(&Arb3);
	deleteArbArray(&Arb6);
	
	//////////////////////////////////////////////////////////////////////
	Arb6 = newArbArray(5,2);
	Temp2 = (void *) newMidPoint(-5, 1, 99, 9);
	rv = setArbElement(Arb6,1, Temp2);
	rv = setArbElement(Arb6,1, Temp2);
	void * Temp4 = (void *) newMidPoint(-5, 1, 99, 9);
	
	printf("Testing:printArbArray\n");
	rv = printArbArray(NULL);
	assert(rv==-1);
	ArbArray Arb7 = newArbArray(9,2);
	rv = printArbArray(Arb7);
	assert(rv==0);
	
	printf("Testing:getElementsUsed\n");
	rv = getElementsUsed(NULL);
	assert(rv==-1);
	rv = getElementsUsed(Arb7);
	assert(rv==0);
	rv = setArbElement(Arb7,0,Temp4);
	assert(rv==0);
	rv = getElementsUsed(Arb7);
	assert(rv==1);

	printf("Testing:getElementsReserved\n");
	rv = getElementsReserved(NULL);
	assert(rv==-1);
	rv = getElementsReserved(Arb7);
	assert(rv==9);

	printf("Testing:removeArbElement\n");
	rv = removeArbElement(NULL, 0);
	assert(rv==-1);
	rv = removeArbElement(Arb7, -1);
	assert(rv==-1);
	rv = removeArbElement(Arb7, 9);
	assert(rv==-1);
	rv = removeArbElement(Arb7, 1);
	assert(rv==-1);
	rv = removeArbElement(Arb7, 0);
	assert(rv==0);
	rv = getElementsUsed(Arb7);
	assert(rv==0);

	printf("Testing:getMP\n");
	MidPoint mp101;
	mp101 = getMP(Arb7,-1);
	assert(mp101==NULL);
	mp101 = getMP(NULL,0);
	assert(mp101==NULL);
	mp101 = getMP(Arb7,0);
	assert(mp101==NULL);
	mp101 = getMP(Arb7,9);
	assert(mp101==NULL);
	rv = setArbElement(Arb7,5,Temp4);
	assert(rv==0);
	mp101 = getMP(Arb7,5);
	assert(mp101!=NULL);
	assert(getMP_id(mp101)==1);
	mp101 = getMP(Arb6,0);
	assert(mp101==NULL);

	printf("Testing:getMPnei1\n");
	rv = getMPnei1(NULL, 0);
	assert(rv==-1);
	rv = getMPnei1(Arb7, 9);
	assert(rv==-1);
	rv = getMPnei1(Arb7, -1);
	assert(rv==-1);
	rv = getMPnei1(Arb7,0);
	assert(rv==-1);
	rv = getMPnei1(Arb7,5);
	assert(rv==99);
	rv = getMPnei1(Arb6,3);
	assert(rv==-1);

	printf("Testing:getMPnei2\n");
	rv = getMPnei2(NULL, 0);
	assert(rv==-1);
	rv = getMPnei2(Arb7, 9);
	assert(rv==-1);
	rv = getMPnei2(Arb7, -1);
	assert(rv==-1);
	rv = getMPnei2(Arb7,0);
	assert(rv==-1);
	rv = getMPnei2(Arb7,5);
	assert(rv==9);
	rv = getMPnei2(Arb6,3);
	assert(rv==-1);

	printf("Testing:getMPOrder\n");
	rv = getMPOrder(NULL, 0);
	assert(rv==-1);
	rv = getMPOrder(Arb7, 9);
	assert(rv==-1);
	rv = getMPOrder(Arb7, -1);
	assert(rv==-1);
	rv = getMPOrder(Arb7,0);
	assert(rv==-1);
	rv = getMPOrder(Arb7,5);
	assert(rv==-5);
	rv = getMPOrder(Arb6,3);
	assert(rv==-1);
	
	deleteMidPoint((MidPoint)Temp2);
	deleteMidPoint((MidPoint)Temp4);
	deleteArbArray(&Arb6);
	deleteArbArray(&Arb7);
	//////////////////////////////////////////////////////////////////////
	Arb7 = newArbArray(9,2);
	Temp4 = (void *) newMidPoint(-5, 1, 99, 9);
	rv = setArbElement(Arb7,0,Temp4);
	rv = setArbElement(Arb7,5,Temp4);
	
	mp101 = getMP(Arb7,5);
	assert(mp101!=NULL);
	ArbArray Arb9 = newArbArray(5,0);
	printf("Testing:addToOrLL\n");
	rv = addToOrLL(NULL,1, &mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, -1, &mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 5, &mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 0, NULL);
	assert(rv==-1);
	rv = addToOrLL(Arb7, 0, &mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 0, &mp101);
	assert(rv==0);
	rv = getElementsUsed(Arb9);
	assert(rv==1);
	OrderMagLL OMLL11 = getOrderLL(Arb9, 0);
	assert(getOMLL_size(OMLL11)==1);
	assert(mp101!=NULL);
  rv = addToOrLL(Arb9, 0, &mp101);
	assert(rv==-2);
	MidPoint mp30 = newMidPoint(-3,121,32,55);
	rv = addToOrLL(Arb9, 0, &mp30);
	assert(rv==0);
	OrderMagLL OMLL12 = getOrderLL(Arb9,0);
	assert(getOMLL_size(OMLL12)==2);
	rv = getElementsUsed(Arb9);
	assert(rv==1);

	printf("Testing:setDefaultArbElement\n");
	rv = setDefaultArbElem(NULL, 1, 3);
	assert(rv==-1);
	rv = setDefaultArbElem(Arb9, -1, 3);
	assert(rv==-1);
	rv = setDefaultArbElem(Arb9, 5, 3);
	assert(rv==-1);
	rv = setDefaultArbElem(Arb7, 1, 3);
	assert(rv==-1);
	rv = setDefaultArbElem(Arb9, 1, 3);
	assert(rv==0);
	rv = getElementsUsed(Arb9);
	assert(rv==2);
	rv = setDefaultArbElem(Arb9, 0, 3);
	assert(rv==0);
	
  rv = getElementsUsed(Arb9);
	assert(rv==2);

  printf("Clearing heap\n");	
	deleteMidPoint((MidPoint)Temp4);
	deleteMidPoint(mp30);
	deleteArbArray(&Arb7);
	deleteArbArray(&Arb9);
  return 0;
}
