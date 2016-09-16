#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "midpoint.h"

int main() {
	int nei1, nei2, order;
	int rv;
	nei1=22;
	nei2=74;
	order=-2;
	MidPoint mp;
  MidPoint mp2;
  MidPoint mp3;
  MidPoint mp4;
  MidPoint mp5;
  MidPoint mp6;

	printf("Testing:newMidPoint\n");
	mp = newMidPoint(order, -1, nei1, nei2);
	assert(mp==NULL);
	mp = newMidPoint(order, 0, nei1, nei1);
	assert(mp==NULL);
	mp = newMidPoint(order, 0, -1, nei2);
	assert(mp==NULL);
	mp = newMidPoint(order, 0, nei1, -1);
	assert(mp==NULL);
	mp=newMidPoint(order,0,5,5);
	assert(mp==NULL);
	mp = newMidPoint(order, 0, nei1, nei2);
	assert(mp!=NULL);

	printf("Testing:deleteMidPoint\n");
	rv = deleteMidPoint(NULL);
	assert(rv==-1);
	rv = deleteMidPoint(mp);
	assert(rv==0);

	mp = newMidPoint(order, 0, nei1, nei2);
	assert(mp!=NULL);
	printf("Testing:printMP\n");
	printf("Should print Order of -2, nei1 22, nei2 74, ID 0\n");
	rv = printMP(NULL);
	assert(rv==-1);
	rv = printMP(mp);
	assert(rv==0);

	printf("Testing:getMP_order\n"); //added condition
	rv = getMP_order(NULL);
	assert(rv==-1);
	rv = getMP_order(mp);
  assert(rv==-2);
	printf("Output should be -2\n%d\n",rv);

	printf("Testing:getMP_id\n");
	rv = getMP_id(NULL);
	assert(rv==-1);
	rv = getMP_id(mp);
	assert(rv==0);
	printf("Output should be 0\n%d\n",rv);

	printf("Testing:setMP_id\n");
	rv = setMP_id(mp,3);
	assert(rv==0);
	rv = setMP_id(NULL,2);
	assert(rv==-1);
	rv = setMP_id(mp,-1);
	assert(rv==-1);
	rv = getMP_id(mp);
	assert(rv==3);

	printf("Testing:CompareNeiMidPoint\n");
	mp2=newMidPoint(2,2,15,nei1);
	mp3=newMidPoint(12,4,nei1,122);
	mp4=newMidPoint(14,1,11,nei2);
	mp5=newMidPoint(16,5,nei2,19);
	mp6=newMidPoint(18,12,21,212);
	rv=CompareNeiMidPoint(NULL,mp);
	assert(rv==-1);
	rv=CompareNeiMidPoint(mp,NULL);
	assert(rv==-1);
	rv=CompareNeiMidPoint(mp,mp2);
	assert(rv==1);
	rv=CompareNeiMidPoint(mp,mp3);
	assert(rv==1);
	rv=CompareNeiMidPoint(mp,mp4);
	assert(rv==1);
	rv=CompareNeiMidPoint(mp,mp5);
	assert(rv==1);
	rv=CompareNeiMidPoint(mp,mp6);
	assert(rv==0);

  deleteMidPoint(mp);
  deleteMidPoint(mp2);
  deleteMidPoint(mp3);
  deleteMidPoint(mp4);
  deleteMidPoint(mp5);
  deleteMidPoint(mp6);
  return 0;
}
