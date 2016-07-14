#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "linklist2.h"

int main(){ //couldn't test blankLLlist
	
	printf("Starting\n");

	printf("Testing: newLLNode\n");
	LLNode2 LLrv=newLLNode2(-5);
	assert(LLrv!=NULL);
	LLrv=newLLNode2(4);
	assert(LLrv!=NULL);
	
	printf("Testing: newLinkList\n");
	linklist2 LL = newLinkList2(-1);
	assert(LL!=NULL);
	LL = newLinkList2(0);
	assert(LL!=NULL);


	int rv;
	printf("Testing:addLLNode2\n");
	rv = addLLNode2(NULL,1);
	assert(rv==-1);
	rv = addLLNode2(NULL,-5);
	assert(rv==-1);
	rv = addLLNode2(LL, 0);
	assert(rv!=-1); 
	rv = addLLNode2(LL, -1);
	assert(rv!=-1);
	rv = addLLNode2(LL, 34);
	assert(rv==0);
	printf("Expect 0 0 -1 and 34: ");
	printLL2(LL);
	rv = addLLNode2(LL, 34);
	assert(rv!=-1); 
	rv = addLLNode2(LL,54);
	rv=addLLNode2(LL,64);
	assert(rv==0);
	printf("Expect 0 0 -1 34, 34, 54, 64:");
	printLL2(LL);
	getLLlength2(LL);

	printf("Testing:removeLLNode2\n");
	rv = removeLLNode2(NULL, 1);
	assert(rv==-1);
	rv=removeLLNode2(LL, -5);
	assert(rv==-1);
	rv = removeLLNode2(LL, 34);
	assert(rv==-1);
	printf("Expect 0 0 -1 34 34 54 and 64: ");
	printLL2(LL);
	rv = removeLLNode2(LL,0);
	assert(rv==-1);
  rv = removeLLNode2(LL,1);
	printf("Expect 0 -1 34 34 54 and 64: ");
	printLL2(LL);
 
  printf("Testing:removeLLNodes2\n");
  rv = removeLLNodes2(LL,3,4);
	assert(rv==0);
  printf("Expect 0 -1 54 and 64: ");
	printLL2(LL);
 
  	
  printf("Testing:deleteLL\n");
	rv=deleteLL2(NULL);
	assert(rv==-1);
	rv = deleteLL2(&LL);
	assert(rv==0);


	printf("Testing: deleteLLNode2");
	LLNode2 nrv = newLLNode2(25);
	rv=deleteLLNode2(NULL);
	assert(rv==-1);
	rv=deleteLLNode2(nrv);
	assert(rv==0);
	printf("Expect Nothing: ");
	printLL2(LL);
	linklist2 LL2=newLinkList2(25); //recreating deleted linklist2

	printf("Testing: getLLlength");
	addLLNode2(LL, 74);
	rv=getLLlength2(NULL);
	assert(rv==-1);
	printf("Expect 1: "); 
	rv=getLLlength2(LL2);
	printLL2(LL2);
	assert(rv==1);

	printf("Testing: printLL");
	rv=printLL2(NULL);
	assert(rv==-1);
	printf("Expect 13 and 74: ");
	rv=printLL2(LL2); //previously tested
	assert(rv==0);
	
	addLLNode2(LL, 15);
	printf("Testing: getLLstartID\n");
	rv=getLLstartID2(NULL);
	assert(rv==-1);
	rv=getLLstartID2(LL2);
	assert(rv==25);

	printf("End");

	return 0;
}
