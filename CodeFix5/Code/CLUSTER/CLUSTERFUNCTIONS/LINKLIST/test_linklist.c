#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "linklist.h"

int main(){ //couldn't test blankLLlist
	
	printf("Starting\n");

	printf("Testing: newLLNode\n");
	LLNode LLrv=newLLNode(-5);
	assert(LLrv==NULL);
	LLrv=newLLNode(4);
	assert(LLrv!=NULL);
	
	printf("Testing: newLinkList\n");
	linklist LL = newLinkList(-1);
	assert(LL==NULL);
	LL = newLinkList(0);
	assert(LL!=NULL);


	int rv;
	printf("Testing:addLLNode\n");
	rv = addLLNode(NULL,1);
	assert(rv==-1);
	rv = addLLNode(NULL,-5);
	assert(rv==-1);
	rv = addLLNode(LL, 0);
	assert(rv==-1); //not less than 0?
	rv = addLLNode(LL, -1);
	assert(rv==-1);
	rv = addLLNode(LL, 34);
	assert(rv==0);
	printf("Expect 34: ");
	printLL(LL);
	rv = addLLNode(LL, 34);
	assert(rv==-1); 
	rv = addLLNode(LL,54);
	rv=addLLNode(LL,64);
	assert(rv==0);
	printf("Expect 34, 54, 64:");
	printLL(LL);
	getLLlength(LL);

	printf("Testing:removeLLNode\n");
	rv = removeLLNode(NULL, 1);
	assert(rv==-1);
	rv=removeLLNode(LL, -5);
	assert(rv==-1);
	rv = removeLLNode(LL, 34);
	assert(rv==0);
	printf("Expect 54 and 64: ");
	printLL(LL);
	rv = removeLLNode(LL,0);
	assert(rv==0); //first node is zero
	
	printf("Testing:deleteLL\n");
	rv=deleteLL(NULL);
	assert(rv==-1);
	rv = deleteLL(&LL);
	assert(rv==0);


	printf("Testing: deleteLLNode");
	LLNode nrv = newLLNode(25);
	rv=deleteLLNode(NULL);
	assert(rv==-1);
	rv=deleteLLNode(nrv);
	assert(rv==0);
	printf("Expect Nothing: ");
	printLL(LL);
	linklist LL2=newLinkList(25); //recreating deleted linklist

	printf("Testing: getLLlength");
	addLLNode(LL, 74);
	rv=getLLlength(NULL);
	assert(rv==-1);
	printf("Expect 1: "); 
	rv=getLLlength(LL2);
	printLL(LL2);
	assert(rv==1);

	printf("Testing: printLL");
	rv=printLL(NULL);
	assert(rv==-1);
	printf("Expect 13 and 74: ");
	rv=printLL(LL2); //previously tested
	assert(rv==0);
	
	addLLNode(LL, 15);
	printf("Testing: getLLstartID\n");
	rv=getLLstartID(NULL);
	assert(rv==-1);
	rv=getLLstartID(LL2);
	assert(rv==25);

	printf("End");

	return 0;
}
