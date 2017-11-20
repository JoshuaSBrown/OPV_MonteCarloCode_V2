#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "neighll.h"
#include "../NEIGHNODE/neighnode.h"

int main() {
	int rv;
	double rvd;
	NeighNode Nei;
	NeighNode Nei2;
  NeighNode Nei3;
  NeighLL neighLL;
  Hop h;
	//////////////////////////////////////////////////////////////////////
	printf("Testing:newNeighNode\n");
	Nei = newNeighNode(-1);
	assert(Nei==NULL);
	Nei = newNeighNode(3);
	assert(Nei!=NULL);

	printf("Testing:deleteNeighNode\n");
	rv = deleteNeighNode(NULL);
	assert(rv==-1);
	rv = deleteNeighNode(Nei);
	assert(rv==0);
	Nei=newNeighNode(5);
	setNeighNode_hopstart(Nei,newHop());
	rv=deleteNeighNode(Nei);
	assert(rv==0);

	Nei = newNeighNode(5);
	printf("Testing:printNeighNode\n");
	rv = printNeighNode(NULL);
	assert(rv==-1);
	printf("Neigh ID should be 5,hoplength 0\n");
	rv = printNeighNode(Nei);
	assert(rv==0);

	printf("Testing:getNeighNode_p\n");
	rvd = getNeighNode_p(NULL,1);
	assert(rvd==-1.0);
	rvd = getNeighNode_p(Nei,3);
	assert(rvd==-1.0);
	rvd=getNeighNode_p(Nei,-1);
	assert(rvd==-1.0);
	rvd=getNeighNode_p(Nei,1);
  setNeighNode_hopstart(Nei,newHop());
	rv = setNeighNode_p(Nei,5.3,1);			
	assert(rv==0);
  rvd=getNeighNode_p(Nei,1);
	assert(rvd==5.3);

	Nei2=newNeighNode(12);
	printf("Testing:setNextNeighNode\n");
	rv=setNextNeighNode(&Nei,NULL);
	assert(rv==-1);
	rv=setNextNeighNode(NULL,&Nei2);
	assert(rv==-1);
	rv=setNextNeighNode(&Nei,&Nei2);
	assert(rv==0);

	printf("Testing:getNextNeigh\n");
	Nei2 = getNextNeigh(NULL);
	assert(Nei2==NULL);
	Nei2 = getNextNeigh(Nei);
	assert(Nei2!=NULL);
	assert(getNextNeigh(Nei)==Nei2);

	printf("Testing:setNeighNodeNew_p\n");
	rv = setNeighNodeNew_p(NULL,23.3);
	assert(rv==-1);
	rv = setNeighNodeNew_p(Nei,-12.3);
	assert(rv==-1);
	rv = setNeighNodeNew_p(Nei,32.9);
	assert(rv==0);
	rvd = getNeighNode_p(Nei,1);
	assert(rvd=32.9);

	printf("Testing:getNeighNode_id\n");
	rv=getNeighNode_id(NULL);
	assert(rv==-1);
	rv=getNeighNode_id(Nei2);
	assert(rv==12);

	printf("Testing:setNeighNode_id\n");
	rv = setNeighNode_id(NULL,3);
	assert(rv==-1);
	rv = setNeighNode_id(Nei,-32);
	assert(rv==-1);
	rv = setNeighNode_id(Nei,99);
	assert(rv==0);
	rv = getNeighNode_id(Nei);
	assert(rv==99);
	//deleteNeighNode(Nei);

	printf("Testing:setNeighNode_p\n");
	rv=setNeighNode_p(NULL,13.2,2);
	assert(rv==-1);
	rv=setNeighNode_p(Nei,-1.5,2);
	assert(rv==-1);
	rv=setNeighNode_p(Nei,13.2,0);
	assert(rv==-1);
	rv=setNeighNode_p(Nei,13.2,15);
	assert(rv==-1);
	rv=setNeighNode_p(Nei2,13.2,1);
	assert(rv==-1);
	rv=setNeighNode_p(Nei,13.2,1);
	assert(rv==0);

	printf("Testing:getNeighNode_t\n");
	rv=getNeighNode_t(NULL,2);
	assert(rv==-1);
	rv=getNeighNode_t(Nei,0);
	assert(rv==-1);
	rv=getNeighNode_t(Nei,15);
	assert(rv==-1);

	printf("Testing:setNeighNode_t\n");
	rv=setNeighNode_t(NULL,11.4,1);
	assert(rv==-1);
	rv=setNeighNode_t(Nei,-1.5,1);
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,0);
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,15);
	assert(rv==-1);
	rv=setNeighNode_t(Nei2,11.4,1); //for ->start==NULL?
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,1);
	assert(rv==0);
	rvd=getNeighNode_t(Nei,1);
	assert(rvd==11.4);

	printf("Testing: getNeighNode_hoplength\n");
	rv=getNeighNode_hoplength(NULL);
	assert(rv==-1);
	rv=getNeighNode_hoplength(Nei);
	assert(rv=2);

	printf("Testing:setNeighNode_hopstart\n"); //no getter
	h = newHop();
  rv=setNeighNode_hopstart(NULL,h);
	assert(rv==-1);
	rv=setNeighNode_hopstart(Nei,NULL);
	assert(rv==-1);
  Nei3 = newNeighNode(5);	
  rv=setNeighNode_hopstart(Nei3,h);
	assert(rv==0);

//////////////////////////////////////////////////////////////////////
	printf("Testing:newNeighLL\n");
	neighLL = newNeighLL();
	assert(neighLL!=NULL);

	printf("Testing:deleteNeighLL\n");
	rv = deleteNeighLL(NULL);
	assert(rv==-1);
	rv = deleteNeighLL(neighLL);
	assert(rv==0);
	
	printf("Testing:printNeighLL\n");
	rv = printNeighLL(NULL);
	assert(rv==-1);
	neighLL = newNeighLL();
	printf("Should print size of 0\n");
	rv = printNeighLL(neighLL);
	assert(rv==0);

	printf("Testing:deleteNeighLLAll");
	rv=deleteNeighLLAll(NULL);
	assert(rv==-1);
	rv=deleteNeighLLAll(neighLL);
	assert(rv==0);
	neighLL=newNeighLL();
	setNeighLL_start(neighLL,newNeighNode(15));
	assert(deleteNeighLLAll(neighLL)==0);
  neighLL=newNeighLL();
  setNeighLL_addNeighNode(neighLL, newNeighNode(20));
  setNeighLL_addNeighNode(neighLL, newNeighNode(21));
	assert(getNeighLL_numNeigh(neighLL)==2);
	rv=deleteNeighLLAll(neighLL);
	assert(rv==0);

	printf("Testing:getNeighLL_numNeigh");
	rv=getNeighLL_numNeigh(NULL);
	assert(rv==-1);
	neighLL=newNeighLL();
	rv=getNeighLL_numNeigh(neighLL);
	assert(rv==0);

	printf("Testing:setNeighLL_start"); //no getter
	NeighNode NeiNod =newNeighNode(25);
	rv=setNeighLL_start(NULL,NeiNod);
	assert(rv==-1);
	rv=setNeighLL_start(neighLL,NULL);
	assert(rv==-1);
	rv=setNeighLL_start(neighLL,NeiNod);
	assert(rv==0);
  deleteNeighLLAll(neighLL);
  deleteNeighNode(Nei);
  deleteNeighNode(Nei2);
  deleteNeighNode(Nei3);

  return 0;
}
