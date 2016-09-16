#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "node.h"

int main() {
	
	int rv;
	double rvd;
	Node nd;
	Node nd2;

	printf("Testing:newNode\n");
	nd = newNode(-1);
	assert(nd==NULL);
	nd = newNode(0);
	assert(nd!=NULL);

	printf("Testing:deleteNode\n");
	rv = deleteNode(NULL);
	assert(rv==-1);
	rv = deleteNode(nd);
	assert(rv==0);

	printf("Testing:printNode\n");
	nd = newNode(2);
	assert(nd!=NULL);
	rv = printNode(NULL);
	assert(rv==-1);
	printf("Node ID should be 2 followed by 0's\n");
	rv = printNode(nd);
	assert(rv==0);

	printf("Testing:getNode_id\n");
	rv = getNode_id(NULL);
	assert(rv==-1);
	rv = getNode_id(nd);
	assert(rv==2);

	printf("Testing:getNextNode\n");
	nd2 = getNextNode(NULL);
	assert(nd2==NULL);
	nd2 = getNextNode(nd);
	assert(nd2==NULL);
	Node nd3=getNextNode(nd);
	assert(nd3==nd2);
	
	printf("Testing:setNode_id\n");
	rv = setNode_id(NULL,1);
	assert(rv==-1);
	rv = setNode_id(nd,-2);
	assert(rv==-1);
	rv = setNode_id(nd,88);
	assert(rv==0);
	rv = getNode_id(nd);
	assert(rv==88);

	printf("Testing:getNode_p\n");
	rvd = getNode_p(NULL);
	assert(rvd==-1.0);
	rvd = getNode_p(nd);
	assert(rvd==0.0);

	printf("Testing:setNode_p\n");
	rv = setNode_p(NULL, 1.1);
	assert(rv==-1);
	rv = setNode_p(nd,-2.3);
	assert(rv==-1);
	rv = setNode_p(nd,5.5); //p-value can't be greater than 1? check clusterfunctions
	assert(rv==0);
	rvd = getNode_p(nd);
	assert(rvd=5.5);

	printf("Testing:getFlagFro\n");
	rv = getFlagFro(NULL);
	assert(rv==-1);
	rv = getFlagFro(nd);
	assert(rv==0);
	printf("Testing:getFlagBeh\n");
	rv = getFlagBeh(NULL);
	assert(rv==-1);
	rv = getFlagBeh(nd);
	assert(rv==0);
	printf("Testing:getFlagLef\n");
	rv = getFlagLef(NULL);
	assert(rv==-1);
	rv = getFlagLef(nd);
	assert(rv==0);
	printf("Testing:getFlagRig\n");
	rv = getFlagRig(NULL);
	assert(rv==-1);
	rv = getFlagRig(nd);
	assert(rv==0);
	printf("Testing:getFlagAbo\n");
	rv = getFlagAbo(NULL);
	assert(rv==-1);
	rv = getFlagAbo(nd);
	assert(rv==0);
	printf("Testing:getFlagBel\n");
	rv = getFlagBel(NULL);
	assert(rv==-1);
	rv = getFlagBel(nd);
	assert(rv==0);
	printf("Testing:getFlag\n");
	rv = getFlag(NULL,3);
	assert(rv==-1);
	rv = getFlag(nd,-1);
	assert(rv==-1);
	rv = getFlag(nd,6);
	assert(rv==-1);
	rv = getFlag(nd,0);
	assert(rv==0);

	printf("Testing:setFlagFro\n");
	rv = setFlagFro(NULL);
	assert(rv==-1);
	rv = setFlagFro(nd);
	assert(rv==0);
	rv = getFlagFro(nd);
	assert(rv==1);
	printf("Testing:setFlagBeh\n");
	rv = setFlagBeh(NULL);
	assert(rv==-1);
	rv = setFlagBeh(nd);
	assert(rv==0);
	rv = getFlagBeh(nd);
	assert(rv==1);
	printf("Testing:setFlagLef\n");
	rv = setFlagLef(NULL);
	assert(rv==-1);
	rv = setFlagLef(nd);
	assert(rv==0);
	rv = getFlagLef(nd);
	assert(rv==1);
	printf("Testing:setFlagRig\n");
	rv = setFlagRig(NULL);
	assert(rv==-1);
	rv = setFlagRig(nd);
	assert(rv==0);
	rv = getFlagRig(nd);
	assert(rv==1);
	printf("Testing:setFlagAbo\n");
	rv = setFlagAbo(NULL);
	assert(rv==-1);
	rv = setFlagAbo(nd);
	assert(rv==0);
	rv = getFlagAbo(nd);
	assert(rv==1);
	printf("Testing:setFlagBel\n");
	rv = setFlagBel(NULL);
	assert(rv==-1);
	rv = setFlagBel(nd);
	assert(rv==0);
	rv = getFlagBel(nd);
	assert(rv==1);
	rv = deleteNode(nd);

	//assert(rv==0); //why needed?
	nd = newNode(5);
	assert(rv==0);
	rv = getNode_id(nd);
	assert(rv==5);
	printf("Testing:setFlag\n");
	rv = getFlagRig(nd);
	assert(rv==0);
	rv = setFlag(nd, 3);
	assert(rv==0);
	rv = getFlagRig(nd);
	assert(rv==1);
	rv=setFlag(NULL,3);
	assert(rv==-1);
	rv=setFlag(nd,-1);
	assert(rv==-1);
	rv=setFlag(nd,6);
	assert(rv==-1);
	deleteNode(nd);
	deleteNode(nd2);
	deleteNode(nd3);

  return 0;
}
