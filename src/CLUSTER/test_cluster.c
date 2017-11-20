#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "cluster.h"
#include "../ELECTRODE/electrode.h"
#include "../NEIGHLL/neighll.h"
#include "../NEIGHNODE/neighnode.h"
#include "../NODE/node.h"
#include "../MIDPOINT/midpoint.h"

int main() {
	int rv;
	double rvd;
	NeighLL neighLL2;
  ClusterLL ClLL;
	ClusterLL ClLL2;
	ClusterLL ClLL3;
  ClusterLL ClLL6;
  ClusterLL ClLL5;
  MidPoint mp2;

  

	//////////////////////////////////////////////////////////////////////
	printf("Testing:newClusterLL\n");
	ClLL = newClusterLL(-1);
	assert(ClLL==NULL);
	ClLL = newClusterLL(1);
	assert(ClLL!=NULL);
	
	printf("Testing:deleteClusterLL\n");
	rv = deleteClusterLL(NULL);
	assert(rv==-1);
	rv = deleteClusterLL(ClLL);
	assert(rv==0);

	printf("Testing:deleteClusterLLNodes\n");
	ClLL = newClusterLL(3);
	assert(ClLL!=NULL);
	rv = deleteClusterLLNodes(ClLL);
	assert(rv==0);

	printf("Testing:deleteAllClusterLL\n");
	ClLL = newClusterLL(4);
	rv = deleteAllClusterLL(NULL);
	assert(rv==-1);
	rv = deleteAllClusterLL(&ClLL);
	assert(rv==0);	

	printf("Testing:printNodesClusterLL\n");
	rv = printNodesClusterLL(NULL);
	assert(rv==-1);
	ClLL = newClusterLL(5);
	printf("Should print no Nodes but the cluster id should be 5\n");
	rv = printNodesClusterLL(ClLL);
	assert(rv==0);
	
	printf("Testing:getNextClusterLL\n");
	ClLL2 = getNextClusterLL(NULL);
	assert(ClLL2==NULL);
	ClLL2 = newClusterLL(84);
	deleteClusterLL(ClLL2);
	ClLL2 = getNextClusterLL(ClLL);
	assert(ClLL2==NULL);
	deleteClusterLL(ClLL2);

	printf("Testing:printClusterLL\n");
	rv = printClusterLL(NULL);
	assert(rv==-1);
	rv = printClusterLL(ClLL);
	assert(rv==0); //start here

	
	printf("Testing:getCluster_NeighLL");
	neighLL2=getCluster_NeiLL(NULL);
	assert(neighLL2==NULL);
	neighLL2=getCluster_NeiLL(ClLL);
	assert(neighLL2==0); //have to create and set neighbor of ClLL

	printf("Testing:setCluster_NeiLL");

	printf("Testing:getCluster_id\n");
	rv = getCluster_id(NULL);
	assert(rv==-1);
	rv = getCluster_id(ClLL);
	assert(rv==5);

	printf("Testing:getCluster_numNodes\n");
	rv = getCluster_numNodes(NULL);
	assert(rv==-1);
	rv = getCluster_numNodes(ClLL);
	assert(rv==0);

	printf("Testing:setCluster_elecXFid");

	printf("Testing:setCluster_elecXBid");

	printf("Testing:setCluster_elecYLid");

	printf("Testing:setCluster_elecYRid");

	printf("Testing:setCluster_elecZBid");

	printf("Testing:setCluster_elecZAid");

	printf("Testing:getCluster_NeighLL");

	printf("Testing:setCluster_NeiLL");

	printf("Testing:setCluster_elecXid\n");
	rv = setCluster_elecXid(&ClLL, -1);
	assert(rv==-1);
	rv = setCluster_elecXid(&ClLL, 2);
	assert(rv==-1);
	rv = setCluster_elecXid(NULL, 1);
	assert(rv==-1);
	rv = setCluster_elecXid(&ClLL, 1);
	assert(rv==0);
	rv = setCluster_elecXid(&ClLL, 0);
	assert(rv==0);
	rv = setCluster_elecXid(&ClLL, 1);
	assert(rv==-1);

	printf("Testing:setCluster_elecYid\n");
	rv = setCluster_elecYid(&ClLL, -1);
	assert(rv==-1);
	rv = setCluster_elecYid(&ClLL, 2);
	assert(rv==-1);
	rv = setCluster_elecYid(NULL, 1);
	assert(rv==-1);
	rv = setCluster_elecYid(&ClLL, 1);
	assert(rv==0);
	rv = setCluster_elecYid(&ClLL, 0);
	assert(rv==0);
	rv = setCluster_elecYid(&ClLL, 1);
	assert(rv==-1);
	
	printf("Testing:setCluster_elecZid\n");
	rv = setCluster_elecZid(&ClLL, -1);
	assert(rv==-1);
	rv = setCluster_elecZid(&ClLL, 2);
	assert(rv==-1);
	rv = setCluster_elecZid(NULL, 1);
	assert(rv==-1);
	rv = setCluster_elecZid(&ClLL, 1);
	assert(rv==0);
	rv = setCluster_elecZid(&ClLL, 0);
	assert(rv==0);
	rv = setCluster_elecZid(&ClLL, 1);
	assert(rv==-1);

	ClLL2 = newClusterLL(66);
	int ElecXF = 1;
	printf("Testing:setCluster_elecXsum\n");
	rv = setCluster_elecXsum(ClLL, 1.1, ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(ClLL, -1.1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(NULL, 1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(ClLL2,1,ElecXF);
	assert(rv==-1);
  Electrode elxf = newElectrode();
  rv = setCluster_elecXF(&ClLL2,elxf);
	assert(rv==0);
  rv = setCluster_elecXsum(ClLL2, 1.1, ElecXF);  
  assert(rv==0); 
  
	int ElecYR = 1;
	printf("Testing:setCluster_elecYsum\n");
	rv = setCluster_elecYsum(ClLL, 2.1, ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(ClLL, -2, ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(NULL, 1, ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(ClLL2,3, ElecYR);
	assert(rv==-1);
	ClLL6 = newClusterLL(63);
  Electrode elyr = newElectrode();
  rv = setCluster_elecYR(&ClLL6,elyr);
  assert(rv==0);
  rv = setCluster_elecYsum(ClLL6, 2.1, ElecYR);  
  assert(rv==0); 

	int ElecZA = 1;
	printf("Testing:setCluster_elecZsum\n");
	rv = setCluster_elecZsum(ClLL, 2.9, ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(ClLL, -1.1, ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(NULL, 1, ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(ClLL2,9.1, ElecZA);
	assert(rv==-1);
	ClLL5 = newClusterLL(69);
  Electrode elza = newElectrode();
  rv = setCluster_elecZA(&ClLL5,elza);
  assert(rv==0);
  rv = setCluster_elecZsum(ClLL5, 1.1, ElecZA);  
  assert(rv==0); 

	printf("Testing:addToCluster_elecXsum\n");
	rv = addToCluster_elecXsum(ClLL,-2,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(NULL,1,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(ClLL,1,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(ClLL2,1,ElecXF);
	assert(rv==0);

	printf("Testing:addToCluster_elecYsum\n");
	rv = addToCluster_elecYsum(ClLL,-2, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(NULL,1, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(ClLL2,1, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(ClLL,1, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(ClLL6,1, ElecYR);
	assert(rv==0);

	printf("Testing:addToCluster_elecZsum\n");
	rv = addToCluster_elecZsum(ClLL,-2, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(NULL,1, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(ClLL2,1, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(ClLL,1, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(ClLL5,1, ElecZA);
	assert(rv==0);

	printf("Testing:getCluster_elecXsum\n");
	rvd = getCluster_elecXsum(NULL,ElecXF);
	assert(rvd==-1);
	rvd = getCluster_elecXsum(ClLL2,ElecXF);
	assert(rvd==2.1);

	printf("Testing:getCluster_elecYsum\n");
	rvd = getCluster_elecYsum(NULL,ElecYR);
	assert(rvd==-1);
	rvd = getCluster_elecYsum(ClLL2,ElecYR);
	assert(rvd==-1);
	rvd = getCluster_elecYsum(ClLL6,ElecYR);
	assert(rvd==3.1);

	printf("Testing:getCluster_elecZsum\n");
	rvd = getCluster_elecZsum(NULL,ElecZA);
	assert(rvd==-1);
	rvd = getCluster_elecZsum(ClLL2,ElecZA);
	assert(rvd==-1);
	rvd = getCluster_elecZsum(ClLL5,ElecZA);
	assert(rvd==2.1);

	printf("Testing:setCluster_id\n");
	rv = setCluster_id(ClLL,-1);
	assert(rv==-1);
	rv = setCluster_id(NULL,2);
	assert(rv==-1);
	rv = setCluster_id(ClLL,80);
	assert(rv==0);
	rv = getCluster_id(ClLL);
	assert(rv==80);

	printf("Testing:getCluster_Sum\n");
	rvd = getCluster_Sum(NULL);
	assert(rvd==-1.0);
	rvd = getCluster_Sum(ClLL);
	assert(rvd==0);

	printf("Testing:addToClusterSum\n");
	rv = addToClusterSum(NULL,3.3);
	assert(rv==-1);
	rv = addToClusterSum(ClLL,-13);
	assert(rv==-1);
	rv = addToClusterSum(ClLL,32);
	assert(rv==0);
	rvd = getCluster_Sum(ClLL);
	assert(rvd==32);
	
	printf("Testing:setNextClusterLL\n");
	rv = setNextClusterLL(NULL, ClLL2);
	assert(rv==-1);
	assert(ClLL2!=NULL);
	rv = setNextClusterLL(ClLL, ClLL2);
	assert(rv==0);
	ClLL3 = getNextClusterLL(ClLL);
	assert(ClLL3!=NULL);
	assert(getCluster_id(ClLL2)==getCluster_id(ClLL3));
	rv = setNextClusterLL(ClLL, NULL);
	assert(rv==0);
	ClLL3 = getNextClusterLL(ClLL);
	assert(ClLL3==NULL);
	deleteAllClusterLL(&ClLL);
	deleteAllClusterLL(&ClLL2);
	deleteAllClusterLL(&ClLL6);
	deleteAllClusterLL(&ClLL5);
  deleteElectrode(&elxf);	
  deleteElectrode(&elyr);	
  deleteElectrode(&elza);	
	//////////////////////////////////////////////////////////////////////
	ClLL = newClusterLL(5);
	rv = setCluster_elecXid(&ClLL, 0);
	rv = setCluster_elecYid(&ClLL, 0);
	rv = setCluster_elecZid(&ClLL, 0);
	rv = setCluster_elecXid(&ClLL, 1);
	rv = setCluster_elecYid(&ClLL, 1);
	rv = setCluster_elecZid(&ClLL, 1);
	rv = setCluster_elecXsum(ClLL, 1.1, ElecXF);
	rv = setCluster_elecYsum(ClLL, 2.1, ElecYR);
	rv = setCluster_elecZsum(ClLL, 2.9, ElecZA);
	rv = addToCluster_elecXsum(ClLL,1,ElecXF);
	rv = addToCluster_elecYsum(ClLL,1,ElecYR);
	rv = addToCluster_elecZsum(ClLL,1, ElecZA);
	rv = setCluster_id(ClLL,80);
	rv = addToClusterSum(ClLL,32);
	ClLL2 = newClusterLL(84);
	
	Node nd4;
	printf("Testing:getStartNode\n");
	nd4 = getStartNode(ClLL);
	assert(nd4==NULL);
	nd4 = getStartNode(NULL);
	assert(nd4==NULL);

	NeighNode neighN;
	printf("Testing:getStartNeigh\n");
	neighN = getStartNeigh(ClLL);
	assert(neighN==NULL);
	neighN = getStartNeigh(NULL);
	assert(neighN==NULL);

	printf("Testing:getCluster_numNodes\n");
	rv = getCluster_numNodes(NULL);
	assert(rv==-1);
	rv = getCluster_numNodes(ClLL);
	assert(rv==0);

	printf("Testing:setCluster_elecXsum\n");
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(32);
  Electrode elxf2 = newElectrode();
  setCluster_elecXF(&ClLL,elxf2);
  rv = setCluster_elecXsum(ClLL, 3.3,ElecXF);
	assert(rv==0);
	rv = setCluster_elecXsum(ClLL, -1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(NULL,2.5,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(ClLL2,2.3,ElecXF);
	assert(rv==-1);

	printf("Testing:setCluster_elecYsum\n");
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(82);
  Electrode elyr2 = newElectrode();
  setCluster_elecYR(&ClLL,elyr2);
  rv = setCluster_elecYsum(ClLL, 3.3,ElecYR);
	assert(rv==0);
	rv = setCluster_elecYsum(ClLL, -1,ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(NULL,2.5,ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(ClLL2,2.3,ElecYR);
	assert(rv==-1);
	
	printf("Testing:setCluster_elecZsum\n");
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(102);
  Electrode elza2 = newElectrode();
  setCluster_elecZA(&ClLL,elza2);
	rv = setCluster_elecZsum(ClLL, 3.3,ElecZA);
	assert(rv==0);
	rv = setCluster_elecZsum(ClLL, -1,ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(NULL,2.5,ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(ClLL2,2.3,ElecZA);
	assert(rv==-1);

	printf("Testing:getCluster_p\n");
	rvd = getCluster_p(ClLL);
	assert(rvd==0.0);
	rvd = getCluster_p(NULL);
	assert(rvd==-1.0);

	printf("Testing:setCluster_p\n");
	rv = setCluster_p(ClLL,1.1);
	assert(rv==0);
	rvd = getCluster_p(ClLL);
	assert(rvd==1.1);
	rv = setCluster_p(NULL,2.2);
	assert(rv==-1);
	rv = setCluster_p(ClLL,-1.3);
	assert(rv==-1);
	rvd = getCluster_p(ClLL);
	assert(rvd==1.1);

	printf("Testing:getCluster_elecXid\n");
  deleteElectrode(&elxf2);
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(32);
  elxf2 = newElectrode();
  Electrode elxf3 = newElectrode();
  setCluster_elecXB(&ClLL,elxf3);
  setCluster_elecXF(&ClLL,elxf2);
	rv = getCluster_elecXid(ClLL);
	assert(rv==2);
	rv = getCluster_elecXid(ClLL2);
	assert(rv==-1);
	rv = getCluster_elecXid(NULL);
	assert(rv==-1);

	printf("Testing:getCluster_elecYid\n");
  deleteElectrode(&elyr2);
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(82);
  elyr2 = newElectrode();
  Electrode elyr3 = newElectrode();
  setCluster_elecYL(&ClLL,elyr3);
  setCluster_elecYR(&ClLL,elyr2);
	rv = getCluster_elecYid(ClLL);
	assert(rv==2);
	rv = getCluster_elecYid(ClLL2);
	assert(rv==-1);
	rv = getCluster_elecYid(NULL);
	assert(rv==-1);

	printf("Testing:getCluster_elecZid\n");
  deleteElectrode(&elza2);
  deleteAllClusterLL(&ClLL);
  ClLL = newClusterLL(102);
  elza2 = newElectrode();
  Electrode elza3 = newElectrode();
  setCluster_elecZB(&ClLL,elza3);
  setCluster_elecZA(&ClLL,elza2);
	rv = getCluster_elecZid(ClLL);
	assert(rv==2);
	rv = getCluster_elecZid(ClLL2);
	assert(rv==-1);
	rv = getCluster_elecZid(NULL);
	assert(rv==-1);

	printf("Testing:getCluster_time\n");
	rvd = getCluster_time(NULL);
	assert(rvd==-1.0);
	rvd = getCluster_time(ClLL);
	assert(rvd==0.0);

	deleteNode(nd4);	
	deleteAllClusterLL(&ClLL);
	deleteAllClusterLL(&ClLL2);
  deleteElectrode(&elxf2);	
  deleteElectrode(&elxf3);	
  deleteElectrode(&elyr2);	
  deleteElectrode(&elyr3);	
  deleteElectrode(&elza2);	
  deleteElectrode(&elza3);	
	//////////////////////////////////////////////////////////////////////
	printf("Testing:addNodeEndClusterLL\n");
	ClusterLL ClLL9 = newClusterLL(77);
	Node nd11;
	rv = addNodeEndClusterLL(NULL,2);
	assert(rv==-1);
	rv = addNodeEndClusterLL(ClLL9,-1);
	assert(rv==-1);
	rv = addNodeEndClusterLL(ClLL9,73);
	assert(rv==0);
	rv = getCluster_numNodes(ClLL9);
	assert(rv==1);
	nd11 = getStartNode(ClLL9);
	assert(getNode_id(nd11)==73);
	assert(getNextNode(nd11)==NULL);
	rv = addNodeEndClusterLL(ClLL9,82);
	assert(rv==0);
	assert(getCluster_numNodes(ClLL9)==2);
	assert(getNextNode(nd11)!=NULL);
	nd11 = getNextNode(nd11);
	assert(getNode_id(nd11)==82);
	deleteAllClusterLL(&ClLL9);

  deleteElectrode(&elxf2);
	//////////////////////////////////////////////////////////////////////
	ClLL = newClusterLL(5);
	rv = setCluster_elecXid(&ClLL, 0);
	rv = setCluster_elecYid(&ClLL, 0);
	rv = setCluster_elecZid(&ClLL, 0);
	rv = setCluster_elecXid(&ClLL, 1);
	rv = setCluster_elecYid(&ClLL, 1);
	rv = setCluster_elecZid(&ClLL, 1);
	rv = setCluster_elecXsum(ClLL, 1.1, ElecXF);
	rv = setCluster_elecYsum(ClLL, 2.1, ElecYR);
	rv = setCluster_elecZsum(ClLL, 2.9, ElecZA);
	rv = addToCluster_elecXsum(ClLL,1,ElecXF);
	rv = addToCluster_elecYsum(ClLL,1,ElecYR);
	rv = addToCluster_elecZsum(ClLL,1, ElecZA);
	rv = setCluster_id(ClLL,80);
	rv = addToClusterSum(ClLL,32);
	
	printf("Testing:addNodeNewCluster\n");
	mp2 = newMidPoint(-2, 1, 3, 298);
	Node nd6;
	ClusterLL ClLL22;
	ClLL22=NULL;
	ClusterLL ClLL4;
	rv = addNodeNewCluster(&ClLL22,mp2);
	assert(rv==-1);
	rv = addNodeNewCluster(&ClLL,NULL);
	assert(rv==-1);
	rv = addNodeNewCluster(&ClLL,mp2);
	assert(rv==0);
	rv = getCluster_numNodes(ClLL);
	assert(rv==0);
	nd6 = getStartNode(ClLL);
	assert(nd6==NULL);
	ClLL4 = getNextClusterLL(ClLL);
	rv = getCluster_numNodes(ClLL4);
	assert(rv==2);
	nd6 = getStartNode(ClLL4);
	assert(getNode_id(nd6)==3);
	nd6 = getNextNode(nd6);
	assert(getNode_id(nd6)==298);
	assert(getNextClusterLL(ClLL4)==NULL);
	deleteMidPoint(mp2);
	deleteAllClusterLL(&ClLL);

	//////////////////////////////////////////////////////////////////////
	printf("Testing:addNodeToCluster\n");
  	ClusterLL ClLL12 = newClusterLL(200);
	assert(getCluster_id(ClLL12)==200);
	ClusterLL ClLL13 ;
	MidPoint mp22 = newMidPoint(-3,33,5,4);
	Node nd13;
	rv = addNodeToCluster( NULL, mp22);
	assert(rv==-1);
	rv = addNodeToCluster( ClLL12, NULL);
	assert(rv==-1);
	rv = addNodeToCluster( ClLL12, mp22);
	assert(rv==0);
	rv = getCluster_numNodes(ClLL12);
	assert(rv==2);
	nd13 = getStartNode(ClLL12);
	assert(getNode_id(nd13)==5);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==4);
	deleteMidPoint(mp22);
	
	MidPoint mp23 = newMidPoint(-3, 22, 8, 11);
	rv = addNodeToCluster( ClLL12, mp23);
	assert(rv==0);
	rv = getCluster_numNodes(ClLL12);
	assert(rv==2);
	ClLL13 = getNextClusterLL(ClLL12);
	assert(ClLL13!=NULL);
	rv = getCluster_id(ClLL13);
	assert(rv==201);
	rv = getCluster_numNodes(ClLL13);
	assert(rv==2);
	deleteMidPoint(mp23);

	MidPoint mp24 = newMidPoint(-3,33,5,4);
	rv = addNodeToCluster(ClLL12, mp24);
	assert(rv==-1);
	assert(getCluster_numNodes(ClLL12)==2);
	deleteMidPoint(mp24);

	MidPoint mp25 = newMidPoint(-3,33,5,22);
	rv = addNodeToCluster(ClLL12, mp25);
	assert(rv==0);
	assert(getCluster_numNodes(ClLL12)==3);
	nd13 = getStartNode(ClLL12);
	nd13 = getNextNode(nd13);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==22);
	deleteMidPoint(mp25);

	MidPoint mp26 = newMidPoint(-3,33,8,28);
	rv = addNodeToCluster(ClLL12, mp26);
	assert(rv==0);
	ClLL13 = getNextClusterLL(ClLL12);
	assert(getCluster_numNodes(ClLL13)==3);
	nd13 = getStartNode(ClLL13);
	nd13 = getNextNode(nd13);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==28);
	deleteMidPoint(mp26);	

	MidPoint mp27 = newMidPoint(-3,33,22,28);
	rv = addNodeToCluster(ClLL12,mp27);
	assert(rv==0);
	assert(getCluster_numNodes(ClLL12)==6);
	assert(getNextClusterLL(ClLL12)==NULL);
	deleteMidPoint(mp27);
	nd13 = getStartNode(ClLL12);
	assert(getNode_id(nd13)==5);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==4);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==22);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==8);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==11);
	nd13 = getNextNode(nd13);
	assert(getNode_id(nd13)==28);
	nd13 = getNextNode(nd13);
	assert(nd13==NULL);
	deleteAllClusterLL(&ClLL12);
	
	//////////////////////////////////////////////////////////////////////
	ClLL = newClusterLL(5);
	rv = setCluster_elecXid(&ClLL, 0);
	rv = setCluster_elecYid(&ClLL, 0);
	rv = setCluster_elecZid(&ClLL, 0);
	rv = setCluster_elecXid(&ClLL, 1);
	rv = setCluster_elecYid(&ClLL, 1);
	rv = setCluster_elecZid(&ClLL, 1);
	rv = setCluster_elecXsum(ClLL, 1.1, ElecXF);
	rv = setCluster_elecYsum(ClLL, 2.1, ElecYR);
	rv = setCluster_elecZsum(ClLL, 2.9, ElecZA);
	rv = addToCluster_elecXsum(ClLL,1,ElecXF);
	rv = addToCluster_elecYsum(ClLL,1,ElecYR);
	rv = addToCluster_elecZsum(ClLL,1, ElecZA);
	rv = setCluster_id(ClLL,80);
	rv = addToClusterSum(ClLL,32);
	mp2 = newMidPoint(-2, 1, 3, 298);
	rv = addNodeNewCluster(&ClLL,mp2);
	ClLL4 = getNextClusterLL(ClLL);
	
	printf("Testing:addNeighNodeToCluster\n");
	NeighNode NeighN7;
	ClusterLL ClLL24 = NULL;
	rv = addNeighNodeToCluster(&ClLL24,55);
	assert(rv==-1);
	ClusterLL ClLL8 = newClusterLL(99);
	rv = addNeighNodeToCluster(&ClLL8,-1);
	assert(rv==-1);
	rv = addNeighNodeToCluster(&ClLL8,34);
	assert(rv==0);
	NeighN7 = getStartNeigh(ClLL8);
	assert(getNeighNode_id(NeighN7)==34);
	rv = getCluster_numNeigh(ClLL8);
	assert(rv==1);
	rv = addNeighNodeToCluster(&ClLL8,34);
	printClusterLL(ClLL8);
  assert(rv==-1);
	rv = getCluster_numNeigh(ClLL8);
	assert(rv==1);
	NeighN7 = getStartNeigh(ClLL8);
	assert(NeighN7!=NULL);
	assert(getNextNeigh(NeighN7)==NULL);
	rv = addNeighNodeToCluster(&ClLL8,82);
	assert(rv==0);
	rv = getCluster_numNeigh(ClLL8);
	assert(rv==2);
	NeighN7 = getStartNeigh(ClLL8);
	NeighN7 = getNextNeigh(NeighN7);
	assert(getNeighNode_id(NeighN7)==82);
	assert(getNextNeigh(NeighN7)==NULL);
	
	printf("Testing:addPvalClusterNode\n");
	rv = addPvalClusterNode(NULL,1,3.4);
	assert(rv==-1);
	rv = addPvalClusterNode(ClLL,-1,34);
	assert(rv==-1);
	rv = addPvalClusterNode(ClLL,1,-89);
	assert(rv==-1);
	rv = addPvalClusterNode(ClLL,1,32);
	assert(rv==-1);
	rv = addPvalClusterNode(ClLL4,1,2.2);
	assert(rv==-1);
	rv = addPvalClusterNode(ClLL4,3, 5.5);
	assert(rv==0);
	rvd = getCluster_time(ClLL4);
	assert(rvd==1/5.5);
	nd6 = getStartNode(ClLL4);
	assert(getNode_p(nd6)==5.5);

	printf("Testing:addPvalClusterNeighNode\n");
	ClusterLL ClLL54 = newClusterLL(92);
	rv = addPvalClusterNeighNode(NULL, 34, 3.3);
	assert(rv==-1);
	rv = addPvalClusterNeighNode(ClLL8,-2, 3.3);
	assert(rv==-1);
	rv = addPvalClusterNeighNode(ClLL8,3, 3.3);
	assert(rv==-1);
	rv = addPvalClusterNeighNode(ClLL8, 34, -1);
	assert(rv==-1);
	rv = addPvalClusterNeighNode(ClLL54, 34, 3.3);
	assert(rv==-1);
	rv = addPvalClusterNeighNode(ClLL8, 34, 2.2);
	assert(rv==0);
	NeighN7 = getStartNeigh(ClLL8);
	assert(NeighN7!=NULL);
	assert(getNeighNode_p(NeighN7,1)==2.2);
	
	deleteMidPoint(mp2);
	deleteAllClusterLL(&ClLL);
	deleteAllClusterLL(&ClLL8);
	deleteAllClusterLL(&ClLL54);
  return 0;
}
