#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "cluster.h"

int main() {
	int nei1, nei2, order;
	int rv;
	double rvd;
	void * Temp;
	void * Temp2;
	nei1=22;
	nei2=74;
	order=-2;
	MidPoint mp;
	OrderMagLL OMLL;
	Node nd;
	Node nd2;
	Hop h;
	NeighNode Nei;
	NeighNode Nei2;
	Electrode el;
	NeighLL neighLL;
	ClusterLL ClLL;
	ClusterLL ClLL2;
	ClusterLL ClLL3;

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
	assert(rv==0);
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
	mp6=newMidPOint(18,12,21,212);
	rv=CompareNeiMidPoint(NULL,mp);
	assert(rv==NULL);
	rv=CompareNeiMidPoint(mp,NULL);
	assert(rv==NULL);
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

///////////////////////////////////////////////////////////////////////////////
//Order of Magnitude Link List

	printf("Testing:newOrLL\n");
	OMLL = newOrLL(1);
	assert(OMLL!=NULL);

	printf("Testing:checkNewOrLL");
	rv=checkNewOrLL(NULL,mp);
	assert(rv==-1);
	rv=checkNewOrLL(OMLL,NULL);
	assert(rv==-1);
	OMLL->start=newMidPoint(3,17,nei1,nei2);
	deleteMidPoint(OMLL->start);
	rv=checkNewOrLL(OMLL,OMLL->start);
	assert(rv==-1);
	mp=newMidPoint(1,15,nei1,nei2);
	rv=checkNewOrLL(OMLL,mp);
	assert(rv==0);

	printf("Testing:deleteOrLL\n");
	rv = deleteOrLL(NULL);
	assert(rv==-1);
	rv = deleteOrLL(OMLL);
	assert(rv==-1);
	rv=deleteOrLL(OMLL);
	assert(rv==0);

	printf("Testing:deleteAllOrLL\n"); 
	OMLL = newOrLL(order);
	assert(OMLL!=NULL);
	rv = deleteAllOrLL(NULL);
	assert(rv==-1);
	rv = deleteAllOrLL(OMLL);
	assert(rv==0);
	OMLL=newOrLL(-3);
	OMLL->start=newMidPoint(order,6,nei1,nei2);
	OMLL->start->next=newMidPoint(order,7,nei2,nei1);
	rv=deleteAllorLL(OMLL);
	assert(rv==0); 
	OMLL=newOrLL(-4);
	OMLL->start=newMidPoint(order,8,nei1,nei2);
	OMLL->start->next=newMidPoint(order,9,nei1,nei2);
	OMLL->start->next->next=newMidPoint(order,10,15,12);
	rv=deleteAllorLL(OMLL);
	assert(rv==0);

	OMLL = newOrLL(5);
	printf("Testing:printOrLL\n");
	rv = printOrLL(NULL);
	assert(rv==-1);
	printf("Should print OrderMag of 5 and Size of LL of 0\n");
	rv = printOrLL( OMLL );
	assert(rv==0);
	OMLL->start=newMidPoint(5,5,nei1,nei2);
	printf("Should print OrderMag of 5 and Size of 1\n");
	rv=printOrLL(OMLL);
	assert(rv==0);

	printf("Testing:getOMLL_size\n");
	rv=getOMLL_size(NULL);
	assert(rv==-1);
	rv=getOMLL_size(OMLL);
	assert(rv==1);

	printf("getOMLL_order\n");
	rv=getOMLL_order(NULL);
	assert(rv==-1);
	rv=getOMLL_order(OMLL);
	assert(rv==5);

	printf("Testing:getOMLLstartMP\n");
	mp = getOMLLstartMP(NULL);
	assert(mp==NULL);
	mp = getOMLLstartMP(OMLL);
	assert(mp!=NULL); //
	deleteAllOrLL(OMLL);

///////////////////////////////////////////////////////////////////////////////

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
	Node nd2 = getNextNode(NULL);
	assert(nd2==NULL);
	nd2 = getNextNode(nd);
	assert(nd2==NULL);
	nd2=nd->next;
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
	//////////////////////////////////////////////////////////////////////
	//Tools for accessing Hop
	h=newHop();
	printf("Testing:deleteAllHop\n");
	rv=deleteAllHop(NULL);
	assert(rv==-1);
	rv=deleteAllHop(h);
	assert(rv==0);
	h=newHop();
	Hop h2=newHop();
	h2 = h->next;
	rv=deleteAllHop(h);
	assert(rv==0);

	h=newHop();
	printf("Testing:printHop\n");
	rv=printHop(NULL);
	assert(rv==-1);
	printf("Expecting p and t = 0.0");
	rv=printHop(h);
	assert(rv==0);

	printf("Testing:setHop_t\n");
	rv==setHop_t(NULL,5);
	assert(rv==-1);
	rv==setHop_t(h,-5);
	assert(rv==-1);
	rv==setHop_t(h,5);
	assert(rv==0);
	
	printf("Testing:setHop_p\n");
	rv==setHop_p(NULL,51);
	assert(rv==-1);
	rv==setHop_p(h,51);
	assert(rv==-1);
	rv==setHop_p(h,51);
	assert(rv==0);
	printf("Expecting:t=5,p=51");
	printHop(h);

	printf("Testing:setHop_next\n");
	h2=newHop();
	rv=setHop_next(NULL,h2);
	assert(rv==-1);
	rv=setHop_next(h1,NULL);
	assert(rv==-1);
	rv=setHop_next(h1,h2);
	assert(rv==0);
	//////////////////////////////////////////////////////////////////////
	printf("Testing:newNeighNode\n");
	NeighNode Nei = newNeighNode(-1);
	assert(Nei==NULL);
	Nei = newNeighNode(3);
	assert(Nei!=NULL);

	printf("Testing:deleteNeighNode\n");
	rv = deleteNeighNode(NULL);
	assert(rv==-1);
	rv = deleteNeighNode(Nei);
	assert(rv==0);
	Nei=newNeighNode(5);
	Nei->start=newHop();
	Nei->hoplength=1;
	rv=deleteNeighNode(Nei);
	assert(rv==0);

	Nei = newNeighNode(5);
	Nei->hoplength=2;
	printf("Testing:printNeighNode\n");
	rv = printNeighNode(NULL);
	assert(rv==-1);
	printf("Neigh ID should be 5,hoplength 2\n");
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
	Nei->start->p=5.3;			
	rvd=getNeighNode_p(Nei,1);
	assert(rvd==5.3);

	NeighNode Nei2=newNeighNode(12);
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
	assert(Nei2==NULL);
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
	rv==setNeighNode_p(Nei,-1.5,2);
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
	rv==setNeighNode_t(Nei,-1.5,1);
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,0);
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,15);
	assert(rv==-1);
	rv=setNeighNode_t(Nei2,11.4,1); //for ->start==NULL?
	assert(rv==-1);
	rv=setNeighNode_t(Nei,11.4,1);
	assert(rv==0);
	rv=getNeighNode_t(Nei,2);
	assert(rv==11.4);

	printf("Testing: getNeighNode_hoplength\n");
	rv=getNeighNode_hoplength(NULL);
	assert(rv==-1);
	rv==getNeighNode_hoplength(Nei);
	assert(rv=2);

	printf("Testing:setNeighNode_hopstart\n"); //no getter
	rv=setNeighNode_hopstart(NULL,h);
	assert(rv==-1);
	rv=setNeighNode_hopstart(Nei,NULL);
	assert(rv==-1);
	rv=setNeighNode_hopstart(Nei,h);
	assert(rv==0);

	printf("Testing:setNeighNode_hoplength\n");
	rv=setNeighNode_hoplength(NULL,5);
	assert(rv==-1);
	rv=setNeighNode_hoplength(Nei,0);
	assert(rv==-1);
	rv=setNeighnode_hoplength(Nei,5);
	assert(rv==0);
	assert(getNeighNode_hoplength(Nei)==5);

	//////////////////////////////////////////////////////////////////////

	printf("Testing:newElectrode\n");
	el = newElectrode();
	assert(el!=NULL);

	printf("Testing:deleteElectrode\n");
	rv = deleteElectrode(NULL);
	assert(rv==-1);
	rv = deleteElectrode(&el);
	assert(rv==0);

	printf("Testing:printElectrode\n");
	rv = printElectrode(NULL);
	assert(rv==-1);
	el = newElectrode();
	printf("Should print p of 0.0\n");
	rv = printElectrode(el);
	assert(rv==0);

	printf("Testing:getElectrode_alpha\n");
	assert(getElectrode_alpha(NULL)==-1);
	assert(getElectrode_alpha(el)==0);

	printf("Testing:setElectrode_alpha\n");
	rv = setElectrode_alpha(NULL, 34);
	assert(rv==-1);
	rv = setElectrode_alpha(el, 12);
	assert(rv==0);
	assert(getElectrode_alpha(el)==12);

	printf("Testing:getElectrode_Charges\n");
	rv = getElectrode_Charges(NULL);
	assert(rv==-1);
	rv = getElectrode_Charges(el);
	assert(rv==0);

	printf("Testing:setElectrode_Charges\n");
	rv = setElectrode_Charges(NULL, 1);
	assert(rv==-1);
	rv = setElectrode_Charges(el, -1);
	assert(rv==-1);
	//assert(getElectrode_Charges(el)==0);
	rv = setElectrode_Charges(el, 32);
	assert(rv==0);
	assert(getElectrode_Charges(el)==32);

	printf("Testing:Electrode_addCharge\n");
	rv = Electrode_addCharge(NULL);
	assert(rv==-1);
	rv = Electrode_addCharge(el);
	assert(rv==0);
	assert(getElectrode_Charges(el)==33);

	printf("Testing:Electrode_minusCharge\n");
	rv = Electrode_minusCharge(NULL);
	assert(rv==-1);
	rv = Electrode_minusCharge(el);
	assert(rv==0);
	assert(getElectrode_Charges(el)==32);

	printf("Testing:getElectrode_FermiEnergy\n");
	rvd = getElectrode_FermiEnergy(el);
	assert(rvd==-1);
	rvd = getElectrode_FermiEnergy(el);
	assert(rvd==0.0);

	printf("Testing:setElectrode_FermiEnergy\n");
	rv = setElectrode_FermiEnergy(NULL, 2.3);
	assert(rv==-1);
	rv = setElectrode_FermiEnergy(el, -1.1);
	assert(rv==0);
	assert(getElectrode_FermiEnergy(el)==-1.1);

	printf("Testing:getElectrode_Sum\n");
	rv = getElectrode_Sum(NULL);
	assert(rv==-1);
	rv = getElectrode_Sum(el);
	assert(rv==-1);

	printf("Testing:setElectrode_Sum\n");
	rv=setElectrode_Sum(NULL,15.4);
	assert(rv==-1);
	rv=setElectrode_Sum(el,-1);
	assert(rv==-1);
	rv=setElectrode_Sum(el,15.4);
	assert(rv==0);
	assert(getElectrode_Sum(el)==15.4);

	printf("Testing:getElectrode_HopRates\n");
	assert(getElectrode_HopRates(NULL)==NULL);
	assert(getElectrode_HopRates(el)==0.0);

	printf("Testing:setElectrode_HopRates\n"); //can be negative?
	assert(setElectrode_HopRates(NULL, (int) -5.5)==NULL);
	assert(setElectrode_HopRates(el, NULL==-1);
	assert(setElectrode_HopRates(el, (int) -5.5 == 0); 
	assert(getElectrode_HopRates(el)==-5.5);

	printf("Testing:getElectrode_AdjacentSites\n");
	assert(getElectrode_AdjacentSites(NULL)==NULL);
	assert(getElectrode_AdjacentSites(el)==0.0);

	printf("Testing:setElectrode_AdjacentSites\n"); 
	assert(setElectrode_AdjacentSites(NULL, (int)5)==NULL);
	assert(setElectrode_AdjacentSites(el,NULL==NULL);
	assert(setElectrode_AdjacentSites(el, (int)5==0);
	assert(getElectrode_AdjacentSites(el)==5);
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
	neighLL->start=newNeighNode(15);
	assert(deleteNeighLL(neighLL)==-1);
	neighLL->start->next=newNeighNode(20);
	assert(getNeighLL_numNeigh==2);
	rv=deleteNeighLL(neighLL);
	assert(rv==0);

	printf("Testing:getNeighLL_numNeigh");
	rv=getNeighLL_numNeigh(NULL);
	assert(rv==-1);
	neighLL=newNeighLL();
	rv=getNeighLL_numNeigh(neighLL);
	assert(rv==0);

	printf("Testing:setNeighLL_start"); //no getter
	neighNode=newNeighNode(25);
	rv=setNeighLL_start(NULL,neighNode);
	assert(rv==-1);
	rv=setNeighLL_start(newNeighLL,NULL);
	assert(rv==-1);
	rv=setNeighLL_start(newNeighLL,neighNode);
	assert(rv==0);

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
	rv=getCluster_NeiLL(NULL);
	assert(rv=NULL);
	rv=getCluster_NeiLL(ClLL);
	assert(rv= //have to create and set neighbor of ClLL

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
	assert(rv==0);
	rv = setCluster_elecXsum(ClLL, -1.1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(NULL, 1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(ClLL2,1,ElecXF);
	assert(rv==-1);

	int ElecYR = 1;
	printf("Testing:setCluster_elecYsum\n");
	rv = setCluster_elecYsum(ClLL, 2.1, ElecYR);
	assert(rv==0);
	rv = setCluster_elecYsum(ClLL, -2, ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(NULL, 1, ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(ClLL2,3, ElecYR);
	assert(rv==-1);

	int ElecZA = 1;
	printf("Testing:setCluster_elecZsum\n");
	rv = setCluster_elecZsum(ClLL, 2.9, ElecZA);
	assert(rv==0);
	rv = setCluster_elecZsum(ClLL, -1.1, ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(NULL, 1, ElecZA);
	assert(rv==-1);
	rv = setCluster_elecZsum(ClLL2,9.1, ElecZA);
	assert(rv==-1);

	printf("Testing:addToCluster_elecXsum\n");
	rv = addToCluster_elecXsum(ClLL,-2,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(NULL,1,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(ClLL2,1,ElecXF);
	assert(rv==-1);
	rv = addToCluster_elecXsum(ClLL,1,ElecXF);
	assert(rv==0);

	printf("Testing:addToCluster_elecYsum\n");
	rv = addToCluster_elecYsum(ClLL,-2, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(NULL,1, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(ClLL2,1, ElecYR);
	assert(rv==-1);
	rv = addToCluster_elecYsum(ClLL,1, ElecYR);
	assert(rv==0);

	printf("Testing:addToCluster_elecZsum\n");
	rv = addToCluster_elecZsum(ClLL,-2, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(NULL,1, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(ClLL2,1, ElecZA);
	assert(rv==-1);
	rv = addToCluster_elecZsum(ClLL,1, ElecZA);
	assert(rv==0);

	printf("Testing:getCluster_elecXsum\n");
	rvd = getCluster_elecXsum(NULL,ElecXF);
	assert(rvd==-1);
	rvd = getCluster_elecXsum(ClLL2,ElecXF);
	assert(rvd==-1);
	rvd = getCluster_elecXsum(ClLL,ElecXF);
	assert(rvd==2.1);

	printf("Testing:getCluster_elecYsum\n");
	rvd = getCluster_elecYsum(NULL,ElecYR);
	assert(rvd==-1);
	rvd = getCluster_elecYsum(ClLL2,ElecYR);
	assert(rvd==-1);
	rvd = getCluster_elecYsum(ClLL,ElecYR);
	assert(rvd==3.1);

	printf("Testing:getCluster_elecZsum\n");
	rvd = getCluster_elecZsum(NULL,ElecZA);
	assert(rvd==-1);
	rvd = getCluster_elecZsum(ClLL2,ElecZA);
	assert(rvd==-1);
	rvd = getCluster_elecZsum(ClLL,ElecZA);
	assert(rvd==3.9);

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
	rv = setCluster_elecXsum(ClLL, 3.3,ElecXF);
	assert(rv==0);
	rv = setCluster_elecXsum(ClLL, -1,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(NULL,2.5,ElecXF);
	assert(rv==-1);
	rv = setCluster_elecXsum(ClLL2,2.3,ElecXF);
	assert(rv==-1);

	printf("Testing:setCluster_elecYsum\n");
	rv = setCluster_elecYsum(ClLL, 3.3,ElecYR);
	assert(rv==0);
	rv = setCluster_elecYsum(ClLL, -1,ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(NULL,2.5,ElecYR);
	assert(rv==-1);
	rv = setCluster_elecYsum(ClLL2,2.3,ElecYR);
	assert(rv==-1);
	
	printf("Testing:setCluster_elecZsum\n");
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
	rv = getCluster_elecXid(ClLL);
	assert(rv==2);
	rv = getCluster_elecXid(ClLL2);
	assert(rv==-1);
	rv = getCluster_elecXid(NULL);
	assert(rv==-1);

	printf("Testing:getCluster_elecYid\n");
	rv = getCluster_elecYid(ClLL);
	assert(rv==2);
	rv = getCluster_elecYid(ClLL2);
	assert(rv==-1);
	rv = getCluster_elecYid(NULL);
	assert(rv==-1);

	printf("Testing:getCluster_elecZid\n");
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
	MidPoint mp2 = newMidPoint(-2, 1, 3, 298);
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

	printf("Testing:getNextMP\n");
	mp101 = getNextMP(NULL);
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
	rv = addToOrLL(NULL,1, mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, -1, mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 5, mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 0, NULL);
	assert(rv==-1);
	rv = addToOrLL(Arb7, 0, mp101);
	assert(rv==-1);
	rv = addToOrLL(Arb9, 0, mp101);
	assert(rv==0);
	rv = getElementsUsed(Arb9);
	assert(rv==1);
	OrderMagLL OMLL11 = getOrderLL(Arb9, 0);
	assert(getOMLL_size(OMLL11)==1);
	rv = addToOrLL(Arb9, 0, mp101);
	assert(rv==-1);
	MidPoint mp30 = newMidPoint(-3,121,32,55);
	rv = addToOrLL(Arb9, 0, mp30);
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
	
	deleteMidPoint((MidPoint)Temp4);
	deleteMidPoint(mp30);
	deleteArbArray(&Arb7);
	deleteArbArray(&Arb9);

	//////////////////////////////////////////////////////////////////////
	printf("Testing:setNeighNodeNew_p\n");
	NeighNode Nei99 = newNeighNode(99);
	setNeighNodeNew_p(Nei99, 0.2);
	printNeighNode(Nei99);
	
	printf("Testing:setNeighNode_t\n");
	setNeighNode_t(Nei99,32,1);
	printNeighNode(Nei99);
	
	printf("Testing:deleteNeighNode\n");
	rv = deleteNeighNode(Nei99);
	assert(rv==0);

}
