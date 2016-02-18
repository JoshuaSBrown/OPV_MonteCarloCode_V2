#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>

#include "cluster.h"
//#include "../../../MEM/mem.h"

struct _MidPoint{
	int id;
	int orderMag; //log(HopRate)
	//Neigbors between the mid points
	int nei1;
	int nei2;
	MidPoint next;
};

struct _OrderMagLL{
	int orderMag;
	int size;
	MidPoint start;
};

struct _Node{
	int id;
	double p; 
	Node next;
	int flag[6];
};

struct _NeighNode{
	int id;
	int hoplength;
	NeighNode next;
	Hop start;
};

struct _Hop{
	double p;
	double t;
	Hop next;
};

struct _Electrode{
	double sum;
	double alpha;
	int Charges;
	double FermiEnergy;
	void * HopRates;
	void * AdjacentSites;
};

struct _NeighLL{
	int numNeigh;
	NeighNode start;
};

struct _ClusterLL{
	int numNodes;
	int id;
	double clusterp;
	//time = tcluster = tprob*1/20  where tprob 
	//is the time to escape the cluster
	double time;
	Node startNode;
	NeighLL Neigh;
	int elecXF_id;
	int elecXB_id;
	int elecYR_id;
	int elecYL_id;
	int elecZB_id;
	int elecZA_id;
	Electrode elecXF;
	Electrode elecXB;
	Electrode elecYR;
	Electrode elecYL;
	Electrode elecZA;
	Electrode elecZB;
	ClusterLL next;
	double sum;
};

struct _ArbArray{
	int type;
	int reserved;
	int used;
	void * Num[0];
};

////////////////////////////////////////////////////////////////////////////////
//Tools for accessing Mid point array and Mid points
MidPoint newMidPoint(int order,int Mid_ID, int nei1, int nei2){

	if( nei1<0 || nei2<0 || Mid_ID<0 || nei1==nei2){
		return NULL;
	}

	MidPoint mp = (MidPoint) malloc(sizeof(struct _MidPoint));

	//Can not be combined with above if statement if the Midpoint is 
	//successfully created but there is a problem with one of the
	//parameters then it would need to be deleted. 
	if(mp==NULL ){
		return NULL;
	}

	mp->orderMag=order;
	mp->nei1=nei1;
	mp->nei2=nei2;
	mp->next=NULL;
	mp->id=Mid_ID;
	return mp;
}

int deleteMidPoint(MidPoint mp){
	if(mp==NULL){
		return -1;
	}
	free(mp);
	return 0;
}

int printMP(const_MidPoint mp){
	if(mp==NULL){
		return -1;
	}
	printf("Id: %d \t OrderMag %d \t nei1: %d \t nei2: %d\n",mp->id,mp->orderMag, mp->nei1, mp->nei2);
	return 0;
}

int getMP_order(const_MidPoint mp){
	if(mp==NULL){ //added condition
		return -1;
	}
	return mp->orderMag;
}

int getMP_id(const_MidPoint mp){
	if(mp==NULL){
		return -1;
	}
	return mp->id;
}

int setMP_id(MidPoint mp, int ID){
	if(mp==NULL||ID<0){
		return -1;
	}
	mp->id=ID;
	return 0;
}

int CompareNeiMidPoint(const_MidPoint mp1, MidPoint mp2) {

	if(mp1==NULL || mp2==NULL){
		return -1;
	}

	if(mp1->nei1==mp2->nei1){
		return 1;
	}
	if(mp1->nei2==mp2->nei2){
		return 1;
	}
	if(mp1->nei1==mp2->nei2){
		return 1;
	}
	if(mp1->nei2==mp2->nei1){
		return 1;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////
//Tools for accessing Order of Magnitude Link List
OrderMagLL newOrLL(int orderMag){
	OrderMagLL OMLL= (OrderMagLL) malloc(sizeof(struct _OrderMagLL));
	if(OMLL==NULL){
		return NULL;
	}
	OMLL->orderMag=orderMag;
	OMLL->size=0;
	OMLL->start=NULL;
	return OMLL;
}

int checkNewOrLL(OrderMagLL OMLL, MidPoint mp){ //use when adding newMidPoint to OMLL
	if(OMLL==NULL || mp==NULL){
		return -1; 
	}
	while(mp != NULL){
		if(OMLL->orderMag != mp->orderMag ){
			return -1;
		}
		mp=mp->next;
	}
	return 0;
}

int deleteOrLL(OrderMagLL OMLL){ //used only if no midpoints
	if(OMLL==NULL || OMLL->start!=NULL){
		return -1;
	}	
	free(OMLL);
	return 0;
}

int deleteAllOrLL(OrderMagLL OMLL){
	if(OMLL==NULL){
		return -1;
	}
	//Delete All the MidPoints in the Link List first
	MidPoint tempmp;
	MidPoint tempmp2;
	int rv;

	if (OMLL->start!=NULL) {
		tempmp=OMLL->start;
		if (tempmp->next!=NULL){
			tempmp2=tempmp->next;
			while(tempmp2!=NULL){
				rv = deleteMidPoint(tempmp);
				if(rv==-1){
					return -1;
				}
				tempmp=tempmp2;
				tempmp2=tempmp2->next;
			}
		}
		rv = deleteMidPoint(tempmp);
		if (rv==-1){
			return -1;
		}
	}
	free(OMLL);
	return 0;
}

int printOrLL(const_OrderMagLL OMLL) {

	if (OMLL==NULL){
		return -1;
	}
	MidPoint tempmp;	
	tempmp=OMLL->start;
	printf("\nOrderMag %d Size of LL: %d\n",OMLL->orderMag, OMLL->size);
	if (tempmp) {
		printMP(tempmp);
		while(tempmp->next){
			tempmp=tempmp->next;
			printMP(tempmp);

		}
	}
	return 0;
}

int getOMLL_size(const_OrderMagLL OMLL){
	if(OMLL==NULL){
		return -1;
	}

	return OMLL->size;
}

int getOMLL_order(const_OrderMagLL OMLL){
	if(OMLL==NULL){
		return -1;
	}

	return OMLL->orderMag;
}

MidPoint getOMLLstartMP(const_OrderMagLL OMLL){
	if(OMLL==NULL){
		return NULL;
	}
	return OMLL->start;
}

////////////////////////////////////////////////////////////////////////
//Tools for accessing Nodes
Node newNode(int N_ID){

	if(N_ID<0){
		return NULL;
	}

	Node Nod = (Node) malloc(sizeof(struct _Node));

	if(Nod==NULL){
		return NULL;
	}
	Nod->id=N_ID;
	Nod->p=0;
	Nod->flag[1]=0;
	Nod->flag[0]=0;
	Nod->flag[2]=0;
	Nod->flag[3]=0;
	Nod->flag[5]=0;
	Nod->flag[4]=0;
	Nod->next=NULL;
	return Nod;
}

int deleteNode(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	free(Nod);
	return 0;
}

int printNode(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	printf("Node id %d\n",Nod->id);
	printf("Pval %g\n",Nod->p);
	printf("Fro \t Beh \t Lef \t Rig \t Abo \t Bel\n");
	printf("%d \t %d \t %d \t %d \t %d \t %d\n",\
			Nod->flag[1], Nod->flag[0], Nod->flag[2], Nod->flag[3],\
			Nod->flag[5], Nod->flag[4]);
	return 0;
}

int getNode_id(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->id;
}

Node getNextNode(const_Node Nod){
	if(Nod==NULL){
		return NULL;
	}
	return Nod->next;
}

int setNode_id(Node Nod, int ID){
	if(Nod==NULL || ID<0){
		return -1;
	}
	Nod->id=ID;
	return 0;
}

double getNode_p(const_Node Nod){
	if(Nod==NULL){
		return -1.0;
	}
	return Nod->p;
}

int setNode_p(Node Nod, double pval){
	if(Nod==NULL || pval<0.0){
		return -1;
	}
	Nod->p=pval;
	return 0;
}

int getFlagFro(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[1];
}

int getFlagBeh(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[0];
}

int getFlagLef(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[2];
}

int getFlagRig(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[3];
}

int getFlagAbo(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[5];
}

int getFlagBel(const_Node Nod){
	if(Nod==NULL){
		return -1;
	}
	return Nod->flag[4];
}

int getFlag(const_Node Nod, int index){
	if(Nod==NULL|| index<0 || index>5){
		return -1;
	}
	return Nod->flag[index];
}

int setFlagFro(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[1]=1;
	return 0;
}

int setFlagBeh(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[0]=1;
	return 0;
}

int setFlagLef(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[2]=1;
	return 0;
}

int setFlagRig(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[3]=1;
	return 0;
}

int setFlagAbo(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[5]=1;
	return 0;
}
int setFlagBel(Node Nod){
	if(Nod==NULL){
		return -1;
	}
	Nod->flag[4]=1;
	return 0;
}

int setFlag(Node Nod, int index){
	if(Nod==NULL || index<0 || index>5){
		return -1;
	}
	Nod->flag[index]=1;
	return 0;
}

////////////////////////////////////////////////////////////////////////
//Tools for accessing Neigh Nodes
NeighNode newNeighNode(int N_ID){
	if(N_ID<0){
		return NULL;
	}
	NeighNode NeighNod = (NeighNode) malloc(sizeof(struct _NeighNode)); 

	if(NeighNod==NULL){
		return NULL;
	}
	NeighNod->next=NULL;
	NeighNod->id=N_ID;
	NeighNod->start=NULL;
	NeighNod->hoplength=0;
	return NeighNod;
}

int deleteNeighNode(NeighNode NeighNod){
	if(NeighNod==NULL){
		return -1;
	}

	if(NeighNod->start!=NULL){
		if(NeighNod->hoplength>0){
			printf("Hop length %d\n",NeighNod->hoplength);
			deleteAllHop(NeighNod->start);
			NeighNod->start=NULL;
			NeighNod->hoplength=0;
		}
	}

	free(NeighNod);
	return 0;
}

int printNeighNode(const_NeighNode NeighNod){
	if(NeighNod==NULL){
		return -1;
	}
	printf("NeighNode id %d HopLength %d\n",NeighNod->id,NeighNod->hoplength);

	Hop h = NeighNod->start;
	while(h!=NULL){
		printHop(h);
		h = h->next;
	}
	return 0;
}

int getNeighNode_id(const_NeighNode NeighNod){
	if(NeighNod==NULL){
		return -1;
	}
	return NeighNod->id;
}

int setNextNeighNode(NeighNode * NeighNod,NeighNode * NeighNod2){
	if(NeighNod==NULL || NeighNod2==NULL){
		return -1;
	}
	
	if(*NeighNod==NULL || *NeighNod2==NULL){ //don't know how to test
		return -1;
	}

	(*NeighNod)->next=(*NeighNod2);

	return 0;
}

NeighNode getNextNeigh(const_NeighNode NeighNod){
	if(NeighNod==NULL || NeighNod->next==NULL){ //second condition (added)?
		return NULL;
	}
	return NeighNod->next;
}

int setNeighNodeNew_p(NeighNode NeighNod, double pval){
	if (NeighNod==NULL || pval<0.0){
		return -1;
	}
	if (NeighNod->start==NULL){
		NeighNod->start = newHop();
		NeighNod->start->p=pval;
		NeighNod->hoplength++;
	}else{
		Hop h = NeighNod->start;
		while(h->next!=NULL){
			h = h->next;
		}
		h->next = newHop();
		h = h->next;
		h->p = pval;
		NeighNod->hoplength++;
	}
	return 0;
}

int setNeighNode_p(NeighNode NeighNod, double pval, int Elem){

	if (NeighNod==NULL || pval<0.0){
		return -1;
	}
	if (NeighNod->start==NULL || Elem>NeighNod->hoplength || Elem<1){
		return -1;
	}

	int inc=1;
	Hop h = NeighNod->start;
	while (inc<Elem){
		inc++;
		h = h->next;
	}
	h->p = pval;

	return 0;
}

int setNeighNode_t(NeighNode NeighNod, double time, int Elem){

	if (NeighNod==NULL || time<0.0){
		return -1;
	}

	if (NeighNod->start==NULL || Elem>NeighNod->hoplength || Elem<1){			return -1;
	}

	int inc=1;
	Hop h = NeighNod->start;
	printHop(h);
	while (inc<Elem){
		inc++;
		h = h->next;
		printHop(h);
	}
	h->t = time;

	return 0;
}

double getNeighNode_p(const_NeighNode NeighNod, int Elem){

	if (NeighNod==NULL){
		printf("NeighNode is NULL\n");
		return -1;
	}

	if (NeighNod->start==NULL || Elem>NeighNod->hoplength || Elem<1){
		printf("There is no hop array\n");
		return -1;
	}

	int inc=1;
	Hop h = NeighNod->start;

	while (inc<Elem){
		inc++;
		h = h->next;
	}

	return h->p;

}

int getNeighNode_hoplength(NeighNode NeighNod){
	if(NeighNod==NULL){
		return -1;
	}

	return NeighNod->hoplength;
}

double getNeighNode_t(NeighNode NeighNod, int Elem){

	if (NeighNod==NULL){
		return -1;
	}

	if (NeighNod->start==NULL || Elem>NeighNod->hoplength || Elem<1){
		return -1;
	}

	int inc=1;
	Hop h = NeighNod->start;
	while (inc<Elem){
		inc++;
		h = h->next;
	}
	return h->t;
}

int setNeighNode_id(NeighNode NeighNod, int ID){
	if(NeighNod==NULL || ID<0){
		return -1;
	}
	NeighNod->id=ID;
	return 0;
}

int setNeighNode_hoplength(NeighNode NeighNod, int hl){
	if(NeighNod==NULL || hl<0){
		return -1;
	}
	NeighNod->hoplength = hl;

	return 0;
}

int setNeighNode_hopstart(NeighNode NeighNod, Hop h){
	if(NeighNod==NULL || h==NULL){
		return -1;
	}
	NeighNod->start = h;
	return 0;
}
//////////////////////////////////////////////////////////////
Electrode newElectrode(void){

	Electrode el = (Electrode) malloc(sizeof(struct _Electrode));

	if(el==NULL){
		return NULL;
	}

	el->sum=0.0;
	el->alpha=0.0;
	el->Charges=0.0;
	el->FermiEnergy=0.0;
	el->HopRates = NULL;
	el->AdjacentSites = NULL;
	return el;
}

int deleteElectrode(Electrode * el){
	if(el==NULL){
		return -1;
	}
	if((*el)==NULL){
		return -1;
	}

	free(*el);
	*el=NULL;
	return 0;
}

int setElectrode_alpha(Electrode el, double alpha){
	if(el==NULL || alpha<=0){
		return -1;
	}

	el->alpha = alpha;
	return 0;
}

double getElectrode_alpha(Electrode el){
	return el->alpha;
}

int getElectrode_Charges(Electrode el){
	if(el==NULL){
		return -1;
	}

	return el->Charges;
}

int setElectrode_Charges(Electrode el, double NumCharge){
	if(el==NULL || NumCharge<0){
		return -1;
	}

	el->Charges = NumCharge;
	return 0;
}

int Electrode_addCharge(Electrode el){
	if(el==NULL){
		return -1;
	}
	el->Charges++;
	return 0;
}

int Electrode_minusCharge(Electrode el){
	if(el==NULL || el->Charges==0){
		return -1;
	}

	el->Charges--;
	return 0;
}


int setElectrode_FermiEnergy(Electrode el, double FermiE){
	if(el==NULL){
		return -1;
	}

	el->FermiEnergy = FermiE;
	return 0;
}

double getElectrode_FermiEnergy(Electrode el){
	if(el==NULL){
		return -1;
	}
	return el->FermiEnergy;
}

int setElectrode_HopRates(Electrode el, void * HopS){
	if(el==NULL || HopS==NULL ){
		return -1;
	}
	el->HopRates=HopS;
	return 0;
}

void * getElectrode_HopRates(Electrode el){ 
	if(el==NULL){
		return NULL;
	}
	return el->HopRates;
}

int setElectrode_AdjacentSites(Electrode el, void * AdjacentSites){
	if(el==NULL || AdjacentSites==NULL){
		return -1;
	}
	el->AdjacentSites = AdjacentSites;
	return 0;
}

void * getElectrode_AdjacentSites(Electrode el){
	if(el==NULL){
		return NULL;
	}
	return el->AdjacentSites;
}

int setElectrode_Sum(Electrode el, double sum){
	if(el==NULL || sum<0){
		return -1;
	}
	el->sum=sum;
	return 0;
}

double getElectrode_Sum(Electrode el){
	return el->sum;
}

int printElectrode(const_Electrode el){
	if(el==NULL){
		return -1;
	}
	printf("Electrode sum %g\n",el->sum);
	return 0;
}
////////////////////////////////////////////////////////////////////////
NeighLL newNeighLL(void){

	NeighLL neiLL = (NeighLL) malloc(sizeof(struct _NeighLL));

	if(neiLL==NULL){
		return NULL;
	}

	neiLL->numNeigh=0;
	neiLL->start=NULL;
	return neiLL;
}

int deleteNeighLL(NeighLL neighLL){
	if(neighLL==NULL){
		return -1;
	}

	free(neighLL);
	return 0;
}

int deleteNeighLLAll(NeighLL neighLL){
	if(neighLL==NULL){
		return -1;
	}

	NeighNode NeighNod;
	NeighNode tempNeighNod;
	NeighNod = neighLL->start;

	if(NeighNod!=NULL){
		while(NeighNod->next!=NULL){
			tempNeighNod=NeighNod->next;
			NeighNod->next=NULL;
			deleteNeighNode(NeighNod);
			NeighNod=tempNeighNod;
		}
		NeighNod->next=NULL;
		deleteNeighNode(NeighNod);
		neighLL->start=NULL;
	}

	deleteNeighLL(neighLL);
	return 0;
}

int printNeighLL(const_NeighLL neighLL){

	if(neighLL==NULL){
		return -1;
	}

	printf("NeighLL size %d\n",neighLL->numNeigh);

	NeighNode NeighNod;
	NeighNod = neighLL->start;

	if(NeighNod!=NULL){
		while(NeighNod->next!=NULL){
			printNeighNode(NeighNod);
			NeighNod = NeighNod->next;
		}
		printNeighNode(NeighNod);
	}
	return 0;
}

int getNeighLL_numNeigh( NeighLL Nei){
	if(Nei==NULL){
		return -1;
	}
	return Nei->numNeigh;
}

int setNeighLL_start(NeighLL Nei, NeighNode NeighNod){
	if(Nei==NULL || NeighNod==NULL){
		return -1;
	}

	Nei->start = NeighNod;
	return 0;
}

int setNeighLL_numNeigh(NeighLL Nei, int NumNei){
	if(Nei==NULL || NumNei<0){
		return -1;
	}

	Nei->numNeigh = NumNei;
	return 0;
}

////////////////////////////////////////////////////////////////////////
//Tools for accessing Hop
Hop newHop(void){

	Hop h = (Hop) malloc(sizeof(struct _Hop));

	if(h==NULL){
		return NULL;
	}

	h->t=0.0;
	h->p=0.0;
	h->next=NULL;
	return h;
}

int deleteAllHop(Hop h){//come back to deleteNeighNode

	if(h==NULL){
		return -1;
	}
	Hop temp;
	while( h->next!=NULL){
		temp = h->next;
		h->next=NULL; //why necessary?
		free(h);
		h = temp;
	}
	h->next=NULL; //necessary condition, why restated here?
	free(h);

	return 0;
}

int printHop(Hop h){

	if(h==NULL){
		return -1;
	}

	printf("Hop time %g pval %g\n",h->t,h->p);
	return 0;

}

int setHop_t(Hop h, double t){
	if(h==NULL || t<0){
		return -1;
	}
	h->t = t;
	return 0;
}

int setHop_p(Hop h, double p){
	if(h==NULL || p<0){
		return -1;
	}
	h->p = p;
	return 0;
}

int setHop_next(Hop h, Hop h2){
	if(h==NULL || h2==NULL){
		return -1;
	}
	h->next = h2;
	return 0;
}	
////////////////////////////////////////////////////////////////////////
//Tools for accessing Clusters and Cluster Link List
ClusterLL newClusterLL(int ID){

	if(ID<0 ){
		return NULL;
	}

	ClusterLL clLL = (ClusterLL) malloc(sizeof(struct _ClusterLL));

	if(clLL==NULL){
		return NULL;
	}
	clLL->id=ID;
	clLL->numNodes=0;
	clLL->clusterp=0.0;
	clLL->time=0;
	clLL->startNode=NULL;
	clLL->Neigh=NULL;

	clLL->elecXF_id = 0;
	clLL->elecXB_id = 0;
	clLL->elecYL_id = 0;
	clLL->elecYR_id = 0;
	clLL->elecZB_id = 0;
	clLL->elecZA_id = 0;

	clLL->elecXF=NULL;
	clLL->elecXB=NULL;
	clLL->elecYL=NULL;
	clLL->elecYR=NULL;
	clLL->elecZB=NULL;
	clLL->elecZA=NULL;
	clLL->next=NULL;
	clLL->sum=0;
	return clLL;
}


int deleteClusterLL(ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}
	free(clLL);
	return 0;
}

int deleteClusterLLNodes(ClusterLL clLL){

	if(clLL==NULL ){
		return -1;
	}

	printf("Cluster id %d Number of nodes %d\n",clLL->id,clLL->numNodes);
	Node tempNode;
	Node tempNode2;

	tempNode=clLL->startNode;
	if(tempNode!=NULL){
		if (tempNode->next!=NULL){
			tempNode2=tempNode->next;
			while(tempNode2!=NULL){
				tempNode->next=NULL;
				deleteNode(tempNode);
				tempNode=tempNode2;
				tempNode2=tempNode2->next;
			}
		}
		tempNode->next=NULL;
		deleteNode(tempNode);
	}else{
		free(clLL);
		return 0;
	}

	clLL->startNode=NULL;
	free(clLL);
	return 0;
}

int deleteAllClusterLL(ClusterLL *clLL){

	if((clLL)==NULL){
		return -1;
	}
	if((*clLL)==NULL){
		return -1;
	}

	//Delete All the Nodes in the cluster Link List 
	ClusterLL tempclLL;
	tempclLL=(*clLL)->next;
	Node tempNode;
	Node tempNode2;
	NeighLL neighLL;
	NeighNode tempNeighNode;
	NeighNode tempNeighNode2;

	while(tempclLL!=NULL) {

		printf("Cluster id %d number of Nodes %d\n",(*clLL)->id,(*clLL)->numNodes);
		//Deleting Nodes
		if ((*clLL)->startNode!=NULL) {
			tempNode=(*clLL)->startNode;
			if (tempNode->next!=NULL){
				tempNode2=tempNode->next;
				while(tempNode2!=NULL){
					tempNode->next=NULL;
					deleteNode(tempNode);						
					tempNode=tempNode2;
					tempNode2=tempNode2->next;
				}
			}
			tempNode->next=NULL;
			deleteNode(tempNode);
			(*clLL)->startNode=NULL;
		}

		//Deleting Neighboring Nodes
		neighLL = (*clLL)->Neigh;
		if (neighLL !=NULL){
			printf("Number of Neighbors %d",neighLL->numNeigh);
			if(neighLL->start!=NULL){
				tempNeighNode=neighLL->start;
				if(tempNeighNode->next!=NULL){
					tempNeighNode2=tempNeighNode->next;
					while(tempNeighNode2!=NULL){
						tempNeighNode->next=NULL;
						deleteNeighNode(tempNeighNode);
						tempNeighNode=tempNeighNode2;
						tempNeighNode2=tempNeighNode2->next;
					}
					tempNeighNode->next=NULL;
					deleteNeighNode(tempNeighNode);
				}
				neighLL->start=NULL;
			}
			//Delete neighLL
			deleteNeighLL(neighLL);
			(*clLL)->Neigh=NULL;
		}

		/*
		//Deleting Electrode
		if((*clLL)->elecXF!=NULL){
			deleteElectrode(&((*clLL)->elecXF));
			(*clLL)->elecXF=NULL;
		}
		if((*clLL)->elecXB!=NULL){
			deleteElectrode(&((*clLL)->elecXB));
			(*clLL)->elecXB=NULL;
		}
		if((*clLL)->elecYL!=NULL){
			deleteElectrode(&((*clLL)->elecYL));
			(*clLL)->elecYL=NULL;
		}
		if((*clLL)->elecYR!=NULL){
			deleteElectrode(&((*clLL)->elecYR));
			(*clLL)->elecYR=NULL;
		}
		if((*clLL)->elecZA!=NULL){
			deleteElectrode(&((*clLL)->elecZA));
			(*clLL)->elecZA=NULL;
		}
		if((*clLL)->elecZB!=NULL){
			deleteElectrode(&((*clLL)->elecZB));
			(*clLL)->elecZB=NULL;
		}
*/
		(*clLL)->next=NULL;
		free((*clLL));
		*clLL=NULL;
		(*clLL)=tempclLL;
		tempclLL=tempclLL->next;

	}
	if((*clLL)!=NULL) {
		//Deleting Nodes
		if((*clLL)->startNode!=NULL){
			tempNode=(*clLL)->startNode;
			if (tempNode->next!=NULL){
				tempNode2=tempNode->next;
				while(tempNode2!=NULL){
					tempNode->next=NULL;
					deleteNode(tempNode);
					tempNode=tempNode2;
					tempNode2=tempNode2->next;
				}
			}
			tempNode->next=NULL;
			deleteNode(tempNode);
			(*clLL)->startNode=NULL;
		}

		//Deleting Neighboring Nodes
		neighLL = (*clLL)->Neigh;
		if (neighLL !=NULL){
			if(neighLL->start!=NULL){
				tempNeighNode=neighLL->start;
				if(tempNeighNode->next!=NULL){
					tempNeighNode2=tempNeighNode->next;
					while(tempNeighNode2!=NULL){
						tempNeighNode->next=NULL;
						deleteNeighNode(tempNeighNode);
						tempNeighNode=tempNeighNode2;
						tempNeighNode2=tempNeighNode2->next;
					}
				}
				tempNeighNode->next=NULL;
				printf("Deleting Neigh Nodes part 2b\n");
				deleteNeighNode(tempNeighNode);
			}
			//Delete neighLL
			neighLL->start=NULL;
			deleteNeighLL(neighLL);
			(*clLL)->Neigh=NULL;
		}
	}

	/*
	//Deleting Electrode
	if((*clLL)->elecXF!=NULL){
		deleteElectrode(&((*clLL)->elecXF));
		(*clLL)->elecXF=NULL;
	}
	if((*clLL)->elecXB!=NULL){
		deleteElectrode(&((*clLL)->elecXB));
		(*clLL)->elecXB=NULL;
	}
	if((*clLL)->elecYR!=NULL){
		deleteElectrode(&((*clLL)->elecYR));
		(*clLL)->elecYR=NULL;
	}
	if((*clLL)->elecYL!=NULL){
		deleteElectrode(&((*clLL)->elecYL));
		(*clLL)->elecYL=NULL;
	}
	if((*clLL)->elecZA!=NULL){
		deleteElectrode(&((*clLL)->elecZA));
		(*clLL)->elecZA=NULL;
	}
	if((*clLL)->elecZB!=NULL){
		deleteElectrode(&((*clLL)->elecZB));
		(*clLL)->elecZB=NULL;
	}
*/
	//Freeing (*clLL)
	free((*clLL));
	*clLL=NULL;
	
	return 0;
}

ClusterLL getNextClusterLL(const_ClusterLL clLL){
	if(clLL==NULL){
		return NULL;
	}
	return clLL->next;
}

int printNodesClusterLL(const_ClusterLL clLL){

	if(clLL==NULL){
		return -1;
	}
	Node tempNode;

	printf("\nNodes in Cluster %d\n",clLL->id);
	tempNode=clLL->startNode;
	while(tempNode!=NULL){
		printf("%d \t",tempNode->id);
		tempNode=tempNode->next;
	}
	printf("\n");
	return 0;
}

int printClusterLL(const_ClusterLL clLL){

	if(clLL==NULL){
		return -1;
	}

	Node tempNode;
	NeighNode tempNeighNode;
	const_ClusterLL tempclLL;
	NeighLL neighLL;

	tempclLL=clLL;
	//Cycle through ClusterLL
	while(tempclLL!=NULL){
		//Cycle through Nodes in given LL
		printf("\nNodes within Cluster id %d Number of Nodes %d\n",tempclLL->id,tempclLL->numNodes);
		printf("Time %g\n",tempclLL->time);
		if(tempclLL->startNode==NULL){
			printf("No Nodes\n");
		} else {
			tempNode = tempclLL->startNode;
			while(tempNode!=NULL){
				//printf("%d\n",tempNode->id);	
				printNode(tempNode);
				tempNode=tempNode->next;
			}
		}

		if(tempclLL->Neigh==NULL){
			printf("No Neighboring node list\n");
		}else{

			neighLL = tempclLL->Neigh;
			printf("\nNodes Neighboring Cluster Number of Neighbors %d\n",neighLL->numNeigh);
			if(neighLL->start==NULL){
				printf("No Nodes within Neighbor list\n");
			}else{

				tempNeighNode = neighLL->start;
				while(tempNeighNode!=NULL){
					printNeighNode(tempNeighNode);
					//printf("%d\n",tempNeighNode->id);
					tempNeighNode=tempNeighNode->next;

				}
			}
		}

		if(tempclLL->elecXF!=NULL){
			printElectrode(tempclLL->elecXF);
		}
		if(tempclLL->elecXB!=NULL){
			printElectrode(tempclLL->elecXB);
		}
		if(tempclLL->elecYL!=NULL){
			printElectrode(tempclLL->elecYL);
		}
		if(tempclLL->elecYR!=NULL){
			printElectrode(tempclLL->elecYR);
		}
		if(tempclLL->elecZA!=NULL){
			printElectrode(tempclLL->elecZA);
		}
		if(tempclLL->elecZB!=NULL){
			printElectrode(tempclLL->elecZB);
		}
		tempclLL=getNextClusterLL(tempclLL);
	}

	return 0;
}

int setCluster_NeighLL(ClusterLL clLL, NeighLL NeiLL){
	if(clLL==NULL || NeiLL==NULL){
		return -1;
	}
	clLL->Neigh = NeiLL;
	return 0;
}

NeighLL getCluster_NeiLL(const_ClusterLL clLL){
	if(clLL==NULL){
		return NULL;
	}

	return clLL->Neigh;

}

int getCluster_id(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}
	return clLL->id;
}

int getCluster_numNodes(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}
	return clLL->numNodes;
}

int setCluster_elecXFid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecXF_id = Elec;
	return 0;

}

int setCluster_elecXBid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecXB_id = Elec;	
	
	return 0;

}

int setCluster_elecYLid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecYL_id = Elec;
	return 0;

}

int setCluster_elecYRid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecYR_id = Elec;
	return 0;

}

int setCluster_elecZBid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecZB_id = Elec;
	return 0;

}

int setCluster_elecZAid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( Elec<0 || Elec>1){
		return -1;
	}

	(*clLL)->elecZA_id = Elec;
	return 0;

}

int setCluster_elecXid(ClusterLL * clLL, int Elec){
	if(clLL==NULL || Elec<0 || Elec>1 ){
		return -1;
	}

	if(*clLL==NULL){
		return -1;
	}

	if(Elec==0){
		if( (*clLL)->elecXB_id!=0){
			//If the above statement is triggered
			//It means the electrode has already
			//been recorded
			//Elec value 0 back side 
			//Elec value 1 front side
			return -1;
		}else{
			(*clLL)->elecXB_id = 1;
			return 0;
		}
	}else{
		if((*clLL)->elecXF_id!=0){
			return -1;
		}else{
			(*clLL)->elecXF_id = 1;
			return 0;
		}
	}
}

int setCluster_elecYid(ClusterLL * clLL, int Elec){
	if(clLL==NULL || Elec<0 || Elec>1 ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}

	if(Elec==0){
		if((*clLL)->elecYL_id!=0){
			return -1;
		}else{
			(*clLL)->elecYL_id = 1;
			return 0;
		}
	}else{

		if((*clLL)->elecYR_id!=0){
			//If the above statement is triggered
			//Electrode has already been stored
			//Elec value 0 left side 
			//Elec value 1 right side
			return -1;
		}else{
			(*clLL)->elecYR_id = 1;
			return 0;
		}
	}

}

int setCluster_elecZid(ClusterLL * clLL, int Elec){
	if(clLL==NULL || Elec<0 || Elec>1 ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}

	if(Elec==0){
		if((*clLL)->elecZB_id!=0){
			//If the above statement is triggered
			//It meand the bottom electrode has already 
			//been recorded
			//Elec value 0 below side 
			//Elec value 1 above side
			return -1;
		}else{
			(*clLL)->elecZB_id = 1;
			return 0;
		}

	}else{
		if((*clLL)->elecZA_id!=0){
			return -1;
		}else{
			(*clLL)->elecZA_id = 1;
			return 0;
		}
	}
}

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecXF_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecXF_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecXF = Elec;
	return 0;
}

int setCluster_elecXB(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecXB_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecXB_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecXB = Elec;
	return 0;
}
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecYR_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecYR_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecYR = Elec;
	return 0;
}

int setCluster_elecYL(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecYL_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecYL_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecYL = Elec;
	return 0;
}

int setCluster_elecZA(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecZA_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecZA_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecZA = Elec;
	return 0;
}

int setCluster_elecZB(ClusterLL * clLL, Electrode Elec ){
	
	if(clLL==NULL ){
		return -1;
	}
	if(*clLL==NULL){
		return -1;
	}
	if( (*clLL)->elecZB_id!=1 && Elec==NULL){
		return -1;
	}
	if( (*clLL)->elecZB_id!=0 && Elec!=NULL){
		return -1;
	}

	(*clLL)->elecZB = Elec;
	return 0;
}

int setCluster_elecXsum(ClusterLL clLL, double sum, int Elec){
	if(clLL==NULL){
		printf("ERROR clLL is NULL\n");
		return -1;
	}
	if(sum<0.0 ){
		printf("Sum is %g should be greater than 0\n",sum);
		return -1;
	}
	if(Elec>1 || Elec<0){
		printf("Value of Elec should be either 0 or 1 it is %d\n",Elec);
		return -1;
	}

	if(Elec==0){
		if(clLL->elecXB==NULL){
			printf("ERROR elecXB does not exist\n");
			return -1;
		}else{
			clLL->elecXB->sum=sum;
		}
	}else{
		if(clLL->elecXF==NULL){
			printf("ERROR elecXF does not exist\n");
			return -1;
		}else{
			clLL->elecXF->sum=sum;
		}
	}
	return 0;
}

int setCluster_elecYsum(ClusterLL clLL, double sum, int Elec){
	if(clLL==NULL || sum<0.0){
		return -1;
	}
	if(Elec>1 || Elec<0){
		return -1;
	}

	if(Elec==0){
		if(clLL->elecYL==NULL){
			return -1;
		}else{
			clLL->elecYL->sum = sum;
		}
	}else{
		if(clLL->elecYR==NULL){
			return -1;
		}else{
			clLL->elecYR->sum = sum;
		}
	}
	return 0;
}

int setCluster_elecZsum(ClusterLL clLL, double sum, int Elec){
	if(clLL==NULL || sum<0.0 ){
		return -1;
	}
	if(Elec>1 || Elec<0){
		return -1;
	}

	if(Elec==0){
		if(clLL->elecZB==NULL){
			return -1;
		}else{
			clLL->elecZB->sum = sum;
		}
	}else{
		if(clLL->elecZA==NULL){
			return -1;
		}else{
			clLL->elecZA->sum = sum;
		}
	}
	return 0;
}

int addToCluster_elecXsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL || val<0.0 ){
		return -1;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecXB==NULL){
			return -1;
		}else{
			clLL->elecXB->sum+=val;
		}
	}else{
		if(clLL->elecXF==NULL){
			return -1;
		}else{
			clLL->elecXF->sum+=val;
		}
	}
	return 0;
}

int addToCluster_elecYsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL || val<0.0 ){
		return -1;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecYL==NULL){
			return -1;
		}else{
			clLL->elecYL->sum+=val;
		}
	}else{
		if(clLL->elecYR==NULL){
			return -1;
		}else{
			clLL->elecYR->sum+=val;
		}
	}
	return 0;
}

int addToCluster_elecZsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL || val<0.0 ){
		return -1;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecZB==NULL){
			return -1;
		}else{
			clLL->elecZB->sum+=val;
		}
	}else{
		if(clLL->elecZA==NULL){
			return -1;
		}else{
			clLL->elecZA->sum+=val;
		}
	}
	return 0;
}

int getCluster_elecidXB(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecXB_id;
}

int getCluster_elecidXF(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecXF_id;
}

int getCluster_elecidYL(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecYL_id;
}

int getCluster_elecidYR(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecYR_id;
}

int getCluster_elecidZB(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecZB_id;
}

int getCluster_elecidZA(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	return clLL->elecZA_id;
}

double getCluster_elecXsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
		return -1.0;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecXB==NULL){
			return -1;
		}else{
			return clLL->elecXB->sum;
		}
	}else{
		if(clLL->elecXF==NULL){
			return -1;
		}else{
			return clLL->elecXF->sum;
		}
	}
}

double getCluster_elecYsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
		return -1.0;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecYL==NULL){
			return -1;
		}else{
			return clLL->elecYL->sum;
		}
	}else{
		if(clLL->elecYR==NULL){
			return -1;
		}else{
			return clLL->elecYR->sum;
		}
	}
}

double getCluster_elecZsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
		return -1.0;
	}
	if(Elec<0 || Elec>1){
		return -1;
	}
	if(Elec==0){
		if(clLL->elecZB==NULL){
			return -1;
		}else{
			return clLL->elecZB->sum;
		}
	}else{
		if(clLL->elecZA==NULL){
			return -1;
		}else{
			return clLL->elecZA->sum;
		}
	}
}

int setCluster_id(ClusterLL clLL, int ID){
	if (clLL==NULL || ID<0){
		return -1;
	}
	clLL->id=ID;
	return 0;
}

double getCluster_Sum(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1.0;
	}
	return clLL->sum;
}

int addToClusterSum(ClusterLL clLL, double s){
	if(clLL==NULL || s<=0){
		return -1;
	}
	clLL->sum=clLL->sum+s;
	return 0;
}

int setClusterSum(ClusterLL clLL, double s){
	if(clLL==NULL || s<0){
		return -1;
	}
	clLL->sum=s;
	return 0;
}

int setNextClusterLL(ClusterLL clLL, ClusterLL Nex){
	if(clLL==NULL){
		return -1;
	}

	clLL->next=Nex;
	return 0;
}

Node getStartNode(const_ClusterLL clLL){	
	if(clLL==NULL){
		return NULL;
	}
	return clLL->startNode;
}

Node getCluster_Node( const_ClusterLL clLL, int id){
	if(clLL==NULL || id<0){
		printf("Incorrect id of node or cluster does not exist\n");
		return NULL;
	}

	Node tempNode = clLL->startNode;

	while(tempNode!=NULL){

		if (tempNode->id==id){
			return tempNode;
		}	

		tempNode=tempNode->next;
	}

	printf("Could not find Node with given id in cluster\n");
	return NULL;
}

NeighNode getStartNeigh(const_ClusterLL clLL){
	if(clLL==NULL || clLL->Neigh==NULL){
		return NULL;
	}
	return clLL->Neigh->start;
}

int getCluster_numNeigh(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}
	if(clLL->Neigh==NULL){
		return -1;
	}

	return clLL->Neigh->numNeigh;
}

double getCluster_p(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1.0;
	}
	return clLL->clusterp;
}

int setCluster_p(ClusterLL clLL, double p){
	if(clLL==NULL || p<0.0){
		return -1;
	}
	clLL->clusterp=p;
	return 0;
}

int checkCluster_elec(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1;
	}

	if (clLL->elecXF==NULL && clLL->elecXB==NULL &&\
			clLL->elecYL==NULL && clLL->elecYR==NULL  &&\
			clLL->elecZA==NULL && clLL->elecZB==NULL ){
		return 0;
	}

	return 1;

}

int getCluster_elecXid(const_ClusterLL clLL){
	if(clLL==NULL || (clLL->elecXF==NULL && clLL->elecXB==NULL)){
		return -1;
	}

	if(clLL->elecXF!=NULL && clLL->elecXB!=NULL){
		return 2;
	}else if(clLL->elecXF!=NULL){
		return 1;
	}else {
		return 0;
	}
}

int getCluster_elecYid(const_ClusterLL clLL){
	if(clLL==NULL || (clLL->elecYL==NULL && clLL->elecYR==NULL)){
		return -1;
	}

	if(clLL->elecYR!=NULL && clLL->elecYL!=NULL){
		return 2;
	}else if(clLL->elecYR!=NULL){
		return 1;
	}else{
		return 0;
	}
}

int getCluster_elecZid(const_ClusterLL clLL){
	if(clLL==NULL || (clLL->elecZB==NULL && clLL->elecZA==NULL)){
		return -1;
	}

	if(clLL->elecZB!=NULL && clLL->elecZA!=NULL){
		return 2;
	}else if(clLL->elecZA!=NULL){
		return 1;
	}else{
		return 0;
	}
}

int setCluster_time(ClusterLL clLL, double t){
	if(clLL==NULL){
		return -1;
	}

	clLL->time=t;
	return 0;
}

double getCluster_time(const_ClusterLL clLL){
	if(clLL==NULL){
		return -1.0;
	}
	return clLL->time;
}

////////////////////////////////////////////////////////////////////////////////////
int addNodeEndClusterLL(ClusterLL clLL, int nei){

	if(clLL==NULL || nei<0){
		return -1;
	}

	Node tempNode;
	tempNode=clLL->startNode;

	if (tempNode==NULL){
		clLL->numNodes++;
		clLL->startNode=newNode(nei);
		return 0;
	}

	while (tempNode->next!=NULL){
		tempNode=tempNode->next;
	}
	clLL->numNodes++;
	tempNode->next=newNode(nei);
	return 0;
}

int addClusterLLNode(ClusterLL *clLL, Node nd){

	if(clLL==NULL || nd==NULL){
		return -1;
	}

	Node tempNode;
	tempNode = (*clLL)->startNode;

	if(tempNode==NULL){
		(*clLL)->startNode = nd;
	}else{
		while(tempNode->next!=NULL){
			tempNode = tempNode->next;
		}
		tempNode->next = nd;
	}

	return 0;
}

int addNodeNewCluster(ClusterLL * clLL, MidPoint mp){

	if( (*clLL)==NULL){
		printf("ERROR clLL is NULL\n");
		return -1;
	}
	if( mp==NULL ){
		printf("ERROR mp is NULL\n");
		return -1;
	}
	if((*clLL)->next!=NULL){
		printf("Neighboring Cluster is NULL\n");
		return -1;
	}

	ClusterLL clLL2;
	int ClusterId = (*clLL)->id+1;
	(*clLL)->next= newClusterLL(ClusterId);
	clLL2=(*clLL)->next;
	clLL2->numNodes=clLL2->numNodes+2;

	clLL2->startNode = newNode(mp->nei1);
	Node tempNode = clLL2->startNode;
	tempNode->next= newNode(mp->nei2);
	return 0;
}

int addNodeToCluster( ClusterLL clLL, MidPoint mp){

	if(clLL==NULL || mp==NULL){
		return -1;
	}

	//Check first cluster. If first cluster does not exist 
	//means no other clusters exist so we must create one
	if (clLL->startNode==NULL) {

		//New cluster should have an initial size of 0
		clLL->numNodes =clLL->numNodes+2;
		clLL->startNode = newNode(mp->nei1);
		Node tempNode = clLL->startNode;
		tempNode->next = newNode(mp->nei2);
		return 0;
	}

	ClusterLL tempclLL=clLL;
	Node tempNode= clLL->startNode;

	int ClusterId = clLL->id;
	int unfoundnei;

	while(tempNode->id!=mp->nei1 && tempNode->id!=mp->nei2){

		//Reach end of nodes move to next link list
		if(tempNode->next==NULL){
			//If next ClusterLL is NULL then must create a new one
			if(tempclLL->next==NULL){
				addNodeNewCluster(&tempclLL, mp);
				return 0;

			} else {
				//Move to the next ClusterLL
				tempclLL = tempclLL->next;
				ClusterId = tempclLL->id;
				tempNode =  tempclLL->startNode;
			}
		}else{
			tempNode=tempNode->next;
		}
	}

	//If the code gets to this point it means that one of the
	//nodes already exists in the cluster. we need to make 
	//sure that the other node does not exist in a separate cluster
	if (tempNode->id==mp->nei1){
		unfoundnei=mp->nei2;
	}else{
		unfoundnei=mp->nei1;
	}

	ClusterLL foundclLL=tempclLL;

	Node endNodeFoundLL=tempNode;
	//Search within the rest of the current cluster
	while(tempNode->next!=NULL){
		if(tempNode->next->id==unfoundnei){
			//This means both sites are within the same
			//cluster already so nothing is done
			return -1;
		}
		tempNode=tempNode->next;
		endNodeFoundLL=tempNode;
	}

	//Now we are going to cycle through the rest of the 
	//clusters 
	ClusterLL PrevClusterLL = tempclLL;

	tempclLL=tempclLL->next;
	//Make sure the next cluster exists if doesn't add to
	//the site to the formally found node;
	if(tempclLL==NULL){
		addNodeEndClusterLL(foundclLL, unfoundnei);
		return 0;
	}
	ClusterId=tempclLL->id;
	tempNode= tempclLL->startNode;

	while(tempNode->id!=unfoundnei){
		if(tempNode->next==NULL){
			//If next ClusterLL is NULL then must 
			//add second node to the found cluster
			if(tempclLL->next==NULL){
				addNodeEndClusterLL(foundclLL, unfoundnei);
				return 0;

			} else {
				//Move to the next ClusterLL
				PrevClusterLL = tempclLL;
				tempclLL = tempclLL->next;
				ClusterId = tempclLL->id;
				tempNode = tempclLL->startNode;
			}
		}else{
			tempNode=tempNode->next;
		}
	}

	//If we get to this point it means that we found
	//the second node on a second ClusterLL. This
	//means that the two different ClusterLL must
	//be consolidated. To establish a protocal we will 
	//say that the higher id will be merged
	//with the lower cluster id so that the higher
	//one no longer exists. 
	//Attaching list to end of found node
	endNodeFoundLL->next=  tempclLL->startNode;
	//Increasing size of cluster to account for new nodes
	printf("\nJoining Clusters %d size: %d and %d Size %d\n",foundclLL->id,foundclLL->numNodes,\
			tempclLL->id, tempclLL->numNodes);

	//Merging can only be done if the electrodes 
	//and neighLL have not yet been assigned;
	if(tempclLL->Neigh!=NULL || tempclLL->elecXF!=NULL || tempclLL->elecXB!=NULL ||\
			tempclLL->elecYL!=NULL || tempclLL->elecYR!=NULL ||\
			tempclLL->elecZA!=NULL || tempclLL->elecZB!=NULL ||\
			foundclLL->Neigh!=NULL || foundclLL->elecXF!=NULL || foundclLL->elecXB!=NULL ||\
			foundclLL->elecYL!=NULL || foundclLL->elecYR!=NULL ||\
			foundclLL->elecZA!=NULL || foundclLL->elecZB!=NULL){
		return -1;
	}

	tempNode = tempclLL->startNode;
	printf("Nodes in %d\n",tempclLL->id);
	while(tempNode!=NULL){
		tempNode=tempNode->next;
	}

	tempNode= foundclLL->startNode;
	printf("Nodes in %d\n",foundclLL->id);

	do{
		tempNode=tempNode->next;
	} while (tempNode->next!=NULL);

	foundclLL->numNodes=foundclLL->numNodes+tempclLL->numNodes;
	endNodeFoundLL->next= tempclLL->startNode;
	printf("Final Cluster %d size: %d\n",foundclLL->id,foundclLL->numNodes);

	//By passing second ClusterLL
	PrevClusterLL->next=tempclLL->next;
	//Deleting the cluster LL without deleting nodes etc
	deleteClusterLL(tempclLL);

	tempNode= foundclLL->startNode;
	printf("Nodes in Cluster %d\n",foundclLL->id);
	while(tempNode!=NULL){
		tempNode=tempNode->next;
	}

	return 0;
}

int addNeighNodeToCluster( ClusterLL* clLL, int Neigh_ID){

	if((*clLL)==NULL || Neigh_ID<0){
		return -1;
	}

	if ((*clLL)->Neigh==NULL){
		(*clLL)->Neigh = newNeighLL();
		(*clLL)->Neigh->start = newNeighNode(Neigh_ID);
		(*clLL)->Neigh->numNeigh++;
		return 0;
	}

	NeighNode tempNeigh =  (*clLL)->Neigh->start;

	while(tempNeigh->id!=Neigh_ID){

		if(tempNeigh->next==NULL){
			//Add node at the end of the link list
			tempNeigh->next = newNeighNode(Neigh_ID);
			(*clLL)->Neigh->numNeigh++;
			return 0;
		}
		tempNeigh = tempNeigh->next;
	};

	//Means that the Node was already added to the list
	return -1;
}

int addPvalClusterNode(ClusterLL clLL, int Node_ID, double pval){

	if(clLL==NULL || Node_ID<0 || pval<0.0){
		return -1;
	}
	//Search through Cluster till find Node_ID
	int N_ID;
	Node tempNode = getStartNode(clLL);
	if (tempNode==NULL){
		return -1;
	}
	N_ID = getNode_id(tempNode);
	while( N_ID!=Node_ID){
		if(tempNode->next==NULL){
			return -1;
		}
		tempNode=tempNode->next;
		N_ID = getNode_id(tempNode);
	}

	//Used to calculate the total time for a charge
	//to exist hop anywhere within the cluster
	clLL->time=1/pval+clLL->time;
	setNode_p(tempNode, pval+getNode_p(tempNode));
	return 0;
}

int addPvalClusterNeighNode(ClusterLL clLL, int Node_ID, double pval){

	if(clLL==NULL || Node_ID<0 || pval<0.0 ){
		return -1;
	}

	if(clLL->Neigh==NULL){
		return -1;
	}

	//Search through Cluster till find Node_ID
	int Neigh_id;
	NeighNode tempNeighNode =  getStartNeigh(clLL);
	if(tempNeighNode==NULL){
		return -1;
	}
	Neigh_id = getNeighNode_id(tempNeighNode);
	while( Neigh_id!=Node_ID){
		if (tempNeighNode->next==NULL){
			return -1;
		}
		tempNeighNode = tempNeighNode->next;
		Neigh_id = getNeighNode_id(tempNeighNode);
	}

	setNeighNodeNew_p(tempNeighNode, pval);
	return 0;
}

/////////////////////////////////////////////////////////////////////////
//Tools for accessing ArbArray
ArbArray newArbArray(int len, int type){

	if(len<0 || type>=3 || type<0){
		return NULL;
	}

	ArbArray Arb;	

	//Each element contains an Order of Magnitude Link List
	if(type==0) {
		Arb = (ArbArray) malloc(sizeof(struct _ArbArray)+sizeof(OrderMagLL)*len);
		//Each element contains a Cluster Link List
	} else if (type==1) {
		Arb = (ArbArray) malloc(sizeof(struct _ArbArray)+sizeof(ClusterLL)*len);
		//Each element contains a mid point
	} else {
		Arb = (ArbArray) malloc(sizeof(struct _ArbArray)+sizeof(MidPoint)*len);
	}

	if(Arb==NULL){
		return NULL;
	}
	Arb->reserved=len;
	Arb->used=0;
	Arb->type=type;
	int i;
	for(i=0;i<len;i++){
		Arb->Num[i]=NULL;
	}

	return Arb;
}

int ArbArrayCheck(const_ArbArray Arb){
	if(Arb==NULL){
		return -1;
	}
	printf("Arb Type Check %d\n",Arb->type);
	return 0;
}

int deleteAllMidPointArray(ArbArray * Arb){

	if(Arb==NULL){
		return -1;
	}
	if((*Arb)==NULL){
		return -1;
		
	}
	
	if ((*Arb)->type==2){
		
		int i;
		for(i=0;i<(*Arb)->reserved;i++){
			if((*Arb)->Num[i]!=NULL){
				deleteMidPoint((*Arb)->Num[i]);
			}
		}
		free(*Arb);
		*Arb=NULL;
		return 0;
	}
	printf("ERROR not a midpoint array!\n");
	return -1;
}


int deleteArbArray(ArbArray * Arb){

	if(Arb==NULL){
		return -1;
	}
	if((*Arb)==NULL){
		return -1;
	}
	int rv;
	printf("Deleting Arb Array\n");
	printf("Reserved %d\n",(*Arb)->reserved);
	printf("type %d\n",(*Arb)->type);
	//OrderMagLL
	if((*Arb)->type==0){
		int i;
		for(i=0;i<(*Arb)->reserved;i++){
			printf("Cycling through reserved\n");
			if ((*Arb)->Num[i]!=NULL){
				printf("Non NULL element found so deleting\n");
				deleteOrLL((OrderMagLL) (*Arb)->Num[i]);
				(*Arb)->Num[i]=NULL;
			}
		}
		free((*Arb));
		*Arb=NULL;
		return 0;
		//ClusterLL
	} else if ((*Arb)->type==1) {
		//Must Cycle through all the elements in the array and get rid of 
		//all the clusters first
		int i;
		ClusterLL clLL;
		for(i=0;i<(*Arb)->reserved;i++){
			if ((*Arb)->Num[i]!=NULL){
				printf("Deleting Array\n");
				clLL = (*Arb)->Num[i];
				rv = deleteAllClusterLL(&clLL );
				(*Arb)->Num[i]=NULL;
			}
		}
		free(*Arb);
		*Arb = NULL;
		return 0;
	} else if ((*Arb)->type==2){
		
		//WARNING we can not delete the midpoints here
		//because more than one arbitrary array may use them
		//They need to be delted externally
		/*
		int i;
		for(i=0;i<Arb->reserved;i++){
			if(Arb->Num[i]!=NULL){
				deleteMidPoint(Arb->Num[i]);
			}
		}
		*/
		free(*Arb);
		*Arb=NULL;
		return 0;
	}

	return -1;
}

void * getArbElement(const_ArbArray Arb, int element){
	if(Arb==NULL || element<0 || element>=Arb->reserved){
		return NULL;
	}
	return Arb->Num[element];
}

int NullArbElement(ArbArray Arb, int element){

	if(element<0){
		return -1;
	}

	if(Arb->Num[element]!=NULL){
		Arb->used--;
	}
	Arb->Num[element]=NULL;

	return 0;

}

int setArbElement(ArbArray Arb, int element, void * ptr){

	if(Arb==NULL || element<0 || element>=Arb->reserved || ptr==NULL){
		return -1;
	}

	if(Arb->Num[element]!=NULL){
		Arb->Num[element]=ptr;
	}else{
		Arb->used++;
		Arb->Num[element]=ptr;
	}
	return 0;
}

int printArbArray(const_ArbArray Arb, ...) {
	if(Arb==NULL){
		return -1;
	}

	if(Arb->type==0) {
		int i;
		printf("Printing Array of Link Lists\n");
		printf("Total number of Linked Lists: %d\n",Arb->used);
		for(i=0;i<Arb->used;i++) {
			printOrLL((OrderMagLL) Arb->Num[i]);
		}	
	}else if(Arb->type==1) {
		va_list argument;
		va_start (argument, Arb);
		int orderLow = va_arg(argument, int);
		int element;
		ClusterLL tempclLL;
		for(element=0;element<Arb->reserved;element++){
			printf("\nOrder of Magnitude %d\n",element+orderLow);
			tempclLL = (ClusterLL) getArbElement(Arb, element);
			if (tempclLL==NULL){
				printf("Cluster does not exist.\n");
			} else {
				printClusterLL(tempclLL);
			}
		}
		va_end(argument);
	}else if(Arb->type==2) {
		int i;
		for(i=0;i<Arb->used;i++){
			printMP((MidPoint) Arb->Num[i]);
		}
	}

	return 0;
}

int getElementsUsed(const_ArbArray Arb){
	if(Arb==NULL){
		return -1;
	}
	return Arb->used;
}

int getElementsReserved(const_ArbArray Arb){
	if(Arb==NULL){
		return -1;
	}
	return Arb->reserved;
}

int removeArbElement(ArbArray Arb, int element){
	if(Arb==NULL || element<0 || element>=Arb->reserved){
		return -1;
	}

	if(Arb->Num[element]!=NULL){
		Arb->used--;
		Arb->Num[element]=NULL;
		return 0;
	}
	return -1;
}

MidPoint getMP(const_ArbArray Arb,int element) {
	if(Arb==NULL ){
		return NULL;
	}
	if (Arb->type!=2 || element<0 || element>=Arb->reserved){
		return NULL;
	}
	MidPoint mp = (MidPoint) Arb->Num[element];
	return mp;
}

MidPoint getNextMP(const_MidPoint mp){
	if(mp==NULL){
		return NULL;
	}

	return mp->next;
}

int getMPnei1(const_ArbArray Arb, int Mid_ID){
	if(Arb==NULL ){
		return -1;
	}
	if( Mid_ID<0 || Mid_ID>=Arb->reserved || Arb->type!=2){
		return -1;
	}
	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if( mp!=NULL){
		return mp->nei1;
	}
	return -1;
}

int getMPnei2(const_ArbArray Arb, int Mid_ID){
	if(Arb==NULL){
		return -1;
	}
	if(Mid_ID<0 || Mid_ID>=Arb->reserved || Arb->type!=2){
		return -1;
	}


	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if (mp!=NULL){
		return mp->nei2;
	}
	return -1;
}

int getMPOrder(const_ArbArray Arb, int Mid_ID){
	if(Arb==NULL){
		return -1;
	}
	if(Mid_ID<0 || Mid_ID>=Arb->reserved || Arb->type!=2){
		return -1;
	}
	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if(mp!=NULL){
		return mp->orderMag;
	}
	return -1;
}

int addToOrLL(ArbArray Arb, int element , MidPoint mp){
	//Ensure that input parameters are error free
	if(Arb==NULL || element<0 || element>=Arb->reserved || mp==NULL){
		return -1;
	}
	if(Arb->type!=0){
		return -1;
	}
	//increment the size of the order of magnitude link list
	OrderMagLL OMLL = (OrderMagLL) Arb->Num[element];
	if(OMLL==NULL){
		Arb->used++;
		int Order = getMP_order(mp);
		printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
		printf("Order %d\n",Order);
		Arb->Num[element]=(void *) newOrLL(Order);
		OMLL = (OrderMagLL) Arb->Num[element];	
	}
	//Create a temporary midpoint
	MidPoint tempmp;
	tempmp=OMLL->start;

	if (tempmp==NULL) {
		OMLL->start=mp;
		OMLL->size++;
		return 0;

	}
	//Cycle through the mid points to find the end of the link list and to
	//ensure that the mid point has not already been added. 

	if(tempmp->id==mp->id){
		return -1;
	}

	while(tempmp->next){
		if(tempmp->id==mp->id){
			return -1;
		}
		tempmp=tempmp->next;
	}
	tempmp->next=mp;
	OMLL->size++;

	return 0;
}

OrderMagLL getOrderLL(const_ArbArray Arb, int element){
	if(Arb==NULL|| element>=Arb->reserved || element<0 || Arb->type!=0){
		return NULL;
	}
	return (OrderMagLL) Arb->Num[element];
}

int setDefaultArbElem(ArbArray Arb, int element, int orderMag){
	if(Arb==NULL || element<0 || element>=Arb->reserved){
		return -1;
	}
	if(Arb->type!=0){
		return -1;
	}

	//Setting default values
	if(Arb->Num[element]==NULL){
		Arb->used++;
	}else{
		//If not NULL need to remove
		//the old version
		deleteOrLL(Arb->Num[element]);
	}
	Arb->Num[element]=newOrLL(orderMag);
	return 0;
}
