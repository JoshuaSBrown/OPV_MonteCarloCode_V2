#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clustersitenode.h"
#include "../SITENODE/sitenode.h"
#include "../CLUSTER/cluster.h"
#include "../CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../ARBARRAY/arbarray.h"
#include "../NODE/node.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"
#include "../ERROR/error.h"

int printSNfun(SiteNode sn, int (*printfunc)(void *)) {

  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in printSNfun\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
	if(printfunc==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR printfunc is NULL in printSNfun\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
  #endif

	printf("\ninitE: %d\n",getInitE(sn));
	printf("dwellStatus: %d\n", getDwelStat(sn));
	printf("visitFreq: %g\n",getVisFreq(sn));
	printf("visit: %g\n",getVis(sn));
	printf("energy: %g\n",getEnergy(sn));
	if(getType(sn)==0){
		//Should be calling printPoint
		printfunc( getPoint(sn));
	}else if(getType(sn)==1){
		//Should be calling printClusterLL twice
		printfunc(getPoint(sn));
		printfunc(getClusterList(sn));
	}
	printf("\n");
	return 0;
}
/*
int printSNarray(SNarray snA, int (*printfunc)(void *)) {

	if(snA==NULL || printfunc==NULL){
		return -1;
	}

	printf("\nlength: %d\n",getAlen(snA));
	printf("width: %d\n",getAwid(snA));
	printf("height: %d\n",getAhei(snA));
	printf("total: %d\n",getAtotal(snA));
	printf("**************\n");

	SiteNode sn;
	int i;
	for(i=0;i<getAtotal(snA);i++) {
		printf("\nNode %d",i);

		sn = getSNwithInd(snA,i);
		if(getType(sn)==0){
			printSN(sn, &printPoint);
		} else {
			printSN(sn, printfunc);
		}
	}
	return 0;
}
*/
int PrintFile_xyz( int orderLow, SNarray snA, ArbArray *ArClLL, char * FileName){

	if(ArClLL==NULL || snA==NULL){
		return -1;
	}

	if( (*ArClLL)==NULL){
		return -1;
	}

	printf("Printing %s.xyz\n",FileName);
	//int length;
	int element;
	int NodeID;
	int i, j, k;
	int order;
	int TotSize=0;
	ClusterLL ClLL;
	Node Nod;
	SiteNode sn;
	//length = strlen(FileName);

	//strcat(FileName,".xyz");
	char buf[256];
	snprintf(buf,sizeof buf,"%s%s",FileName,"Cluster.xyz");

	FILE *Out;
	if(( Out = fopen( buf,"w"))==NULL){
		printf("Error! unable to open .xyz file\n");
	}else {
		//Calculate the number of Nodes total
		for( element=0; element<getElementsReserved(*ArClLL);element++){
			ClLL = (ClusterLL) getArbElement((*ArClLL), element);
			while(ClLL!=NULL){
				printf("Number of Nodes %d\n",getCluster_numNodes(ClLL));
				TotSize+=getCluster_numNodes(ClLL);
				ClLL=getNextClusterLL(ClLL);
			}
		}
		fprintf(Out,"%d\n\n",TotSize);
		for( element=0; element<getElementsReserved((*ArClLL));element++){
			ClLL = (ClusterLL) getArbElement((*ArClLL), element);
			order = orderLow+element;
			while(ClLL!=NULL){
				//Cycle through Nodes in given LL
				if(getStartNode(ClLL)==NULL){
					printf("No Nodes\n");
				} else {
					Nod = getStartNode(ClLL);
					while(Nod!=NULL){
						NodeID=getNode_id(Nod);
						getLoc(&i, &j, &k, NodeID, snA);
						sn=getSN(snA,i,j,k);
						fprintf(Out,"C\t%f\t%f\t%f\t%f\t%f\t%f\n",\
								(double)i,(double)j,(double)k,(double)order,(double)getCluster_id(ClLL),getEnergy(sn));
						Nod=getNextNode(Nod);
					}
				}
				ClLL=getNextClusterLL(ClLL);
			}

		}
		fclose(Out);
	}
	return 0;
}

int PrintNeighFile_xyz( int orderLow, SNarray snA, ArbArray *ArClLL,  char * FileName){

	if(ArClLL==NULL || snA==NULL){
		return -1;
	}

	if(*ArClLL==NULL){
		return -1;
	}

	printf("Printing %s.xyz\n",FileName);
	//int length;
	int element;
	int NodeID;
	int i, j, k;
	int order;
	int TotSize=0;
	ClusterLL ClLL;
	NeighNode NeighNod;
	SiteNode sn;
	//length = strlen(FileName);

	//strcat(FileName,".xyz");
  char buf[256];
	snprintf(buf, sizeof buf,"%s%s",FileName,"NeighCluster.xyz");
	FILE * Out;
	if(( Out = fopen( buf,"w"))==NULL){
		printf("Error! unable to open .xyz file\n");
	}else {
		//Calculate the number of Nodes total
		for( element=0; element<getElementsReserved(*ArClLL);element++){
			ClLL = (ClusterLL) getArbElement(*ArClLL, element);
			while(ClLL!=NULL){
				TotSize+=getCluster_numNeigh(ClLL);
				ClLL=getNextClusterLL(ClLL);
			}
		}

		fprintf(Out,"%d\n\n",TotSize);
		for( element=0; element<getElementsReserved(*ArClLL);element++){
			ClLL = (ClusterLL) getArbElement(*ArClLL, element);
			order = orderLow+element;
			while(ClLL!=NULL){
				//Cycle through Nodes in given LL
				if(getStartNeigh(ClLL)==NULL){
					printf("No Nodes\n");
				} else {
					NeighNod =  getStartNeigh(ClLL);
					while(NeighNod!=NULL){
						NodeID=getNeighNode_id(NeighNod);
						getLoc(&i, &j, &k, NodeID, snA);
						sn=getSN(snA,i,j,k);
						fprintf(Out,"N\t%f\t%f\t%f\t%f\t%f\t%f\n",\
								(double)i,(double)j,(double)k,(double)order,(double)getCluster_id(ClLL),getEnergy(sn));
						NeighNod=getNextNeigh(NeighNod);
					}
				}
				ClLL=getNextClusterLL(ClLL);
			}

		}
		fclose(Out);
	}
	return 0;
}
/*
double getsum(const_SiteNode sn){

	if(sn==NULL){
		return -1;
	}

	if ( getType(sn)==0) {
		Point pt = (Point) getPoint(sn);
		return getsumPt(sum);
	} else {
		//Grabs the sum from the neighbors
		//Tells the time it would take to hop
		//off on average
		ClusterLL ClLL = (ClusterLL) sn->NeighList;
		return getClusterSum(NeighClLL);
	}

}
*/

int OccAllNeiCluster(SNarray snA,int i1,int j1,int k1) {

	if(snA==NULL){
		return -1;
	}
	//Need to cycle through the neighbor cluster link list
	//And determine if any of the sites are unoccupied
	SiteNode sn = getSN(snA,i1,j1,k1);
	SiteNode tempsn;
	printf("Type %d\n",getType(sn));
	if ( getType(sn)!=1 || getClusterList(sn)==NULL){
		return -1;
	}

	int NeighID;
	int i, j, k;

	ClusterLL ClLL = (ClusterLL) getClusterList(sn);
	NeighNode NeighNod = (NeighNode) getStartNeigh(ClLL);

	if( NeighNod==NULL){
		return -1;
	}

	while (NeighNod!=NULL){
		NeighID=getNeighNode_id(NeighNod);
		getLoc(&i, &j, &k, NeighID, snA);
		tempsn=getSN(snA, i, j, k);
		if(getDwelStat(tempsn)==-1){
			//If any of the sites has a value of -1
			//it means there is a Neigh site that is unoccupied
			printf("Site that is unoccupied %d\n",NeighID);
			return 0;
		}
		NeighNod=getNextNeigh(NeighNod);
	}
	//All sites are occupied
	printf("All sites are occupied\n");
	return 1;
}


int OccAllCluster(SNarray snA,int i1,int j1,int k1) {
  #ifdef _ERROR_CHECKING_ON_	
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if( i1<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i1 is less than 0 in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
    
  }
  if( i1>=getAlen(snA) ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i1 is greater or equal to getAlen(snA) in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(j1<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j1 is less than 0 in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(j1>=getAwid(snA)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j1 is greater or equal to getAwid(snA) in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
	if(k1<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k1 is less than 0 in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(k1>=getAhei(snA) ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k1 is greater or equal to getAhei(snA) in OccAllCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//Need to cycle through the neighbor cluster link list
	//And determine if any of the sites are unoccupied
	SiteNode sn = getSN(snA,i1,j1,k1);
	SiteNode tempsn;
	int NodeID;
	int i, j, k;

	ClusterLL ClLL = (ClusterLL) getClusterList(sn);

	if (ClLL==NULL){
		return -1;
	}

	Node Nod =  getStartNode(ClLL);
	while (Nod!=NULL){
		NodeID=getNode_id(Nod);
		//printf("Checking Node %d\n",NodeID);
		getLoc(&i, &j, &k, NodeID, snA);
		tempsn=getSN(snA, i, j, k);
		if(getDwelStat(tempsn)==-1){
			//If any of the sites does has a value of -1
			//it means there is a site within the cluster
			//that is unoccupied
			printf("Node not occupied %d\n",NodeID);
			return 0;
		}
		Nod=getNextNode(Nod);
	}
	//All sites are occupied
	return 1;
}

/*
int getNeighClusterPvalHigh(SNarray snA,int i1,int j1,int k1, double * pvalHigh,long double * sum) {

	if(snA==NULL){
		return -1;
	}
	//Need to cycle through the neighbor cluster link list
	//And determine if any of the sites are unoccupied
	SiteNode sn = getSN(snA,i1,j1,k1);
	SiteNode tempsn;

	if(sn->type!=1 || sn->NeighList==NULL){
		return -1;
	}

	int NeighID;
	int init=0;
	*sum=0.0;
	int i, j, k;

	ClusterLL NeighClLL = (ClusterLL) getNeighList(sn);
	//printf("Cluster id %d\n",getCluster_id(NeighClLL));
	NeighNode NeighNod = (NeighNode) getStart(NeighClLL);
	if(NeighNod==NULL){
		printf("WARNING: Node is NULL!\n");
	}

	while (NeighNod!=NULL){
		NeighID=getNeighNode_id(NeighNod);
		getLoc(&i, &j, &k, NeighID, snA);
		tempsn=getSN(snA, i, j, k);
		printf("DwelStat Node in cluster %d\n",getDwelStat(tempsn));
		//-1 if site is unoccupied
		if(getDwelStat(tempsn)==-1){
			*sum+=(long double) getNeighNode_p(NeighNod);
			printf("Neigh ID %d pval %g sum %Lg \n",NeighID,getNeighNode_p(NeighNod),*sum);
			//printf("Node id %d sum %Lg\n",NeighID,*sum);
			if (init==0){
				init=1;
				*pvalHigh=getNeighNode_p(NeighNod);
				printf("Neigh ID %d pvalHigh %g\n",NeighID,*pvalHigh);
			}else{
				if(*pvalHigh<getNeighNode_p(NeighNod)){
					*pvalHigh=getNeighNode_p(NeighNod);
					printf("Neigh ID %d pvalHigh %g\n",NeighID,*pvalHigh);
				}
			}
		}
		NeighNod=getNextNeigh(NeighNod);
	}
	if(getElec(NeighClLL)!=0){
		printf("Value Electrode %g\n",getElecP(NeighClLL));
		*sum+=getElecP(NeighClLL);
		if(*pvalHigh>getElecP(NeighClLL)){
			*pvalHigh=getElecP(NeighClLL);
		}
	}

	return 0;
}

void getClusterPvalLow(SNarray snA,int i1,int j1,int k1, double * pvalLow,long double * sum2) {
	assert(snA!=NULL);
	//Need to cycle through the neighbor cluster link list
	//And determine if any of the sites are unoccupied
	SiteNode sn = getSN(snA,i1,j1,k1);
	SiteNode tempsn;
	assert(sn->type==1);

	int NodeID;
	int init=0;
	int i, j, k;
	printf("Finding pvalLow\n");
	ClusterLL ClLL = (ClusterLL) sn->DataStru;
	//Check if the cluster neighbors the electrode
	Node Nod = (Node) getStart(ClLL);
	assert(Nod!=NULL);
	while (Nod!=NULL){
		NodeID=getNode_id(Nod);
		getLoc(&i, &j, &k, NodeID, snA);
		tempsn=getSN(snA, i, j, k);

		//-1 means the site is unoccupied
		if(getDwelStat(tempsn)==-1){
			printf("pval %g\n",getNode_p(Nod));
			*sum2+=(long double) getNode_p(Nod);
			if (init==0){
				init=1;

				*pvalLow=getNode_p(Nod);
			}else{
				if(*pvalLow<getNode_p(Nod)){
					*pvalLow=getNode_p(Nod);
				}
			}

		}else{
			printf("Site occupied ID %d\n",NodeID);
		}
		Nod=getNextNode(Nod);
	}
}
*/
int HopOffCluster(SNarray snA, int ID, double position2, int * newID, double * timeOff){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    fprintf(stderr,"ERROR snA is NULL in HopOffCluster\n");
    return return_error_val();
  }
  if(ID<0){
    fprintf(stderr,"ERROR ID is less than 0 in HopOffCluster\n");
    return return_error_val();
  }
  if(position2<0){
    fprintf(stderr,"ERROR position2 is less than 0 in HopOffCluster\n");
    return return_error_val();
  }
  if(position2>1.001){
    fprintf(stderr,"ERROR position2 is greater than 1.001 in HopOffCluster\n");
		return return_error_val();
	}
  #endif

	SiteNode sn = getSNwithInd(snA,ID);
	ClusterLL ClLL = (ClusterLL) getClusterList(sn);
	int hops;
	int elem;
	NeighNode NeighNod = getStartNeigh(ClLL);
	//Grab the nodes neighboring the cluster
  printNeighNodesClusterLL(ClLL);
	while(NeighNod!=NULL){

		//Also must cycle through the different hops
		hops = getNeighNode_hoplength(NeighNod);

		for(elem=1;elem<(hops+1);elem++){
			
			printf("Value rand %g ID of node %d Pval node %g elem %d hops+1 %d\n",position2,getNeighNode_id(NeighNod),getNeighNode_p(NeighNod,elem),elem,(hops+1));
			if(position2<getNeighNode_p(NeighNod, elem)){
				*newID = getNeighNode_id(NeighNod);
				*timeOff = getNeighNode_t(NeighNod, elem);
				if(*timeOff<=0){
          printf("1 timeOff in HopOffCluster 0 or negative %g\n",*timeOff);
          exit(1);
        }else{
          printf("2 timeOff HopOffCluster %g\n",*timeOff);
        }
        return 0;
			}
		}

		NeighNod = getNextNeigh(NeighNod);
	}
	//If no neighboring site has been returned at this point
	//there is a problem
	return -1;

}

int HopWithinCluster(SNarray snA, int ID, double position2, int * newID){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in HopWithinCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
	if( position2<0 ){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR position2 is less than 0 in HopWithinCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
	if(position2>1.001){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR position2 is greater than 1.001 in HopWithinCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
  #endif
	
	SiteNode sn = getSNwithInd(snA,ID);
	ClusterLL ClLL = (ClusterLL) getClusterList(sn);

	Node Nod = getStartNode(ClLL);
	//grab the nodes that are within the cluster

	while (Nod!=NULL){

		if(position2<getNode_p(Nod)){
			//return the ID of the chosen node
			*newID = getNode_id(Nod);
			return 0;
		}
		Nod = getNextNode(Nod);
	}

	//If code reaches this point something
	//is wrong
	return -1;
}

int getClusterIDGivenSiteNodeID(SNarray snA, int SiteID){
  #ifdef _ERROR_CHECKING_ON_
  if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA in getClusterIDGivenSiteNodeID is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(SiteID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR SiteID in getClusterIDGivenSiteNodeID is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  SiteNode sn = getSNwithInd(snA, SiteID);
  return getClusterIDGivenSiteNode(sn);
}

int getClusterIDGivenSiteNode(SiteNode sn){

  #ifdef _ERROR_CHECKING_ON_
  if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn in getClusterIDGivenSiteNode is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  ClusterLL ClLL = (ClusterLL) getClusterList(sn);
  printf("Cluster ID %d\n",getCluster_id(ClLL));
  printf("Site Energy %g\n",getEnergy(sn));
  return getCluster_id(ClLL);
}

int mergeClusterLLGivenSiteNodeIDs( ClusterLL ClLLAll, SNarray snA,matrix MasterM,\
                                    int SiteID1, int SiteID2,ParameterFrame PF){

  #ifdef _ERROR_CHECKING_ON_
  if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA in mergeClusterLLGivenSiteNodeIDs is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(SiteID1<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR SiteID1 in mergeClusterLLGivenSiteNodeID is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(SiteID2<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR SiteID2 in mergeClusterLLGivenSiteNodeID is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  //If we get to this point it means that we found
  //the second node on a second ClusterLL. This
  //means that the two different ClusterLL must
  //be consolidated. To establish a protocal we will 
  //say that the higher id will be merged
  //with the lower cluster id so that the higher
  //one no longer exists. 
  //Attaching list to end of found node
  int ClusterID1;
  int ClusterID2;
  ClusterID1           = getClusterIDGivenSiteNodeID(snA,SiteID1);
  ClusterID2           = getClusterIDGivenSiteNodeID(snA,SiteID2);
  ClusterLL clLL1      = getClusterGivenClusterID(ClLLAll, ClusterID1);
  ClusterLL clLL2      = getClusterGivenClusterID(ClLLAll, ClusterID2);
  Node startNode       = getStartNode(clLL2);
  addClusterLLNode(&clLL1, startNode);
  //Increasing size of cluster to account for new nodes
  int numNodes1;
  int numNodes2;

  SiteNode tempSN = NULL;

  numNodes1 = getCluster_numNodes(clLL1);
  numNodes2 = getCluster_numNodes(clLL2);
  printf("MergeStep1\n");
  setCluster_NumNodes(clLL1,numNodes1+numNodes2);
  //We want to fill the hole in ClLLAll because ClLL2 is 
  //being removed from the link list
  printf("Cluster id1 %d Cluster id2 %d\n",ClusterID1,ClusterID2);
  removeClusterLLfromClusterLL(&ClLLAll, clLL2);
  printf("MergeStep2\n");

  deleteClusterLL_NeighNodes(&clLL1);
  deleteClusterLL_NeighNodes(&clLL2);
  
  printf("MergeStep3\n");
  //Remove the dependence of the SiteNode on the obsolete 
  //cluster and replace it with the correct one
  if(ClusterID1<ClusterID2){
  printf("MergeStep4\n");
    tempSN = getSNwithInd(snA,SiteID2);
    setDataStruc(&tempSN,1,(void *)clLL1);
    //Determine orientation of nodes with respect to each other
  printf("MergeStep5\n");
    DetermineNodeOrientationSingleCluster(&clLL1 , snA );
    //Calculate the Neighboring Nodes for the Cluster
  printf("MergeStep6\n");
    CalculateNeighNodesForSingleClusterLL(clLL1,snA,PFget_Px(PF),PFget_Py(PF),PFget_Pz(PF));
    CalculateSumAndPGivenSingleClusterLL(snA, clLL1, MasterM, PFget_Attempts(PF),\
        PFget_Px(PF),PFget_Py(PF),PFget_Pz(PF));
  printf("MergeStep7\n");

  }else{
  printf("MergeStep8\n");
    tempSN = getSNwithInd(snA,SiteID1);
    setDataStruc(&tempSN,1,(void *)clLL2);
    //Determine orientation of nodes with respect to each other
  printf("MergeStep9\n");
    DetermineNodeOrientationSingleCluster(&clLL2 , snA );
    //Calculate the Neighboring Nodes for the Cluster
  printf("MergeStep10\n");
    CalculateNeighNodesForSingleClusterLL(clLL2,snA,PFget_Px(PF),PFget_Py(PF),PFget_Pz(PF));
  printf("MergeStep11\n");
    CalculateSumAndPGivenSingleClusterLL(snA, clLL2, MasterM, PFget_Attempts(PF),\
        PFget_Px(PF),PFget_Py(PF),PFget_Pz(PF));
  }

  return 0;
}

int DetermineNodeOrientationSingleCluster(ClusterLL * TempClLL ,\
                                          const_SNarray snA){

  #ifdef _ERROR_CHECKING_ON_
  printf("HERE\n");
  if (TempClLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR TempClLL is NULL in DetermineNodeOrientationSingleCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  printf("HERE2\n");
  if (*TempClLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *TempClLL is NULL in DetermineNodeOrientationSingleCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if ( snA == NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in DetermineNodeOrientationSingleCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  //This is done by first checking all the hop rates off the cluster
  //Take all hops that are within an order of magnitude

  int Node_ID;
  int Node_IDFro, Node_IDBeh, Node_IDLef, Node_IDRig, Node_IDAbo, Node_IDBel;
  int i, j, k;
  Node tempNode, tempNode2;

  tempNode = getStartNode(*TempClLL);

  printf("DetermineOrientation Clustersitenode\n");

  while(tempNode!=NULL){
    //Cycle Nodes
    Node_ID=getNode_id(tempNode);
    printf("NodeID %d\n",Node_ID);
    getLoc( &i, &j, &k, Node_ID, snA);

    //Must be careful when considering periodicity in 
    //the x direction sites may act like a cluster if 
    //periodic and they are clumped at the edges of 
    //the simulation box. As soon as the electrode is
    //attached the cluster no longer acts like a
    //cluster
    Node_IDFro=getIndFroP(snA, i, j, k);
    Node_IDBeh=getIndBehP(snA, i, j, k);
    Node_IDLef=getIndLefP(snA, i, j, k);
    Node_IDRig=getIndRigP(snA, i, j, k);
    Node_IDAbo=getIndAboP(snA, i, j, k);
    Node_IDBel=getIndBelP(snA, i, j, k);
    tempNode2=getNextNode(tempNode);
    //Cycle through Cluster and find which neighboring sites are within the cluster
    while(tempNode2!=NULL){

      //If the charge is hopping to the end electrode we do not want
      //to exclude it. Clusters do not exist if it is easy for the 
      //charge to hop to an electrode at either end, but if it is
      //periodic in the X direction the cluster does need to be included

      //In the code below determing which nodes with respect to x, y and z
      //or tempNode are already in the cluster. 
      //printf("Fro %d Beh %d NodeID %d\n",Node_IDFro,Node_IDBeh,getNode_id(tempNode2));
      if(Node_IDFro==getNode_id(tempNode2) ){
        setFlagFro(tempNode);
        setFlagBeh(tempNode2);
      } else if (Node_IDBeh==getNode_id(tempNode2) ){
        setFlagBeh(tempNode);
        setFlagFro(tempNode2);
      } else if (Node_IDLef==getNode_id(tempNode2) ){
        setFlagLef(tempNode);
        setFlagRig(tempNode2);
      } else if (Node_IDRig==getNode_id(tempNode2) ){
        setFlagRig(tempNode);
        setFlagLef(tempNode2);
      } else if (Node_IDAbo==getNode_id(tempNode2) ){
        setFlagAbo(tempNode);
        setFlagBel(tempNode2);
      } else if (Node_IDBel==getNode_id(tempNode2) ){
        setFlagBel(tempNode);
        setFlagAbo(tempNode2);
      }

      tempNode2=getNextNode(tempNode2);
    }//Compare nodes within cluster
    tempNode = getNextNode(tempNode);
  }
  return 0;
}

int checkSNconnectedSameCluster(const_SiteNode sn1,
                                const_SiteNode sn2){

  #ifdef _ERROR_CHECKING_ON_
  if(sn1==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn1 is NULL in checkSNconnectedSameCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(sn2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn2 is NULL in checkSNconnectedSameCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(getClusterList(sn2)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn2->Cluster is NULL in checkSNconnectedSameCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }if(getClusterList(sn1)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn1->Cluster is NULL in checkSNconnectedSameCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int id1;
  int id2;
  id1 = getCluster_id((ClusterLL)getClusterList(sn1));
  id2 = getCluster_id((ClusterLL)getClusterList(sn2));

  if(id1==id2){
    return 1;
  }
  return 0;
}

int HopOnOffCluster(SNarray snA, int ID, double position){
	
	if(snA==NULL || position<0 || position>1.001){
		return -1;
	}
	
	SiteNode sn = getSNwithInd(snA,ID);
	if(sn==NULL){
		printf("ERROR sn found to be NULL at index %d\n",ID);
		exit(1);
	}	
	ClusterLL ClLL = (ClusterLL) getClusterList(sn);
	
	if(ClLL==NULL){
		printf("ERROR Cluster is NULL\n");
		exit(1);
	}

	double tcluster = 0.0;
	//double tprob = 0.0;
	//double SsumTotal = 0.0;
	double ProbStay = 0.0;
	
	tcluster = getCluster_time(ClLL);
	//tprob = tcluster*20;
	//SsumTotal = (1/tprob+1/tcluster);
	//ProbStay = (1/tcluster);
	//ProbStay = ProbStay/SsumTotal;
	
	//printf("snA total %d\t",getAtotal(snA));
	//printf("ID %d \t",ID);
	printf("tcluster %g\t",tcluster);
	//printf(" SsumTotal %g \t",SsumTotal);
	//printf(" tprop %g\t", tprob);	
	//printf("ProbStay %g\n",ProbStay);
	//printf("position %g\n",position);

	if( position < ProbStay){
		//This means the charge will stay
		//within the cluster
		return 0;
	}else{
		//This means the charge will try to jump
		//to a site off the cluster
		return 1;
	}

/*
	int NeighID;
	int i, j, k;
	double pval=0.0;
	long double sumCheck=0.0;
	SiteNode sn = getSN(snA,i,j,k);
	SiteNode tempsn;
	ClusterLL NeighClLL = (ClusterLL) getNeighList(sn);
	NeighNode NeighNod = (NeighNode) getStart(NeighClLL);
	printf("Cluster id %d\n",getCluster_id(NeighClLL));

	while (NeighNod!=NULL){
		NeighID=getNeighNode_id(NeighNod);
		getLoc(&i, &j, &k, NeighID, snA);
		tempsn=getSN(snA, i, j, k);
		if(getDwelStat(tempsn)==-1){
			pval+=getNeighNode_p(NeighNod)/((double) sum);
			sumCheck+=(long double) getNeighNode_p(NeighNod);
			printf("Neigh ID %d sumCheck %Lg sum %Lg\n",NeighID,sumCheck, sum);
			//printf("Value of sum %Lg sumCheck %g\n",sum,sumCheck);
			printf("pval %g\n",pval);
			assert(sumCheck<=sum);
			assert(pval<=1.01);
			if(position<=pval){
				*NewID=NeighID;
				return 1;
			}
		}
		NeighNod  = getNextNeigh(NeighNod);
	}


	//Check electrodes
	if(getElec(NeighClLL)!=0){
		pval+=getElecP(NeighClLL)/sum;
		assert(pval<1.01);
		printf("pval Elec%g\n",pval);
		if (getElec(NeighClLL)==1){
			if(position<=pval){
				//Hopped to Right electrode
				return 2;
			}
		} else{
			if(position<=pval){
				//Hopped to Left Electrode
				return 3;
			}
		}
	}

	assert(pval<1.01);
	//Check to see if stays on Cluster
	pval+=getClusterP(NeighClLL)/sum;
	assert(pval<1.01);
	printf("pval Stay%g\n",pval);
	if(position<=pval){
		return 4;
	}
	return -1;
	*/
}

/*
double getp(SiteNode sn, int num, int choice) {

	if(sn==NULL || choice<0.0 || choice>2){
		return -1;
	}

	if(choice==0){

		if(sn->type!=0 || num<0 || num>5){
			return -1;
		}

		Point pt = (Point) sn->DataStru;
		return pt->p[num];
	}else{

		if(sn->type!=1){
			return -1;
		}

		ClusterLL ClLL = (ClusterLL) sn->DataStru;
		Node tempNode = (Node) getStart(ClLL);
		int Node_ID;

		if(tempNode==NULL){
			return -1;
		}

		Node_ID = getNode_id(tempNode);
		while( Node_ID!=num){

			if(getNextNode(tempNode)==NULL){
				return -1;
			}

			tempNode=getNextNode(tempNode);
			Node_ID=getNode_id(tempNode);
		}
		return getNode_p(tempNode);
	}

}

//choice = 0 set P value of Point
//choice = 1 set P value of Node in cluster (num = SiteID)
//choice = 2 set P value of Neighbor of Cluster (num = SiteID)
//choice = 3 set P value of Electrode (num = 0)
//choice = 4 set P value of Cluster as whole (num = 0)
int setp(SiteNode sn,int num, double var, int choice) {

	if(sn==NULL || choice<0.0 || choice>4){
		return -1;
	}

	if((choice>=1 && choice<=4) && sn->type!=1){
		return -2;
	}

	if (choice==0) {

		if(sn->type!=0){
			return -2
		}
		if(num<0.0 && num>5.0){
			return -3;
		}

		Point pt = (Point) sn->DataStru;
		pt->p[num]=var;

	}else if(choice==1){

		if(num<0){
			return -3;
		}
		ClusterLL ClLL = (ClusterLL) sn->DataStru;
		Node tempNode = (Node) getStart(ClLL);
		int Node_ID;

		if(tempNode==NULL){
			return -4
		}

		Node_ID = getNode_id(tempNode);
		while( Node_ID!=num){
			if(tempNode==NULL){
				return -4
			}
			tempNode=getNextNode(tempNode);
			Node_ID=getNode_id(tempNode);
		}
		setNode_p(tempNode,var);

	}else if(choice==2){
		assert(sn->NeighList!=NULL);
		ClusterLL NeighClLL = (ClusterLL) sn->DataStru;
		NeighNode tempNeighNode = (NeighNode) getStart(NeighClLL);
		int NeighNode_ID;

		if(tempNeighNode==NULL){
			return -4;
		}

		assert(tempNeighNode!=NULL);
		NeighNode_ID= getNeighNode_id(tempNeighNode);
		while(NeighNode_ID!=num){
			if(tempNeighNode==NULL){
				return -4;
			}
			tempNeighNode=getNextNeigh(tempNeighNode);
			NeighNode_ID=getNeighNode_id(tempNeighNode);
		}
		setNeighNode_p(tempNeighNode,var);
	} else if(choice==3) {
		if (num!=0){
			return -3;
		}
		ClusterLL ClLL = (ClusterLL) sn->NeighList;
		assert(getElec(ClLL)!=0);
		setElecP(ClLL, var);
	} else {
		if (num!=0){
			return -3;
		}
		assert(sn->NeighList!=NULL);
		ClusterLL ClLL = (ClusterLL) sn->NeighList;
		setClusterP(ClLL,var);
	}

	return 0;
}
*/
