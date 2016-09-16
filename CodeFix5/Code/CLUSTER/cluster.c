#include <stdio.h>
#include <stdlib.h>

#include "cluster.h"
#include "../ERROR/error.h"

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

////////////////////////////////////////////////////////////////////////
//Tools for accessing Clusters and Cluster Link List
ClusterLL newClusterLL(int ID){

	if(ID<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR id IN newClusterLL is less than 0\n");
    #endif
		return NULL;
	}

	ClusterLL clLL = (ClusterLL) malloc(sizeof(struct _ClusterLL));

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newClusterLL\n");
    #endif
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

int appendClusterLL(ClusterLL clLL, ClusterLL clLLnew){

  #ifdef _ERROR_CHECKING_ON_
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in appendClusterLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(clLLnew==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLLnew is NULL in appendClusterLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  ClusterLL tempClLL;
  tempClLL = clLL;
  while(tempClLL->next!=NULL){
    tempClLL = getNextClusterLL(tempClLL);
  }

  tempClLL->next = clLLnew;
  return 0;

}

int deleteClusterLL(ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in deleteCluster is NULL\n");
    #endif
		return -1;
	}
	free(clLL);
	return 0;
}

int deleteClusterLL_NeighNodes(ClusterLL * clLL){
  
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in deleteCluterLL\n");
    #endif
    return -1;
  }

  NeighLL neighLL = (*clLL)->Neigh;
  if (neighLL !=NULL){
    printf("Number of Neighbors %d",getNeighLL_numNeigh(neighLL));
    deleteNeighLLAll(neighLL); 
    //Delete neighLL
    (*clLL)->Neigh=NULL;
  }
  return 0;
}

int removeClusterLLfromClusterLL(ClusterLL * ClLLAll, ClusterLL clLL){

  if(ClLLAll==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ClLLAll in removeClusterLLfromClusterLL is NULL\n");
    #endif
    return -1;
  }
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in removeClusterLLfromClusterLL is NULL\n");
    #endif
    return -1;
  }
  
  if((*ClLLAll)->id==clLL->id){
    (*ClLLAll)=clLL->next;
    deleteClusterLL(clLL);
    return 0;    
  }

  ClusterLL ClLLstart = *ClLLAll;
  while(ClLLstart->next!=NULL){
    
    if(ClLLstart->next->id==clLL->id){
      ClLLstart->next=clLL->next;
      deleteClusterLL(clLL);
      return 0;
    }
    ClLLstart=ClLLstart->next;
  }

  return -1;
}

int deleteClusterLLNodes(ClusterLL clLL){

	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in deleteClusterLLNodes is NULL\n");
    #endif
		return -1;
	}

	printf("Cluster id %d Number of nodes %d\n",clLL->id,clLL->numNodes);
	Node tempNode;

	tempNode=clLL->startNode;
	if(tempNode!=NULL){
	  deleteNodeAll(&tempNode);
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in deleteAllClusterLL is NULL\n");
    #endif
		return -1;
	}
	if((*clLL)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL in deleteAllClusterLL is NULL\n");
    #endif
		return -1;
	}

	//Delete All the Nodes in the cluster Link List 
	ClusterLL tempclLL;
	tempclLL=(*clLL)->next;
	Node tempNode;
	NeighLL neighLL;

	while(tempclLL!=NULL) {

		printf("Cluster id %d number of Nodes %d\n",(*clLL)->id,(*clLL)->numNodes);
		//Deleting Nodes
		if ((*clLL)->startNode!=NULL) {
			tempNode=(*clLL)->startNode;
			deleteNodeAll(&tempNode);
      (*clLL)->startNode=NULL;
		}

		//Deleting Neighboring Nodes
		neighLL = (*clLL)->Neigh;
		if (neighLL !=NULL){
			printf("Number of Neighbors %d",getNeighLL_numNeigh(neighLL));
			//Delete neighLL
			deleteNeighLLAll(neighLL);
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
			deleteNodeAll(&tempNode);
			(*clLL)->startNode=NULL;
		}

		//Deleting Neighboring Nodes
		neighLL = (*clLL)->Neigh;
		if (neighLL !=NULL){
			//Delete neighLL
			deleteNeighLLAll(neighLL);
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getNextClusterLL\n");
    #endif
		return NULL;
	}
	return clLL->next;
}

int printNodesClusterLL(const_ClusterLL clLL){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in printNodesClusterLL is NULL\n");
    #endif
		return -1;
	}
	Node tempNode;

	printf("\nNodes in Cluster %d\n",clLL->id);
	tempNode=clLL->startNode;
	while(tempNode!=NULL){
		printf("%d \t",getNode_id(tempNode));
		tempNode=getNextNode(tempNode);
	}
	printf("\n");
	return 0;
}

int printClusterLL(const_ClusterLL clLL){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in printClusterLL\n");
    #endif
		return -1;
	}

	Node tempNode;
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
				printNode(tempNode);
				tempNode=getNextNode(tempNode);
			}
		}

		if(tempclLL->Neigh==NULL){
			printf("No Neighboring node list\n");
		}else{
			neighLL = tempclLL->Neigh;
      printNeighLL(neighLL);
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

ClusterLL getClusterGivenClusterID(ClusterLL clLL, int ClusterID){

  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in getClusterGivenClusterID is NULL\n");
    #endif
    return NULL;
  }
  if(ClusterID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ClusterID is less than 0 in getClusterGivenClusterID\n");
    #endif
    return NULL;
  }

  ClusterLL Clustertemp = clLL;
  while(Clustertemp!=NULL){
    if(Clustertemp->id==ClusterID){
      return Clustertemp;
    }
    Clustertemp=Clustertemp->next;
  }
  return NULL;
}

int setCluster_NeighLL(ClusterLL clLL, NeighLL NeiLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_NeighLL\n");
    #endif
		return -1;
	}
	if(NeiLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeiLL is NULL in setCluster_NeighLL\n");
    #endif
		return -1;
	}
	clLL->Neigh = NeiLL;
	return 0;
}

int setCluster_NumNodes(ClusterLL clLL, int NumNodes){

  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_NumNodes\n");
    #endif
    return -1;
  }
  if(NumNodes<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NumNodes in setCluster_NumNodes is less than 0\n");
    #endif
    return -1;
  }

  clLL->numNodes=NumNodes;
  return 0;
}

NeighLL getCluster_NeiLL(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_NeiLL\n");
    #endif
		return NULL;
	}

	return clLL->Neigh;

}

int getCluster_id(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_id\n");
    #endif
		return -1;
	}
	return clLL->id;
}

int getCluster_numNodes(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_numNodes\n");
    #endif
		return -1;
	}
	return clLL->numNodes;
}

Node getCluster_LastNode(const_ClusterLL clLL){

  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in getCluster_LastNode is NULL\n");
    #endif
    return NULL;
  }

  Node LastNode = getNode_lastNode(clLL->startNode);

  return LastNode;
}

int setCluster_elecXFid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXFid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecXFid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecXFid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecXFid\n");
    #endif
		return -1;
	}

	(*clLL)->elecXF_id = Elec;
	return 0;

}

int setCluster_elecXBid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXBid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecXBid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecXBid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecXBid\n");
    #endif
		return -1;
	}

	(*clLL)->elecXB_id = Elec;	
	
	return 0;

}

int setCluster_elecYLid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYLid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecYLid\n");
    #endif
		return -1;
	}
	if( Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecYLid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecYLid\n");
    #endif
		return -1;
	}

	(*clLL)->elecYL_id = Elec;
	return 0;

}

int setCluster_elecYRid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYRid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecYRid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecYRid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecYRid\n");
    #endif
		return -1;
	}

	(*clLL)->elecYR_id = Elec;
	return 0;

}

int setCluster_elecZBid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZBid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecZBid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecZBid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 0 in setCluster_elecZBid\n");
    #endif
		return -1;
	}

	(*clLL)->elecZB_id = Elec;
	return 0;

}

int setCluster_elecZAid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZAid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecZAid\n");
    #endif
		return -1;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecZAid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecZAid\n");
    #endif
		return -1;
	}

	(*clLL)->elecZA_id = Elec;
	return 0;

}

int setCluster_elecXid(ClusterLL * clLL, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXid\n");
    #endif
		return -1;
	}
	if( Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecXid\n");
    #endif
		return -1;
	}
	if( Elec>1 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecXid\n");
    #endif
		return -1;
	}

	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecXid\n");
    #endif
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
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYid\n");
    #endif
		return -1;
	}
	if(Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecYid\n");
    #endif
		return -1;
	}
	if( Elec>1 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecYid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecYid\n");
    #endif
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
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZid\n");
    #endif
		return -1;
	}
	if( Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecZid\n");
    #endif
		return -1;
	}
	if( Elec>1 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecZid\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecZid\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXF\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecXF\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXB\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecXB\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYR\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecYR\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYL\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecYL\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZA\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecZA\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZB\n");
    #endif
		return -1;
	}
	if(*clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in setCluster_elecZB\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(sum<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sum is less than 0 in setCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecXsum\n");
    #endif
		return -1;
	}

	if(Elec==0){
		if(clLL->elecXB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR elecXB does not exist yet you are tyring set the sum value using setCluster_elecXsum\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecXB,sum);
		}
	}else{
		if(clLL->elecXF==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR elecXF does not exist yet you are tyring set the sum value using setCluster_elecXsum\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecXF,sum);
		}
	}
	return 0;
}

int setCluster_elecYsum(ClusterLL clLL, double sum, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(sum<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sum is less than 0 in setCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecYsum\n");
    #endif
		return -1;
	}

	if(Elec==0){
		if(clLL->elecYL==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to set the sum o elecYL using setCluster_elecYsum yet elecYL does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecYL,sum);
		}
	}else{
		if(clLL->elecYR==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to set the sum o elecYR using setCluster_elecYsum yet elecYR does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecYR,sum);
		}
	}
	return 0;
}

int setCluster_elecZsum(ClusterLL clLL, double sum, int Elec){
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(sum<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sum is less than 0 in setCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in setCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in setCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecZB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to set sum for elecZB using setCluster_elecZsum yet elecZB is not defined\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecZB, sum);
		}
	}else{
		if(clLL->elecZA==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to set sum for elecZA using setCluster_elecZsum yet elecZA is not defined\n");
      #endif
			return -1;
		}else{
			setElectrode_Sum(clLL->elecZA,sum);
		}
	}
	return 0;
}

int addToCluster_elecXsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addToCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(val<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR val is less than 0 in addToCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in addToCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in addToCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecXB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to add to elecXB sum using addToCluster_elecXsum yet elecXB does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecXB,val);
		}
	}else{
		if(clLL->elecXF==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to add to elecXF sum using addToCluster_elecXsum yet elecXF does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecXF,val);
		}
	}
	return 0;
}

int addToCluster_elecYsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addToCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(val<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR val is less than 0 in addToCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in addToCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in addToCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecYL==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to sum for elecYL using addToCluster_elecYsum yet elecYL does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecYL,val);
		}
	}else{
		if(clLL->elecYR==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to sum for elecYR using addToCluster_elecYsum yet elecYR does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecYR,val);
		}
	}
	return 0;
}

int addToCluster_elecZsum(ClusterLL clLL, double val, int Elec){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addToCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(val<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR val is less than 0 in addToCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in addToCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in addToCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecZB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to sum for elecZB using addToCluster_elecZsum yet elecZB does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecZB,val);
		}
	}else{
		if(clLL->elecZA==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to add to sum for elecZA using addToCluster_elecZsum yet elecZA does not exist\n");
      #endif
			return -1;
		}else{
			setElectrode_AddToSum(clLL->elecZA,val);
		}
	}
	return 0;
}

int checkCluster_elecXdefined(const_ClusterLL clLL){
  #ifdef _ERROR_CHECKING_ON_
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in checkCluster_elecXid\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int rv = 0;

  if(clLL->elecXF!=NULL){
    rv = 1;
  }
  if(clLL->elecXB!=NULL){
    rv+=2;
  }

  return rv;
}

int checkCluster_elecYdefined(const_ClusterLL clLL){
  #ifdef _ERROR_CHECKING_ON_
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in checkCluster_elecYid\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int rv = 0;

  if(clLL->elecYR!=NULL){
    rv = 1;
  }
  if(clLL->elecYL!=NULL){
    rv+=2;
  }

  return rv;
}

int checkCluster_elecZdefined(const_ClusterLL clLL){
  #ifdef _ERROR_CHECKING_ON_
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in checkCluster_elecZid\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int rv = 0;

  if(clLL->elecZA!=NULL){
    rv = 1;
  }
  if(clLL->elecZB!=NULL){
    rv+=2;
  }

  return rv;
}

int getCluster_elecidXB(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCLuster_elecidXB\n");
    #endif
		return -1;
	}

	return clLL->elecXB_id;
}

int getCluster_elecidXF(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecidXF\n");
    #endif
		return -1;
	}

	return clLL->elecXF_id;
}

int getCluster_elecidYL(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecidYL\n");
    #endif
		return -1;
	}

	return clLL->elecYL_id;
}

int getCluster_elecidYR(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecidYR\n");
    #endif
		return -1;
	}

	return clLL->elecYR_id;
}

int getCluster_elecidZB(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecidZB\n");
    #endif
		return -1;
	}

	return clLL->elecZB_id;
}

int getCluster_elecidZA(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getClsuter_elecidZA\n");
    #endif
		return -1;
	}

	return clLL->elecZA_id;
}

double getCluster_elecXsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecXsum\n");
    #endif
		return -1.0;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in getCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in getCluster_elecXsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecXB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to access sum of elecXB using getCluster_elecXsum yet elecXB is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecXB);
		}
	}else{
		if(clLL->elecXF==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to access sum of elecXF using getCluster_elecXsum yet elecXF is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecXF);
		}
	}
}

double getCluster_elecYsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL for getCluster_elecYsum\n");
    #endif
		return -1.0;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 1 in getCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in getCluster_elecYsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecYL==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to get the sum for elecYL in getCluster_elecYsum but elecYL is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecYL);
		}
	}else{
		if(clLL->elecYR==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to get the sum for elecYR in getCluster_elecYsum but elecYR is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecYR);
		}
	}
}

double getCluster_elecZsum(const_ClusterLL clLL, int Elec){
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL for getCluster_elecZsum\n");
    #endif
		return -1.0;
	}
	if( Elec>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is greater than 0 in getCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elec is less than 0 in getCluster_elecZsum\n");
    #endif
		return -1;
	}
	if(Elec==0){
		if(clLL->elecZB==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to get sum of elecZB in getCluster_elecZum but elecZB is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecZB);
		}
	}else{
		if(clLL->elecZA==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR trying to get sum of elecZA in getCluster_elecZum but elecZA is not defined for clLL\n");
      #endif
			return -1;
		}else{
			return getElectrode_Sum(clLL->elecZA);
		}
	}
}

int setCluster_id(ClusterLL clLL, int ID){
	if (clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL for setCluster_id\n");
    #endif
		return -1;
	}
	if (ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ID is less than 0 for setCluster_id\n");
    #endif
		return -1;
	}
	clLL->id=ID;
	return 0;
}

double getCluster_Sum(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL for getCluster_Sum\n");
    #endif
		return -1.0;
	}
	return clLL->sum;
}

int addToClusterSum(ClusterLL clLL, double s){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL for addToClusterSum\n");
    #endif
		return -1;
	}
	if( s<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR s is less than or equal to 0 in addToClusterSum\n");
    #endif
		return -1;
	}
	clLL->sum=clLL->sum+s;
	return 0;
}

int setClusterSum(ClusterLL clLL, double s){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setClusterSum\n");
    #endif
		return -1;
	}
	if(s<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR s is less than 0 in setClusterSum\n");
    #endif
		return -1;
	}
	clLL->sum=s;
	return 0;
}

int setNextClusterLL(ClusterLL clLL, ClusterLL Nex){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setNextClusterLL\n");
    #endif
		return -1;
	}

	clLL->next=Nex;
	return 0;
}

Node getStartNode(const_ClusterLL clLL){	
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getStartNode\n");
    #endif
		return NULL;
	}
	return clLL->startNode;
}

Node getCluster_Node( const_ClusterLL clLL, int id){
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_Node\n");
    #endif
		return NULL;
	}
	if(id<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR id is less than 0 in getCluster_Node\n");
    #endif
		return NULL;
	}

	Node tempNode = clLL->startNode;

	while(tempNode!=NULL){

		if (getNode_id(tempNode)==id){
			return tempNode;
		}	

		tempNode=getNextNode(tempNode);
	}

	printf("Could not find Node with given id in cluster\n");
	return NULL;
}

NeighNode getStartNeigh(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getStartNeigh\n");
    #endif
		return NULL;
	}
	if(clLL->Neigh==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Neigh link list is not defined for clLL in getStartNeigh\n");
    #endif
		return NULL;
	}
	return getNeighLL_start(clLL->Neigh);
}

int getCluster_numNeigh(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getClsuter_numNeigh\n");
    #endif
		return -1;
	}
	if(clLL->Neigh==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Neigh link list is not defined for clLL in getCluster_numNeigh\n");
    #endif
		return -1;
	}

	return getNeighLL_numNeigh(clLL->Neigh);
}

double getCluster_p(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_p\n");
    #endif
		return -1.0;
	}
	return clLL->clusterp;
}

int setCluster_p(ClusterLL clLL, double p){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_p\n");
    #endif
		return -1;
	}
	if( p<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR p is less than 0 in setCluster_p\n");
    #endif
		return -1;
	}
	clLL->clusterp=p;
	return 0;
}

int checkCluster_elec(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in checkCluster_elec\n");
    #endif
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
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecXid\n");
    #endif
		return -1;
	}
	if((clLL->elecXF==NULL && clLL->elecXB==NULL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR No electrodes are defined cannot get electrode id in getCluster_elecXid\n");
    #endif
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
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecYid\n");
    #endif
		return -1;
	}
	if((clLL->elecYL==NULL && clLL->elecYR==NULL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR No electrodes are defined cannot get electrode id in getCluster_elecYid\n");
    #endif
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
	if(clLL==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_elecZid\n");
    #endif
		return -1;
	}
	if((clLL->elecZB==NULL && clLL->elecZA==NULL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR No electrodes are defined cannot get electrode id in getCluster_elecZid\n");
    #endif
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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in setCluster_time\n");
    #endif
		return -1;
	}

	clLL->time=t;
	return 0;
}

double getCluster_time(const_ClusterLL clLL){
	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in getCluster_time\n");
    #endif
		return -1.0;
	}
	return clLL->time;
}

////////////////////////////////////////////////////////////////////////////////////
int addNodeEndClusterLL(ClusterLL clLL, int nei){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addNodeEndClusterLL\n");
    #endif
		return -1;
	}
	if(nei<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nei is less than 0 in addNodeEndClusterLL\n");
    #endif
		return -1;
	}

	Node tempNode;
	tempNode=clLL->startNode;

	if (tempNode==NULL){
		clLL->startNode=newNode(nei);
	}else{
    tempNode = getNode_lastNode(tempNode);
    setNextNode(&tempNode,newNode(nei));
	}
	clLL->numNodes++;
  return 0;
}

int addClusterLLNode(ClusterLL *clLL, Node nd){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addClusterLLNode\n");
    #endif
		return -1;
	}
	if( nd==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nd is NULL in addClusterLLNode\n");
    #endif
		return -1;
	}

	Node tempNode;
	tempNode = (*clLL)->startNode;

	if(tempNode==NULL){
		(*clLL)->startNode = nd;
	}else{
    tempNode = getNode_lastNode(tempNode);
		setNextNode(&tempNode,nd);
	}
	(*clLL)->numNodes++;

	return 0;
}

int addNodeNewCluster(ClusterLL * clLL, MidPoint mp){

	if( (*clLL)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addNodeNewCluster\n");
    #endif
		return -1;
	}
	if( mp==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in addNodeNewCluster\n");
    #endif
		return -1;
	}
	if((*clLL)->next!=NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in addNodeNewCluster\n");
    #endif
		return -1;
	}

	ClusterLL clLL2;
	int ClusterId = (*clLL)->id+1;
	(*clLL)->next= newClusterLL(ClusterId);
	clLL2=(*clLL)->next;
	clLL2->numNodes=clLL2->numNodes+2;

	clLL2->startNode = newNode(getMP_nei1(mp));
	Node tempNode = clLL2->startNode;
	setNextNode(&tempNode,newNode(getMP_nei2(mp)));
	return 0;
}

int addNodesToClusterGivenSites(ClusterLL clLL, int siteID1, int siteID2){

  #ifdef _ERROR_CHECKING_ON_
  if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL in addNodesToClusterGivenSites is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(siteID1<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR siteID1 in addNodesToClusterGivenSites is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(siteID2<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR siteID2 in addNodesToClusterGivenSites is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(siteID1==siteID2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR siteID2==siteID1 in addNodesToClusterGivenSites\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  clLL->numNodes = clLL->numNodes+2;
  clLL->startNode = newNode(siteID1);
  Node tempNode = clLL->startNode;
  setNextNode(&tempNode,newNode(siteID2));
  return 0;

}

int addNodeToCluster( ClusterLL clLL, MidPoint mp){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addNodeToCluster\n");
    #endif
		return -1;
	}
	if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in addNodeToCluster\n");
    #endif
		return -1;
	}

	//Check first cluster. If first cluster does not exist 
	//means no other clusters exist so we must create one
	if (clLL->startNode==NULL) {

		//New cluster should have an initial size of 0
		clLL->numNodes =clLL->numNodes+2;
		clLL->startNode = newNode(getMP_nei1(mp));
		Node tempNode = clLL->startNode;
		setNextNode(&tempNode,newNode(getMP_nei2(mp)));
		return 0;
	}

	ClusterLL tempclLL=clLL;
	Node tempNode= clLL->startNode;

	//int ClusterId = clLL->id;
	int unfoundnei;

	while(getNode_id(tempNode)!=getMP_nei1(mp) && getNode_id(tempNode)!=getMP_nei2(mp)){

		//Reach end of nodes move to next link list
		if(getNextNode(tempNode)==NULL){
			//If next ClusterLL is NULL then must create a new one
			if(tempclLL->next==NULL){
				addNodeNewCluster(&tempclLL, mp);
				return 0;

			} else {
				//Move to the next ClusterLL
				tempclLL = tempclLL->next;
				//ClusterId = tempclLL->id;
				tempNode =  tempclLL->startNode;
			}
		}else{
			tempNode=getNextNode(tempNode);
		}
	}

	//If the code gets to this point it means that one of the
	//nodes already exists in the cluster. we need to make 
	//sure that the other node does not exist in a separate cluster
	if (getNode_id(tempNode)==getMP_nei1(mp)){
		unfoundnei=getMP_nei2(mp);
	}else{
		unfoundnei=getMP_nei1(mp);
	}

	ClusterLL foundclLL=tempclLL;

	Node endNodeFoundLL=tempNode;
	//Search within the rest of the current cluster
	while(getNextNode(tempNode)!=NULL){
		if(getNode_id(getNextNode(tempNode))==unfoundnei){
			//This means both sites are within the same
			//cluster already so nothing is done
			return -1;
		}
		tempNode=getNextNode(tempNode);
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
	//ClusterId=tempclLL->id;
	tempNode= tempclLL->startNode;

	while(getNode_id(tempNode)!=unfoundnei){
		if(getNextNode(tempNode)==NULL){
			//If next ClusterLL is NULL then must 
			//add second node to the found cluster
			if(tempclLL->next==NULL){
				addNodeEndClusterLL(foundclLL, unfoundnei);
				return 0;

			} else {
				//Move to the next ClusterLL
				PrevClusterLL = tempclLL;
				tempclLL = tempclLL->next;
				//ClusterId = tempclLL->id;
				tempNode = tempclLL->startNode;
			}
		}else{
			tempNode=getNextNode(tempNode);
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
	setNextNode(&endNodeFoundLL,tempclLL->startNode);
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

	foundclLL->numNodes=foundclLL->numNodes+tempclLL->numNodes;
	printf("Final Cluster %d size: %d\n",foundclLL->id,foundclLL->numNodes);

	//By passing second ClusterLL
	PrevClusterLL->next=tempclLL->next;
	//Deleting the cluster LL without deleting nodes etc
	deleteClusterLL(tempclLL);

	return 0;
}

int addNeighNodeToCluster( ClusterLL* clLL, int Neigh_ID){

	if((*clLL)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *clLL is NULL in addNeighNodeToCluster\n");
    #endif
		return -1;
	}
	if(Neigh_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Neigh_ID is less than 0 in addNeighNodeToCluster\n");
    #endif
		return -1;
	}

	if ((*clLL)->Neigh==NULL){
		(*clLL)->Neigh = newNeighLL();
		setNeighLL_start((*clLL)->Neigh,newNeighNode(Neigh_ID));
		return 0;
	}

  int rv;
  NeighNode Nei = newNeighNode(Neigh_ID);
  rv = setNeighLL_addNeighNode((*clLL)->Neigh,Nei);
	
  //Means that the Node was already added to the list
  if(rv==-2){
    deleteNeighNode(Nei); 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR cannot add NeighNode to clLL in addNieghNodeToCluster because a NeighNode with the Neigh_ID is already defined in the Cluster link list\n");
    #endif
    return -1;
  }

	return 0;
}

int addPvalClusterNode(ClusterLL clLL, int Node_ID, double pval){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addPvalClusterNode\n");
    #endif
		return -1;
	}
	if(Node_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Node_ID is less than 0 in addPvalClusterNode\n");
    #endif
		return -1;
	}
	if(pval<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pval is less than 0 in addPvalClusterNode pval must be positive or 0 because it is a probability\n");
    #endif
		return -1;
	}
	//Search through Cluster till find Node_ID
	int N_ID;
	Node tempNode = getStartNode(clLL);
	if (tempNode==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR the starting node for clLL in addPvalClusterNode is not defined this means we can not find teh node with the same id as Node_ID\n");
    #endif
		return -1;
	}
	N_ID = getNode_id(tempNode);
	while( N_ID!=Node_ID){
		if(getNextNode(tempNode)==NULL){
			return -1;
		}
		tempNode=getNextNode(tempNode);
		N_ID = getNode_id(tempNode);
	}

	//Used to calculate the total time for a charge
	//to exist hop anywhere within the cluster
	clLL->time=1/pval+clLL->time;
	setNode_p(tempNode, pval+getNode_p(tempNode));
	return 0;
}

int addPvalClusterNeighNode(ClusterLL clLL, int Node_ID, double pval){

	if(clLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR clLL is NULL in addPvalClusterNeighNode\n");
    #endif
		return -1;
	}
	if( Node_ID<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Node_ID is less than 0 in addPvalClusterNeighNode\n");
    #endif
		return -1;
	}
	if(pval<0.0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pval is less than 0 in addPvalClusterNeighNode pval is a probability it must be 0 or positive\n");
    #endif
		return -1;
	}

	if(clLL->Neigh==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR there is no NeighLL attached to clLL in addPvalClusterNeighNode\n");
    #endif
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
		if (getNextNeigh(tempNeighNode)==NULL){
			return -1;
		}
		tempNeighNode = getNextNeigh(tempNeighNode);
		Neigh_id = getNeighNode_id(tempNeighNode);
	}

	setNeighNodeNew_p(tempNeighNode, pval);
	return 0;
}

