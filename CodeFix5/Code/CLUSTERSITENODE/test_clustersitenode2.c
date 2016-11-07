#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "clustersitenode.h"
#include "../SITENODE/sitenode.h"
#include "../CLUSTER/cluster.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"
#include "../CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../ERROR/error.h"

int main(void){


	printf("Beginning Testing\n");

	int rv;
	int OrderL;
	int OrderH;
	int MidPtsTotal;
	int PeriodicX = 1;
	int PeriodicY = 1;
	int PeriodicZ = 1;
	int TotalOrders;
	int attempts = 15;
	int XElecOn = 1;
	int YElecOn = 1;
	int ZElecOn = 1;
	int SLength = 6;
  int SWidth = 7;
  int SHeight = 8;
  double AttemptToHop = 10E13;
	double gamma = 2E9;
	double KT = 1;
	double reOrg = 1;
	double SiteDist = 1E-9;

  int x, y, z;

	SNarray snA = newSNarray(SLength,SWidth,SHeight);
  SiteNode tempSN;
	//Create Clusters
	SiteNode sn = getSN(snA,0,1,1);
	setEnergy(sn, -4);
	sn = getSN(snA, 1,1,1);
	setEnergy(sn, -4);
  printf("Index of two of the sites that should be tagged as a cluster\n");
  int Site1 = getIndex(snA,0,1,1);
  int Site2 = getIndex(snA,1,1,1);

  printf("Index %d\n",Site1);
  printf("Index %d\n",Site2);

	sn = getSN(snA, 4,4,3);
	setEnergy(sn,-4);
	sn = getSN(snA, 3,4,3);
	setEnergy(sn,-4);
	sn = getSN(snA, 3,3,3);
	setEnergy(sn,-4);
	sn = getSN(snA, 3,2,3);
	setEnergy(sn,-4);
	sn = getSN(snA, 4,2,3);
	setEnergy(sn,-4);
 
  int Site3 = getIndex(snA,4,4,3);
  int Site4 = getIndex(snA,3,4,3);
  int Site5 = getIndex(snA,3,3,3);
  int Site6 = getIndex(snA,3,2,3);
  int Site7 = getIndex(snA,4,2,3);
  
  matrix MasterM = CalculateAllHops(snA, 0.001, 0, 0, \
           KT, reOrg, SiteDist  , AttemptToHop, gamma,\
           PeriodicX, PeriodicY , PeriodicZ );
  
  //Creating a cluster
  ClusterLL ClLLAll = newClusterLL(1);
  rv = addNodesToClusterGivenSites(ClLLAll, Site1, Site2);
  printf("Printing location of nodes added to cluster\n");
  getLoc(&x,&y,&z,Site1,snA);
  printf("Node %d x: %d y: %d z: %d\n",Site1,x,y,z);
  getLoc(&x,&y,&z,Site2,snA);
  printf("Node %d x: %d y: %d z: %d\n",Site2,x,y,z);

  printf("Testing: DetermineNodeOrientationSingleCluster\n");
  ClusterLL clLL = NULL;
  rv = DetermineNodeOrientationSingleCluster(&clLL,snA);
  assert(rv==-1);
  rv = DetermineNodeOrientationSingleCluster(NULL,snA);
  assert(rv==-1);
  rv = DetermineNodeOrientationSingleCluster(&ClLLAll,NULL);
  assert(rv==-1);
  rv = DetermineNodeOrientationSingleCluster(&ClLLAll,snA);
  assert(rv==0);
  printClusterLL(ClLLAll);

  printf("Testing: CalculateNeighNodesForSingleClusterLL\n");
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,PeriodicX,PeriodicY, -1);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,PeriodicX,PeriodicY, 2);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,PeriodicX,2, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,PeriodicX,-1, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,2,PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,-1,PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,NULL,PeriodicX,PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(NULL,snA,PeriodicX,PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,snA,PeriodicX,PeriodicY, PeriodicZ);
  assert(rv==0);
  printClusterLL(ClLLAll);
  
  printf("Testing: CalculateSumAndPGivenSingleClusterLL\n");
  rv = CalculateSumAndPGivenSingleClusterLL(NULL, ClLLAll,\
         MasterM,attempts, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, NULL,\
     MasterM,attempts, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, NULL\
    ,attempts, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,0, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, -1, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, 2, PeriodicY, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, PeriodicX, -1, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, PeriodicX, 2, PeriodicZ);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, PeriodicX, PeriodicY, -1);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, PeriodicX, PeriodicY, 2);
  assert(rv==-1);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLLAll, MasterM,attempts, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==0);
   
  printClusterLL(ClLLAll);
  /*
  //Do book keeping connect sites (0,1,1) and (1,1,1) with cluster
  tempSN =  getSNwithInd(snA,Site1);
  setDataStruc(&tempSN,1,(void *)ClLLAll);
  tempSN =  getSNwithInd(snA,Site2);
  setDataStruc(&tempSN,1,(void *)ClLLAll);
*/
  //Now we are going to go ahead and add the second cluster and then
  //see what happens if we append a site to the cluster
/*  ClusterLL ClLL2 = newClusterLL(2);
  rv = getCluster_id(ClLL2);
  assert(rv==2);
  rv = addNodesToClusterGivenSites(ClLL2, Site3, Site4);
  assert(rv==0);
  rv = DetermineNodeOrientationSingleCluster(&ClLL2,snA);
  assert(rv==0);
  rv = CalculateNeighNodesForSingleClusterLL(ClLL2,snA,PeriodicX,PeriodicY, PeriodicZ);
  assert(rv==0);
  rv = CalculateSumAndPGivenSingleClusterLL(snA, ClLL2, MasterM,attempts, PeriodicX, PeriodicY, PeriodicZ);
  assert(rv==0);

  rv = appendClusterLL(ClLLAll,ClLL2);
  assert(rv==0);*/
  /* Get the sitenodes to point to the correct cluster */
/*  tempSN = getSNwithInd(snA,Site3);
  setDataStruc(&tempSN,1,(void *)ClLL2);
  tempSN = getSNwithInd(snA,Site4);
  setDataStruc(&tempSN,1,(void *)ClLL2);

  printf("Cluster id %d\n",getCluster_id(ClLL2));
  rv = getClusterIDGivenSiteNode(getSNwithInd(snA,Site3));
  assert(rv==2);
  rv = getClusterIDGivenSiteNode(getSNwithInd(snA,Site4));
  assert(rv==2);

  //Clear the Cluster from ClLL2 pointer
  ClLL2 = NULL;
  //Grab it now assuming we just have the site
  int ClusterID2;
  ClusterID2 = getClusterIDGivenSiteNodeID(snA,Site4);
  
  printf("********************************\n\n");
  printClusterLL(ClLLAll);
  //assert(ClusterID2 == 2);
  //ClLL2 = getClusterGivenClusterID(ClLLAll,ClusterID2);
  //assert(ClLL2!=NULL);

  //printf("Testing: mergeClusterLLGivenSiteNodeIDs\n");
  
  //printf("Creating ParameterFrame\n");
  ParameterFrame PF = newParamFrame();
  
  //Just need these functions for when calling merge
  PFset_Px(PF,PeriodicX);
  PFset_Py(PF,PeriodicY);
  PFset_Pz(PF,PeriodicZ);

  //rv = mergeClusterLLGivenSiteNodeIDs( ClLLAll, 
  //                                     snA,
  //                                     MasterM, 
  //                                     PF);
  
  deleteParamFrame(&PF);*/
  deleteAllClusterLL(&ClLLAll);
  deleteMatrix(&MasterM);	
 	deleteSNarray(&snA);
	return 0;
}
