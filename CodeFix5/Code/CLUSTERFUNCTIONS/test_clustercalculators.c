#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "clustercalculators.h"
#include "../CLUSTER/cluster.h"
#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"

int main() {
	
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	int attempts = 15;
	int rv;
	
	double SiteDistance = 1;
    double AttemptToHop = 1E-13;
	double gamma = 2;

	/* To test we will begin by creating a world that is 
     * non periodic */
	int PeriodicX = 0;
    int PeriodicY = 0;
    int PeriodicZ = 0;

	SNarray snA = newSNarray(5,4,4);
	SiteNode sn;
	int SiteID;
    int NumNodes;
	Node first;
	/* Create Low energy points for cluster */
	sn = getSN(snA,1,2,1);
	setEnergy(sn,-3);
	sn = getSN(snA,2,2,1);
	setEnergy(sn,-3);
	sn = getSN(snA,3,2,1);
	setEnergy(sn,-3);
	sn = getSN(snA,1,1,1);
	setEnergy(sn,-3);
	sn = getSN(snA,2,1,1);
	setEnergy(sn,-3);

	matrix MasterM = CalculateAllHops(snA,
									  0,
									  0,
									  0,
									  1,
									  1,
									  SiteDistance,
									  AttemptToHop,
									  gamma,
									  PeriodicX,
									  PeriodicY,
									  PeriodicZ);

	/* Create a cluster and add nodes to the
     * cluster */
	ClusterLL clLL = newCluster(0);
	SiteID = getIndex(snA,1,2,1);
	addNodeEndClusterLL(ClLL,SiteID);
	SiteID = getIndex(snA,2,2,1);
	addNodeEndClusterLL(ClLL,SiteID);
	SiteID = getIndex(snA,3,2,1);
	addNodeEndClusterLL(ClLL,SiteID);
	SiteID = getIndex(snA,1,1,1);
	addNodeEndClusterLL(ClLL,SiteID);
	SiteID = getIndex(snA,2,1,1);
	addNodeEndClusterLL(ClLL,SiteID);

	/* Get the number of Nodes that are within the cluster */
	NumNodes = getCluster_numNodes(ClLL);
	assert(NumNodes==5); 
	/* Create the matrix to store hopping option information */	
	matrix mtxHopOpt = newMatrix(NumNodes,3);

	/* Get first Node in the cluster */
	first = getStartNode(ClLL);

	CountOptions(first,
                 &mtxHopOpt,
                 snA);

	printMatrix(mtxHopOpt);

	/* We need to calculate the sum values for each of the 
     * SiteNodes in the snA, the sum values contain the sum
     * of the rates off the site */

	deleteMatrix(&mtxHopOpt);
	deleteClusterLL(&ClLL);
	deleteMatrix(&MasterM);
	deleteSNarray(&snA);
	return 0;
}
