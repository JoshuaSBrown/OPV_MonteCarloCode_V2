#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "clusterfunctions.h"
#include "../CLUSTER/cluster.h"
#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"
#include "../CLUSTERSITENODE/clustersitenode.h"
#include "../CHARGE/charge.h"
#include "../CHARGESITENODE/chargesitenode.h"

int main() {

	printf("Beginning Testing\n");

	int GlobalClusterID = 0;
	int ClusterID = -1;
	int rv;
	
	int Threshold = 10;
	int PeriodicX = 1;
	int PeriodicY = 1;
	int PeriodicZ = 1;
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

	ParameterFrame PF = newParamFrame();
	PFset_Px(PF,PeriodicX);
	PFset_Py(PF,PeriodicY);
	PFset_Pz(PF,PeriodicZ);
	PFset_Attempts(PF,attempts);
	PFset_ClusterAlgTrigger(PF,Threshold);
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
	sn = getSN(snA, 4,1,3);
    setEnergy(sn,-4);
	sn = getSN(snA, 4,1,2);
    setEnergy(sn,-4);

    int Site3 = getIndex(snA,4,4,3);
    int Site4 = getIndex(snA,3,4,3);
    int Site5 = getIndex(snA,3,3,3);
    int Site6 = getIndex(snA,3,2,3);
    int Site7 = getIndex(snA,4,2,3);
    int Site8 = getIndex(snA,4,1,3);
    int Site9 = getIndex(snA,4,1,2);

	matrix MasterM = CalculateAllHops(snA, 0.001, 0, 0, \
            KT, reOrg, SiteDist  , AttemptToHop, gamma,\
            PeriodicX, PeriodicY , PeriodicZ );

	ClusterLL ClLLAll;
	ArbArray  ArClLL = newArbArray(1,1);

	ClLLAll = newClusterLL(GlobalClusterID);
	GlobalClusterID++;
	rv = addNodesToClusterGivenSites(ClLLAll,Site1,Site2);
	assert(rv==0);
	rv = DetermineNodeOrientationSingleCluster(&ClLLAll,snA);
	assert(rv==0);
	rv = CalculateNeighNodesForSingleClusterLL(ClLLAll,
                                               snA,
                                               PeriodicX,
                                               PeriodicY,
                                               PeriodicZ);
	assert(rv==0);
	rv = CalculateSumAndPGivenSingleClusterLL(snA,ClLLAll,
                                              MasterM,
                                              PFget_Attempts(PF),
                                              PeriodicX,
                                              PeriodicY,
                                              PeriodicZ);
	assert(rv==0);
	tempSN = getSNwithInd(snA,Site1);
	setDataStruc(&tempSN,1,(void *)ClLLAll);
	tempSN = getSNwithInd(snA,Site2);
	setDataStruc(&tempSN,1,(void *)ClLLAll);

	setArbElement(ArClLL,0,(void *)ClLLAll);
	/* Create a charge */
	Charge one = newCharge();
	/* Initialize charge path so it has a memory of four nodes */
	initChargePath(&one,4);	
	setCx(one,3);
	setCy(one,3);
	setCz(one,3);
	/* Updating teh charge path or memory */
	ClusterID = -1;
	for(int i=0;i<10;i++){
		updatePath(snA,one,Site5,ClusterID);
		updatePath(snA,one,Site6,ClusterID);
	}
	updatePath(snA,one,Site5,ClusterID);

	matrix FutureSites = newMatrix(1,1);
	setE(FutureSites,1,1,Site6);

	/* Checking for Clusters */
	ClusterChargePath(one,
                      0,
					  FutureSites,
                      snA,
                      MasterM,
                      PF,
                      &ArClLL,
				      &GlobalClusterID);
                                 
	printClusterLL(ClLLAll);
	printCharge(one);

	for(int i=0;i<10;i++){
		updatePath(snA,one,Site7,ClusterID);
		updatePath(snA,one,Site6,1);
	}

	updatePath(snA,one,Site6,1);
	updatePath(snA,one,Site7,ClusterID);
	
	setE(FutureSites,1,1,Site7);

	ClusterChargePath(one,
                      0,
					  FutureSites,
                      snA,
                      MasterM,
                      PF,
                      &ArClLL,
				      &GlobalClusterID);
	
	printClusterLL(ClLLAll);
	printCharge(one);

	/* Create second cluster */
	for(int i=0;i<10;i++){
		updatePath(snA,one,Site8,ClusterID);
		updatePath(snA,one,Site9,ClusterID);
	}

	updatePath(snA,one,Site8,ClusterID);
	updatePath(snA,one,Site9,ClusterID);

	printCharge(one);
	
	setE(FutureSites,1,1,Site8);

	ClusterChargePath(one,
                      0,
					  FutureSites,
                      snA,
                      MasterM,
                      PF,
                      &ArClLL,
				      &GlobalClusterID);
	
	printClusterLL(ClLLAll);
	printCharge(one);

	/* Should now join second and 1st cluster */

	for(int i=0;i<10;i++){
		updatePath(snA,one,Site8,2);
		updatePath(snA,one,Site7,1);
	}

	updatePath(snA,one,Site8,2);
	updatePath(snA,one,Site7,1);

	printCharge(one);

	setE(FutureSites,1,1,Site8);

	ClusterChargePath(one,
                      0,
					  FutureSites,
                      snA,
                      MasterM,
                      PF,
                      &ArClLL,
                      &GlobalClusterID);


	printClusterLL(ClLLAll);
	printCharge(one);

	deleteArbArray(&ArClLL);
	deleteParamFrame(&PF);
	deleteMatrix(&FutureSites);
	deleteMatrix(&MasterM);
	deleteSNarray(&snA);
	deleteCharge(one);
	return 0;
}
