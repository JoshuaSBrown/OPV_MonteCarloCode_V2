#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CLUSTERFUNCTIONS/clusterfunctions.h"
#include "clustersitenode.h"

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
	double AttemptToHop = 10E13;
	double gamma = 2E9;
	double KT = 1;
	double reOrg = 1;
	double SiteDist = 1E-9;

	SNarray snA = newSNarray(6,7,8);

	//Create Clusters
	SiteNode sn = getSN(snA,0,1,1);
	setEnergy(sn, -4);
	sn = getSN(snA, 1,1,1);
	setEnergy(sn, -4);

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

	//0 bias, 1 kT, 1 reOrgEnergy
	matrix MasterM = CalculateAllHops(snA, 0, 1, 1, KT, reOrg, SiteDist, AttemptToHop, gamma, PeriodicX, PeriodicY, PeriodicZ);

	ArbArray mpA = MPsort( &OrderL, &OrderH, &MidPtsTotal, MasterM, snA, PeriodicX, PeriodicY, PeriodicZ);

	TotalOrders = OrderH-OrderL+1;

	ArbArray ArClLL2 = SortOrderMag(TotalOrders, OrderL, mpA);

	ArbArray ArClLL = ClusterSort( TotalOrders, OrderL, ArClLL2);

	FilterCluster( TotalOrders, OrderL, MasterM, &ArClLL, snA, PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);

	CalculateNeighNodes(TotalOrders, OrderL, &ArClLL, snA, PeriodicX, PeriodicY, PeriodicZ);

	CalculateSumAndP(TotalOrders, snA, &ArClLL, MasterM, attempts, PeriodicX, PeriodicY, PeriodicZ);
	printf("\n*********************Starting PrintCheck****************\n");
	PrintCheck( TotalOrders, OrderL, ArClLL, snA, MasterM);

	printf("Elements Reserved %d\n",getElementsReserved( ArClLL));

	char c[] = "Coordinates\0";
	PrintFile_xyz(OrderL, snA, &ArClLL, &c[0]);

	char d[] = "CoordinatesNeigh\0";
	PrintNeighFile_xyz(OrderL, snA, &ArClLL, &d[0]);

	ConnectClusterSN(TotalOrders, snA, ArClLL);
	
	printf("Testing:OccAllNeiCLuster\n");
	//Check within both clusters ensure that 
	//none of the neighbor nodes are occupied
	//
	rv = OccAllNeiCluster(snA, 1, 1, 1);
	assert(rv==0);
	rv = OccAllNeiCluster(snA, 4, 4, 3);
	assert(rv==0);

	//Correctly identifies when a site is not
	//part of the cluster
	rv = OccAllNeiCluster(snA, 5,5,5);
	assert(rv==-1);


	//Grab the nodes around the first cluster 
	//and Make them occupied
	//sn = getSN(snA,0,1,1);
	//sn = getSN(snA, 1,1,1);

	sn = getSN(snA,0,0,1);
	setDwelStat(sn, 1);
	sn = getSN(snA,0,2,1);
	setDwelStat(sn, 2);
	sn = getSN(snA,0,1,0);
	setDwelStat(sn, 3);
	sn = getSN(snA,0,1,2);
	setDwelStat(sn, 4);
	sn = getSN(snA,5,1,1);
	setDwelStat(sn, 5);
	
	sn = getSN(snA, 2,1,1);
	setDwelStat(sn, 6);
	sn = getSN(snA, 1,0,1);
	setDwelStat(sn, 7);
	sn = getSN(snA, 1,2,1);
	setDwelStat(sn, 8);
	sn = getSN(snA, 1,1,0);
	setDwelStat(sn, 9);
	sn = getSN(snA, 1,1,2);
	setDwelStat(sn, 10);
	
	rv = OccAllNeiCluster(snA, 1, 1, 1);
	assert(rv==1);
	
	sn = getSN(snA, 1,1,2);
	setDwelStat(sn, -1);
	
	rv = OccAllNeiCluster(snA, 1, 1, 1);
	assert(rv==0);

	/*
	printf("Testing:OccAllCluster\n");
	rv = OccAllCluster(snA, 1,1,1);
	assert(rv==0);
	rv = OccAllCluster(snA, 0,1,1);
	assert(rv==0);
	rv = OccAllCluster(snA, 3,3,3);
	assert(rv==0);
	rv = OccAllCluster(snA, 5,5,5);
	assert(rv==-1);
	
	sn = getSN(snA, 0,1,1);
	setDwelStat(sn, 11);
	sn = getSN(snA, 1,1,1);
	setDwelStat(sn, 12);
	rv = OccAllCluster(snA, 1,1,1);
	assert(rv==1);
	*/
	deleteArbArray(&ArClLL);
	deleteArbArray(&ArClLL2);
	deleteMatrix(MasterM);
	deleteAllMidPointArray(&mpA);
	deleteSNarray(snA);
	return 0;
}
