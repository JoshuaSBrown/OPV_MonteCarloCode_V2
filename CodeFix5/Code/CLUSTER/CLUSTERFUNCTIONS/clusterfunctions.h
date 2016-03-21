#ifndef _CLUSTERFUNCTIONS_H_
#define _CLUSTERFUNCTIONS_H_

#include "SITENODE/sitenode.h"
#include "DATASTRUCT/cluster.h"
#include "MATRIX/matrix.h"
#include "../../PARAMETERS/read.h"

/* Calculates and returns the hops of every jump in the
	 form of a matrix
*/
matrix CalculateAllHops(const_SNarray snA, const double electricEnergyX, \
												const double electricEnergyY, const double electricEnergyZ,\
												const double KT, const double reOrgEnergy, const double SiteDistance,\
												const double AttemptToHop, const double gamma,\
												const int PeriodicX, const int PeriodicY, const int PeriodicZ);

/* Sorts all the hops into Midpoints and splits the midpoints
	 into an array.
*/
ArbArray MPsort(int * orderL, int * orderH, int * MidPtsTotal, const_matrix MasterM, const_SNarray snA, \
								const int PeriodicX, const int PeriodicY, const int PeriodicZ);

/* Sorts the midpoints into  a linked list with an Array. Each element
	 of the array starts a link list. Every node in the link list is
	 of the same order of magnitude.
*/
ArbArray SortOrderMag(const int TotalOrders,const int orderLow, const_ArbArray mpA);

/* Sorts the nodes into clusters based on their order of magnitude
	 Each element of the arb array points to a clusterll that is of 
	 a specefic order of magnitude. Any other ClusterLL connected
	 is of the same order of magnitude. 
	 Each ClusterLL is a group of nodes with the same order of magnitude
	 connecting them.
	 The first cluster of every ArbArray is given an id of 1 when clusters are
	 joined they take the lowest id
*/
ArbArray ClusterSort(const int TotalOrders,const int orderLow, const_ArbArray ArLL);

/* This function removes clusters that a charge could easily hop 
	 out of and leaves behind the clusters that are hard to escape
	 Basically ClArLL array is filtered and then returned as a 
	 separate ArbArray
*/
int FilterCluster(const int TotalOrders,const int orderLow,const_matrix MasterM,\
									ArbArray * ClArLL, const_SNarray snA, const int PeriodicX,\
									const int PeriodicY, const int PeriodicZ, const int XElecOn,\
									const int YElecOn, const int ZElecOn);

/* Function checks to see if a given cluster forms a percolation
	 pathway across the sample if this is indeed the case
	 the cluster is not really a cluster and needs to be deleted
*/
int FeelPercolation(ClusterLL ClLL, const_SNarray snA, int PeriodicY, int PeriodicZ);

/* Prints out all the cluster and the order of magnitude
	 associated with each of the nodes
	 if completes succesfully returns 0 else -1
*/
int PrintCheck(int TotalOrders, int orderLow, const_ArbArray ClArLL, const_SNarray snA, const_matrix MasterM);

/* Creates a separate array that contains a list of all the
	 neigbors to the clusters previously found.
*/
int CalculateNeighNodes(int TotalOrders, int orderLow, ArbArray * ClArLL, SNarray snA, int PeriodicX, int PeriodicY, int PeriodicZ);

/* Calculates the sum and the p values of the sites in the
	 clusters.
*/
int CalculateSumAndP(const int TotalOrders, const_SNarray snA, ArbArray * ClArLL, const_matrix MasterM,const int attempts,const int PeriodicX,const int PeriodicY,const int PeriodicZ);

/* Function calculates the number of hopping options a given site
	 has. The columns are:
	 col 1 - hops within the cluster
	 col 2 - hops off the cluster
	 col 3 - ID of the site
*/
int CountOptions(Node tempNode, matrix * mtxHopOpt, const_SNarray snA);

/* Calculates the using numerical technique the probability a charge
	 will hop to a given site within the cluster irrespective of the
	 time it will take.
*/
matrix CalculateProb(const_ClusterLL TempClLL, matrix mtxHopOpt, const_SNarray snA, const int attempts);

/* Calculates the probability a charge will hop to a neighbor of the
	 cluster with respect to time the matrix is normalized
*/
matrix CalculateProbNeighDwell(const int countNeighOpts, matrix mtxDwellTime, \
		                            const_matrix mtxProb, \
																Node tempNode, const_matrix MasterM, \
																const_SNarray snA, const double rateN,\
																int PeriodicX, int PeriodicY, int PeriodicZ);

/* Calculates the time a charge will spend on a given site also calculates
	 the time a charge will take to exit the cluster
*/
matrix CalculateDwellTimeAndRateN(ClusterLL * TempClLL, matrix * mtxTimes,\
																	const_matrix mtxProbNeigh, const_matrix MasterM,\
																	const_SNarray snA, double * rateN,\
																	const int PeriodicX, const int PeriodicY, const int PeriodicZ);

/* Finally we calculate the pval for the Nodes in the cluster
	 and the Nodes outside of the cluster.
*/
int CalculatePvalNeigh(ClusterLL * TempClLL, const_matrix mtxTimes,\
												const_matrix mtxProbNeighDwell);

/* And we calculate the pval for the Nodes within the cluster
*/
int CalculatePvalNodes(ClusterLL * TempClLL, matrix mtxProb, matrix mtxDwellTime);

/* This function works to connect the appropriaty link list to the SiteNodes
	 that are within a cluster
*/
int ConnectClusterSN(int TotalOrders, SNarray snA, ArbArray ClArLL);

/* Connects the appropriate clusters with the correct electrodes
	 is this really necessary if we know it is next to the electrode
	 do we really need to connect to it?
*/
int ConnectClusterElec( ArbArray * ClArLL,\
												Electrode elXB, Electrode elXF,\
												Electrode elYL, Electrode elYR,\
												Electrode elZB, Electrode elZA);

/* Find Cluster function is designed to:
	 1. Calculate the hop rates between all sites
	 2. If the hop rates are of the same order of magnitude
	    a mid point is created and classfied based on the
			order of magnitude.
	 3. Mid points of the same order of magnitude are then 
	    stored in a link list. 
	 4. Mid points within each link list are compared to 
	    determine if they share sites. 
	 5. Mid points of the same order of magnitude that 
	  	share sites are considered part of a cluster. The
			sites that are part of a cluster are stored in a 
			separate link list. 
*/
int FindCluster( int * OrderL, SNarray snA, double electricEnergyX, \
								 double electricEnergyY, double ElectricEnergyZ, \
								 ArbArray * ClArLL1,\
								 double kT, ParameterFrame PF);


/* This function will save all the cluster information so that it does
	 not have to be recalculated
*/
int SaveCluster(char * FileName, int OrderL, SNarray  snA, double electricEnergyX,\
								double electricEnergyY, double electricEnergyZ, ArbArray ClArLL1,\
								double kT, ParameterFrame PF, Electrode elXb, Electrode elXf,\
								Electrode elYl, Electrode elYr, Electrode elZb, Electrode elZa);

/* This function will load all the information needed to simulate clusters if it
	 was previously calculated
*/
int LoadCluster(char * FileName, int * OrderL, SNarray * snA, double * electricEnergyX,\
								double * electricEnergyY, double * electricEnergyZ, ArbArray * ClArLL1,\
								double * kT, ParameterFrame * PF);


int LoadCluster_Data( char * FileName, int * OrderL, SNarray * snA, double electricEnergyX,\
											double electricEnergyY, double electricEnergyZ,\
											ArbArray *ClArLL1, double kT);

int LoadCluster_Only(char * FileName, int * OrderL, SNarray * snA, double electricEnergyX,\
										 double electricEnergyY, double electricEnergyZ,\
										 ArbArray *ClArLL, double kT);

#endif
