#ifndef _CLUSTERCALCULATORS_H_
#define _CLUSTERCALCULATORS_H_

#include "../SITENODE/sitenode.h"
#include "../CLUSTER/cluster.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"
#include "../ARBARRAY/arbarray.h"
#include "../CHARGE/charge.h"

/* Given a node within the cluster this function will count the 
 * number of hops a given site has that are off a cluster as
 * well as the number of sites that are within the cluster 
 * mtxHopOpt Column:
 * 1 - Number of hops within cluster
 * 2 - Number of Hops out the cluster
 * 3 - ID of site 
 */
int CountOptions(Node tempNode,
                 matrix * mtxHopOpt,
                 const_SNarray snA);

/* This function is designed to calculate the probability a charge
 * will hop to a site within the cluster. It works by initially
 * assiging a probability of 1 to every site within the cluster. 
 * Information is then exchanged between nearest neighbors and 
 * the probability a charge will occupy a site within the cluster
 * is updated and normalized. This process is repeated for the 
 * specified number of attempts. 
 *
 * Needed before calling function 
 *   1. Determined which nodes are part of the cluster
 *   2. Join nodes to TempClLL
 * 
 * Values returned 
 *   1. mtxProb
*/
matrix CalculateProb(const_ClusterLL TempClLL,
                     matrix mtxHopOpt,
                     const_SNarray snA,
                     const int attempts);

/* Function will calculate the escape rate from each of the
 * the sites that is within the cluster. The output mtxTescape
 * has a number of rows determined by the number of nodes in the 
 * cluster and a column width of 2. The first column is the 
 * escape time. The second column is the id of the node. 
 *
 * Inputs
 *   1. mtxNeighRate
 * Output
 *   1. mtxTescape
 */
matrix CalculateTescape(ClusterLL * TempClLL,
                        matrix mtxNeighRate);

/* Will calculate the cluster course grained dwell time
 * To do this the following inputs are needed:
 *   1. mtxHopOffSiteProb
 *   2. mtxTescape
 *
 * Output
 *   1. (double) ClusterDwell
 */
double CalculateClusterDwellTime(matrix mtxHopOffSiteProb,
                                 matrix mtxTescape);

matrix CalculateNeighPval(ClusterLL * TempClLL,
                          const_matrix mtxNeighProb);

matrix CalculateHopOffSiteProb(ClusterLL * TempClLL,
                               SNarray snA,
                               matrix mtxDwellProb,
                               matrix mtxNeighRate,
                               matrix mtxProb);

matrix CalculateNeighProb(ClusterLL * TempClLL,
                          SNarray snA,
                          matrix mtxDwellProb,
                          matrix mtxNeighRate,
                          matrix mtxProb);

matrix CalculateNodeDwellProb(ClusterLL * TempClLL,
                              SNarray snA);

matrix CalculateNeighRateProb(ClusterLL * TempClLL,
                              matrix mtxProb,
                              matrix MasterM,
                              SNarray snA,
                              int PeriodicX,
                              int PeriodicY,
                              int PeriodicZ);

int CalculatePvalNodes(ClusterLL * TempClLL,
                       matrix mtxProb,
                       matrix mtxDwellTime);

#endif
