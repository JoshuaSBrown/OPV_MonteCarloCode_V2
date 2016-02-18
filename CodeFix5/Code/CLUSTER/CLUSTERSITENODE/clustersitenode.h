#ifndef _CLUSTERSITENODE_
#define _CLUSTERSITENODE_

#include "../CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"

/* Prints the contents of a site node the correct
	 print function pointer must be passed as the second
	 argument. Designed to use the pointers to both
	 printPoint and printClusterLL
	 Will return -1 if sn or the printfunc is NULL
	 Will return 0 if successful
 */
int printSNfun(SiteNode sn, int (*printfunc)(void *));

/* Prints the contents of a whole SiteNode array
	 need to pass the function for printing contents
	 of cluster though
	 Will return -1 if snA is NULL or printfunc is NULL
	 Will return 0 if successful
 */
//int printSNarray(SNarray snA, int (*printfunc)(void *));

/* Prints the coordinates of the clusters their id's and
	 the order of magnitude of the hops within the cluster
 */
int PrintFile_xyz( int orderLow, SNarray snA, ArbArray *ArClLL, char * FileName);

/* Prints the coordinates of the clusters their id's and
	 the order of magnitude of the hops within the cluster
 */
int PrintNeighFile_xyz( int orderLow, SNarray snA, ArbArray *ArClLL, char * FileName);

//double getsum(const_SiteNode sn);

/* Checks if any of the sites neighboring a cluster are
	 unoccupied.
	 The point of this function is to not waste time. If
	 a charge exists in a cluster and all the neighboring 
	 sites have charges already on them the charge has no 
	 chance of escaping and thus it is simpler to just consider
	 movements within the cluster
	 Returns a 0 if at least one site is unoccupied
	 Returns a -1 is not part of a cluster. 
	 If this is the case should simply consider
	 the site as a point not a part of a cluster
	 returns a 1 if all the neighboring sites are occupied
 */
int OccAllNeiCluster(SNarray snA, int i1, int j1, int k1);

/* Checks within the Cluster to see if charges can move around
	 within the cluster
	 Returns a -1 if malformed input or if site is not connected
	 to a cluster
	 0 if a site exists within the cluster that is unoccupied
	 1 if there are no available sites within the cluster
 */
int OccAllCluster(SNarray snA, int i1, int j1, int k1);
//
///* Function cycles through all the sites neighboring the cluster
//	 that site i1, j1, k1 is a part of. Looks at all the p values
//	 and finds the highest one.
//	 It also calcultes the sum of the pvalues. If there is an
//	 electrode neighboring the cluster is also takes into account
//	 the pvalue of the electrode
// */
//void getNeighClusterPvalHigh(SNarray snA, int i1, int j1, int k1, double * pvalHigh,long double * sum);
//
///* Function cycles through all the sites within the cluster
//	 that i1, j1, k1 is a part of. Looks at all the p values
//	 and finds the lowest one.
//	 Also calculates the sum of all the pvalues does not
//	 account for electrodes because electrodes are not allowed
//	 to exist in the cluster.
// */
//void getClusterPvalLow(SNarray snA, int i1, int j1, int k1, double * pvalLow,long double * sum2);

/* Function finds the new site neighboring the cluster that
	 the charge is to move to, as well as the time it takes
	 to hop
*/
int HopOffCluster(SNarray snA, int ID, double position2, int * newID, double * timeOff);

/* If a charge hops within cluster then this function 
	 determines if the charge moves to a separate site 
	 within the cluster:
	 0 - means moved to a new site within the cluster
	 -1 - means malformed input
*/
int HopWithinCluster(SNarray snA, int ID, double position2, int * newID);

/* Function determines if a charge occupying the site with
	 index ID will hop on or off the cluster it is a part of
	 depending on the value of (position). Position is just
	 a random number between 0 and 1.
	 1 - a neighboring site
	 0 - stays on the cluster;
	 -1 - fails malformed input
*/
int HopOnOffCluster(SNarray snA, int ID, double position);
//
//
///* Gets the p values of the given sitenode
//	 choice = 0 assumes you want to grab it from a
//	 point datastructure
//	 choice = 1 assumes you want to grab from a
//	 cluster LL
// */
//double getp(SiteNode sn, int num, int choice);
//
///* Same as get p, with the exception you are
//	 setting p as opposed to getting it.
//	 num descripes which p value you are grabbing
//	 Choice describes what action is being taken
//
//	 choice = 0 setting the pval of a point
//	 num must be between 0 - 5
//
//	 choice = 1 setting the pval of a site in
//	 a cluster
//	 num must be the id of a site within the cluster
//
//	 choice = 2 setting the pval of a site Neighboring
//	 a cluster
//	 num must be the id of a site neighboring the
//	 cluster
//
//	 choice = 3 setting the pval of an electrode
//	 neighboring the cluster
//	 num must be set to 0
//
//	 choice = 4 setting the pval of a cluster as
//	 a whole.
//	 num must be set to 0
//
//	 In all cases p is the value being assigned
// */
//void setp(SiteNode sn, int num, double var, int choice);
//
#endif
