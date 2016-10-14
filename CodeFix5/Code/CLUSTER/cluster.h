#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <stdarg.h>
#include <stdlib.h>

#include "../MIDPOINT/midpoint.h"
#include "../ELECTRODE/electrode.h"
#include "../NEIGHLL/neighll.h"
#include "../NEIGHNODE/neighnode.h"
#include "../NODE/node.h"

/* The ClusterLL is a link list that points to
	 all the nodes in the cluster
	 It also points to another link list containing
	 all the nieghboring sites
	 And it points to the electrode datastructure
*/
typedef struct _ClusterLL * ClusterLL;
typedef struct _ClusterLL const * const_ClusterLL;

/////////////////////////////////////////////////////////
//Tools for accessing Clusters and Cluster Link Lists
/* Creates a Cluster Link list which contains all the 
	 clusters in the material also contains either
	 Nodes or Neighnodes specefied by type:
	 type 0 - Equivalent to Nodes
	 type 1 - Equivalent to NeighNodes
*/
ClusterLL newClusterLL(int ID);

/* This function appends a cluster to the end of another
 * ClusterLl
 */
int appendClusterLL(ClusterLL clLL, ClusterLL clLLnew);

/* Just deletes the ClusterLL datastructure leaves the 
	 nodes etc
*/
int deleteClusterLL( ClusterLL clLL);

/* Just deletes the nieghbor nodes 
 */
int deleteClusterLL_NeighNodes(ClusterLL * clLL);

/* Removes the second cluster from the first clusterLL
 * assuming that the second cluster is linked somewhere
 * within the first one. Returns 
 *  0 if success
 * -1 if failed
 */
int removeClusterLLfromClusterLL(ClusterLL * ClLLAll, ClusterLL clLL);

/* This function deletes all the nodes linked to clLL
	 However it will not work if clLL->next points to 
	 another link list. This is a safety measure to ensure
	 that we still have access to the other link lists 
	 attached to clLL. 
	 it does not delete any of the other ClusterLL that
	 are attached to clLL. It only deletes itself. 
*/
int deleteClusterLLNodes( ClusterLL clLL);

/* Deletes all the Link Lists attached to ClusterLL clLL
	 Will delete:
	 o All neighLL 
	 	- All neighNodes in neighLL
	 o The Electrodes if they exist
	 o All the Nodes within clLL
	 o all the ClusterLL attached to clLL and subsequently
	 all the electrodes Nodes and neighLL and NeighNodes. 
*/
int deleteAllClusterLL( ClusterLL * clLL);

/* Gets the next ClusterLL that is attached to the 
	 original clLL
*/
ClusterLL getNextClusterLL(const_ClusterLL clLL);

/* Prints out all the nodes in the cluster passed
	 does not print out clusters that are attached 
	 to the link list.
*/
int printNodesClusterLL(const_ClusterLL clLL);

int printNeighNodesClusterLL(const_ClusterLL clLL);

/* Prints all the Clusters in a given Cluster LL
	 and their respective nodes or site ids
	 If clLL uses Nodes or NeighNodes it should
	 function correctly
*/
int printClusterLL(const_ClusterLL clLL);

/* This function will search through clLL until 
 * it finds a cluster with the matching ClusterID
 * if unable to find it, the function will return 
 * NULL
 */
ClusterLL getClusterGivenClusterID(ClusterLL clLL, int ClusterID);

int setCluster_NeighLL(ClusterLL clLL, NeighLL NeiLL);

/* Gran the NeighLL from the clsuter clLL
*/
NeighLL getCluster_NeiLL(const_ClusterLL clLL);

/* Set the number of Nodes in the Cluster returns
 *  0 if success
 * -1 if fails
 */
int setCluster_NumNodes(ClusterLL clLL, int NumNodes);

/* Get the cluster ID and returns as an integer
*/
int getCluster_id(const_ClusterLL clLL);

/* Gets the size of the cluster and returns it
*/
int getCluster_numNodes(const_ClusterLL clLL);

/* Get the last Node in the Cluster
 */
Node getCluster_LastNode(const_ClusterLL clLL);

int setCluster_elecXBid(ClusterLL * clLL, int Elec);
int setCluster_elecXFid(ClusterLL * clLL, int Elec);
int setCluster_elecYLid(ClusterLL * clLL, int Elec);
int setCluster_elecYRid(ClusterLL * clLL, int Elec);
int setCluster_elecZBid(ClusterLL * clLL, int Elec);
int setCluster_elecZAid(ClusterLL * clLL, int Elec);

/* sets the electrode attribute for the cluster link list
	 Electrode objectect Elec is connected to the cluster
	 will return -1 if the cluster was not previously identified
	 as being next to the electrode
*/
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYL(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZA(ClusterLL * clLL, Electrode Elec);

/* Sets the hop probability attribute for electrodes
	 along the X axis
*/
int setCluster_elecXsum(ClusterLL clLL, double val, int Elec);

/* Sets the hop probability attribute for electrodes
	 along the Y axis
*/
int setCluster_elecYsum(ClusterLL clLL, double val, int Elec);

/* Sets the hop probability attribute for electrodes
	 along the z axis
*/
int setCluster_elecZsum(ClusterLL clLL, double val, int Elec);

/* Adds to sum already stored in the cluster clLL
*/
int addToCluster_elecXsum(ClusterLL clLL, double val, int Elec);
int addToCluster_elecYsum(ClusterLL clLL, double val, int Elec);
int addToCluster_elecZsum(ClusterLL clLL, double val, int Elec);

/* Functions will check to see if respective electrodes
 * are defined.
 */
/* If front electrode is defined will return 1
 * If just back electrode is defined will return 2
 * If both electrodes are defined will return 3
 */
int checkCluster_elecXdefined(const_ClusterLL clLL);
/* If right electrode is defined will return 1
 * If just left electrode is defined will return 2
 * If both electrodes are defined will return 3
 */
int checkCluster_elecYdefined(const_ClusterLL clLL);
/* If above electrode is defined will return 1
 * If just below electrode is defined will return 2
 * If both electrodes are defined will return 3
 */
int checkCluster_elecZdefined(const_ClusterLL clLL);

int getCluster_elecidXB(const_ClusterLL clLL);
int getCluster_elecidXF(const_ClusterLL clLL);
int getCluster_elecidYL(const_ClusterLL clLL);
int getCluster_elecidYR(const_ClusterLL clLL);
int getCluster_elecidZB(const_ClusterLL clLL);
int getCluster_elecidZA(const_ClusterLL clLL);


/* Gets the hop probability to the electrode
*/
double getCluster_elecXsum(const_ClusterLL clLL, int Elec);
double getCluster_elecYsum(const_ClusterLL clLL, int Elec);
double getCluster_elecZsum(const_ClusterLL clLL, int Elec);

/* Sets the Cluster id
*/
int setCluster_id(ClusterLL clLL, int ID);

/* Add to the value of sum already stored in the
	 cluster.
*/
int addToClusterSum(ClusterLL clLL, double s);

/* Sets the sum attribute of the cluster
*/
int setClusterSum(ClusterLL clLL, double s);

/* Gets the sum attribute of the cluster
*/
double getCluster_Sum(const_ClusterLL clLL);

/* Sets the next cluster to Nex
*/
int setNextClusterLL(ClusterLL clLL, ClusterLL Nex);

/* Gets the Node that starts the Link List
*/
Node getStartNode(const_ClusterLL clLL);

/* Returns the Node with identity id
*/
Node getCluster_Node(const_ClusterLL clLL, int id);

/* Gets the NeighNode that starts the Link List
*/
NeighNode getStartNeigh(const_ClusterLL clLL);

/* Gets the Number of neighboiring nodes to cluster
	 clLL
*/
int getCluster_numNeigh(const_ClusterLL clLL);

/* gets the size of the cluster the number of nodes in it
*/
int getCluster_numNodes(const_ClusterLL clLL);

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);
/* If front electrode is defined will return 1
 * If just back electrode is defined will return 2
 * If both electrodes are defined will return 3
 */
int checkCluster_elecXdefined(const_ClusterLL clLL);

int getCluster_elecidXB(const_ClusterLL clLL);
int getCluster_elecidXF(const_ClusterLL clLL);
int getCluster_elecidYL(const_ClusterLL clLL);
int getCluster_elecidYR(const_ClusterLL clLL);
int getCluster_elecidZB(const_ClusterLL clLL);
int getCluster_elecidZA(const_ClusterLL clLL);


/* Gets the hop probability to the electrode
*/
double getCluster_elecXsum(const_ClusterLL clLL, int Elec);
double getCluster_elecYsum(const_ClusterLL clLL, int Elec);
double getCluster_elecZsum(const_ClusterLL clLL, int Elec);

/* Sets the Cluster id
*/
int setCluster_id(ClusterLL clLL, int ID);

/* Add to the value of sum already stored in the
	 cluster.
*/
int addToClusterSum(ClusterLL clLL, double s);

/* Sets the sum attribute of the cluster
*/
int setClusterSum(ClusterLL clLL, double s);

/* Gets the sum attribute of the cluster
*/
double getCluster_Sum(const_ClusterLL clLL);

/* Sets the next cluster to Nex
*/
int setNextClusterLL(ClusterLL clLL, ClusterLL Nex);

/* Gets the Node that starts the Link List
*/
Node getStartNode(const_ClusterLL clLL);

/* Returns the Node with identity id
*/
Node getCluster_Node(const_ClusterLL clLL, int id);

/* Gets the NeighNode that starts the Link List
*/
NeighNode getStartNeigh(const_ClusterLL clLL);

/* Gets the Number of neighboiring nodes to cluster
	 clLL
*/
int getCluster_numNeigh(const_ClusterLL clLL);

/* gets the size of the cluster the number of nodes in it
*/
int getCluster_numNodes(const_ClusterLL clLL);

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);
/* If front electrode is defined will return 1
 * If just back electrode is defined will return 2
 * If both electrodes are defined will return 3
 */
int checkCluster_elecXdefined(const_ClusterLL clLL);

int getCluster_elecidXB(const_ClusterLL clLL);
int getCluster_elecidXF(const_ClusterLL clLL);
int getCluster_elecidYL(const_ClusterLL clLL);
int getCluster_elecidYR(const_ClusterLL clLL);
int getCluster_elecidZB(const_ClusterLL clLL);
int getCluster_elecidZA(const_ClusterLL clLL);


/* Gets the hop probability to the electrode
*/
double getCluster_elecXsum(const_ClusterLL clLL, int Elec);
double getCluster_elecYsum(const_ClusterLL clLL, int Elec);
double getCluster_elecZsum(const_ClusterLL clLL, int Elec);

/* Sets the Cluster id
*/
int setCluster_id(ClusterLL clLL, int ID);

/* Add to the value of sum already stored in the
	 cluster.
*/
int addToClusterSum(ClusterLL clLL, double s);

/* Sets the sum attribute of the cluster
*/
int setClusterSum(ClusterLL clLL, double s);

/* Gets the sum attribute of the cluster
*/
double getCluster_Sum(const_ClusterLL clLL);

/* Sets the next cluster to Nex
*/
int setNextClusterLL(ClusterLL clLL, ClusterLL Nex);

/* Gets the Node that starts the Link List
*/
Node getStartNode(const_ClusterLL clLL);

/* Returns the Node with identity id
*/
Node getCluster_Node(const_ClusterLL clLL, int id);

/* Gets the NeighNode that starts the Link List
*/
NeighNode getStartNeigh(const_ClusterLL clLL);

/* Gets the Number of neighboiring nodes to cluster
	 clLL
*/
int getCluster_numNeigh(const_ClusterLL clLL);

/* gets the size of the cluster the number of nodes in it
*/
int getCluster_numNodes(const_ClusterLL clLL);

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);

int getCluster_elecidXB(const_ClusterLL clLL);
int getCluster_elecidXF(const_ClusterLL clLL);
int getCluster_elecidYL(const_ClusterLL clLL);
int getCluster_elecidYR(const_ClusterLL clLL);
int getCluster_elecidZB(const_ClusterLL clLL);
int getCluster_elecidZA(const_ClusterLL clLL);


/* Gets the hop probability to the electrode
*/
double getCluster_elecXsum(const_ClusterLL clLL, int Elec);
double getCluster_elecYsum(const_ClusterLL clLL, int Elec);
double getCluster_elecZsum(const_ClusterLL clLL, int Elec);

/* Sets the Cluster id
*/
int setCluster_id(ClusterLL clLL, int ID);

/* Add to the value of sum already stored in the
	 cluster.
*/
int addToClusterSum(ClusterLL clLL, double s);

/* Sets the sum attribute of the cluster
*/
int setClusterSum(ClusterLL clLL, double s);

/* Gets the sum attribute of the cluster
*/
double getCluster_Sum(const_ClusterLL clLL);

/* Sets the next cluster to Nex
*/
int setNextClusterLL(ClusterLL clLL, ClusterLL Nex);

/* Gets the Node that starts the Link List
*/
Node getStartNode(const_ClusterLL clLL);

/* Returns the Node with identity id
*/
Node getCluster_Node(const_ClusterLL clLL, int id);

/* Gets the NeighNode that starts the Link List
*/
NeighNode getStartNeigh(const_ClusterLL clLL);

/* Gets the Number of neighboiring nodes to cluster
	 clLL
*/
int getCluster_numNeigh(const_ClusterLL clLL);

/* gets the size of the cluster the number of nodes in it
*/
int getCluster_numNodes(const_ClusterLL clLL);

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYL(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZA(ClusterLL * clLL, Electrode Elec);

int setCluster_elecXF(ClusterLL * clLL, Electrode Elec);
int setCluster_elecXB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYR(ClusterLL * clLL, Electrode Elec);
int setCluster_elecYL(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZB(ClusterLL * clLL, Electrode Elec);
int setCluster_elecZA(ClusterLL * clLL, Electrode Elec);

/* Sets the id of the electrodes
	 0 - back side 1 - front side
	 0 - left side 1 - right side
	 0 - bottom side 1 - above side
*/
int setCluster_elecXid(ClusterLL * clLL, int Elec);
int setCluster_elecYid(ClusterLL * clLL, int Elec);
int setCluster_elecZid(ClusterLL * clLL, int Elec);

/* Sets the pvalue of the electrode
*/
int setCluster_elecXp(ClusterLL clLL, double p, int Elec);
int setCluster_elecYp(ClusterLL clLL, double p, int Elec);
int setCluster_elecZp(ClusterLL clLL, double p, int Elec);

/* Gets the p value of the cluster as a whole
*/
double getCluster_p(const_ClusterLL clLL);

/* Sets the p value of a cluster as a whole
*/
int setCluster_p(ClusterLL clLL, double p);

/* Check if the Cluster clLL is attached to
	 any electrodes return -1 if input is 
	 NULL
	 0 if not attached to any electrodes
	 1 if attached to at least 1 electrode
*/
int checkCluster_elec(const_ClusterLL clLL);

/* get the id of the electrode
	 will return:
	 2 - if both electrodes are touching
	 the cluster
	 1 - if either the Forward, Left, Below
	 electrode are touching the cluster
	 0 - if the either the Backward, Right,
	 and Top electrode are touching the 
	 cluster
	 -1 if no electrodes are touching
*/
int getCluster_elecXid(const_ClusterLL clLL);
int getCluster_elecYid(const_ClusterLL clLL);
int getCluster_elecZid(const_ClusterLL clLL);

/* Sets the time of the Cluster
*/
int setCluster_time(ClusterLL clLL, double t);

/* Grabs the time of the cluster
*/
double getCluster_time(const_ClusterLL clLL);

/* Adds a Node at the end of a ClusterLL provided that the 
	 site ID is provided via nei, the node is actually created
	 in this function
*/
int addNodeEndClusterLL(ClusterLL clLL, int nei);

/* Adds an actual node to the cluster link list, the node
	 was created outside of the cluster link list. This function
	 simply links the node with the cluste link list. 
*/
int addClusterLLNode(ClusterLL *clLL, Node nd);

/* This function adds the pval to the correct Node that is
	 within the cluster
	 returns a -1 if the node was not found in the link list
*/
int addPvalClusterNode( ClusterLL clLL, int Node_ID, double pval);

/* This function adds Nodes to a cluster when a mid point is
	 passed to it. The mid point is between two nodes. 
	 If both nodes are not found amidst any of the clusters a
	 new cluster is made and both are added to the same cluster 
	 function returns 0
	 If one node is found in a cluster but the other is not
	 the second node is added to the first cluster 
	 function returns 0
	 If both nodes are found in the same cluster nothing is done
	 function returns -1
	 If Each node is in a separate cluster the two clusters are
	 merged taking the id of the lower cluster.
	 function returns 0
*/
int addNodeToCluster( ClusterLL clLL, MidPoint mp);

/* This function adds NeighNodes to a Cluster link list an 
	 integer is returned based on the execuation.
	 0 - means NeighNode was added at the beginning of the 
	 list
	 0 - means was added at the end of the list
	 0 - means that the node was already found in the list
	 and nothing was changed
	 -1 - means incorrect inputs
*/
int addNeighNodeToCluster(ClusterLL * clLL, int Neigh_ID);

/* Creates a new cluster and adds the first two nodes
	 nei1 and nei2 (Both these represent site IDs)
	 The cluster_id is assigned as 1 greater than the last
	 non NULL cluster.
*/
int addNodeNewCluster(ClusterLL * clLL, MidPoint mp);

/* Adds the siteID1 and siteID2 nodes to cluster clLL
 */
int addNodesToClusterGivenSites(ClusterLL clLL, int siteID1, int siteID2);

/* This function works by adding the pval to the correct Node 
	 in the NeichclLL link list. If the (Node_ID) returns a -1
	 it means the Node was not found in the link list.
	 If it is found it add a new Hop structure and add the pval
	 to the new Hop structure. 
*/
int addPvalClusterNeighNode( ClusterLL NeighclLL, int Node_ID, double pval);
#endif
