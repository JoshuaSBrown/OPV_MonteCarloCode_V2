#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <stdarg.h>
#include <stdlib.h>

/* MidPoint structure are points between 
	 sites in a lattice. They contain the information
	 regarding the order of magnitude of the hopping 
	 rates between the two sites. 
	 Information is only stored in a mid point when both
	 hop rates from site1 to site2 and from site2 to site1
	 are on the same order of magnitude.
*/
typedef struct _MidPoint * MidPoint;
typedef struct _MidPoint const * const_MidPoint;

/* Multiple Link list of this type are used to sort the
	 mid pts such that all the mid points of a given order
	 exists in a given link list. 
	 Mid_ID ranges between 0 and < total reserved
*/
typedef struct _OrderMagLL * OrderMagLL;
typedef struct _OrderMagLL const * const_OrderMagLL;

typedef struct _Node * Node;
typedef struct _Node const * const_Node;

typedef struct _NeighNode * NeighNode;
typedef struct _NeighNode const * const_NeighNode;

/* Each NeighborNode can have multiple hops that
	 end on it. The data structure Hop is used
	 to keep track of how many hops are allowed
	 as well as the time it takes when a hop is 
	 made and the probability that the specific 
	 hop is made.
*/
typedef struct _Hop * Hop;
typedef struct _Hop const * const_Hop;

typedef struct _Electrode * Electrode;
typedef struct _Electrode const * const_Electrode;

typedef struct _NeighLL * NeighLL;
typedef struct _NeighLL const * const_NeighLL;

/* The ClusterLL is a link list that points to
	 all the nodes in the cluster
	 It also points to another link list containing
	 all the nieghboring sites
	 And it points to the electrode datastructure
*/
typedef struct _ClusterLL * ClusterLL;
typedef struct _ClusterLL const * const_ClusterLL;

/* This type is an array where each element points
	 to either a MidPoint an OrderMagLL or a ClusterLL
	 The type of data held in the array is recorded as
	 a type attribute in the array. i.e.
	 Arb->type=0 - Stores OrderMagLL
	 Arb->type=1 - Stores ClusterLL
	 Arb->type=2 - Stores MidPoints
*/
typedef struct _ArbArray * ArbArray;
typedef struct _ArbArray const * const_ArbArray;
/////////////////////////////////////////////////////////
//Tools for accessing Mid Points
/* newMidPoint creates a mid point between two neigboring
	 sites. It stores the order of magnitude of the hoprate
	 between the two sites as well.
	 Function will return NULL if it is unable to allocate
	 memory for the point or if the other parameters are 
	 inappropriate
*/
MidPoint newMidPoint(int order, int Mid_ID, int nei1, int nei2);

/* Delete MidPoint simply deletes the Mid point passed
	 Will return -1 if mp is NULL
	 Will return 0 if successful
*/
int deleteMidPoint(MidPoint mp);

/* For printing the contents of a mid point
	 Will return 0 if successful and -1 if
	 mp is NULL
 */
int printMP(const_MidPoint mp);

/* Grabs the order of magnitude of the midpoint
*/
int getMP_order(const_MidPoint mp);

/* Grab the Mid Point Id return as an integer
	 or returns -1 if mp is NULL
*/
int getMP_id(const_MidPoint mp);

/* Set the Mid Point id will return 0 if success
	 and -1 if failed
*/
int setMP_id(MidPoint mp, int ID);

/* This function compares the neighbors of two mid points
	 if the mid points share neighboring sites than the 
	 function returns a 1. If they do not share sites then
	 a 0 is returned. 
*/
int CompareNeiMidPoint(const_MidPoint mp1, MidPoint mp2);

/////////////////////////////////////////////////////////
//Tools for accessing Order Magnitude Link List
/* This function creates an order of magnitude link list
	 starting at the mid point mp. The order of magnitude
	 is also defined and stored in the orderMag attribute
*/
OrderMagLL newOrLL(int orderMag);

/* This function only deletes the OrderMagnitude data
	 type. It does not delete the MidPoints referenced
	 by the link List. 
*/
int deleteOrLL(OrderMagLL OMLL);

/* This function deletes the Midpoints held in the link
	 list as well as the actual OMLL.
*/
int deleteAllOrLL(OrderMagLL OMLL);

/* Print mid points from link list
*/
int printOrLL( const_OrderMagLL OMLL);

/* Grab size of OMLL
*/
int getOMLL_size(const_OrderMagLL OMLL);

/* Grans the value of the order
*/
int getOMLL_order(const_OrderMagLL OMLL);

/* Function grabs the starting Midpoint of the OMLL
*/
MidPoint getOMLLstartMP(const_OrderMagLL OMLL);

/////////////////////////////////////////////////////////
//Tools for accessing Nodes
/* Function creates a Node. The Node acts like a
	 Link list it contains the id of a node and points to
	 another Node that is also in the same cluster
*/
Node newNode(int N_ID);

/* Deletes a node
*/
int deleteNode(Node Nod);

/* Prints contents of node attributes
*/
int printNode(const_Node Nod);

/* get the node id and return as an integer
*/
int getNode_id(const_Node Nod);

/* get the next node
*/
Node getNextNode(const_Node Nod);

/* Set the Node id
*/
int setNode_id(Node Nod, int ID);

/* Set Node's probability value
*/
int setNode_p(Node Nod, double pval);

/* Get Node's probability value
*/
double getNode_p(const_Node Nod);

/* Set flags
	 If a node is within the 
	 cluster it is set to 1
	 if it is outside it is
	 left as a 0
*/
int setFlagFro(Node Nod);
int setFlagBeh(Node Nod);
int setFlagLef(Node Nod);
int setFlagRig(Node Nod);
int setFlagAbo(Node Nod);
int setFlagBel(Node Nod);

int getFlagFro(const_Node Nod);
int getFlagBeh(const_Node Nod);
int getFlagLef(const_Node Nod);
int getFlagRig(const_Node Nod);
int getFlagAbo(const_Node Nod);
int getFlagBel(const_Node Nod);

int getFlag(const_Node Nod, int index);
int setFlag(Node Nod, int index);
/////////////////////////////////////////////////////////
//Tools for accessing Neighbor Nodes
NeighNode newNeighNode(int N_ID);

/* Deletes a node
*/
int deleteNeighNode(NeighNode NeighNod);

/* Prints contents of node
*/
int printNeighNode(const_NeighNode NeighNod);

/* Get the NeighNode ID and set as an integer
*/
int getNeighNode_id(const_NeighNode NeighNod);

int setNextNeighNode(NeighNode * NeighNod, NeighNode * NeighNod2);

/* Gets the next Neighbor in line
*/
NeighNode getNextNeigh(const_NeighNode NeighNod);

/* Sets the probability of moving to this neighboring node
	 for a new hop
*/
int setNeighNodeNew_p(NeighNode NeighNod, double pval);

/* Set the probability of moving to this neighboing node
	 taking the hop route specified by Elem
*/
int setNeighNode_p(NeighNode NeighNod, double pval, int Elem);

/* Sets the time of moving to this neighboring node
	 given that it takes the hop route specified by 
	 Elem.
*/
int setNeighNode_t(NeighNode NeighNod, double time, int Elem);

/* Gets the value of probability of hopping to NeighNod from
	 the hop route speciefied by Elem.
*/
double getNeighNode_p(const_NeighNode NeighNod, int Elem);

/* Gets the length of the number of hops to a given neighboring
	 site off the cluster
*/
int getNeighNode_hoplength(NeighNode NeighNod);

/* Gets the value of the time of hopping to NeighNod from
	 the hop route specified by Elem.
*/
double getNeighNode_t(NeighNode NeighNod, int Elem);

/* Set the NeighNode Id
*/
int setNeighNode_id(NeighNode NeighNod, int ID);

/* Set the NeighNode hoplength
*/
int setNeighNode_hoplength(NeighNode NeighNod, int hl);

/* Sets the starting hop from the NeighNod
*/
int setNeighNode_hopstart(NeighNode NeighNod, Hop h);

//////////////////////////////////////////////////////////
/* Creates a new Electrode
	 ID defined where the electrode is
	 0 - (-(x,y,z))
	 1 - (+(x,y,z))
*/
Electrode newElectrode(void);

/* Deletes the electrode
*/
int deleteElectrode(Electrode *el);

/* Sets the tunneling constant for the electrode
*/
int setElectrode_alpha(Electrode el, double alpha);

/* Grab the tunneling constant for hopping 
	 from the electrode to the medium
*/
double getElectrode_alpha(Electrode el);

/* Sets the fermi energy of the electrode
*/
int setElectrode_FermiEnergy(Electrode el, double FermiE);

/* gets the fermi energy of the electrode
*/
double getElectrode_FermiEnergy(Electrode el);

/* Get number of charges in electrode
*/
int getElectrode_Charges(Electrode el);

/* Sets the number of charges on the electrode
*/
int setElectrode_Charges(Electrode el, double NumCharge);

/* Increases the number of charges on the electrode by 1
*/
int Electrode_addCharge(Electrode el);

/* Decreases the number of charges on the electrode by 1
	 if there are 0 and try to reduce will return -1
*/
int Electrode_minusCharge(Electrode el);

/* connects to a structure that has the pvals of hops to 
	 the site neighboring the electrode.
*/
int setElectrode_HopRates(Electrode el, void * HopS);

/* Grabs the datastructure containing hop rates off
	 the electrode
*/
void * getElectrode_HopRates(Electrode el);

/* Sets the data structure containing the adjacent site
	 information
*/
int setElectrode_AdjacentSites(Electrode el, void * snA);

/* Gets the data structure with the adjacent site information
*/
void * getElectrode_AdjacentSites(Electrode el);

/* Sets the sum of the electrode used when calculating
	 the transition time off the electrode
*/
int setElectrode_Sum(Electrode el, double sum);

/* Grabs the sum of the electrode
*/
double getElectrode_Sum(Electrode el);

/* Prints contents of electrode id and p value
*/
int printElectrode(const_Electrode el);

///////////////////////////////////////////////////////////
/* Creates a NeighLL datastructure contains the size, id and
	 the neighnodes attached to a cluster
*/
NeighLL newNeighLL(void);

/* Deletes the neighLL
*/
int deleteNeighLL(NeighLL neighLL);

/* Deletes the neighLL and all the NeighNodes in the neighLL
*/
int deleteNeighLLALL(NeighLL neighLL);

/* Prints contents of neighLL
*/
int printNeighLL(const_NeighLL neighLL);

int setNeighLL_start(NeighLL Nei, NeighNode NeighNod);

int setNeighLL_numNeigh(NeighLL Nei, int NumNei);

/////////////////////////////////////////////////////////
//Tools for accessing Hop data structure (no getters puroposefully)
Hop newHop(void);

//Deletes all Hops linked to the hop h
int deleteAllHop(Hop h);

int printHop(Hop h);

int setHop_t(Hop h, double t);

int setHop_p(Hop h, double p);

int setHop_next(Hop h, Hop h2);
/////////////////////////////////////////////////////////
//Tools for accessing Clusters and Cluster Link Lists
/* Creates a Cluster Link list which contains all the 
	 clusters in the material also contains either
	 Nodes or Neighnodes specefied by type:
	 type 0 - Equivalent to Nodes
	 type 1 - Equivalent to NeighNodes
*/
ClusterLL newClusterLL(int ID);

/* Just deletes the ClusterLL datastructure leaves the 
	 nodes etc
*/
int deleteClusterLL( ClusterLL clLL);

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

/* Prints all the Clusters in a given Cluster LL
	 and their respective nodes or site ids
	 If clLL uses Nodes or NeighNodes it should
	 function correctly
*/
int printClusterLL(const_ClusterLL clLL);

int setCluster_NeighLL(ClusterLL clLL, NeighLL NeiLL);

/* Gran the NeighLL from the clsuter clLL
*/
NeighLL getCluster_NeiLL(const_ClusterLL clLL);

/* Get the cluster ID and returns as an integer
*/
int getCluster_id(const_ClusterLL clLL);

/* Gets the size of the cluster and returns it
*/
int getCluster_numNodes(const_ClusterLL clLL);

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

/* This function works by adding the pval to the correct Node 
	 in the NeichclLL link list. If the (Node_ID) returns a -1
	 it means the Node was not found in the link list.
	 If it is found it add a new Hop structure and add the pval
	 to the new Hop structure. 
*/
int addPvalClusterNeighNode( ClusterLL NeighclLL, int Node_ID, double pval);

/////////////////////////////////////////////////////////
//Tools for accessing Arbirtrary array

/* Creates a new Arbitrary Array the size of the array is 
	 input as an integer and stored as a reserved attribute
	 the number of elements used in the array is defined by
	 the used attribute.
	 The elements range between 0 and < the total reserved
	 The array is defined as a specefic type. Where the 
	 elements are a certain data structure type i.e.
	 type = 0 - Order of Magnitude Link List Array
	 type = 1 - Cluster Link List Array
	 type = 2 - Mid Point Array
*/
ArbArray newArbArray(int len, int type);

/* The purpose of this function is to ensure that stack
	 overflow did not occur will return a -1 if the type 
	 of ArbArray has been changed to a number that is not
	 0, or 1 (Only checks to ensure ArbArray is an Array
	 of Link list or Order of magnitude link list not 
	 a mid point array which would be of type 2). 
*/
int ArbArrayCheck(const_ArbArray Arb);

int deleteAllMidPointArray(ArbArray * Arb);

/* This function deletes the arbitray array according to
	 it's type if it is of type 0 it only deletes the Order
	 of magnitude data structure and does not delete the 
	 Mid Points. 
	 If it is of type 1 it deletes all the nodes and the
	 link list.
	 If it is of type 2 it does not delete the mid points
	 only the Arbitrary array itself.
	 Returns a 0 if succesful and -1 if failed
*/
int deleteArbArray(ArbArray * Arb);

/* Prints the contents of the arbitrary array according
	 to it's type, a second parameter should be submitted
	 if the array is of type 1. The lowest order of 
	 magnitude should be submitted as a second argument:
	 printArbArray(ArbArray Arb, int orderLow);
*/
int printArbArray(const_ArbArray Arb, ...);

/*Sets the element of Arb to NULL
*/
int NullArbElement(ArbArray Arb, int element);

/* Sets an element of the ArbArray, increments the number
	 of elements used if it is NULL
*/
int setArbElement(ArbArray Arb, int element, void * ptr);

/* Grabs the elements of the ArbArray does not remove it
	 from the array though.
*/
void * getArbElement(const_ArbArray Arb, int element);

/* Returns the number of elements used
*/
int getElementsUsed(const_ArbArray Arb);

/* Returns the number of elements resereved
*/
int getElementsReserved(const_ArbArray Arb);

/* Removes the elements from the array and makes it NULL
*/
int removeArbElement(ArbArray Arb, int element);

/* Grab neighbor one of midpoint Mid_ID specified in the
	 array mpA
	 will return the nei1 id or -1 if failed
*/
int getMPnei1(const_ArbArray Arb, int Mid_ID);

/* Grab neighbor two of midpoint Mid_ID specified in the 
	 array mpA
	 will return the nei1 id or -1 if failed
*/
int getMPnei2(const_ArbArray Arb, int Mid_ID);

/* Get the mid point from mid point array mpA with id 
	 Mid_ID
*/
MidPoint getMP(const_ArbArray Arb, int Mid_ID);

/* Get the next midpoint in the link list
*/
MidPoint getNextMP(const_MidPoint mp);

/* Get the order of magnitude of the hopping rate for mid
	 Mid_ID from midpoint array mpA.
*/
int getMPOrder(const_ArbArray Arb, int Mid_ID);

/* Add a Mid Point to the Order of Magnitude Link list 
	 array specefied as Arb->Num[element]. This link list 
	 contains elements that point to the order of mangitude
	 link lists. 
	 If sucessful return 0
	 If Mid point was previously added returns -1
*/
int addToOrLL(ArbArray Arb, int element, MidPoint mp);

/* This function gets the "OrderMagLL" and returns the 
	 pointer.
*/
OrderMagLL getOrderLL(const_ArbArray Arb, int element);

/* This function sets the default values for the order of 
	 magnitude LL in Arb->Num[element] it will only work
	 if the Arbitrary array used is of type 1
	 If a new LL is added and the elements was not NULL than
	 the size will be incremented.
*/
int setDefaultArbElem(ArbArray Arb, int element, int orderMag);

#endif
