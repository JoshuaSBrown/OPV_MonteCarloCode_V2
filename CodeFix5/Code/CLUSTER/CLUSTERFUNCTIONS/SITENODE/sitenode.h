#ifndef _SITENODE_H_
#define _SITENODE_H_

typedef struct _SiteNode * SiteNode;
typedef struct _SiteNode const * const_SiteNode;

typedef struct _SNarray * SNarray;
typedef struct _SNarray const * const_SNarray;

typedef struct _Point * Point;
typedef struct _Point const * const_Point;

/* Creates a new Site node default properties of 0
	 Except for the dwellstatus which are initially set to -1
	 This symbolizes that the site is not occupied.
	 SiteNode can be of type:
	 0 - single points
	 1 - part of a cluster

	 Will return NULL if was unable to allocate memory
*/
SiteNode newSN(void);

/* Creates a point. This is the default value contains
	 the p value and sum of a point if it is not 
	 considered to be in a cluster
	 
	 Will return NULL if unable to allocate memory
*/
Point newPoint(void);

/* Creates a new Site node array length i width j and height k
	 Will return NULL if:
	 	o Unable to allocate memory
		o The width, length or height is less than 0
*/	 
SNarray newSNarray(int len, int wid, int hei);

/* Deletes a Point
	 Will return -1 pt is found to be NULL
	 Will return 0 if successful 
*/
int deletePoint(Point pt);

/* Deletes inidividual sitenode
	 Will return -1 if sn is NULL
	 Will return 0 if successful
*/
int deleteSN(SiteNode sn);

/* Deletes sitenode array and all sitenodes in the array
	 Will return -1 if snA is found to be NULL
	 Will return 0 if successful
*/
int deleteSNarray(SNarray snA);

/* Prints the contents of a point the sum and p values
	 Will return -1 if pt is NULL
	 Will return 0 if successful
*/
int printPoint(void* vpt);

/* This function prints the id of each site followed by
	 the frequency the site was visited by a charge, whether
	 the site was visited or not and the energy of the site
	 in eV:
	 col 1	col 2		col 3 col 4 	col 5				col 6
	 i			j				k			Freq		(0 or 1)		Energy [eV]
	
	 e.g.

	 3			4				1			231			1						-5.02
	 1			3				9			0				0						-4.90

	 Site (3,4,1) was visited 231 times, the 1 indicates 
	 charges hopped to the site and the actual energy of the
	 site is -5.02 eV.
	 Site (1,3,9) was not visited thus 0's for col 4 and 5
	 and had an energy of -4.90 eV

	 The output is stored in .xyz file format and by default
	 is labeled VisitFreq.xyz

	 Will return -1 if snA is NULL
	 Will return 0 if successful
*/
int printVisitFreq(const_SNarray snA, char * FileName);

/* Grab a site node from a site node array structure
 	 by specifiying the correct indices.
	 Will return NULL if snA does not exist or if the
	 call to getindex returns a negative value
	 Else it will return a ptr to a sitenode
	 i, j and k must be from 0-(length/width/height)-1
 */
SiteNode getSN(const_SNarray snA, int i, int j, int k);

/* Given the index will grab the correct sitenode
*/
SiteNode getSNwithInd(const_SNarray snA, int Ind);

/* Sets the ptr to DataStruct can either point
	 to a Cluster Link list or to a point
	 ty must be a 0 or 1
	 0 - is to point to a Point
	 1 - is to point to a ClusterLinkList
	 Warning this function does not delete old Points
	 when it replaces them
*/
int setDataStruc(SiteNode sn, int ty, void * ptr);

/* Prints contents of Sitenode
*/
int printSN(SiteNode sn);

/* Gets the pointer to the datastructure point
*/
void * getPoint(const_SiteNode sn);

/* Gets the pointer to the datastructe clusterLL
*/
void * getClusterList(const_SiteNode sn);

/* Sets the ptr to the link list that contains
	 all the sites nieghboring the cluster
*/
//int setNeighList(SiteNode sn, void * ptr);

/* Gets the pointer to the Neighbor link list
*/
//void * getNeighList(const_SiteNode sn);

// Grabs the dwell status of the node
int getDwelStat(const_SiteNode sn);

/* The Dwell Status describes whether a
	 charge is occupying the site or not
	 if it is -1 it is unoccupied else it
	 will have a positive number or 0 which
	 describes which charge is occupying the 
	 site.
	 A -1 is returned if a value is submitted
	 which is less than -1 or the sn is NULL
*/
int setDwelStat(SiteNode sn, int stat);

/* Increment VisitFrequency by 1*/
void incVisFreq(SiteNode sn);

int setVisFreq(SiteNode sn, int freq);

double getVisFreq(const_SiteNode sn);

/* Describes whether site has been visited
	 or not can either be a 0 or a 1
*/
int setVis(SiteNode sn, double vis);

double getVis(const_SiteNode sn);

int setInitE(SiteNode sn, int E);

int getInitE(const_SiteNode sn);

/* Sets the Energy of the Sitenode 
	 Energy must be real cannot be nan
	 The site node must not be NULL
*/
int setEnergy(SiteNode sn, double Energy);

double getEnergy(const_SiteNode sn);

int setTime(SiteNode sn, double time);

int addTime(SiteNode * sn, double time);

double getTime(SiteNode sn);

int setsum(SiteNode sn, double s);

double getsum(SiteNode sn);

int setSN_p(SiteNode sn, int i, double val);

double getSN_p(SiteNode sn, int i);

double getsumPt( const_Point pt);

/* Gets the type
	 0 - means its a point
	 1 - means it is a Cluster ptr
*/
int getType(const_SiteNode sn);

/* Grabs the i j and k locations of the site given the index
	 Will return a -1 if snA is NULL or if the Index is negative
	 Will return 0 if successful
*/
int getLoc( int * i, int * j, int * k, int Index, const_SNarray snA);

/* Grabs the index of the site node array i,j and k
   must be 0<=x<snA->(length , width, height) as is
	 appropriate
	 Will return a -1 if the index exceeds the width
	 height or length of the snA or if the snA is NULL
	 Will return the Index otherwise
	 */
int getIndex( const_SNarray snA, int i, int j, int k);

/* Grabs the index of the site i-1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if i-1 is outside of the boundaries of snA
	 If less than 0 (i equivalent to x)
*/
int getIndexFrontP( const_SNarray snA, int Index);
int getIndexFront( const_SNarray snA, int Index);

/* Grabs the index of the site i+1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if i+1 is outside of the boundaries of snA
	 If greater than snA->length (i equivalent to x)
*/
int getIndexBehindP( const_SNarray snA, int Index);
int getIndexBehind( const_SNarray snA, int Index);

/* Grabs the index of the site j-1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if j-1 is outside of the boundaries of snA
	 If less than 0 (j equivalent to y)
*/
int getIndexLeftP( const_SNarray snA, int Index);
int getIndexLeft( const_SNarray snA, int Index);

/* Grabs the index of the site j+1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if j+1 is outside of the boundaries of snA
	 If greater than snA->width (j equivalent to y)
*/
int getIndexRightP( const_SNarray snA, int Index);
int getIndexRight( const_SNarray snA, int Index);

/* Grabs the index of the site k-1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if k-1 is outside of the boundaries of snA
	 If less than 0 (k equivalent to z) 
*/
int getIndexBelowP( const_SNarray snA, int Index);
int getIndexBelow( const_SNarray snA, int Index);

/* Grabs the index of the site k+1 of site (i,j,k) which
	 is given in this case as the Index
	 returns -1 if k+1 is outside of the boundaries of snA
	 If greater than snA->height (k equivalent to z)
*/
int getIndexAboveP( const_SNarray snA, int Index);
int getIndexAbove( const_SNarray snA, int Index);


/* Grabs the index of the site j-1 of site (i,j,k) which
	 is given in this case as the Index
	 Considers that sample is periodic in the y direction
	 (j equivalent to y)
*/
int getIndexLeftPeriodic( const_SNarray snA, int Index);

/* Grabs the index of the site j+1 of site (i,j,k) which
	 is given in this case as the Index
	 Considers that the sample is periodic in the y
	 direction (j equivalent to y)
*/
int getIndexRightPeriodic( const_SNarray snA, int Index); 
int getIndexBelowPeriodic( const_SNarray snA, int Index); 
int getIndexAbovePeriodic( const_SNarray snA, int Index); 


/* The following 8 functions do the same thing as the 8
	 above the difference being (i,j,k) are already konwn
	 so it is not necessary to find them each time a
	 function is called.
	 If i,j and k are outside the lattice a -1 will be
	 returned
	 If after grabbing the index next to (i,j,k) and
	 the site is outside the lattice for the non-
	 periodic commands a -1 will be returned. 
*/
int getIndFro( const_SNarray snA, int i, int j, int k);

int getIndBeh( const_SNarray snA, int i, int j, int k);

int getIndLef( const_SNarray snA, int i, int j, int k);

int getIndRig( const_SNarray snA, int i, int j, int k);

int getIndBel( const_SNarray snA, int i, int j, int k);

int getIndAbo( const_SNarray snA, int i, int j, int k);

int getIndFroP( const_SNarray snA, int i, int j, int k);

int getIndBehP( const_SNarray snA, int i, int j, int k);

int getIndLefP( const_SNarray snA, int i, int j, int k);

int getIndRigP( const_SNarray snA, int i, int j, int k);

int getIndBelP( const_SNarray snA, int i, int j, int k);

int getIndAboP( const_SNarray snA, int i, int j, int k);

/* Finds the index assuming periodic in the j direction (y-axis)
*/
int getIndexPeriodicX( const_SNarray snA, int i, int j, int k);

/* Finds the index assuming periodic in the j direction (y-axis)
*/
int getIndexPeriodicY( const_SNarray snA, int i, int j, int k);

/* Finds the index assuming periodic in the k direction (z-axis)
*/
int getIndexPeriodicZ( const_SNarray snA, int i, int j, int k);

/* Gets index if i j and k are all considered periodic
*/
int getIndexPeriodic( const_SNarray snA, int i, int j, int k);

int getAlen( const_SNarray snA);

int getAwid( const_SNarray snA);

int getAhei( const_SNarray snA);

int getAtotal( const_SNarray snA);

/* Calculates the total number of unoccupied sites along the yz-plane,
	 xz-plane and xy-plane;
   at index l (equivalent to some i,j, or k) along the x-axis,y-axis
	 or z-axis.
	 l can be between 0<=l<snA->length, 0<=l<snA->width, 0<=l<snA->height 
	 Will return a -1 otherwise or if snA is NULL
*/
int getUnOccYZplane(const_SNarray snA,const int l);
int getUnOccXZplane(const_SNarray snA,const int l);
int getUnOccXYplane(const_SNarray snA,const int l);

// If site in ahead of (i,j,k) in x direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccXpos(const_SNarray snA, int i, int j,int k);

// If site in behind (i,j,k) in x direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccXneg(const_SNarray snA,int i,int j,int  k);

// If site in ahead of (i,j,k) in Y direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccYpos(const_SNarray snA, int i, int j, int k);

// If site in behind (i,j,k) in y direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccYneg(const_SNarray snA, int i, int j, int k);

// If site in ahead of (i,j,k) in Y direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccYposPeriodic(const_SNarray snA,int i,int j,int k);

// If site in behind (i,j,k) in y direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccYnegPeriodic(const_SNarray snA,int i,int j,int k);

// If site in ahead of (i,j,k) in z direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccZpos(const_SNarray snA,int i,int j,int k);

// If site in behind (i,j,k) in z direction is unoccupied 
// unoccupied returns 0 Occupied returns 1 
int OccZneg(const_SNarray snA,int i, int j,int k);

/* Checks to see if there are any unoccupied sites around site 
	 (i,j,k) if there are unoccupied sites returns 0 if there are
	 no unoccupied sites returns 1
	 */
int OccAllNei(const_SNarray snA,int i,int j,int k);

int OccAllNeiPeriodicY(const_SNarray snA,int i,int j,int k);

/* Sets default values to the Site Node array
	 */
int setDefaultSNa(SNarray snA);

int printSNarray( const_SNarray snA);

int printSNarray_Detailed(const_SNarray snA);
#endif
