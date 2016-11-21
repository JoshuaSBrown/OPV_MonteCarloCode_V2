#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/types.h>
#include <math.h>
#include <assert.h>

#include "clustercalculators.h"
#include "../ERROR/error.h"
#include "../MIDPOINT/midpoint.h"
#include "../MONTECARLO/montecarlo.h"
#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"
#include "../NEIGHLL/neighll.h"
#include "../NEIGHNODE/neighnode.h"
#include "../LINKLIST/linklist.h"
#include "../CLUSTER/cluster.h"
#include "../ARBARRAY/arbarray.h"
#include "../PARAMETERS/read.h"
#include "../CHARGE/charge.h"
#include "../MATRIX_LINKLIST/matrix_linklist.h"
#include "../CLUSTERSITENODE/clustersitenode.h"
#include "../CHARGESITENODE/chargesitenode.h"

/* This functino calculates the number of hops a site
 * has that are within the cluster 
 * Output
 *   mtxHopOpt
 */
int CountOptions( Node tempNode, 
                  matrix * mtxHopOpt, 
                  const_SNarray snA){
    //Function counts the number of hopping options
    //a given site has within the cluster if it
    //wants to keep hopping within the cluster and
    //off the cluster
    //matrix mtxHopOpt stores the Id of the site
    //and the corresponding number of options
    //col 1 - hops within cluster
    //col 2 - hops off the cluster
    //col 3 - ID of site
	#ifdef _ERROR_CHECKING_ON_
    if ((mtxHopOpt)==NULL){ 
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxHopOpt is NULL\n");
		#endif
		return return_error_val();
    } 
    if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL\n");
		#endif
        return return_error_val();
    }
	#endif 

    int countOpt;
    int countOpt2;
    int inc = 1;
    int Node_ID;
    int i, j, k;
    while (tempNode!=NULL){

        countOpt=0;
        countOpt2=0;

        Node_ID = getNode_id(tempNode);
        getLoc( &i, &j, &k, Node_ID, snA);
        if(getFlagFro(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
        }
        if(getFlagBeh(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
        }
        if(getFlagLef(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
		}
        if(getFlagRig(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
        }
        if(getFlagBel(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
        }
        if(getFlagAbo(tempNode)==1){
            countOpt++;
        }else{
            countOpt2++;
        }

        setE((*mtxHopOpt),inc,1,countOpt);
        setE((*mtxHopOpt),inc,2,countOpt2);
        setE((*mtxHopOpt),inc,3,Node_ID);
        inc++;
        tempNode = getNextNode(tempNode);
    }

    return 0;
}


/* This function calculate mtxProb */
matrix CalculateProb(const_ClusterLL TempClLL, 
                     matrix mtxHopOpt, 
                     const_SNarray snA,
		             const int attempts ){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(mtxHopOpt==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxHopOpt is NULL in CalculateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in CalculateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(attempts<2){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR attempts is less than 2 in CalculateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif
	//This function can be used for periodic or non periodic conditions it does not
	//do anything with neighbors outside of the clusters.
	int inc;
	int Node_ID;
	int numNodes;
	int Node_IDFro, Node_IDBeh;
	int Node_IDLef, Node_IDRig;
	int Node_IDBel, Node_IDAbo;
	int i, j, k;
	int Row;
	int Row2;
	double val;

	//The second column of mtxProb will contain the IDs of the respective sites
	//We need to reinitialize the 2nd column to the ID's of the CNTs
	numNodes = getCluster_numNodes(TempClLL);
	printf("getCluster_numNodes(TempClLL) %d\n",numNodes);
	matrix mtxProb = newMatrixSet( getCluster_numNodes(TempClLL),2, (1/((double)get
					Cluster_numNodes(TempClLL))));
	matrix mtxProbNew = newMatrix(6,1);
	Node tempNode;

	//Initilize Node Ids in the mtxProb
	tempNode = getStartNode(TempClLL);

	printf("\n****************Calculating Prob Matrix********************\n");
	inc = 1;

	while(tempNode!=NULL){
		Node_ID = getNode_id(tempNode);
		setE(mtxProb,inc,2,(double) Node_ID);
		tempNode = getNextNode(tempNode);
		inc++;
	}
	for(int attempt=1;attempt<(attempts*numNodes);attempt++){
		tempNode= getStartNode(TempClLL);
		inc = 1;

		printf("attempt %d\n",attempt);

		while (tempNode!=NULL){

			//1 hop behind              index - 0
			//2 hop infront             index - 1
			//3 hop left                index - 2
			//4 hop right               index - 3
			//5 hop below               index - 4
			//6 hop above               index - 5
			Node_ID = getNode_id(tempNode);

			getLoc( &i, &j, &k, Node_ID, snA);

			if(getFlagFro(tempNode)==1){
				//Node in front is within the cluster
				Node_IDFro = getIndFroP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDFro, 3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDFro,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,2,1,val);
			}
			if(getFlagBeh(tempNode)==1){
				//Node behind is within the cluster
				Node_IDBeh = getIndBehP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDBeh,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDBeh,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,1,1,val);
			}

			if(getFlagLef(tempNode)==1){
				//Node behind is within the cluster
				Node_IDLef = getIndLefP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDLef,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDLef,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,3,1,val);
			}
			if(getFlagRig(tempNode)==1){
				//Node behind is within the cluster
				Node_IDRig = getIndRigP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDRig,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDRig,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,4,1,val);
			}
			if(getFlagBel(tempNode)==1){
				//Node behind is within the cluster
				Node_IDBel = getIndBelP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDBel,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDBel,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,5,1,val);
			}

			if(getFlagAbo(tempNode)==1){
				//Node behind is within the cluster
				Node_IDAbo = getIndAboP(snA,i,j,k);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDAbo,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDAbo,2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,6,1,val);
			}

			val = (SumOfCol(mtxProbNew,1)+getE(mtxProb,inc,1))/2;

			setE(mtxProb,inc,1,val);
			inc++;
			tempNode = getNextNode(tempNode);
			setAll(mtxProbNew, 0.0);
		}

		val = SumOfCol(mtxProb,1);
		//Normalize the matrix
		DivideEachElementCol(&mtxProb,1, val);

	}

	printMatrix(mtxProb);
	deleteMatrix(&mtxProbNew);

	return mtxProb;
}

/* This function is used to calculate the time at which 
 * a charge located on a site will jump off the site and
 * out of the cluster.
 * Consider the following scenario
 *
 *                 n3
 *                  |
 *            n2 - s3 - n4
 *             |    | 
 *       n1 - s1 - s2 - n5
 *             |    |
 *            n7   n6
 * 
 * To calculate the time a charge will take to jump from a site
 * We need
 *  o The rates at which it will hop from a site to its neighbors
 *
 * Tescapesite1 = 1/(RateHopToNeigh1from1+
 *                   RateHopToNeigh2from1+
 *                   RateHopToNeigh7from1)
 * 
 * For Tescapesite2
 *
 * Tescapesite2 = 1/(RateHopToNeigh5from2+
 *                   RateHopToNeigh6from2)
 *
 * Step1 Determine how many sites there are
 * Step2 Create the matrix
 * Step3 Calculate the sum of the rates off a site
 * Step4 Calculate the escape times 
 * 
 */
matrix CalculateTescape(ClusterLL * TempClLL, 
                        matrix mtxNeighRate){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateTescape\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(mtxNeighRate==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxNeighRate is NULL\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif
	/**********Step 1**********/
	int NumNodes =  getCluster_numNodes(*TempClLL);

	/**********Step 2**********/
	matrix mtxRate = newMatrixSet(NumNodes,2,-1.0);

	double Rate;
	int Node_ID;

	/**********Step 3**********/
	/* Cycle through the rates */
	for(int i=1;i<=getRows(mtxNeighRate);i++){
		
		/* Determine the Node_ID associated with the
         * rate */
		Node_ID = getE(mtxNeighRate,2);
		for(int j=1;j<=getRows(mtxRate);j++){
			
			/* Determine if the Node_ID has been stored */
			if(getE(mtxRate(j,2))==-1.0){
                /* If not stored store the node Id and the rate */
				Rate = getE(mtxNeighRate,i,1);
				setE(mtxRate,j,2,Node_ID);
				setE(mtxRate,j,1,Rate);
			}else if(getE(mtxRate(j,2))==Node_ID){
                /* If already stored then add the new rate to what
                 * has already been stored */
				Rate = getE(mtxRate,j,1);
				Rate += getE(mtxNeighRate,i,1);
				setE(mtxRate,j,1,Rate);
			}
		}

	}		

	matrix mtxTescape = mtxRate;
	/**********Step 4**********/
    /* Convert the sum of the rates to the escape time
     * for each site */
	for(int i=1;i<=getRows(mtxTescape);i++){
		Rate = getE(mtxTescape,i,1);
		setE(mtxTescape,i,1,1/Rate);
	}	

	return mtxTescape;
}

/* This function will calculate the dwell time of the cluster as a 
 * whole. The cluster dwell time should describe the rate at which
 * charges exit the cluster. It should based on the literature follow
 * an exponential behavior.
 *
 * EscapeRate = exp(-ClusterDwellTime*t ) 
 * 
 * Where t is the time and the ClusterDwellTime acts as a constant
 * This function simply calculates the constant and returns it.
 */
double CalculateClusterDwellTime(matrix mtxHopOffSiteProb,
                                 matrix mtxTescape){

	#ifdef _ERROR_CHECKING_ON_
	if(mtxHopOffSiteProb==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxHopOffSiteProb is NULL in "
                       "CalculateClusterDwellTime\n");
		#endif
		return (double) return_error_val();
	}
	if(mtxTescape==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxTescape is NULL in "
                       "CalculateClusterDwellTime\n");
		#endif
		return (double) return_error_val();
	}
	#endif

	double Tescape;
	double HopOffProb;
	double ClusterDwell;
	int Node_ID;

	ClusterDwell = 0.0;

	for(int i=1;i<=getRows(mtxHopOffSiteProb),i++){
		Node_ID = getE(mtxHopOffSiteProb,i,2);
		for(int j=1;j<=getRows(mtxTescape);j++){
			if(Node_ID==getE(mtxTescape,j,2)){
				HopOffProb = getE(mtxHopOffProb,i,1);
				Tescape   = getE(mtxTescape,j,1);
				ClusterDwell += HopOffProb*Tescape;
			}
		}
	}

	return ClusterDwell;
}

/* This function will take the probabilities for a cluster associated
 * with the likelyhood that a charge will hop to the neighboring site.
 * The probabilities will be assigned to the neighbor nodes from the
 * mtxNeighProb matrix. They will however follow the format such that
 * the probabilities are added. e.g. consider 3 neighboring sites with
 * the following probabilities:
 *
 * Neigh Site1 0.25
 * Neigh Site2 0.40
 * Neigh site3 0.35
 *
 * The pvals will add to 1
 *
 * Neigh Site1 0.25
 * Neigh Site2 0.65
 * Neigh Site3 1.00
 */
matrix CalculateNeighPval(ClusterLL * TempClLL,const_matrix mtxNeighProb){

	//Calculate Pval & time for NeighNodes 
	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"TempClLL is NULL in CalculateNeighPval\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(mtxNeighProb==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"mtxNeighProb is NULL in CalculateNeighPval\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(getRows(mtxNeighProb)!=getCluster_numNeigh(*TempClLL)){
		#ifdef _ERROR_
		fprintf(stderr,"mtxNeighProb does not have the same number "
                       "of rows as in CalculateNeighPval\n");
		#endif
		return_error_val()
		return NULL;
	}
	#endif

	/* Look through the Neighbors associated with TempClLL */
	NeighNode NeighNod = getStartNeigh(*TempClLL);
	/* pval will start at 0.0 and as we add the pval to 
     * each neighboring node it will be increased by the
     * probabilitiy */
	double Neighpval = 0.0;
	int Neigh_ID;
	int Neigh_ID2;

    /* Calculate the total */
	while(NeighNod!=NULL){
		Neigh_ID = getNeighNode_id(NeighNod);
		for(int i=1;i<=getRows(mtxNeighProb);i++){
			Neigh_ID2 = getE(mtxNeighProb,i,2);
			if(Neigh_ID==Neigh_ID2){
				Neighpval += getE(mtxNeighProb,i,1);
                setNeighNodeNew_p(NeighNod,pval);
			}
		}
		NeighNod = getNextNeigh(NeighNod);
	}

	/* Sanity check ensure that pval is 1.0 at the end */
	if(Neighpval>1.01 || Neighpval < 0.99){
		printf("ERROR in CalculateNeighPval pval is not unity at end.\n");
        exit(1);
	}
	return mtxNeighProb;
}
/* This function is used to calculate the probability
 * that a charge located within the cluster will jump 
 * from the site . 
 * Consider the following scenario
 *
 *                 n3
 *                  |
 *            n2 - s3 - n4
 *             |    | 
 *       n1 - s1 - s2 - n5
 *             |    |
 *            n7   n6
 * 
 * To calculate the probability a charge will jump from a site
 * We need
 *  o The probability a site is on site1 for a given iteration
 *  o The probability a site is on site1 for a given second
 *  o The rate at which it will hop from site1 to a neighbor
 *
 * Psite1 = ProbDwell1*ProbSite1*(RateHopToNeigh1from1+
                                  RateHopToNeigh2from1+
                                  RateHopToNeigh7from1)/AllRates
 * 
 * For Psite2
 *
 * Psite2 = ProbDwell2*ProbSite2*(RateHopToNeigh5from2+
 *                                RateHopToNeigh6from2)/AllRates
 *
 * Step1 Determine how many sites there are
 * Step2 Create the matrix
 * Step3 Calculate the probabilities 
 * 
 * Assumptions:
 */

matrix CalculateHopOffSiteProb(ClusterLL * TempClLL, 
                               SNarray snA,
                               matrix mtxDwellProb,
                               matrix mtxNeighRate,
                               matrix mtxProb){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateHopOffSiteProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in CalculateHopOffSiteProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(getRows(mtxProb)!=getRows(mtxDwellProb)){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxProb rows does not equal mtxDwellProb CalculateHopOffSiteProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif

	/**********Step 1**********/
	int NumberSite = getRows(mtxProb);
		
	/**********Step 2**********/
	matrix mtxHopOffSiteProb = newMatrixSet(NumberSite,2,-1.0);

	/* Increment through the matrix */
    int inc;
	int Node_ID;
	int Node_ID2;
	int flag;

	double Total            /* Store the total so we can normalize at   */
					        /* the end.                                 */
	double HopOffSiteProb;  /* What we are calculating                  */
	double DwellProb;       /* Dwell probability each node has it's own */
	double RateProb;        /* Each rate off the cluster divided by all */
                            /* the other rates                          */
	double RateProbPerSite; /* Sum of the rates for a given site        */
	double SiteProb;        /* Prob a charge will be on a site each node*/
                            /* has it's own                             */

	Node tempNode = getStartNode(*TempClLL);
    inc           = 1;
	Total         = 0.0;
	/********Step 3*********/
	while(tempNode){
        Node_ID  = getNode_id(tempNode);
		sn       = getSNwithInd(snA,Node_ID);
	
		/* mtxProb has the same number of rows as mtxDwell we will 
         * cycle through them and find the appropriate DwellRateProb
         * and SiteProb for Node_ID */
		/* flag is used to break out of the loop early if both values
         * are found */
        flag = 0;
		for(int i=1;i<=getRows(mtxProb);i++){
			if(Node_ID==getE(mtxProb,i,2)){
				SiteProb = getE(mtxProb,i,1);
				flag++;
				if(flag==2){
					break;
				}
			}
			if(Node_ID==getE(mtxDwellProb,i,2)){
				DwellProb = getE(mtxDwellProb,i,1);
				flag++;
				if(flag==2){
					break;
				}
			}
		}

		RateProbPerSite = 0;
		/* Now we will cycle through each of the rates and calculate 
         * the probability of jumping off of a site */
		for(int i=1;i<=getRows(mtxNeighRate);i++){
			
			RateProb = getE(mtxNeighRate,i,1);
			Node_ID2 = getE(mtxNeighRate,i,2);
			
			/* Here is a rate that describes transfer off of Node_ID */
			if(Node_ID==Node_ID2){
				RateProbPerSite += RateProb;
			}		
		}

		/* Finally calculating the probability the site will jump
         * off of Node_ID to a neighboring site */
		HopOffSiteProb = RateProbPerSite*SiteProb*DwellProb;	
		/* Adding the total probabilities together so we can 
         * normalize at the end */
		Total += HopOffSiteProb;	
		setE(mtxHopOffSiteProb,inc,2,Node_ID);
		setE(mtxHopOffSiteProb,inc,1,HopOffSiteProb);
		inc++;

		tempNode = getNextNode(tempNode);
	}	

	/* Now we need to normalize the mtxNeighProb */
	DivideEachElementCol(mtxHopOffSiteProb,1,Total);

	return mtxHopOffSiteProb;
}

/* This function is used to calculate the probability
 * that a charge located within the cluster will jump 
 * to a neighboring site. 
 * Consider the following scenario
 *
 *                 n3
 *                  |
 *            n2 - s3 - n4
 *             |    | 
 *       n1 - s1 - s2 - n5
 *             |    |
 *            n7   n6
 * 
 * To calculate the probability a charge will jump to neighbor
 * 1 we need
 *  o The probability a site is on site1 for a given iteration
 *  o The probability a site is on site1 for a given second
 *  o The rate at which it will hop to neigh 1 compared with all other neighbors
 *
 * Pnei1 = ProbDwell1*ProbSite1*RateHopToNeigh1from1/AllRates
 * 
 * For Pnei2
 *
 * Pnei2 = ProbDwell1*ProbSite1*RateHopToNeigh2from1/AllRates +
 *         ProbDwell2*ProbSite2*RateHopToNeigh2from3/AllRates
 *
 * Step1 Determine how many neighbors there are
 * Step2 Create the matrix
 * Step3 Calculate the probabilities 
 * 
 * Assumptions:
 * This function assumes that the neighnodes have been correctly
 * added to the cluster
 */
matrix CalculateNeighProb(ClusterLL * TempClLL, 
                          SNarray snA,
                          matrix mtxDwellProb,
                          matrix mtxNeighRate,
                          matrix mtxProb){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateNeighProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in CalculateNeighProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(getRows(mtxProb)!=getRows(mtxDwellProb)){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxProb rows does not equal mtxDwellProb CalculateNeighProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif
	
	/********Step 1*********/
	int totalNumberNeigh = getCluster_numNeigh(*TempClLL);

	/********Step 2*********/
    /* The first column stores the probabilities 
     * The second column stores the id of the nieghnode 
     * will assign -1.0 initially as a sanity check */
	matrix mtxNeighProv = newMatrixSet(totaNumberNeigh,2,-1.0);

	/* Increment through the matrix */
    int inc;
	int Node_ID;
	int Neigh_ID;
	int flag;

	double Total       /* Store the total so we can normalize at   */
					   /* the end.                                 */
	double NeighProb;  /* What we are calculating                  */
	double DwellProb;  /* Dwell probability each node has it's own */
	double RateProb;   /* Each rate off the cluster divided by all */
                       /* the other rates                          */
	double SiteProb;   /* Prob a charge will be on a site each node*/
                       /* has it's own                             */

	Node tempNode = getStartNode(*TempClLL);
	total         = 0.0;
	/********Step 3*********/
	while(tempNode){
        Node_ID  = getNode_id(tempNode);
		sn       = getSNwithInd(snA,Node_ID);
	
		/* mtxProb has the same number of rows as mtxDwell we will 
         * cycle through them and find the appropriate DwellRateProb
         * and SiteProb for Node_ID */
		/* flag is used to break out of the loop early if both values
         * are found */
        flag = 0;
		for(int i=1;i<=getRows(mtxProb);i++){
			if(Node_ID==getE(mtxProb,i,2)){
				SiteProb = getE(mtxProb,i,1);
				flag++;
				if(flag==2){
					break;
				}
			}
			if(Node_ID==getE(mtxDwellProb,i,2)){
				DwellProb = getE(mtxDwellProb,i,1);
				flag++;
				if(flag==2){
					break;
				}
			}
		}

		/* Now we will cycle through each of the rates and calculate 
         * the probability of jumping to a neighboring site */
		for(int i=1;i<=getRows(mtxNeighRate);i++){
			
			RateProb = getE(mtxNeighRate,i,1);
			Neigh_ID = getE(mtxNeighRate,i,3);
			
			/* Cycle through mtxNeighProb and ensure we have not
             * already accounted for hops to this neighbor */
			for(int j=1;j<=getRows(mtxNeighProb)){
				if(getE(mtxNeighProb,j,2)==-1){
                    /* Neighbor has not been accounted for */
					NeighProb = 0.0;
					setE(mtxNeighProb,j,2,Neigh_ID);
                    inc = j;
					break;
				}else if(getE(mtxNeighProb,j,2)==Neigh_ID){
                    /* Neighbor has been stored */
                    NeighProb = getE(mtxNeighProb,j,1);
                    inc = j;
					break;
				}		
			}
			
			/* Updating NeighProb and total */
            total     += RateProb*SiteProb*DwellProb;
			NeighProb += RateProb*SiteProb*DwellProb;
			setE(mtxNeighProb,j,1);
		}
		
		tempNode = getNextNode(tempNode);
	}	

	/* Now we need to normalize the mtxNeighProb */
	DivideEachElementCol(mtxNeighProb,1,total);

	return mtxNeighProb;
}

/* This function is designed to calculate the probability that 
 * a charge will dwell on a site within the cluster. I.e
 * consider a cluster with three sites each site has it's own
 * dwell time. For site 1
 * 
 * Tdwell1 = 1/(Rate1_1+Rate1_2+Rate1_3+Rate1_4+Rate1_5+Rate1_6)
 * where there are six rates because there are six neighbors for
 * each site
 *
 * Tdwell2 = 1/(Rate2_1+Rate2_2+Rate2_3+Rate2_4+Rate2_5+Rate2_6)
 * Tdwell3 = 1/(Rate3_1+Rate3_2+Rate3_3+Rate3_4+Rate3_5+Rate3_6)
 *
 * To calculate the dwell prob for site 1 
 *
 * TdwellProb1 = Tdwell1/(Tdwell1+Tdwell2+Tdwell3)
 * 
 * The same is true for sites 2 and 3
 * 
 * Step1 Calculate the sum of the dwell times
 * Step2 Create matrix to store probabilities
 * Step3 Store probabilities in the matrix
 */ 
matrix CalculateNodeDwellProb(ClusterLL * TempClLL, SNarray snA ){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateNodeDwellProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in CalculateNodeDwellProb\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif
	
	Node tempNode = getStartNode(*TempClLL);
	int Node_ID;
	double totalDwell = 0.0;
	int totalNumNodes = getCluster_numNodes(*TempClLL);
	/* The site nodes and Node IDs stored in the 
     * Clusters should be the IDs associated with 
     * each of the sites in the lattice. */
    SiteNode sn;

	/********Step 1*********/
	while(tempNode){
	
        Node_ID  = getNode_id(tempNode);
		sn       = getSNwithInd(snA,Node_ID);
		/* Here sum is equivalent to the sum of the rates
         * (Rate1+Rate2+...Rate6)
         * Thus we can get the sum of the dwell times by
         * 
         * totalDwell = 1/(Rate1_1+Rate1_2+...Rate1_6) + 1/(Rate2_1+...) + ...
         */
		totalDwell += 1/getsum(sn);	
		/* Now that we have the sn we need to calculate
         * the total dwell time */	
		tempNode = getNextNode(tempNode);
	}

	/********Step 2*********/
	/* Now we are free to create the matrix */
	matrix mtxProbDwellProb = newMatrix(totaNumNodes,2);

	if(mtxProbDwellProb==NULL){
		fprintf(stderr,"ERROR mtxProbDwellProb is NULL in "
                       "CalculateNeighDwellProb\n");
		exit(1);
	}	
	/********Step 3*********/
	tempNode = getStartNode(*TempClLL);
	Node_ID = getNode_id(tempNode);
    sn = getSNwithInd(snA,Node_ID);

	/* Used to increment the matrix */
	int    inc = 1;
	double dwell;

	/********Step 1*********/
	while(tempNode){
	
		dwell = 1/getsum(sn);
		setE(mtxProbDwellProb,inc,1,dwell/totalDwell);
		setE(mtxProbDwellProb,inc,2,Node_ID);
        inc++;

		tempNode = getNextNode(tempNode);
        Node_ID  = getNode_id(tempNode);
		sn       = getSNwithInd(snA,Node_ID);
	}

	return mtxProbDwellProb;
}

/* This function will calculate the rate probability of a cluster
 * for each possible hop off the cluster i.e.
 * If if I have cluster of two sites there are a total of 10 hops
 * off the cluster each the first element in the rateProb matrix
 * will contain the rate of site 1 to neighboring site 1 dividted
 * by the sum of all the rates off of the cluster
 * rateProb(1,1) = Rate1/(Rate1+Rate2+...Rate10);
 * rateProb(2,1) = Rate2/(Rate1+Rate2+...Rate10);
 *
 * Step 1 Calculate the number of Hops 
 *        Calculate the total of the rates 
 *        totalRates = (Rate1+Rate2+....Rate10)
 * Step 2 Create the matrix to store the probabilities
 * Step 3 Store Probabilties in matrix 
 */
matrix CalculateNeighRateProb(ClusterLL * TempClLL, 
                              matrix mtxProb,
                              matrix MasterM,
                              SNarray snA,
                              int PeriodicX,
                              int PeriodicY,
                              int PeriodicZ){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculateNeighRateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(mtxProb==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxProb is NULL in CalculateNeighRateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in CalculateNeighRateProb.\n");
		#endif
		return_error_val();
		return NULL;
	}
	#endif

	//1 - Order of Magnitude of hop behind site i,j,k
	//2 - Order of Magnitude of hop infornt of site i,j,k
	//3 - Order of Magnitude of hop to site left of site i,j,k
	//4 - Order of Magnitude of hop to site right of site i,j,k
	//5 - Order of Magnitude of hop to site below site i,j,k
	//6 - Order of Magnitude of hop to site above site i,j,k
	//7 - Hop rate to site behind i,j,k
	//8 - Hop rate to site infront of i,j,k
	//9 - Hop rate to site left of site i,j,k
	//10 - Hop rate to site right of i,j,k
	//11 - Hop rate to site below i,j,k
	//12 - Hop rate to site above i,j,k

	/* Grab the first node in the cluster */
	Node tempNode = getStartNode(*TempClLL);

	/* Need to keep track of the total of the Rates off
     * of the cluster we will use this to normalize the
     * mtxNeighRateProb matrix */
	double totalRate = 0.0;

	/* Keep track of the total number of hops off the 
     * cluster */
    int numHopsOffCluster = 0;

	/******* Step 1 ********/
	while(tempNode!=NULL){
 	
		Node_ID = getNode_id(tempNode);	
		getLoc( &i, &j, &k, Node_ID, snA);

		/* Check to see if the neighbor node infront of the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexFront will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexFront(snA, Node_ID)==1 || PeriodicX!=0 ){
			if(getFlagFro(tempNode)!=1){
			/* getFlagFro determines if the site infront of Node_ID is
             * part of the cluster or if it is a neighboring site
             * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,8);
                totalRate += tempVal;
				numHopsOffCluster++;
			}
            /* If the if statement was not triggered it means the site
             * the neighboring site was within the cluster */

		}else{
		/* If the site in front is not part of the system and the system
         * is not periodic it must mean that we have reached an electrode */

			if(getFlagFro(tempNode)!=1){
		
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}

		/* Check to see if the neighbor node behind the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexBehind will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexBehind(snA, Node_ID)==1 || PeriodicX!=0 ){
			if(getFlagBeh(tempNode)!=1){
				/* getFlagBeh determines if the site is behind Node_ID is
				 * part of the cluster or if it is a neighboring site
				 * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,7);
				totalRate += tempVal;
				numHopsOffCluster++;
			}
			/* If the if statement was not triggered it means the site
			 * the neighboring site was within the cluster */

		}else{
			/* If the site behind is not part of the system and the system
			 * is not periodic it must mean that we have reached an electrode */
			if(getFlagBeh(tempNode)!=1){

				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
	
		/* Check to see if the neighbor node to the left of the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexleft will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexLeft(snA, Node_ID)==1 || PeriodicY!=0 ){
			if(getFlagLef(tempNode)!=1){
				/* getFlagLef determines if the site to the left of Node_ID is
				 * part of the cluster or if it is a neighboring site
				 * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,9);
				totalRate += tempVal;
				numHopsOffCluster++;
			}
			/* If the if statement was not triggered it means the site
			 * the neighboring site was within the cluster */

		}else{
			/* If the site to the left is not part of the system and the system
			 * is not periodic it must mean that we have reached an electrode */
			if(getFlagLef(tempNode)!=1){

				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}

		/* Check to see if the neighbor node to the right of the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexRight will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexRight(snA, Node_ID)==1 || PeriodicY!=0 ){
			if(getFlagRig(tempNode)!=1){
				/* getFlagRig determines if the site to the right of Node_ID is
				 * part of the cluster or if it is a neighboring site
				 * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,10);
				totalRate += tempVal;
				numHopsOffCluster++;
			}
			/* If the if statement was not triggered it means the site
			 * the neighboring site was within the cluster */

		}else{
			/* If the site to the right is not part of the system and the system
			 * is not periodic it must mean that we have reached an electrode */
			if(getFlagRig(tempNode)!=1){

				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}

		/* Check to see if the neighbor node above the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexAbove will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexAbove(snA, Node_ID)==1 || PeriodicZ!=0 ){
			if(getFlagAbo(tempNode)!=1){
				/* getFlagAbo determines if the site above Node_ID is
				 * part of the cluster or if it is a neighboring site
				 * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,12);
				totalRate += tempVal;
				numHopsOffCluster++;
			}
			/* If the if statement was not triggered it means the site
			 * the neighboring site was within the cluster */

		}else{
			/* If the site above is not part of the system and the system
			 * is not periodic it must mean that we have reached an electrode */
			if(getFlagAbo(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}

		/* Check to see if the neighbor node below the 
         * current Node_ID is within the bounds of the system
         * If it is checkBoundsIndexBelow will return 1
         * If it is not but the system is periodic then we will
         * continue anyway */
		if(checkBoundsIndexBelow(snA, Node_ID)==1 || PeriodicZ!=0 ){
			if(getFlagBel(tempNode)!=1){
				/* getFlagBel determines if the site above Node_ID is
				 * part of the cluster or if it is a neighboring site
				 * if it is not part of the cluster it will return 0 */
				tempVal = getE(MasterM, Node_ID+1,11);
				totalRate += tempVal;
				numHopsOffCluster++;
			}
			/* If the if statement was not triggered it means the site
			 * the neighboring site was within the cluster */

		}else{
			/* If the site below is not part of the system and the system
			 * is not periodic it must mean that we have reached an electrode */
			if(getFlagBel(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		tempNode = getNextNode(tempNode);
	}

	/******* Step 2 ********/

	/* Now we know how many hops off the cluster there are as well
     * as the sum of the rates. We are now prepared to make the matrix
     * mtxNeighRateProb that will store the probability that a given
     * rate will be used */
    /* The first column stores the rate prob 
     * column two stores the id of the site the charge is hopping from
     * column three stores the id of the neighboring site the charge
     * is hopping too */
	matrix mtxNeighRateProb = newMatrix(numHopsOffCluster,3);
	
	if(mtxNeighRateProb==NULL){
		fprintf(stderr,"ERROR mtxNeighRateProb is NULL in "
                       "CalculateNeighRateProb\n");
		exit(1);
	}

    /* Once again we grab the first node in the cluster */
	Node tempNode = getStartNode(*TempClLL);

	/* We call the same function and if statements listed in the 
     * former while loop but we are now calculating the probabilities */
	int HopNum = 1;

	/******* Step 3 ********/
	while(tempNode!=NULL){
		
		Node_ID = getNode_id(tempNode);	
		getLoc( &i, &j, &k, Node_ID, snA);
		if(checkBoundsIndexFront(snA, Node_ID)==1 || PeriodicX!=0 ){
			if(getFlagFro(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,8);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i+1,j,k));
				HopNum++;
			}
		}else{
			if(getFlagFro(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		if(checkBoundsIndexBehind(snA, Node_ID)==1 || PeriodicX!=0 ){
			if(getFlagBeh(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,7);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i-1,j,k));
				HopNum++;
			}
		}else{
			if(getFlagBeh(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		if(checkBoundsIndexLeft(snA, Node_ID)==1 || PeriodicY!=0 ){
			if(getFlagLef(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,9);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i,j-1,k));
				HopNum++;
			}
		}else{
			if(getFlagLef(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		if(checkBoundsIndexRight(snA, Node_ID)==1 || PeriodicY!=0 ){
			if(getFlagRig(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,10);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i,j+1,k));
				HopNum++;
			}
		}else{
			if(getFlagRig(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		if(checkBoundsIndexAbove(snA, Node_ID)==1 || PeriodicZ!=0 ){
			if(getFlagAbo(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,12);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i,j,k+1));
				HopNum++;
			}
		}else{
			if(getFlagAbo(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		if(checkBoundsIndexBelow(snA, Node_ID)==1 || PeriodicZ!=0 ){
			if(getFlagBel(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,11);
				setE(mtxNeighRateProb,HopNum,1,tempVal/totalRate);
				setE(mtxNeighRateProb,HopNum,2,Node_ID);
				setE(mtxNeighRateProb,HopNum,3,getIndex(snA,i,j,k-1));
				HopNum++;
			}
		}else{
			if(getFlagBel(tempNode)!=1){
				printf("ERROR in calculate mtxNeigRateProb\n");
				exit(1);
			}
		}
		tempNode = getNextNode(tempNode);
	}
	return mtxNeighRateProb;
}

int CalculatePvalNodes(ClusterLL * TempClLL, matrix mtxProb, matrix mtxDwellTime){

	#ifdef _ERROR_CHECKING_ON_
	if(TempClLL==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR TempClLL is NULL in CalculatePvalNodes\n");
		#endif
		return return_error_val();
	}
	if(mtxProb==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxProb is NULL in CalculatePvalNodes\n");
		#endif
		return return_error_val();
	}
	if(mtxDwellTime==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR mtxDwellTime is NULL in CalculatePvalNodes\n");
		#endif
		return return_error_val();
	}
	#endif

	//Calculate Pval & time for Nodes

	//printf("Number of Nodes %d\n",getCluster_numNodes(*TempClLL));
	Node tempNode = getStartNode((*TempClLL));
	int inc=1;
	double val=0.0;
	double total=0;

	while(tempNode!=NULL){

		val = getE(mtxProb,inc,1)*getE(mtxDwellTime,inc,1);
		//printf("Value of val %g\n",val);
		total += val;
		setNode_p(tempNode,val);		
		inc++;
		tempNode = getNextNode(tempNode);
	}

	tempNode=getStartNode((*TempClLL));

	val = 0.0;
	inc = 1;
	//Normalize
	while(tempNode!=NULL){

		printf("getNode_p %g\n",getNode_p(tempNode));
		val += getNode_p(tempNode)/total;
		setNode_p(tempNode,val);	
		tempNode = getNextNode(tempNode);
	}

	return 0;
}

