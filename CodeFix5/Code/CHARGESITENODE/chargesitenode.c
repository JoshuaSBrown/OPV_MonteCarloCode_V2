#include <stdio.h>
#include <stdlib.h>

#include "../CHARGE/charge.h"
#include "../ERROR/error.h"
#include "../SITENODE/sitenode.h"
#include "../MATRIX_LINKLIST/matrix_linklist.h"
#include "../CLUSTER/cluster.h"

int updatePath(SNarray snA, Charge ch, int SiteID, int ClusterID){
  #ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR can not update charge path ");
    fprintf(stderr,"charge is NULL.\n");
    #endif
    return return_error_val();
  }
  if(SiteID<0){
  #ifdef _ERROR_
    fprintf(stderr,"ERROR SiteID is less than 0 cannot updatePath\n");
    #endif
    return return_error_val();
  }
  if(SiteID>=getAtotal(snA)){
  #ifdef _ERROR_
    fprintf(stderr,"ERROR SiteID is greater than the total \n");
    fprintf(stderr,"number of nodes in the snA in  updatePath\n");
    #endif
    return return_error_val();
  }
  #endif

  printf("Entering Update Path\n");
  //Step 1 Determine if site is already recorded path list 
  double siteID = (double) SiteID;
  int position;
  int inc;
  int SiteID2;
  int ClusterID2;
  int numMatches;
  double visits;
  matrix_linklist mtxll;
  //The first row in the path matrix is reserved for the id
  //of the site
  mtxll = getChargePath(ch);

  #ifdef _ERROR_CHECKING_ON_
  if(mtxll == NULL){
  #ifdef _ERROR_
    fprintf(stderr,"ERROR mtxll returned NULL in updatePath\n");
    #endif
    return return_error_val();
  }
  #endif

  position = getMatrixLLLastMatchAtRow(mtxll,siteID,1);

  if(ClusterID==-1){
    //This means that the site of interest is not connected 
    //to a cluster
    if(position==-1){
      double path[3];
      path[0] = siteID;
      path[1] = 1;
      path[2] = -1;

      addLL_MNodeBegin(&mtxll,path,3);
      removeLL_MNodeEnd(&mtxll);
    }else{
      moveMatrixLLNodeToStart(mtxll,position);
      //Increment the number of times the site has been visited
      incMatrixLLElem(mtxll,1,2);
    }
  }else{

    printf("Cluster Exists\n");
    //This means that the site is part of a cluster
    //Returns a matrix with a list of where in mtxll sequence 
    //the matches occur 
    matrix matches = getMatrixLLMatchAtRow(mtxll,ClusterID,3, &numMatches);

    //This means we need to update the sites listed
    //in the ChargePath so they show which clusters
    //they are associated with
    if(numMatches==0){
      //Cycle through the path
      for(inc = 1; inc<=getMatrixLLlength(mtxll);inc++){
        SiteID2 = getMatrixLLNodeElem(mtxll,inc,1);
        //Check if part of a cluster
        if (checkSNconnectedCluster(getSNwithInd(snA,SiteID2))==1){
          //Determine the Cluster ID
          ClusterID2 = getCluster_id(getClusterList(getSNwithInd(snA,SiteID)));
          //Update the Elem in the path to 
          //show it is attached to the ClusterID2
          setMatrixLLElem(mtxll,inc,3,ClusterID2);
        }
      }
    }

    printCharge(ch);

    //If some of the sites the charge has been hopping
    //to belong to the same cluster we will increment
    //these sites
    for(inc=1;inc<=numMatches;inc++){
      incMatrixLLElem(mtxll,getE(matches,inc,1),2);
    }
    //If position returns -1 there is no match for the site
    //this mean that the charge has not recently hopped to
    //this particular site even if it is in the same cluster
    if(position==-1){
      //The sites are unique but more than one site
      //is part of the same cluster, here the number
      //of times the cluster has been visited is used
      //to initialize the new site that is being placed
      //at the front of the mtxll
      if(numMatches>0){
        printf("matches(1,1) %g\n",getE(matches,1,1));
        visits = getMatrixLLNodeElem(mtxll,(int)getE(matches,1,1),2);
      }else{
        visits = 1;
      }
      double path[3];
      path[0] = siteID;
      path[1] = visits;
      path[2] = ClusterID;
      mtxll = mtxll;
      addLL_MNodeBegin(&mtxll,path,3);
      removeLL_MNodeEnd(&mtxll);
    }else{
      moveMatrixLLNodeToStart(mtxll,position);
    }
    if(numMatches>0){
      deleteMatrix(&matches);
    }
  }
  printf("Leaving UpdatePath\n");
  printCharge(ch);
  return 0;
}



