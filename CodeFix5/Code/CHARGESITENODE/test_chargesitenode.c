#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "chargesitenode.h"
#include "../CHARGE/charge.h"
#include "../SITENODE/sitenode.h"
#include "../ERROR/error.h"

int main(void){

  // Used to collect return value
  int rv;
  double rvd;

  // Used to store the site id
  int SiteID;
  int SiteID2;
  int SiteID3;
  int SiteID4;
    
  // Used to store a charge
  Charge ch;

  // Used to specify cluster id
  int ClusterID;

  // specifiying charge array will 
  // have a length of 5
  int len = 5;

  // Specifying the number of nodes
  // that each charge path will 
  // contain
  int numNodes = 4;

  // Creating a charge array
  ChargeArray chA = newChargeA(len);

  // Creating a sn array for the 
  // charges to move around in a 
  // 5x5x5 system
  SNarray snA = newSNarray(5,5,5);

  printf("Testing: updatePath\n");
  // Grabbing first charge in the 
  // array
  ch = getCharge(chA, 0);
  
  // Chosing SiteID 
  SiteID = getIndex(snA, 0, 0, 0);

  // Returns -1 because Path has not
  // yet been initialized
  rv = updatePath(snA,ch, SiteID, -1);
  assert(rv==-1);
  
  // Initializeing the Charge Path
  initChargeArrayPath(chA,numNodes);

  rv = updatePath(NULL,ch, SiteID, -1);
  assert(rv==-1);
  rv = updatePath(snA,NULL,SiteID,-1);
  assert(rv==-1);
  rv = updatePath(snA,ch, -1, -1);
  assert(rv==-1);

  SiteID = getIndex(snA,4,4,4)+1;
  rv = updatePath(snA,ch, SiteID, -1);
  assert(rv==-1);

  SiteID = getIndex(snA,0,0,0);
  rv = updatePath(snA,ch, SiteID, -1);
  assert(rv==0);
  rvd = getChargePathVisits(ch,1);
  assert(rvd==1);
  rvd = getChargePathVisits(ch,2);
  assert(rvd==0);
  rvd = getChargePathVisits(ch,3);
  assert(rvd==0);
  rvd = getChargePathVisits(ch,3);
  assert(rvd==0);
  rvd = getChargePathVisitsForSite(ch, SiteID);
  assert(rvd==1.0); 
   
  SiteID2 = getIndex(snA,1,0,0);
  rv = updatePath(snA,ch, SiteID2, -1);
  assert(rv==0);
  rv = updatePath(snA,ch, SiteID, -1);
  assert(rv==0);
  rvd = getChargePathVisits(ch,1);
  assert(rvd==2);
  rvd = getChargePathVisits(ch,2);
  assert(rvd==1);
  rvd = getChargePathVisits(ch,3);
  assert(rvd==0);
  rvd = getChargePathVisits(ch,4);
  assert(rvd==0);
  rvd = getChargePathVisitsForSite(ch, SiteID);
  assert(rvd==2.0); 
  rvd = getChargePathVisitsForSite(ch, SiteID2);
  assert(rvd==1.0); 

  // Now we will test what happens if 
  // we assume a charge is hopping to 
  // a cluster
  SiteID3 = getIndex(snA,2,0,0);
  SiteID4 = getIndex(snA,3,0,0);
  ClusterID = 1;
  rv = updatePath(snA,ch, SiteID3, ClusterID);
  assert(rv==0);
  rv = updatePath(snA,ch, SiteID4, ClusterID);
  assert(rv==0);
  // Should be 2 because SiteID3 and
  // SiteID4 are on the same cluser
  rvd = getChargePathVisits(ch,1);
  assert(rvd==2);
  rvd = getChargePathVisits(ch,2);
  assert(rvd==2);
  rvd = getChargePathVisits(ch,3);
  assert(rvd==2);
  rvd = getChargePathVisits(ch,4);
  assert(rvd==1);
  
  /* Deleting ChargeArray */
  deleteChargeA(chA);
  /* Delete SNarray */
  deleteSNarray(&snA);
  printf("Testing Completed\n");
  return 0;
}
