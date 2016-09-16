#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "cluster.h"
#include "../ELECTRODE/electrode.h"
#include "../NEIGHLL/neighll.h"
#include "../NEIGHNODE/neighnode.h"
#include "../NODE/node.h"
#include "../MIDPOINT/midpoint.h"

int main() {
	int rv;

	//////////////////////////////////////////////////////////////////////
  printf("Testing: addNodesToClusterGivenSites\n");
  ClusterLL ClLL = newClusterLL(1);
  assert(ClLL!=NULL);
  rv = addNodesToClusterGivenSites(ClLL, -1, 0);
  assert(rv==-1); 
  rv = addNodesToClusterGivenSites(ClLL, 0, -1);
  assert(rv==-1);
  rv = addNodesToClusterGivenSites(ClLL, 0, 0);
  assert(rv==-1);
  rv = addNodesToClusterGivenSites(NULL, 1, 0);
  assert(rv==-1);
  rv = addNodesToClusterGivenSites(ClLL, 1, 0);
  assert(rv==0);

  printf("Testing: appendClusterLL\n");
  ClusterLL ClLL2 = newClusterLL(2);
  assert(ClLL2!=NULL);
  rv = addNodesToClusterGivenSites(ClLL2, 3,4);
  assert(rv==0);

  rv = appendClusterLL(NULL,ClLL2);
  assert(rv==-1);
  rv = appendClusterLL(ClLL,NULL);
  assert(rv==-1);
  rv = appendClusterLL(ClLL,ClLL2);
  assert(rv==0);

  printClusterLL(ClLL);
  deleteAllClusterLL(&ClLL); 
  return 0;
}
