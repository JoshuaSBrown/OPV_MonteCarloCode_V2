
#include "../CHARGE/charge.h"
#include "../SITENODE/sitenode.h"

/* This function will update the charge path structure so that
 * the site the charge most recently hopped to is at the forefront
 * of the path linklist. If the site was already in the linklist
 * it will increment the number of the times the site was visited
 * as well.
 * If the site is part of the same cluster as other sites in the 
 * linklist all sites that are part of the cluster are incremented. 
 */
int updatePath(SNarray snA, Charge ch, int SiteID,int ClusterID);


