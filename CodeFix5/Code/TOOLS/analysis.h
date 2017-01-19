#ifndef _ANALYSIS_H_
#define _ANALYAIA_H_

#include "../MATRIX/matrix.h"
#include "../LINKLIST2/linklist2.h"

// print usage message
void usage(void);

//Read the .path file store the path of each charge in an array
int readPath(char * FileName);

//Gitven the pathway of a single charge find the percolation path
//and the trap sites. They are printed in .perc and .trap files
//respectively. 
int getPercTrap(char * FileName, const_matrix ChargePath,\
                int ChargeID, int y_max, int z_max);

int uniqueID(int x, int y, int z, int wid, int hei);

int ConvertMatrixLinkList(const_matrix ChargePath, linklist2 * ChargePathLL,\
             int wid, int hei);

#endif
