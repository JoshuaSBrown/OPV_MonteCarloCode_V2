#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

#include "../PARAMETERS/read.h"
#include "chargetransport.h"
#include "../CHARGE/charge.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.h"
#include "../CLUSTER/CLUSTERSITENODE/clustersitenode.h"
#include "../FUNCTIONS/functions.h"
#include "../MEM/mem.h"

int HopToElecX(SNarray snA, Electrode elXb, Charge * one, int * future, int EndY, int EndZ){
	//Elecid defines which electrode could be hopping too
	//where:
	//	getAtotal(snA) + 0 - (-x)
	//	getAtotal(snA) + 1 - (+x)
	//	getAtotal(snA) + 2 - (-y)
	// 	getAtotal(snA) + 3 - (+y)
	//	getAtotal(snA) + 4 - (-z)
	//	getAtotal(snA) + 5 - (+z)

	printf("Hop To ElecX\n");
	if(elXb==NULL){
		printf("elXb is NULL\n");
		return -1;
	}
 	if( snA==NULL){
		printf("snA is NULL\n");
		return -1;
	}
	if(*one == NULL ){
		printf("Charge is NULL\n");
		return -1;
	}	
	if( EndY<0 || EndZ<0){
		return -1;
	}

	SNarray snAelec = (SNarray) getElectrode_AdjacentSites(elXb);
	printf("Length %d Width %d Height %d\n",getAlen(snAelec),getAwid(snAelec),getAhei(snAelec));


	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	int xx = getCx(*one);
	int yy = getCy(*one);
	int zz = getCz(*one);

	int factorX;
	int factorY;
	int factorZ;

	int l;
	int i;
	int j;
	int k;

	double position;

	if(xx<0){
		factorX = -xx/SLength+1;
		i = (xx+SLength*factorX)%SLength;
	}else{
		i = (xx+SLength)%SLength;
	}

	if(yy<0){
		factorY = -yy/SWidth+1;
		j = (yy+SWidth*factorY)%SWidth;
	}else{
		j = (yy+SWidth)%SWidth;
	}

	if(zz<0){
		factorZ = -zz/SHeight+1;
		k = (zz+SHeight*factorZ)%SHeight;
	}else{
		k = (zz+SHeight)%SHeight;
	}

	//   codeY = 1 exclude left hop
	//				 = 2 exclude right hop
	//   codeZ = 1 exclude bottom hop
	//				 = 2 exclude top hop


	printf("Determining code\n");
	int codeX = 0;
	int codeY = 0;
	int codeZ = 0;

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;

	if(xx==0){
		//Posibility to jump to (-x) electrode

		//Have to ensure that charge is not at the edge of the system in the y and z directon
		if(yy == 0 && zz == 0){
			//Next to the bottom left corner
			codeY = 1;
			codeZ = 1;
		}else if(yy == 0 && zz == EndZ*SHeight-1){
			//Next to the left and top corner
			codeY = 1;
			codeZ = 2;
		}else if(zz == 0 && yy == EndY*SWidth-1){
			//Next to the right and bottom corner
			codeY = 2;
			codeZ = 1;
		}else if(zz == EndZ*SHeight-1 && EndY*SWidth-1){
			//Next to the right and top corner
			codeY = 2;
			codeZ = 2;
		}else if(yy == 0){
			//Next to the left side
			codeY = 1;
		}else if(zz == 0){
			//Next to the bottom edge
			codeZ = 1;
		}else if(zz == EndZ*SHeight-1){
			codeZ = 2;
		}else if(yy == EndZ*SWidth-1){ 
			codeY = 2;
		}

	}

	//l value
	//0 - backward 	x-
	//1 - forward 	x+
	//2 - left 			y-
	//3 - right 		y+
	//4 - below 		z-
	//5 - above	 		z+

	//Need to grab the correct site
	snAelec = (SNarray) getElectrode_AdjacentSites(elXb);
	SiteNode site = getSN(snAelec,0,j,k);
	printSN(site);
	l = HoppingToSurroundingSites(site,codeX,codeY,codeZ);

	printf("Assinging future id\n");
	if(l==0){
		//Hopped to electrode behind
		*future = getAtotal(snA) + 0;
	}else if(l==1){
		*future = getIndFroP(snA,i,j,k);
	}else if(l==2){
		*future = getIndLefP(snA,i,j,k);
	}else if(l==3){
		*future = getIndRigP(snA,i,j,k);
	}else if(l==4){
		*future = getIndBelP(snA,i,j,k);
	}else{
		*future = getIndAboP(snA,i,j,k);
	}

	return 0;
}

int HopToElecY(SNarray snA, Electrode elYl, Charge * one, int * future, int EndX, int EndZ){
	//Elecid defines which electrode could be hopping too
	//where:
	//	getAtotal(snA) + 0 - (-x)
	//	getAtotal(snA) + 1 - (+x)
	//	getAtotal(snA) + 2 - (-y)
	// 	getAtotal(snA) + 3 - (+y)
	//	getAtotal(snA) + 4 - (-z)
	//	getAtotal(snA) + 5 - (+z)

	if(elYl==NULL){
		return -1;
	}

	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	int xx = getCx(*one);
	int yy = getCy(*one);
	int zz = getCz(*one);

	int factorX;
	int factorY;
	int factorZ;

	int l;
	int i;
	int j;
	int k;

	double position;

	if(xx<0){
		factorX = -xx/SLength+1;
		i = (xx+SLength*factorX)%SLength;
	}else{
		i = (xx+SLength)%SLength;
	}

	if(yy<0){
		factorY = -yy/SWidth+1;
		j = (yy+SWidth*factorY)%SWidth;
	}else{
		j = (yy+SWidth)%SWidth;
	}

	if(zz<0){
		factorZ = -zz/SHeight+1;
		k = (zz+SHeight*factorZ)%SHeight;
	}else{
		k = (zz+SHeight)%SHeight;
	}

	//   codeX = 1 exclude beh hop
	//				 = 2 exclude front hop
	//   codeZ = 1 exclude bottom hop
	//				 = 2 exclude top hop

	int codeX = 0;
	int codeY = 0;
	int codeZ = 0;

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;

	if(yy==0){
		//Posibility to jump to (-x) electrode

		//Have to ensure that charge is not at the edge of the system in the x and z directon
		if(xx == 0 && zz == 0){
			//Next to the bottom back corner
			codeX = 1;
			codeZ = 1;
		}else if(xx == 0 && zz == EndZ*SHeight-1){
			//Next to the back and top corner
			codeX = 1;
			codeZ = 2;
		}else if(zz == 0 && xx == EndX*SLength-1){
			//Next to the front and bottom corner
			codeX = 2;
			codeZ = 1;
		}else if(zz == EndZ*SHeight-1 && EndX*SLength-1){
			//Next to the front and top corner
			codeX = 2;
			codeZ = 2;
		}else if(xx == 0){
			//Next to the back side
			codeX = 1;
		}else if(zz == 0){
			//Next to the bottom edge
			codeZ = 1;
		}else if(zz == EndZ*SHeight-1){
			codeZ = 2;
		}else if(xx == EndX*SLength-1){ 
			codeX = 2;
		}
	}

	//l value
	//0 - backward 	x-
	//1 - forward 	x+
	//2 - left 			y-
	//3 - right 		y+
	//4 - below 		z-
	//5 - above	 		z+

	//Need to grab the correct site
	SNarray snAelec = (SNarray) getElectrode_AdjacentSites(elYl);
	SiteNode site = getSN(snAelec,i,0,k);
	l = HoppingToSurroundingSites(site, codeX,codeY,codeZ);

	if(l==0){
		*future = getIndBehP(snA,i,j,k);
	}else if(l==1){
		*future = getIndFroP(snA,i,j,k);
	}else if(l==2){
		//Hopping to electrode to the left
		*future = getAtotal(snA)+2;
	}else if(l==3){
		*future = getIndRigP(snA,i,j,k);
	}else if(l==4){
		*future = getIndBelP(snA,i,j,k);
	}else{
		*future = getIndAboP(snA,i,j,k);
	}

	return 0;
}

int HopToElecZ(SNarray snA, Electrode elZb, Charge * one, int * future, int EndX, int EndY){
	//Elecid defines which electrode could be hopping too
	//where:
	//	getAtotal(snA) + 0 - (-x)
	//	getAtotal(snA) + 1 - (+x)
	//	getAtotal(snA) + 2 - (-y)
	// 	getAtotal(snA) + 3 - (+y)
	//	getAtotal(snA) + 4 - (-z)
	//	getAtotal(snA) + 5 - (+z)

	if(elZb==NULL){
		return -1;
	}

	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	int xx = getCx(*one);
	int yy = getCy(*one);
	int zz = getCz(*one);

	int factorX;
	int factorY;
	int factorZ;

	int l;
	int i;
	int j;
	int k;	

	double position;

	if(xx<0){
		factorX = -xx/SLength+1;
		i = (xx+SLength*factorX)%SLength;
	}else{
		i = (xx+SLength)%SLength;
	}

	if(yy<0){
		factorY = -yy/SWidth+1;
		j = (yy+SWidth*factorY)%SWidth;
	}else{
		j = (yy+SWidth)%SWidth;
	}

	if(zz<0){
		factorZ = -zz/SHeight+1;
		k = (zz+SHeight*factorZ)%SHeight;
	}else{
		k = (zz+SHeight)%SHeight;
	}

	//   codeX = 1 exclude beh hop
	//				 = 2 exclude front hop
	//   codeY = 1 exclude left hop
	//				 = 2 exclude right hop

	int codeX = 0;
	int codeY = 0;
	int codeZ = 0;

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;

	if(zz==0){
		//Posibility to jump to (-z) electrode

		//Have to ensure that charge is not at the edge of the system in the x and z directon
		if(xx == 0 && yy == 0){
			//Next to the left back corner
			codeX = 1;
			codeY = 1;
		}else if(xx == 0 && yy == EndY*SWidth-1){
			//Next to the back and right corner
			codeX = 1;
			codeY = 2;
		}else if(yy == 0 && xx == EndX*SLength-1){
			//Next to the front and left corner
			codeX = 2;
			codeY = 1;
		}else if(yy == EndY*SWidth-1 && EndX*SLength-1){
			//Next to the front and right corner
			codeX = 2;
			codeY = 2;
		}else if(xx == 0){
			//Next to the back side
			codeX = 1;
		}else if(yy == 0){
			//Next to the left edge
			codeY = 1;
		}else if(yy == EndY*SWidth-1){
			codeY = 2;
		}else if(xx == EndX*SLength-1){ 
			codeX = 2;
		}
	}

	//l value
	//0 - backward 	x-
	//1 - forward 	x+
	//2 - left 			y-
	//3 - right 		y+
	//4 - below 		z-
	//5 - above	 		z+

	//Need to grab the correct site
	SNarray snAelec = (SNarray) getElectrode_AdjacentSites(elZb);
	SiteNode site = getSN(snAelec,i,j,0);
	printSN(site);
	l = HoppingToSurroundingSites(site,codeX,codeY,codeZ);

	if(l==0){
		*future = getIndBehP(snA,i,j,k);
	}else if(l==1){
		*future = getIndFroP(snA,i,j,k);
	}else if(l==2){
		*future = getIndLefP(snA,i,j,k);
	}else if(l==3){
		*future = getIndRigP(snA,i,j,k);
	}else if(l==4){
		//Hopping to Electrode below
		*future = getAtotal(snA)+4;
	}else{
		*future = getIndAboP(snA,i,j,k);
	}

	return 0;
}

int ElecHopOffX(Electrode el, int * future, SNarray snA){

	int row;
	int col;
	double position;

	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;
	matrix mtx = (matrix) getElectrode_HopRates(el);
	//Cycle through the pvals until the pval
	//is equal or above the random number

	if(mtx==NULL){
		printf("ERROR matrix mtx does not exist in ElecHopOFFX\n");
		exit(1);
	}

	for(row=1;row<=getRows(mtx);row++){
		for(col=1;col<=getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				break;
			}
		}
		if( getE(mtx,row,col)>=position){
			break;
		}
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to future

	*future = getIndex(snA,0,row-1,col-1); 

	if(*future==-1){

		printf("Value or row %d Value of col %d future %d position %g\n",row,col,*future,position);
		printf("ERROR Hopping off Electrode and out of system\n");
		exit(1);
	}

	return 0;
}

int ElecHopOffY(Electrode el, int * future, SNarray snA){

	int row;
	int col;
	double position;

	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;
	matrix mtx = (matrix) getElectrode_HopRates(el);
	//Cycle through the pvals until the pval
	//is equal or above the random number

	for(row=1;row<getRows(mtx);row++){
		for(col=1;col<getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				break;
			}
		}
		if( getE(mtx,row,col)>=position){
			break;
		}
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to futre

	*future = getIndex(snA,row-1,0,col-1); 

	printf("Hopping off Y\n");
	if(*future==-1){
		printf("ERROR Hopping off Electrode and out of system\n");
		exit(1);
	}

	return 0;
}

int ElecHopOffZ(Electrode el, int * future, SNarray snA){

	int row;
	int col;
	double position;

	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;
	matrix mtx = (matrix) getElectrode_HopRates(el);
	//Cycle through the pvals until the pval
	//is equal or above the random number

	for(row=1;row<getRows(mtx);row++){
		for(col=1;col<getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				break;
			}
		}
		if( getE(mtx,row,col)>=position){
			break;
		}
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to futre

	*future = getIndex(snA,row-1,col-1,0); 
	printf("Hopping off Z\n");
	if(*future==-1){
		printf("ERROR Hopping off Electrode and out of system\n");
		exit(1);
	}
	return 0;
}

int SiteHop(SNarray snA, Charge * one, SiteNode site, int * future, const int EndX, const int EndY, const int EndZ, const int XElecOn, const int YElecOn, const int ZElecOn, const int PeriodicX, const int PeriodicY, const int PeriodicZ){

	//totalX is used to measure the total current in the X direction as the charges move around
	//totalY is used to measure the total current in the Y direction as the charges move around
	//totalZ is used to measure the total current in the Z direction as the charges move around
	if((PeriodicX==1 && EndX==1)	|| (PeriodicX==0 && EndX!=1) ||\
			(PeriodicY==1 && EndY==1)	|| (PeriodicY==0 && EndY!=1) ||\
			(PeriodicZ==1 && EndZ==1)	|| (PeriodicZ==0 && EndZ!=1)){
		printf("Periodic Condition triger\n");
		exit(1);
	}

	if((EndX==0 && XElecOn==1) ||\
			(EndY==0 && YElecOn==1) ||\
			(EndZ==0 && ZElecOn==1)){
		printf("Electrode triger\n");
		return -1;
	}

	if(snA==NULL || *one == NULL || site ==NULL){
		printf("NULL datastruct triger\n");
		return -1;
	}

	//Initially, the flag is set to 0 
	//means the charge is able to hop
	int flag = 0; 
	double position;
	int rv;
	int l;
	//Below values define location of the charge in the
	//sample when the periodicity is removed
	int x, xx;
	int y, yy;
	int z, zz;
	int SLength;
	int SWidth;
	int SHeight;

	int factorX;
	int factorY;
	int factorZ;

	SLength = getAlen(snA);
	SWidth = getAwid(snA);
	SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	xx = getCx(*one);
	yy = getCy(*one);
	zz = getCz(*one);

	factorX=0;
	factorY=0;
	factorZ=0;

	if (xx<0){
		factorX = -xx/SLength+1;
		x = (xx+SLength*factorX)%SLength;
	}else{
		x = (xx+SLength)%SLength;
	}
	if (yy<0){
		factorY = -yy/SWidth+1;
		y = (yy+SWidth*factorY)%SWidth;
	}else{
		y = (yy+SWidth)%SWidth;
	}
	if (zz<0){
		factorZ = -zz/SHeight+1;
		z = (zz+SHeight*factorZ)%SHeight;
	}else{
		z = (zz+SHeight)%SHeight;
	}

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;

	if(x!=0 && x!=(SLength-1) && y!=0 && y!=(SWidth-1) && z!=0 && z!=(SHeight-1)){
		rv = OccAllNei(snA, x, y, z);
		if (rv==1){
			//All neighboring sites are occupied
			flag=1;
		}

		//Choose a site regardless function:MakeHop will ensure that the site is unoccupied
		//when the charge tries to move

		l=0;
		//Cycle through the pvals until the pval
		//is equal or above the random number
		//printf("Value of position %g\n",position);

		while( getSN_p(site,l)<position){
			//printf("pval SN %g\n",getSN_p(site,l));
			l++;
		}

	}else{
		//If the charge is at a boundary requires special attention
		int codeX;
		int codeY;
		int codeZ;

		if (PeriodicX==1){
			//Periodic in x
			if(EndX==0){
				//In this case the charge will never hit a boundary
				codeX=0;

			}else{
				//Charge might hit boundary depending on how many samples
				//charge has traveled through

				if (XElecOn==1){
					//If there are electrodes at the boundaries will always include the pvals
					codeX=0;
				}else{
					if(x==0){
						//Need to determine if charge is at the edge of the boundaries.
						if(xx/SLength==0){
							//Means it is at the boundary
							codeX=1;
						}else{
							//Means not at the boundary
							codeX=0;
						}
					}else if(x==(SLength-1)){

						if(xx/SLength==(EndX-1) || xx/SLength==EndX){
							//Hop forward is not allowed 
							codeX=2;
						}else{
							codeX=0;
						}
					}else{
						codeX=0;
					}

				}

			}

		}else{
			//Non periodic in x
			if(XElecOn==0){
				//This means there are no electrodes on the yz electrodes
				//on the x axis
				if(x==0){
					//Located on the back side should only jump forward
					codeX=1;
				}else if(x==SLength-1){
					//Located on the front side should only jump backward
					codeX=2;
				}else{
					codeX=0;
				}
			}else{
				//There are electrodes
				codeX=0;
			}
		}


		if (PeriodicY==1){
			//Periodic in y
			if(EndY==0){
				//In this case the charge will never hit a boundary
				codeY=0;

			}else{
				//Charge might hit boundary depending on how many samples
				//charge has traveled through
				//printf("SHould have made it here\n");
				//printf("Value of YElecOn %d\n",YElecOn);
				if(YElecOn==1){
					codeY=0;
					//printf("Value of codeY %d\n",codeY);
				}else{
					if(y==0){
						if(yy/SWidth==0){
							codeY=1;
						}else{
							codeY=0;
						}
					}else if(y==(SWidth-1)){
						if(yy/SWidth==(EndY-1) || yy/SWidth==EndY){
							codeY=2;
						}else{
							codeY=0;
						}
					}else{
						codeY=0;
					}
				}

			}

		}else{
			//Non periodic in y
			if(YElecOn==0){
				//This means there are no electrodes on the yz electrodes
				//on the x axis
				if(y==0){
					//Should only be able to jump to the right
					codeY=1;
				}else if(y==SWidth-1){
					//Should only be able to jump to the left
					codeY=2;
				}else{
					codeY=0;
				}
			}else{
				codeY=0;
			}

		}

		if (PeriodicZ==1){
			//Periodic in z
			if(EndZ==0){
				//In this case the charge will never hit a boundary
				codeZ=0;

			}else{
				//Charge might hit boundary depending on how many samples
				//charge has traveled through
				if(ZElecOn==1){
					codeZ=0;
				}else{
					if(z==0){
						if(zz/SHeight==0){
							codeZ=1;
						}else{
							codeZ=0;
						}
					}else if(z==(SHeight-1)){
						if(zz/SWidth==EndZ-1 || zz/SWidth==EndZ){
							codeZ=2;
						}else{
							codeZ=0;
						}
					}else{
						codeZ=0;
					}
				}
			}
		}else{
			//Non periodic in z

			if(ZElecOn==0){
				//This means there are no electrodes on the yz electrodes
				//on the x axis
				if(z==0){
					codeZ=1;
				}else if(z==SHeight-1){
					codeZ=2;
				}else{
					codeZ=0;
				}
			}else{
				codeZ=0;
			}
		}

		//l value
		//0 - backward 	x-
		//1 - forward 	x+
		//2 - left 			y-
		//3 - right 		y+
		//4 - below 		z-
		//5 - above	 		z+
		l = HoppingToSurroundingSites(site, codeX, codeY, codeZ);
	}

	if(l==0){
		*future = getIndBehP(snA,x,y,z);
	}else if(l==1){
		*future = getIndFroP(snA,x,y,z);
	}else if(l==2){
		*future = getIndLefP(snA,x,y,z);
	}else if(l==3){
		*future = getIndRigP(snA,x,y,z);
	}else if(l==4){
		*future = getIndBelP(snA,x,y,z);
	}else{
		*future = getIndAboP(snA,x,y,z);
	}

	printf("Value assigned to future %d\n",*future);
	return 0;

}

int ClusterHop(SNarray snA, Charge * ch,  double * tim , int *newID){

	if (snA==NULL || ch==NULL ){
		return -1;
	}

	//Flags determine what options a charge has if it is in a cluster
	int flag1;
	int flag2;

	int rv;
	int x;
	int y;
	int z;
	int ID;

	int factorX;
	int factorY;
	int factorZ;

	double position;
	double position2;

	double timeflag;

	SiteNode sn;

	flag1 = 0;
	flag2 = 0;

	factorX = 0;
	factorY = 0;
	factorZ = 0;

	if(getCx(*ch)<0){
		factorX = -getCx(*ch)/getAlen(snA)+1;
		x = (getCx(*ch)+getAlen(snA)*factorX)%getAlen(snA);
	}else{
		x = (getCx(*ch)+getAlen(snA))%getAlen(snA);
	}
	if(getCy(*ch)<0){
		factorY = -getCy(*ch)/getAwid(snA)+1;
		y = (getCy(*ch)+getAwid(snA)*factorY)%getAwid(snA);
	}else{
		y = (getCy(*ch)+getAwid(snA))%getAwid(snA);
	}
	if(getCz(*ch)<0){
		factorZ = -getCz(*ch)/getAhei(snA)+1;
		z = (getCz(*ch)+getAhei(snA)*factorZ)%getAhei(snA);
	}else{
		z = (getCz(*ch)+getAhei(snA))%getAhei(snA);
	}

	//printf("xx %d yy %d zz %d\n",getCx(*ch),getCy(*ch),getCz(*ch));
	//printf("FactorX %d FactorY %d FactorZ %d\n",factorX, factorY,factorZ);
	//printf("cx cy and cz divided by %d, %d and %d : %d %d %d\n",getAlen(snA), getAwid(snA), getAhei(snA), getCx(*ch)/getAlen(snA),getCy(*ch)/getAwid(snA),getCz(*ch)/getAhei(snA));
	//newID does not change unless
	//hop to a different position
	ID = getIndex(snA, x,y,z);
	*newID = ID;

	position = (double)((double)rand()/(double)RAND_MAX);
	position2 = (double)((double)rand()/(double)RAND_MAX);

	*tim = 0;

	flag1 = OccAllNeiCluster(snA,x,y,z);
	flag2 = OccAllCluster(snA, x,y,z);

	//flag1:
	//0 - At least one unoccupied neighbor
	//1 - all neighboring sites are occupied
	//-1 - malformed input or not part of a cluster

	//flag2
	//0 - at least one site within the cluster is 
	//unoccupied
	//1 - no available sites within the cluster
	//-1 - malformed input or not part of a cluster

	if(flag1==0 && flag2==0){
		//Hopping within cluster and to neighbors is allowed

		//Value of rv
		//1 - hopped off cluster
		//0 - stayed within cluster
		printf("ID %d position %g\n",ID,position);
		rv = HopOnOffCluster(snA,ID,position);
		printf("HopOnOffCluster\n");

		if(rv==1){
			timeflag = 2;

			//Hopped off cluster determining which neighboring 
			//site hopped to as well as the time it takes
			HopOffCluster(snA, ID, position2, newID, tim);

		}else if(rv==0){
			timeflag = 0;
			//Stayed within cluster now will determine which
			//site the charge moves to within the cluster
			HopWithinCluster(snA, ID, position2, newID);
		}

	}else if(flag1==0 && flag2==1){
		//Hopping within cluster not allowed
		//Hopping to neighbors is permitted
		sn = getSN(snA, x,y,z);
		//printf("Hopping Within cluster not allowed\n");
		//Value of rv
		//1 - hopped off cluster
		//0 - stayed within cluster
		rv = HopOnOffCluster(snA,ID,position);

		//printf("HopOnOffCluster %d\n",rv);
		if(rv==0){
			timeflag = 0;
			//should use tcluster to increment the time
			//printf("Staying on cluster charge should not change positions\n");
		}else{
			timeflag = 2;
			//should use hop rate associated with the chosen 
			//neighbor hop

			//Hopped off cluster determining which neighboring 
			//site hopped to as well as the time it takes
			HopOffCluster(snA, ID, position2, newID, tim);
		}


	}else if(flag1==1 && flag2==0){
		//Hopping within cluster is allowed
		//Hopping to neighbors is not allowed

		rv = HopOnOffCluster(snA,ID,position);

		//printf("Hopping to neighbors is not allowed\n");
		if(rv==0){
			timeflag = 0;
			//should use tcluster to increment the time

			HopWithinCluster(snA, ID, position2, newID);
			//Determine which site the charge will move to

		}else{
			timeflag = 1;
			//should use tprob to increment the time
		}


	}else if(flag1==1 && flag2==1){
		//Hopping to neighbors or within cluster
		//is not allowed
		timeflag = 0;
		//printf("Hopping is not allowed\n");
		//Should use tcluster to increment the time

	}else{
		//malformed input or not part of a cluster
		return -1;
	}

	if(timeflag==0){
		SiteNode sn = getSNwithInd(snA,ID);
		ClusterLL ClLL = (ClusterLL) getClusterList(sn);
		*tim = getCluster_time(ClLL);
	}else if(timeflag==1){
		SiteNode sn = getSNwithInd(snA,ID);
		ClusterLL ClLL = (ClusterLL) getClusterList(sn);
		*tim = getCluster_time(ClLL)*20;
	}

	return 0;

}


//Can use MakeHop when jumping off electrode but not when jumping off because (newID?)
int MakeHop(SNarray snA, int newID, Charge *ch, int * totalX, int * totalY, int * totalZ,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
		const int XElecOn, const int YElecOn, const int ZElecOn){

	//Now we need to correctly move the charge to
	//the new position
	int x, y, z;
	int x1, y1, z1;
	int XDiff, YDiff, ZDiff;
	int xD1, yD1, zD1;
	int xD2, yD2, zD2;

	int factorX=0;
	int factorY=0;
	int factorZ=0;

	if(getCx(*ch)<0){
		factorX = -getCx(*ch)/getAlen(snA)+1;
		x = (getCx(*ch)+getAlen(snA)*factorX)%getAlen(snA);
	}else{
		x = (getCx(*ch)+getAlen(snA))%getAlen(snA);
	}
	if(getCy(*ch)<0){
		factorY = -getCy(*ch)/getAwid(snA)+1;
		y = (getCy(*ch)+getAwid(snA)*factorY)%getAwid(snA);
	}else{
		y = (getCy(*ch)+getAwid(snA))%getAwid(snA);
	}
	if(getCz(*ch)<0){
		factorZ = -getCz(*ch)/getAhei(snA)+1;
		z = (getCz(*ch)+getAhei(snA)*factorZ)%getAhei(snA);
	}else{
		z = (getCz(*ch)+getAhei(snA))%getAhei(snA);
	}

	//	getAtotal(snA) + 0 - (-x)
	//	getAtotal(snA) + 1 - (+x)
	//	getAtotal(snA) + 2 - (-y)
	// 	getAtotal(snA) + 3 - (+y)
	//	getAtotal(snA) + 4 - (-z)
	//	getAtotal(snA) + 5 - (+z)

	XDiff = 0;
	YDiff = 0;
	ZDiff = 0;

	//Hopping to the back electrode
	if(newID==getAtotal(snA)){
		XDiff = -1;
		//Hopping to the front electrode etc
	}else if(newID==getAtotal(snA)+1){
		XDiff = 1;
	}else if(newID==getAtotal(snA)+2){
		YDiff = -1;
	}else if(newID==getAtotal(snA)+3){
		YDiff = 1;
	}else if(newID==getAtotal(snA)+4){
		ZDiff = -1;
	}else if(newID==getAtotal(snA)+5){
		ZDiff = 1;
	}else{

		//Hopping to a site within the system
		getLoc(&x1, &y1, &z1, newID, snA);
		SiteNode sn2 = getSNwithInd(snA,newID);
		//Check to see if new Position is occupied
		if(getDwelStat(sn2)==-1){

			XDiff = x1-x;
			YDiff = y1-y;
			ZDiff = z1-z;

			//We can go ahead and assume that if the difference
			//is larger than the length, width or height then
			//the sample has crossed a periodic boundary or
			//it has reached an electrode
			if(PeriodicX==1 || XElecOn==1){
				int SLength = getAlen(snA);
				xD1 = (x1+SLength)-x;
				xD2 = x1-(x+SLength);
				if( XDiff*XDiff<=xD1*xD1 && XDiff*XDiff<=xD2*xD2 ){
					XDiff=XDiff;
				}else if( xD1*xD1<=xD2*xD2 ){
					XDiff=xD1;
				}else{
					XDiff=xD2;
				}
			}

			if(PeriodicY==1 || YElecOn){
				int SWidth = getAwid(snA);
				yD1 = (y1+SWidth)-y;
				yD2 = y1-(y+SWidth);
				if( YDiff*YDiff<=yD1*yD1 && YDiff*YDiff<=yD2*yD2 ){
					YDiff=YDiff;
				}else if( yD1*yD1<=yD2*yD2 ){
					YDiff=yD1;
				}else{
					YDiff=yD2;
				}
			}

			if(PeriodicZ==1 || ZElecOn){
				int SHeight = getAwid(snA);
				zD1 = (z1+SHeight)-z;
				zD2 = z1-(z+SHeight);
				if( ZDiff*ZDiff<=zD1*zD1 && ZDiff*ZDiff<=zD2*zD2 ){
					ZDiff=ZDiff;
				}else if( zD1*zD1<=zD2*zD2 ){
					ZDiff=zD1;
				}else{
					ZDiff=zD2;
				}
			}

			*totalX = *totalX+XDiff;
			*totalY = *totalY+YDiff;
			*totalZ = *totalZ+ZDiff;

			setCx(*ch,getCx(*ch)+XDiff);
			setCy(*ch,getCy(*ch)+YDiff);
			setCz(*ch,getCz(*ch)+ZDiff);

			//printf("After hop Cx %d Cy %d Cz %d\n",getCx(*ch),getCy(*ch),getCz(*ch));
			//This means the charge sucessfully hopped
			//to a new location
			return 0;
		}else{
			//This means the site the charge was going 
			//to hop to was occupied so it didn't move
			//anywhere but time passed
			return -1;
		}
	}

	*totalX = *totalX+XDiff;
	*totalY = *totalY+YDiff;
	*totalZ = *totalZ+ZDiff;

	setCx(*ch,getCx(*ch)+XDiff);
	setCy(*ch,getCy(*ch)+YDiff);
	setCz(*ch,getCz(*ch)+ZDiff);
	return 0;

}

int CheckPt_Test(int * CheckPtNum,int CheckFileExist, char * FileNameCheckPtVersion, int FileNameSize,\
		const double Vx, const double Vy, const double Vz, const double Temperature){


	if(CheckFileExist==0){
		*CheckPtNum = CheckPt_Latest(FileNameCheckPtVersion,FileNameSize,Vx,Vy,Vz,Temperature);
		//Just because a .chpt file exist it does not mean it exists for these parameters
		//Will return a 0 if unable to find a checkpt file with the correct Vx Vy Vz and Temperature
		//Or unable to find a checkpt file at all in the directory

		if(*CheckPtNum>0){
			*CheckPtNum = *CheckPtNum+1;
			//Confirmed that a checkpt file exists
			printf("CheckPoint file found for Vx %g Vy %g Vz %g and Temperature %g\n",Vx,Vy,Vz,Temperature);
			printf("File Name %s\n",FileNameCheckPtVersion);
			return 1;

		} else if(*CheckPtNum==0){
			//Confirmed that a checkpt file did not exist for this version
			printf("CheckPoint file does not exist for Vx %g Vy %g Vz %g and Temperature %g\n",Vx,Vy,Vz,Temperature);
			*CheckPtNum = 1;
			return 0;

		} else if(*CheckPtNum<0){

			//This means there was a serios problem either the directory was deleted
			printf("ERROR unable to open CHECKPOINT directory!\n");
			exit(1);
		}
	}else{
		*CheckPtNum = 1;
	}

	return 0;
}

int Pre_randomWalk(const int CheckPtStatus,char * FileNameCheckPtVersion,char * FileName,\
		long double * t, matrix * Sequence,\
		ChargeArray * chA,matrix * FutureSite,ArbArray * ClArLL, SNarray * snA, ParameterFrame PF,\
		double electricEnergyX, double electricEnergyY, double electricEnergyZ, int r,\
		double Vx, double Vy, double Vz, double Temperature, long int * n, int * nc, int * nca,\
		Electrode * elXb, Electrode * elXf, Electrode * elYl,\
		Electrode * elYr, Electrode * elZb, Electrode * elZa){

	//Declaring constants
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	//Declaring Variables from parameter frame
	int SLength;
	int SWidth;
	int SHeight;
	int EndX;
	int EndY;
	int EndZ;
	int XElecOn;
	int YElecOn;
	int ZElecOn;
	int Ntot;
	int NCh;
	double AttemptToHop;
	double reOrgEnergy;
	double gamma;
	double SiteDistance;

	//Declaring local variables
	int clusterfileExist;
	int Num_elXf;
	int Num_elXb;
	int Num_elYl;
	int Num_elYr;
	int Num_elZb;
	int Num_elZa;
	double MarcusJ0;
	double MarcusCoeff;
	double KT;
	int OrderL;

	//Initializing Variables from parameter frame
	SLength = PFget_Len(PF);
	SWidth = PFget_Wid(PF);
	SHeight = PFget_Hei(PF);
	EndX = PFget_EndX(PF);
	EndY = PFget_EndY(PF);
	EndZ = PFget_EndZ(PF);
	XElecOn = PFget_XElecOn(PF);
	YElecOn = PFget_YElecOn(PF);
	ZElecOn = PFget_ZElecOn(PF);
	Ntot = PFget_Ntot(PF);
	NCh = PFget_NCh(PF);
	AttemptToHop = PFget_AttemptToHop(PF);
	reOrgEnergy = PFget_reOrg(PF);
	gamma = PFget_gamma(PF);
	SiteDistance = PFget_SiteDist(PF);

	//Initializing Local variables
	KT = kB*Temperature;
	*snA = newSNarray(SLength, SWidth, SHeight);
	*Sequence = newMatrix(Ntot,1);
	(*FutureSite) = newMatrix(Ntot,1);
	int rv = printMatrix(*FutureSite);

	printf("Value of Ntot %d\n",Ntot);
	if(rv==-1){
		printf("ERROR FutureSite problem!\n");
		exit(1);
	}

	if(FutureSite==NULL){
		printf("ERROR FutureSite NULL\n");
		exit(1);
	}
	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	//Calculating full Marcus Coefficient;
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(-2*gamma*SiteDistance);

	//printf("MarcusCoeff %g\n",MarcusCoeff);
	//exit(1);

	//This means we are starting from scratch
	if(CheckPtStatus==0){

		//Lets first make sure that a .cluster file does not exist before we
		//create new energies
		clusterfileExist = CheckPt_Cluster(Vx, Vy, Vz, Temperature, r);

		if(clusterfileExist==0){
			//The cluster file exists already so lets just load what was saved
			LoadCluster_Data( &FileName[0], &OrderL, snA, electricEnergyX,\
												electricEnergyY, electricEnergyZ, ClArLL, KT);

		}else{
			//Initialize Site Energies
			printf("Initializing Site Energies\n");
			initSite(electricEnergyX, electricEnergyY, electricEnergyZ,\
					KT, *snA, PF);
		}
		//Create Sequence matrix to store the charges and the order 
		//they should be moved which is based on their dwelltime
		//Charge ids start at 0 and go to Ntot-1
		printf("Initializing Sequence of Charges\n");
		for(int i=0; i<Ntot;i++){
			setE(*Sequence,i+1,1,i);
		}

		//Initialize all charges in the Matrix and creates chargearray
		printf("Initializing Charges\n");
		*chA = initCharget0( *Sequence, *snA,  Ntot, NCh, PFget_D(PF),\
				XElecOn, YElecOn, ZElecOn,EndX, EndY, EndZ);

		//t - global time initially 0 when starting from scratch
		//n - number of steps that charges have been injected starts at 1
		//nc - Number of charges initially in the system equal to the Number
		//		initially injected
		//nca - Number of active charges in the system initiallyequal to the number injected
		*t = 0;
		*n = 1;
		*nc = NCh;
		*nca = NCh;

		////////////////////////////////////////////////////////////////////////
		//Initialize Electrodes
		printf("Initializing Electrodes\n");
		initElec(electricEnergyX, electricEnergyY, electricEnergyZ, MarcusCoeff,\
				KT, *snA,elXb, elXf, elYl, elYr, elZb, elZa, PF);
		
		//Initialize where the first few charges will hop to next
		printf("Initializing future hop sites\n");
		initFutureSite( snA, FutureSite, chA, PF,\
				*elXb, *elYl, *elZb );


		if(FutureSite==NULL || Sequence==NULL || chA==NULL){
			return -1;
		}

		if(clusterfileExist==-1){
			//Go ahead and calculate clusters and save
			FindCluster( &OrderL, (*snA), electricEnergyX,\
									 electricEnergyY, electricEnergyZ,\
									 ClArLL,KT,PF);

			ConnectClusterElec( ClArLL,\
												(*elXb), (*elXf), (*elYl), (*elYr),\
												(*elZb), (*elZa) );

			SaveCluster( &FileName[0], OrderL, (*snA), electricEnergyX,\
									 electricEnergyY, electricEnergyZ, (*ClArLL), KT, PF,
									 *elXb, *elXf, *elYl, *elYr, *elZb, *elZa);

			PrintFile_xyz(OrderL, (*snA), ClArLL, &FileName[0]);
			printFileEnergy((*snA), &FileName[0], electricEnergyX,\
					 						electricEnergyY, electricEnergyZ,PF);

			PrintNeighFile_xyz(OrderL, (*snA), ClArLL, &FileName[0]);
			printMatrix(*FutureSite);

		}


	}else if(CheckPtStatus==1) {
		//This means we are starting from a checkpt file that already exists

		//Lets first make sure that a .cluster file does not exist before we
		//recalculate everything
		clusterfileExist = CheckPt_Cluster(Vx, Vy, Vz, Temperature, r);
		
		//Charge Array chA is created in here
		printf("Loading Charge and Site information from .ckpt\n");
		Load_CheckPt_Data( t, snA, chA, Sequence,\
				FutureSite, FileNameCheckPtVersion, n,nc, nca,\
				&Num_elXb, &Num_elXf, &Num_elYl, &Num_elYr,\
				&Num_elZb, &Num_elZa,Vx,Vy,Vz);
		
		if(FutureSite==NULL || Sequence==NULL || chA==NULL){
			printf("ERROR in the load_checkPt_Data function a datastructure\n");
			printf("Has been found to be NULL\n");
			exit(1);
		}

		////////////////////////////////////////////////////////////////////////
		//Initialize Electrodes
		printf("Initializing Electrodes\n");
		initElec(electricEnergyX, electricEnergyY, electricEnergyZ, MarcusCoeff,\
				KT, *snA,elXb, elXf, elYl, elYr, elZb, elZa, PF);


		printf("Updating number of charges on electrodes based on .ckpt\n");
		if(elXb!=NULL){
			setElectrode_Charges(*elXb,Num_elXb);
		}
		if(elXf!=NULL){
			setElectrode_Charges(*elXf,Num_elXf);
		}
		if(elYl!=NULL){
			setElectrode_Charges(*elYl,Num_elYl);
		}
		if(elYr!=NULL){
			setElectrode_Charges(*elYr,Num_elYr);
		}
		if(elZb!=NULL){
			setElectrode_Charges(*elZb,Num_elZb);
		}
		if(elZa!=NULL){
			setElectrode_Charges(*elZa,Num_elZa);
		}

		if(clusterfileExist==0){
			//The .cluster file exists we have already loaded the site information
			//We only want to load cluster information
			LoadCluster_Only( &FileName[0], &OrderL, snA, electricEnergyX,\
												electricEnergyY, electricEnergyZ, ClArLL, KT);
		}
		
		ConnectClusterElec( ClArLL,\
												(*elXb), (*elXf), (*elYl), (*elYr),\
												(*elZb), (*elZa) );
	
	}else{
		printf("ERROR found in Pre_randomWalk value of checkptstatus %d\n",CheckPtStatus);
		exit(1);
	}


	printf("Printing Future site matrix at end of Pre_randomWalk\n");
	rv = printMatrix(*FutureSite);
	if(FutureSite==NULL || rv==-1){
		printf("FutureSite is NULL\n");
		exit(1);
	}

	return 0;

}

int Post_randomWalk(ArbArray ClArLL, SNarray snA, Electrode elXb, Electrode elXf,\
		Electrode elYl, Electrode elYr, Electrode elZb, Electrode elZa){


	if(snA==NULL){
		printf("ERROR snA found to be NULL\n");
		return -1;
	}
	if(ClArLL==NULL){
		printf("ERROR ClArLL found to be NULL\n");
		return -1;
	}

	SNarray snAmini;
	matrix mtxmini;

	printf("Deleting Electrodes\n");
	if (elXb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXb);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXb);
		deleteMatrix(mtxmini);
		deleteElectrode(&elXb);
		elXb = NULL;
	}
	if (elXf!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXf);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXf);
		deleteMatrix(mtxmini);
		deleteElectrode(&elXf);
		elXf = NULL;
	}


	if (elYl!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYl);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYl);
		deleteMatrix(mtxmini);
		deleteElectrode(&elYl);
		elYl = NULL;
	}
	if (elYr!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYr);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYr);
		deleteMatrix(mtxmini);
		deleteElectrode(&elYr);
		elYr = NULL;
	}
	if (elZb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZb);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZb);
		deleteMatrix(mtxmini);
		deleteElectrode(&elZb);
		elZb = NULL;
	}

	if (elZa!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZa);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZa);
		deleteMatrix(mtxmini);
		deleteElectrode(&elZa);
		elZa = NULL;
	}

	printf("Deleteting Cluster Arbitrary Array Link List\n");
	deleteArbArray(&ClArLL);
	printf("Deleting SiteNode array\n");
	deleteSNarray(snA);

	return 0;
}

int initFutureSite( SNarray * snA, matrix * FutureSite,ChargeArray * chA, ParameterFrame PF,\
		Electrode elXb, Electrode elYl, Electrode elZb ){

	printf("Initializing Future sites\n");
	//Declaring Parameter frame variables
	int SLength;
	int SWidth;
	int SHeight;
	int PeriodicX;
	int PeriodicY;
	int PeriodicZ;
	int EndX;
	int EndY;
	int EndZ;
	int XElecOn;
	int YElecOn;
	int ZElecOn;
	int NCh;
	int rv;

	//Declaring local variables
	int loop;
	int factorX;
	int factorY;
	int factorZ;
	int x;
	int y;
	int z;
	int future;

	Charge one;
	SiteNode site;

	//Initiallizing Parameter Frame variables
	SLength = PFget_Len(PF);
	SWidth = PFget_Wid(PF);
	SHeight = PFget_Hei(PF);
	PeriodicX = PFget_Px(PF);
	PeriodicY = PFget_Py(PF);
	PeriodicZ = PFget_Pz(PF);
	EndX = PFget_EndX(PF);
	EndY = PFget_EndY(PF);
	EndZ = PFget_EndZ(PF);
	XElecOn = PFget_XElecOn(PF);
	YElecOn = PFget_YElecOn(PF);
	ZElecOn = PFget_ZElecOn(PF);
	NCh = PFget_NCh(PF);

	for(loop = 0; loop<NCh; loop++){

		one = getCharge(*chA,loop);

		factorX = 0;
		factorY = 0;
		factorZ = 0;

		if(getCx(one)<0){
			factorX = -getCx(one)/SLength+1;
			x = (getCx(one)+SLength*factorX)%SLength;
		}else{
			x = (getCx(one)+SLength)%SLength;
		}
		if(getCy(one)<0){
			factorY = -getCy(one)/SWidth+1;
			y = (getCy(one)+SWidth*factorY)%SWidth;
		}else{
			y = (getCy(one)+SWidth)%SWidth;
		}
		if(getCz(one)<0){
			factorZ = -getCz(one)/SHeight+1;
			z = (getCz(one)+SHeight*factorZ)%SHeight;
		}else{
			z = (getCz(one)+SHeight)%SHeight;
		}

		site = getSN(*snA,x,y,z);
		if( x==0  && XElecOn==1 ){
			HopToElecX(*snA, elXb, &one, &future, EndY, EndZ);
		}else if( y==0 && YElecOn==1){
			HopToElecY(*snA, elYl, &one, &future, EndX,EndZ);
		}else if( z==0 && ZElecOn==1){
			HopToElecZ(*snA, elZb, &one, &future, EndX, EndY);
		}else{
			rv=SiteHop(*snA, &one, site, &future, EndX, EndY,EndZ,\
					XElecOn,YElecOn,ZElecOn, PeriodicX, PeriodicY, PeriodicZ);
		}
			
		setE(*FutureSite,loop+1,1,future);
	}

	printf("InitFuture Site Printing Matrix FutureSite\n\n");
	//printMatrix(*FutureSite);

	return 0;
}

int randomWalk( SNarray snA,int CheckptNum,\
		char * FileName, double ElectricFieldX,\
		double ElectricFieldY, double ElectricFieldZ,\
		Electrode elXb, Electrode elXf, Electrode elYl,\
		Electrode elYr, Electrode elZb, Electrode elZa,\
		ParameterFrame PF,long double t,matrix Sequence,\
		matrix FutureSite,ChargeArray * chA,\
		long int n,int nc,int nca){

	if(FutureSite==NULL){
		printf("ERROR Future site matrix found to be NULL on entering randomWalk\n");
		return -1;
	}


	time_t now;
	time_t later;
	double seconds;
	time(&now);

	//n - Number of time steps that have passed where charges are injected
	//nc - Number of charges in the system
	//nca - Number of charges that are currently active

	double TStep = PFget_TStep(PF);
	int TCount = PFget_TCount(PF);
	int Time_check = PFget_Time_check(PF);
	int Nstep_av = PFget_Nstep_av(PF);
	int NCh = PFget_NCh(PF);
	double D = PFget_D(PF);
	int XElecOn = PFget_XElecOn(PF);
	int YElecOn = PFget_YElecOn(PF);
	int ZElecOn = PFget_ZElecOn(PF);
	int EndX = PFget_EndX(PF);
	int EndY = PFget_EndY(PF);
	int EndZ = PFget_EndZ(PF);
	int PeriodicX = PFget_Px(PF);
	int PeriodicY = PFget_Py(PF);
	int PeriodicZ = PFget_Pz(PF);
	double SiteDistance = PFget_SiteDist(PF);
	int MovieFrames = PFget_MovieFrames(PF);
	
	printf("will crash here if MovieFrames is 0\n");
	assert(MovieFrames!=0);

	if(snA==NULL || TCount<0 || TStep < 0 || Nstep_av < 0 ||\
			NCh<0 || D < 0 || XElecOn<0 || YElecOn < 0 || ZElecOn <0 ||\
			XElecOn > 1 || YElecOn > 1 || ZElecOn > 1 ||\
			EndX<0 || EndY<0 || EndZ<0 || PeriodicX<0 || PeriodicY<0 ||\
			PeriodicZ<0 || PeriodicX>1 || PeriodicY>1 || PeriodicZ>1 ||\
			t<0){
		return -1;
	}

	//TCount - is the number of time increments that charges are injected
	//TStep - is the size of the time step
	//Nstep_av - this is used to determine how often the data is recorded
	//NCh - Is the number of Charges that are injected per increment where
	// 	 		TCount is the last increment
	//Ntot - is the total number of charges that will ever be put in the
	//			 system
	//rN - this is some kind of conversion factor that is used??? Don't really
	//			know what it is at this point

	printf("Starting Random Walk\n");
	//Declaring Variables
	int Ntot = NCh*TCount;
	int flag;
	int SaveCount;
	//Movie start point
	int Movie = (int) ((double)t/((double)TStep*(double)Nstep_av));

	int x, xx, x1;
	int y, yy, y1;
	int z, zz, z1;

	//Position before charge hopps
	int PrevX, PrevY, PrevZ;
	int i;
	int SLength;
	int SWidth;
	int SHeight;
	int CheckX;
	int CheckY;
	int CheckZ;
	int ClusterYes;
	//Lets assume TotalX, TotalY and TotalZ keep track of the
	//current for a given time increment
	int TotalX;
	int TotalY;
	int TotalZ;

	int TotalXtemp;
	int TotalYtemp;
	int TotalZtemp;

	int Xdrain;
	int Xsource;
	int Ydrain;
	int Ysource;
	int Zdrain;
	int Zsource;

	int CurrentInc;
	int TotalCollected;
	int ElecExit;
	int ChargeID;
	int future;
	int JumpFromElec;

	long double SaveTime;
	double ran;
	double q;

	int factorX;
	int factorY;
	int factorZ;

	int TrackX;
	int TrackY;
	int TrackZ;

	long double TimeTrack1;
	int NumAvgVel;
	long double TotalVelX;
	long double TotalVelY;
	long double TotalVelZ;

	double tim;
	Charge one;
	Charge two;
	SiteNode site;

	q = 1.602E-19;

	NumAvgVel = 0;

	SLength = getAlen(snA);
	SWidth = getAwid(snA);
	SHeight = getAhei(snA);

	//Calculate Current 
	matrix Xcurrent = newMatrix(8,1);
	matrix Ycurrent = newMatrix(8,1);
	matrix Zcurrent = newMatrix(8,1);

	//Keeps track of current in X, Y and Z direction
	//during time increment
	TotalX = 0;
	TotalY = 0;
	TotalZ = 0;
	TrackX = 0;
	TrackY = 0;
	TrackZ = 0;
	CurrentInc = 1;
	TimeTrack1 = 0;

	//Create source and drain for all electrodes that are turned on
	matrix Xelec_Drain = newMatrix(8,1);
	matrix Xelec_Source = newMatrix(8,1);
	matrix Yelec_Drain = newMatrix(8,1);
	matrix Yelec_Source = newMatrix(8,1);
	matrix Zelec_Drain = newMatrix(8,1);
	matrix Zelec_Source = newMatrix(8,1);

	//Drift Velocities should be in units of [m/s]
	matrix Xvelocity = newMatrix(8,1);
	matrix Yvelocity = newMatrix(8,1);
	matrix Zvelocity = newMatrix(8,1);

	matrix System = newMatrix(8,1);

	//Keeps track of the time when each of the matrices are incremented
	matrix timeArray = newMatrix(8,1);

	printf("Initialization Complete!\n");
	printf("Value of nca %d\n",nca);
	//Initilize number of charges reaching source
	//and drian to 0
	Xdrain = 0;
	Xsource = 0;
	Ydrain = 0;
	Ysource = 0;
	Zdrain = 0;
	Zsource = 0;

	TotalVelX = 0;
	TotalVelY = 0;
	TotalVelZ = 0;

	//Continue looping while n is less than the TCount
	//or if there are still charges in the system
	SaveCount = 1;
	SaveTime = 0;
	tim = TStep;

	while( n<TCount || nca>0){

		//If no charges have been inserted in the system we will
		//simply increment the time

		printf("nca %d nc %d elXb %d elXf %d\n",nca,nc,getElectrode_Charges(elXb),getElectrode_Charges(elXf));

		if( nc+getElectrode_Charges(elXb) != nca){
			printf("Charges lost somehow!\n");
			exit(1);
		}

		if (nca==0){

			//Should be the Step time
			tim = TStep;
			//t is incremented after the charge hops
			//the global time is increased
			t += (long double) tim;

			SaveTime += (long double) tim;
			//The time the charge has been in the sample is updated
			Plust(one,(long double) tim);

			if(SaveTime >= (long double)TStep){
				//Here we check to see if we record the data
				//Nstep is used to determine how many timesteps pass
				//before recording
				SaveCount++;

				printf("SaveTime %Lg TStep %g SaveCount %d Nstep_av %d Movie %d MovieFrames %d\n",SaveTime,TStep,SaveCount,Nstep_av,Movie,MovieFrames);

				if((SaveCount%Nstep_av)==0 && Movie<MovieFrames){
					printf("Should have entered the printMovie Routine\n");
					SaveCount = 1;
					printMovie(&Movie,FileName,snA);
					SaveTime = fmod(SaveTime, (long double)TStep);
					SaveDataPoint( &CurrentInc, &NumAvgVel, nc, XElecOn, YElecOn, ZElecOn,\
							&Xcurrent, &Ycurrent, &Zcurrent, &timeArray,\
							&Xvelocity,&Yvelocity, &Zvelocity,&System,\
							&Xelec_Drain,&Yelec_Drain, &Zelec_Drain,\
							&Xelec_Source,&Yelec_Source, &Zelec_Source,\
							&TotalX,&TotalY,&TotalZ,&TrackX,&TrackY,\
							&TrackZ, t,&TimeTrack1,TStep,SiteDistance,\
							&TotalVelX,&TotalVelY,&TotalVelZ,\
							&Xdrain, &Ydrain, &Zdrain, &Xsource, &Ysource,&Zsource,\
							ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
							FileName);

				}

				time(&later);
				seconds = difftime(later,now);

				if(seconds>(double)(Time_check*60)){
					time(&now);
					Save_CheckPt(FileName, &CheckptNum, snA, *chA, Sequence, FutureSite, t, PF,n, nc, nca,\
							elXb, elXf, elYl, elYr, elZb, elZa);
				}

				SaveTime = 0;

			}
		}else{

			//Grab the first charge in the sequence
			JumpFromElec = 0;
			ChargeID = getE(Sequence,1,1);
			one = getCharge(*chA, ChargeID);

			//printf("Number of charges in sample %d active charges %d time %Lg\n",nc,nca, t);
			//printf("Grabbing Charge %d Value of nca %d\n",ChargeID, nca);
			printf("ID of charge %d Position of Charge %d %d %d\n",ChargeID,getCx(one),getCy(one),getCz(one));
			if(getCx(one)<-1 || getCx(one)>SLength){
				printf("Charge position less than -1 or greater than SLength %d\n",SLength);
				exit(1);
			}

			//*********************REGARDLESS OF HOP OR NOT*******************************************

			//Attempt to make the charge hop to site 
			//that was previously determined. If site is
			//occupied return -1. If site unoccupied hop and
			//return 0

			factorX=0;
			factorY=0;
			factorZ=0;

			if(getCx(one)<0){
				factorX = -getCx(one)/SLength+1;
				x1 = (getCx(one)+SLength*factorX)%SLength;
			}else{
				x1 = (getCx(one)+SLength)%SLength;
			}
			if(getCy(one)<0){
				factorY = -getCy(one)/SWidth+1;
				y1 = (getCy(one)+SWidth*factorY)%SWidth;
			}else{
				y1 = (getCy(one)+SWidth)%SWidth;
			}
			if(getCz(one)<0){
				factorZ = -getCz(one)/SHeight+1;
				z1 = (getCz(one)+SHeight*factorZ)%SHeight;
			}else{
				z1 = (getCz(one)+SHeight)%SHeight;
			}

			//Coordinates of charge before hop
			PrevX = getCx(one);
			PrevY = getCy(one);
			PrevZ = getCz(one);

			//First we determine if the charge is hopping from an electrode or
			//from somewhere within the system
			if ((PrevX<0 && EndX!=0) || (PrevY<0 && EndY!=0) || (PrevZ<0 && EndZ!=0)){
				//Hopping from an electrode does not mean hopping to an electrode
				JumpFromElec = 1;
			}else{
				//This is the site jumping from
				site = getSN(snA,x1,y1,z1); 

			}
			//We don't automatically increase the distance if we don't know 
			//if the hop is succesful or not i.e. the charge may not move yet
			TotalXtemp = TotalX;
			TotalYtemp = TotalY;
			TotalZtemp = TotalZ;

			//Function accounts for hops within system from electrodes
			//and two electrodes
			//printMatrix(FutureSite);			
	
			flag = MakeHop(snA, getE(FutureSite,ChargeID+1,1),\
					&one, &TotalXtemp, &TotalYtemp, &TotalZtemp,\
					PeriodicX, PeriodicY, PeriodicZ,\
					XElecOn, YElecOn, ZElecOn);

			//Get the time it took to make the hop
			tim = getDwel(one);

			//Cannot simply increase the time by the dwelstat if not
			//all the charges have yet been inserted

			if (((n+1)*(long int)NCh)<=Ntot && tim>TStep){
				//This means at maximum can only increase the time by the 
				//TStep because more charges need to be inserted. 
				tim = TStep;
				//Because the site does not hop we just make it progress
				//in time a little we set the flag to -1
				flag = -1;
			}


			if(flag==0){
				//it did hop

				//Calculate the velocity
				if (TotalX!=TotalXtemp){
					TrackX +=(TotalXtemp-TotalX);
				}
				if (TotalY!=TotalYtemp){
					TrackY +=(TotalYtemp-TotalY);
				}
				if (TotalZ!=TotalZtemp){
					TrackZ +=(TotalZtemp-TotalZ);
				}

				//Reinitialize TotalX
				TotalX = TotalXtemp;
				TotalY = TotalYtemp;
				TotalZ = TotalZtemp;
			}

			//t is incremented after the charge hops
			//the global time is increased
			t += (long double) tim;

			SaveTime += (long double)tim;
			//The time the charge has been in the sample is updated
			Plust(one,(long double) tim);

			// If the global time has reached the next time step 
			// calculate current

			//printf("SaveTime %Lg TStep %g\n",SaveTime, TStep);
			printf("SaveTime %Lg TStep %g SaveCount %d Nstep_av %d Movie %d MovieFrames %d\n",SaveTime,TStep,SaveCount,Nstep_av,Movie,MovieFrames);

			if(SaveTime >= (long double)TStep){
				//Here we check to see if we record the data
				//Nstep is used to determine how many timesteps pass
				//before recording
				SaveCount++;

				if((SaveCount%Nstep_av)==0 && Movie<MovieFrames){
					//Saving Data
					SaveCount = 1;
					printf("Should have entered the printMovie Routine\n");
					printMovie(&Movie,FileName,snA);
					SaveTime = fmod(SaveTime,(long double)TStep);

					SaveDataPoint( &CurrentInc, &NumAvgVel, nc, XElecOn, YElecOn, ZElecOn,\
							&Xcurrent, &Ycurrent, &Zcurrent, &timeArray,\
							&Xvelocity,&Yvelocity, &Zvelocity,&System,\
							&Xelec_Drain,&Yelec_Drain, &Zelec_Drain,\
							&Xelec_Source,&Yelec_Source, &Zelec_Source,\
							&TotalX,&TotalY,&TotalZ,&TrackX,&TrackY,\
							&TrackZ, t,&TimeTrack1,TStep,SiteDistance,\
							&TotalVelX,&TotalVelY,&TotalVelZ,\
							&Xdrain, &Ydrain, &Zdrain, &Xsource, &Ysource,&Zsource,\
							ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
							FileName);

					//printMatrix(System);
					//printMatrix(timeArray);

				}

				time(&later);
				seconds = difftime(later,now);

				if(seconds>(double)(Time_check*60)){
					time(&now);
					Save_CheckPt(FileName, &CheckptNum, snA, *chA, Sequence,FutureSite, t, PF,n, nc, nca,\
							elXb, elXf, elYl, elYr, elZb, elZa);
				}

				SaveTime = 0;
			}

			//For all the charges still active 
			//need to decrease their dwelltime
			UpdateOccTime(&snA,&one,tim, PF);
			for( i = 1; i<nca; i++ ){
				two = getCharge(*chA,getE(Sequence,i,1));
				MinusDwel(two,tim);
				UpdateOccTime(&snA, &two,tim, PF);
			}
			//*********************************************************************//

			//This is used to keep track of whether or not the charge exited at
			//an electrode or not
			ElecExit = 0;

			////////////////////////////////////////////////////////////////////////
			//Sucessful hop 
			if(flag==0){

				//Hopped from electrode
				if(JumpFromElec==1){
					//Determine which electrode hopped from
					if(PrevX == -1 && EndX!=0){
						printf("Hopping back into system current x %d\n",getCx(one));
						ChargeSystem( &nc, elXb, nca, Sequence, *chA);
					}else if(PrevY == -1 && EndY!=0){
						printf("Hopping back into system current y %d\n",getCy(one));
						ChargeSystem( &nc, elYl, nca, Sequence, *chA);
					}else if(PrevZ == -1 && EndZ!=0){
						printf("Hopping back into system current z %d\n",getCz(one));
						ChargeSystem( &nc, elZb, nca, Sequence, *chA);
					}
					//Hopped from site
				}else{
					setDwelStat(site,-1);
					setVis(site,1);
					incVisFreq(site);
				}
			}


			////////////////////////////////////////////////////////////////////////
			//Determining if Charge left system

			//Current position of charge
			xx = getCx(one);
			yy = getCy(one);
			zz = getCz(one);

			//Checik if a charge has arrived at the front or back electrode
			//If it has updates the following parameters:
			//	Velocities
			//	NumAvgVel
			//	Drain
			//	Source
			//	Number of charges on respective electrodes
			//	TimeTrack1 
			//	t
			//printf("Before CheckIfElecHop of t %Lg\n",t);
			//printf("TotalVelX %Lg\n",TotalVelX);

			CheckIfElecHop(xx, yy, zz,\
					&Xsource, &Ysource, &Zsource,\
					&Xdrain, &Ydrain, &Zdrain,\
					&TrackX, &TrackY, &TrackZ,\
					&TotalVelX, &TotalVelY, &TotalVelZ,\
					SLength,SWidth, SHeight,\
					EndX, EndY, EndZ,\
					&NumAvgVel, &TimeTrack1, t, SiteDistance,\
					&nc,&ElecExit,\
					PeriodicX, PeriodicY, PeriodicZ);	

			//printf("CheckIfElecHop of t %Lg\n",t);
			factorX = 0;
			factorY = 0;
			factorZ = 0;

			//Finding position of charge on sites
			if(xx<0){
				factorX = -xx/SLength+1;
				x = (xx+SLength*factorX)%SLength;
			}else{
				x = (xx+SLength)%SLength;
			}
			if(yy<0){
				factorY = -yy/SWidth+1;
				y = (yy+SWidth*factorY)%SWidth;
			}else{
				y = (yy+SWidth)%SWidth;
			}
			if(zz<0){
				factorZ = -zz/SHeight+1;
				z = (zz+SHeight*factorZ)%SHeight;
			}else{
				z = (zz+SHeight)%SHeight;
			}

			////////////////////////////////////////////////////////////////////////
			//Determining Next hop

			//Before future hop has been determined, tim of hop is 0
			tim=0;
			//printf("Middle Value of t %Lg\n",t);

			//Did not hop to electrode
			if(ElecExit==0){

				//Site jumped too
				site = getSN(snA, x,y,z);
				//Set the DwelStat of new site with the id of 
				//the charge that just jumped
				setDwelStat(site,ChargeID);

				//Check to see if new site is part of a cluster
				ClusterYes=getType(site);

				//Determine future hopping site
				if(ClusterYes==1){
					//If site is part of a cluster ensure that the site
					//is not next to an electrode
					ClusterHopCheck(PeriodicX, PeriodicY, PeriodicZ,\
							SLength, SWidth, SHeight,\
							EndX, EndY, EndZ,\
							x, y, z, site,\
							xx, yy, zz,\
							&CheckX, &CheckY, &CheckZ);

					if(CheckX==1 || CheckY==1 || CheckZ==1){
						//This means near the edge of a sample where 
						//there is an electrode cannot use cluster 
						//approximation

						if( x==0  && XElecOn==1 ){
							printf("Here hopping to elec\n");
							HopToElecX(snA, elXb, &one, &future, EndY, EndZ);
							getLoc(&x,&y,&z,future,snA);
							tim = 1/getsum(getSN(snA,x,y,z));
						}else if( y==0 && YElecOn==1){
							HopToElecY(snA, elYl, &one, &future, EndX, EndZ);
							getLoc(&x,&y,&z,future,snA);
							tim = 1/getsum(getSN(snA,x,y,z));
						}else if( z==0 && ZElecOn==1){
							HopToElecZ(snA, elZb, &one, &future, EndX, EndY);
							getLoc(&x,&y,&z,future,snA);
							tim = 1/getsum(getSN(snA,x,y,z));
						}else{
							//Future Hop To Cluster or Site using cluster algorithm
							ClusterHop(snA, &one, &tim, &future);
						}


					}else{

						//Future Hop To Cluster or Site
						SiteHop(snA, &one, site, &future, EndX, EndY,EndZ,\
								XElecOn,YElecOn,ZElecOn, PeriodicX,PeriodicY, PeriodicZ);
						getLoc(&x,&y,&z,future,snA);
						tim = 1/getsum(getSN(snA,x,y,z));
					}
				}else{
					//Future Hop To Cluster or Site
					SiteHop(snA, &one, site, &future, EndX, EndY,EndZ,\
							XElecOn,YElecOn,ZElecOn, PeriodicX,PeriodicY, PeriodicZ);
					getLoc(&x,&y,&z,future,snA);
					tim = 1/getsum(getSN(snA,x,y,z));

				}

				//Having chosen future site recording it
				setE(FutureSite,ChargeID+1,1,future);

				do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
				//The waiting time of the charge is updated based on
				//it's location and a random number
				//printf("Setting dwel site %lg\n",tim);
				setDwel(one, -log((double) ran/RAND_MAX)*tim);

				if(tim>1){
					printf("tim %g\n",tim);
					printf("WARNING setDwel huge\n");
				}

				//The location of the charge in the sequence is updated
				insertDwelltimePos(nca, *chA, &Sequence);

			}else{
				//If the charge jumped to an electrode it will now have to jump
				//from the electrode back into the system
				if(xx == -1 && EndX!=0){
					//Hopping from back Electrode
					//Updating dwell time of charge
					ElecHopOffX(elXb, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					printf("Hopped to X Source x position %d flag %d\n",xx, flag);
					ChargeElectrode(elXb, &one, &Sequence, *chA, &nc, nca, flag);
					//printf("ChargeElectrode value of t %Lg\n",t);
					//printf("Grabbing dwel %g\n",getDwel(one));

				}else if(xx == EndX*SLength && EndX!=0){
					//Hopped to front electrode
					//Setting dwell time to large number removing charge from system
					printf("Hopped to X Drain x position %d\n",xx);
					ChargeClosure(elXf, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);

				}else if(yy == -1 && EndY!=0){
					//Hopping from left Electrode
					ElecHopOffY(elYl, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					printf("Hopped to Y Source y position %d\n",yy);
					ChargeElectrode(elYl, &one, &Sequence, *chA, &nc, nca, flag);
				}else if(yy == EndY*SWidth && EndY!=0 ){
					//Hopped to right electrode
					printf("Hopped to Y Drain y position %d\n",yy);
					ChargeClosure(elYr, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);

				}else if(zz == -1 && EndZ!=0){
					//Hopping from bottom Electrode
					ElecHopOffZ(elZb, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					printf("Hopped to Z Source z position %d\n",zz);
					ChargeElectrode(elZb, &one, &Sequence, *chA, &nc, nca, flag);
				}else if(zz == EndZ*SHeight && EndZ!=0){
					//Hopped to top electrode
					printf("Hopped to Z Drain z position %d\n",zz);
					ChargeClosure(elZa, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);

				}
			}

		}	

		////////////////////////////////////////////////////////////////////////
		//Regardless of whether a charge hopped or not need to 
		//see if new charges need to be injected


		//printf("Bottom Value of t %Lg\n",t);
		if(t >= ((long double)n)*(long double)TStep && ((n+1)*(long int)NCh)<=Ntot){

			//increment n only up to certain point
			if( (long double)n*TStep<=t || ((n+1)*(long int)NCh)<Ntot){
				n++;
			}
			//Initialize new charges that are inserted
			initCharge( nca, n, chA, Sequence, snA,\
					Ntot, NCh, D,\
					XElecOn, YElecOn, ZElecOn,\
					EndX, EndY, EndZ);

			//Add to the number of charges that are 
			//already in the system

			TotalVelX += (long double)TrackX*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			TotalVelY += (long double)TrackY*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			TotalVelZ += (long double)TrackZ*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			TrackX = 0;
			TrackY = 0;
			TrackZ = 0;
			NumAvgVel++;
			TimeTrack1 = t;
			nc = nc+NCh;
			nca = nca+NCh;
		}

	}

	printf("Printing Transport Data\n");

	printf("Length of rows of System %d\n",getRows(System));
	printTransportData(System, timeArray, Xcurrent, Ycurrent, Zcurrent,\
			Xelec_Drain, Yelec_Drain, Zelec_Drain,\
			Xelec_Source, Yelec_Source, Zelec_Source,\
			Xvelocity, Yvelocity, Zvelocity,\
			XElecOn, YElecOn, ZElecOn, FileName,\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ);

	printf("Deleting matrices\n");
	deleteChargeA(*chA);
	deleteMatrix(timeArray);
	deleteMatrix(Xcurrent);
	deleteMatrix(Ycurrent);
	deleteMatrix(Zcurrent);
	deleteMatrix(Xelec_Drain);
	deleteMatrix(Yelec_Drain);
	deleteMatrix(Zelec_Drain);
	deleteMatrix(Xelec_Source);
	deleteMatrix(Yelec_Source);
	deleteMatrix(Zelec_Source);
	deleteMatrix(Xvelocity);
	deleteMatrix(Yvelocity);
	deleteMatrix(Zvelocity);
	deleteMatrix(System);
	deleteMatrix(Sequence);
	deleteMatrix(FutureSite);
	return 0;
}

int UpdateOccTime(SNarray * snA,Charge * ch, double time, ParameterFrame PF){

	if(snA==NULL || ch == NULL){
		return -1;
	}
	if(*snA==NULL || *ch==NULL){
		return -1;
	}

	int i;
	int j;
	int k;
	int XElecOn;
	int YElecOn;
	int ZElecOn;

	i = getCx(*ch);
	j = getCy(*ch);
	k = getCz(*ch);

	XElecOn = PFget_XElecOn(PF);
	YElecOn = PFget_YElecOn(PF);
	ZElecOn = PFget_ZElecOn(PF);

	if(XElecOn==0){
		i=(i+getAlen(*snA))%(getAlen(*snA));	
	}
	if(YElecOn==0){
		j=(j+getAwid(*snA))%(getAwid(*snA));
	}
	if(ZElecOn==0){
		k=(k+getAhei(*snA))%(getAhei(*snA));
	}

	if(i>=0 && j>=0 && k>=0){

		SiteNode sn = getSN(*snA,i,j,k);
		addTime(&sn,time);
	}

	return 0;
}

int CheckPt_exist(char * File, int buffersize){

	if(buffersize<0){
		printf("ERROR Buffer size less than 0!\n");
		return -1;
	}

	if(File[0]!='\0'){
		printf("ERROR File has not been initialized to \\0\n");
		return -1;
	}

	struct stat st;
	char FileName[256];
	char FileEx[20];
	char FileNameKeep[256];
	int len;
	int i;

	int CheckFileExist;
	//Determines if there is a checkpoint file
	//What is the highest number of the checkpointfile
	int num;
	int keep;

	//Initialize file to Null pointer
	File[0] = '\0';
	CheckFileExist = 0;

	if(stat("CHECKPOINT",&st)==0){
		//CHECKPOINT directory does exist
		DIR *d;
		struct dirent *dir;
		d = opendir("CHECKPOINT");
		if(d) {

			printf("CHECKPOINT directory exists!\n");
			num = 0;
			keep = 0;
			while((dir = readdir(d))!=NULL && CheckFileExist==0){

				strcpy(FileName, dir->d_name);
				//printf("%s\n",FileName);
				len = strlen(FileName);
				if(len>4){
					for(i=4;i>=0;i--){
						FileEx[4-i] = FileName[len-i];
					}
					FileEx[4] = '\0';
					if(strcmp(FileEx,"ckpt\0")==0){
						//It is a ckpt file
						CheckFileExist = 1;
						strcpy(FileNameKeep,FileName);
					}
				}
			}
			closedir(d);
		}

	}else{
		//CHECKPOINT directory does not exist
		return -1;
	}

	if(CheckFileExist==1){
		if(buffersize<len){
			printf("ERROR Buffersize is to small!\n");
			return -1;
		}
		printf("At least one .ckpt file was found!\n");
		printf("%s\n",FileNameKeep);
		strncpy(File,FileNameKeep,buffersize-1);
		return 0;

	}else{
		return -1;
	}
}

int CheckPt_Cluster(const double Vx,const double Vy,const double Vz,const double T, int r){

	if(T < 0){
		printf("ERROR Temperature negative\n");
		return -1;
	}

	struct stat st;
	char FileName[256];
	FILE *file;
	FileName[0] = '\0';
 	sprintf(FileName,"CHECKFILE/DataT%gVx%gVy%gVz%gR%d",T,Vx,Vy,Vz,r);

	if(stat("CLUSTERFILE",&st)==0){
		//CHECKPOINT Cluster directory does exist

		if((file=fopen(FileName,"r"))){
			fclose(file);
			printf("Cluster file exists\n");
			return 0;
		}else{
			printf("Cluster file does not exist\n");
			return -1;
		}

	}else{
		//CHECKPOINT Cluster does not exist
		printf("ERROR CHECKPOINT Cluster directory not found\n");
		return -1;
	}

}


int CheckPt_Latest(char * FileNameFull, int buffersize, const double Vx,const double Vy,const double Vz,const double T){

	if(buffersize<0){
		printf("ERROR Buffer size less than 0 Buffer: %d!\n",buffersize);
		return -1;
	}
	if(T < 0){
		printf("ERROR Temperature negative\n");
		return -1;
	}

	struct stat st;
	char FileName[256];
	char FileEx[20];
	char CheckNum[30];
	char CheckVx[30] = "";
	char CheckVy[30] = "";
	char CheckVz[30] = "";
	char CheckT[30] = "";

	double ChVx;
	double ChVy;
	double ChVz;
	double ChT;
	char FileNameKeep[256];
	int len;
	int stop1;
	int stop2;
	int stop3;
	int stop4;
	int stop5;
	int i;

	int CheckFileExist;
	//Determines if there is a checkpoint file
	//What is the highest number of the checkpointfile
	int num;
	int keep;

	//Initialize file to Null pointer
	FileName[0] = '\0';
	FileNameFull[0] = '\0';
	CheckFileExist = 0;

	ChVx = 0.0;
	ChVy = 0.0;
	ChVz = 0.0;
	ChT = 0.0;

	if(stat("CHECKPOINT",&st)==0){
		//CHECKPOINT directory does exist
		DIR *d;
		struct dirent *dir;
		d = opendir("CHECKPOINT");
		if(d) {

			num = 0;
			keep = 0;
			while((dir = readdir(d))!=NULL){

				strcpy(FileName, dir->d_name);
				//printf("%s\n",FileName);
				len = strlen(FileName);
				if(len>4){
					for(i=4;i>=0;i--){
						FileEx[4-i] = FileName[len-i];
					}
					FileEx[4] = '\0';
					if(strcmp(FileEx,"ckpt\0")==0){
						//It is a ckpt file

						//Need to make sure it is the right checkpoint file

						i = 0;
						while( FileName[len-5-i]!='R'){
							i++;
						}
						stop1 = i;
						for(i=0;i<(stop1-1);i++){
							CheckNum[i] = FileName[len-(5+stop1)+(i+1)];
						}

						CheckNum[stop1-1]='\0';
						sscanf(CheckNum,"%d",&num);

						i = 0;
						while( FileName[len-5-stop1-i]!='z'){
							i++;
						}
						stop2 = i;


						for(i=0;i<(stop2-1);i++){
							CheckVz[i] = FileName[len-(5+stop1+stop2)+(i+1)];
						}

						CheckVz[stop2-1]='\0';
						sscanf(CheckVz,"%lf",&ChVz);

						i = 0;
						while( FileName[len-(5+stop1+stop2+i)]!='y'){
							i++;
						}
						stop3 = i;
						
						for(i=0;i<(stop3-2);i++){

							CheckVy[i] = FileName[len-(5+stop1+stop2+stop3)+(i+1)];
						}

						CheckVz[stop3-1]='\0';
						sscanf(CheckVy,"%lf",&ChVy);

						i = 0;
						while( FileName[len-(5+stop1+stop2+stop3+i)]!='x'){
							i++;
						}
						stop4 = i;
						for(i=0;i<(stop4-2);i++){
							CheckVx[i] = FileName[len-(5+stop1+stop2+stop3+stop4)+(i+1)];
						}

						CheckVz[stop4-1]='\0';
						sscanf(CheckVx,"%lf",&ChVx);


						i = 0;
						while( FileName[len-(5+stop1+stop2+stop3+stop4+i)]!='T'){
							i++;
						}

						stop5 = i;
						for(i=0;i<(stop5-2);i++){
							CheckT[i] = FileName[len-(5+stop1+stop2+stop3+stop4+stop5)+(i+1)];
						}

						CheckVz[stop5-1]='\0';
						sscanf(CheckT,"%lf",&ChT);


						if(num>keep && ChT==T && ChVx==Vx && \
								ChVy==Vy && ChVz==Vz){

							CheckFileExist = 1;
							keep = num;
							strcpy(FileNameKeep,FileName);
						}
					}
				}
			}
			closedir(d);
		}

	}else{
		//CHECKPOINT directory does not exist
		printf("ERROR CHECKPOINT directory not found\n");
		return -1;
	}

	if(CheckFileExist==1){
		printf("FileNameKeep %s\n",FileNameKeep);
		if(buffersize<len){
			printf("Value of len %d value of buffersize %d\n",len,buffersize);
			printf("ERROR Buffersize is to small!\n");
			return -1;
		}
		strncpy(FileNameFull,FileNameKeep,buffersize-1);
		return keep;
	}else{
		return 0;
	}
}

int Load_CheckPt(long double * t, SNarray * snA, ChargeArray * chA, matrix * Sequence,\
		matrix * FutureSite, char * FileName, ParameterFrame *PF,long int *n, int * nc, int *nca,\
		int * Num_elXb, int * Num_elXf, int * Num_elYl, int * Num_elYr, int * Num_elZb,\
		int * Num_elZa){


	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	double KT;

	double electricEnergyX;
	double electricEnergyY;
	double electricEnergyZ;
	double electricFieldX;
	double electricFieldY;
	double electricFieldZ;

	double MarcusJ0;
	double MarcusCoeff;

	//int check;
	//unsigned int position;
	long int intvar0;
	int intvar;
	int intvar2;
	double doublevar;

	int Ntot;
	int i;
	int Max;
	int ii, jj, kk;
	int Xdist;
	double Energy;
	double dwellTime;
	double Time;
	int Occupancy;
	int VisFreq;
	int Vis;

	int Future;
	int ChargeID;

	struct stat st = {0};
	char bufRead[256];
	char bufTemp[256];
	char bufCheck[256];
	FILE * CheckIn;
	SiteNode sn;
	Charge ch;

	//Check if Directory exists
	if(stat("CHECKPOINT",&st)==-1){
		printf("CHECKPOINT directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufTemp, sizeof bufTemp,"%s%s","CHECKPOINT/",FileName);
		snprintf(bufCheck, sizeof bufCheck,"%s%s","CHECKPOINT/",FileName);
		//Now we know which .ckpt file is the last
		if((CheckIn=fopen(bufCheck,"r"))==NULL){
			printf("ERROR: unable to open chosen .ckpt file, there is a problem!\n");
			exit(1);
		}else{

			*PF = newParamFrame();

			//Scanning in Parameters
			fgets(bufRead,256,CheckIn);
			//fscanf(CheckIn,"%s",bufRead);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Len(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Wid(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Hei(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Px(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Py(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Pz(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_EndX(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_EndY(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_EndZ(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_XElecOn(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_YElecOn(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_ZElecOn(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphaxb(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphaxf(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphayl(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphayr(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphazb(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_alphaza(*PF,doublevar*1E9);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_RelativePerm(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_vX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_vY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_vZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_XFermiB(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_XFermiF(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_YFermiL(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_YFermiR(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_ZFermiB(*PF, doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_ZFermiA(*PF, doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VoltageX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VoltageY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VoltageZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_VStepX(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_VStepY(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_VStepZ(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_SiteDist(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_D(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_TCount(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_NCh(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TStep(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Nstep_av(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Time_check(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Rcount(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_CutOff(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_lambda(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_SeedProt(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Attempts(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_FracSeed(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_E0(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_sigma(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_FracTrap(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_Etrap(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_Tsigma(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_AttemptToHop(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TempStart(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_TempStep(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TempInc(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_reOrg(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_gamma(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_MovieFrames(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);

			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			//Scanning in SN information
			fscanf(CheckIn,"%lf %ld %d %d",&doublevar,&intvar0, &intvar, &intvar2);
			*t = (long double) doublevar;
			*n = intvar0;
			*nc = intvar;
			*nca = intvar2;
			printf("nc %d nca %d\n",*nc, *nca);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			//Scan in total number of sites
			fscanf(CheckIn,"%d",&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			i = 0;
			Max = intvar;

			*snA = newSNarray(PFget_Len(*PF),PFget_Wid(*PF),PFget_Hei(*PF));

			while(i<Max){


				fscanf(CheckIn,"%d\t%d\t%d\t%lf\t%d\t%d\t%d\t%lf",&ii,&jj,&kk,&Energy,&Occupancy,&VisFreq,&Vis,&Time);
				//printf("%d %d %d %f %d %d %d\n",ii,jj,kk,Energy,Occupancy, VisFreq,Vis);
				i++;

				sn = getSN(*snA,ii,jj,kk);	
				setEnergy(sn,Energy);
				setVisFreq(sn,VisFreq);
				setDwelStat(sn,Occupancy);
				setTime(sn,Time);
				setVis(sn,Vis);

			}

			//Calculating Hop Rates for Sites
			printf("vX %g vY %g vZ %g Len %d Wid %d Hei %d SiteDist %g\n",PFget_vX(*PF),PFget_vY(*PF),PFget_vZ(*PF),PFget_Len(*PF),PFget_Wid(*PF),PFget_Hei(*PF),PFget_SiteDist(*PF));
			electricFieldX = PFget_VoltageX(*PF)/((double)PFget_Len(*PF)*PFget_SiteDist(*PF));
			electricFieldY = PFget_VoltageY(*PF)/((double)PFget_Wid(*PF)*PFget_SiteDist(*PF));
			electricFieldZ = PFget_VoltageZ(*PF)/((double)PFget_Hei(*PF)*PFget_SiteDist(*PF));
			electricEnergyX = PFget_SiteDist(*PF)*electricFieldX;
			electricEnergyY = PFget_SiteDist(*PF)*electricFieldY;
			electricEnergyZ = PFget_SiteDist(*PF)*electricFieldZ;

			printf("electricFieldX %g ElectricEnergyX %g\n",electricFieldX, electricEnergyX);
			KT = kB*PFget_TempStart(*PF);

			MarcusJ0 = pow( PFget_AttemptToHop(*PF)*hbar*pow(4*PFget_reOrg(*PF)*kB*300/M_PI,1/2),1/2);
			MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*PFget_reOrg(*PF)*KT),1/2)*exp(-2*PFget_gamma(*PF)*PFget_SiteDist(*PF));

			printf("MarcusJ0 %g MarcusCoeff %g\n",MarcusJ0,MarcusCoeff);

			initJumPossibility(electricEnergyX, electricEnergyY, electricEnergyZ,\
					MarcusCoeff, KT,PFget_reOrg(*PF), *snA,\
					PFget_Px(*PF), PFget_Py(*PF), PFget_Pz(*PF),\
					PFget_XElecOn(*PF), PFget_YElecOn(*PF),PFget_ZElecOn(*PF));

			printf("past init\n");


			//Now we are going to scan in the charges
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%d",&intvar);
			Ntot = intvar;
			printf("Ntot %d\n",Ntot);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			i = 0;

			printf("Creating New Matrix\n");
			*chA = newChargeA(Ntot);
			*Sequence = newMatrix(Ntot,1);		
			*FutureSite = newMatrix(Ntot,1);
			printf("Created New Matrix\n");

			while(i<Ntot){


				fscanf(CheckIn,"%d\t%d\t%d\t%d\t%d\t%lf\t%d",&ii,&jj,&kk,&ChargeID,&Xdist,&dwellTime,&Future);
				printf("%d %d %d %d %d %g %d\n",ii,jj,kk,ChargeID,Xdist,dwellTime,Future);
				ch = getCharge(*chA,i);
				setCx(ch,ii);
				setCy(ch,jj);
				setCz(ch,kk);
				setDwel(ch,dwellTime);
				i++;
				setE(*Sequence,i,1,ChargeID);
				intvar = setE(*FutureSite,ChargeID+1,1,Future);
			}

			printf("Last value of rv %d\n",intvar);
			printf("Load_CheckPt Printing FutureSite");
			//printMatrix(**FutureSite);

			//Read in the number of Charges on an electrode
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elXb = intvar;
			printf("Xb %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elXf = intvar;
			printf("Xf %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elYl = intvar;
			printf("Yl %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elYr = intvar;
			printf("Yr %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elZb = intvar;
			printf("Zb %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elZa = intvar;
			printf("Za %d\n",intvar);


			fclose(CheckIn);
		}

		//}



}

return 0;
}

int Load_CheckPt_Data(long double * t, SNarray * snA, ChargeArray * chA, matrix * Sequence,\
		matrix * FutureSite, char * FileName,long int * n,  int * nc, int *nca,\
		int * Num_elXb, int * Num_elXf, int * Num_elYl, int * Num_elYr, int * Num_elZb,\
		int * Num_elZa, double VoltageX, double VoltageY, double VoltageZ){


	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	double KT;

	double electricEnergyX;
	double electricEnergyY;
	double electricEnergyZ;
	double electricFieldX;
	double electricFieldY;
	double electricFieldZ;

	double MarcusJ0;
	double MarcusCoeff;

	//int check;
	//unsigned int position;
	long int intvar0;
	int intvar;
	int intvar2;
	double doublevar;

	int SLength;
	int SWidth;
	int SHeight;
	int PeriodicX;
	int PeriodicY;
	int PeriodicZ;
	int XElecOn;
	int YElecOn;
	int ZElecOn;
	double SiteDistance;
	double AttemptToHop;
	double reOrg;
	double TempStart;
	double gamma;
	int MovieFrames;
	int Ntot;
	int i;
	int Max;
	int ii, jj, kk;
	int Xdist;
	double Energy;
	double dwellTime;
	double Time;
	int Occupancy;
	int VisFreq;
	int Vis;

	int Future;
	int ChargeID;

	struct stat st = {0};
	char bufRead[256];
	char bufTemp[256];
	char bufCheck[256];
	FILE * CheckIn;
	SiteNode sn;
	Charge ch;

	//Check if Directory exists
	if(stat("CHECKPOINT",&st)==-1){
		printf("CHECKPOINT directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufTemp, sizeof bufTemp,"%s%s","CHECKPOINT/",FileName);
		snprintf(bufCheck, sizeof bufCheck,"%s%s","CHECKPOINT/",FileName);
		//Now we know which .ckpt file is the last
		if((CheckIn=fopen(bufCheck,"r"))==NULL){
			printf("ERROR: unable to open chosen .ckpt file, there is a problem!\n");
			exit(1);
		}else{

			//Scanning in Parameters
			fgets(bufRead,256,CheckIn);
			//fscanf(CheckIn,"%s",bufRead);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			SLength = intvar;
			printf("SLength %d\n",SLength);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			SWidth = intvar;
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			SHeight = intvar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PeriodicX = intvar;
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PeriodicY = intvar;
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PeriodicZ = intvar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			XElecOn = intvar;
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			YElecOn = intvar;
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			ZElecOn = intvar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			SiteDistance = doublevar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			AttemptToHop = doublevar;
			printf("AttemptToHop %g\n",AttemptToHop);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			TempStart = doublevar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			reOrg = doublevar;
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			gamma = doublevar*1E9;
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			MovieFrames = intvar;
			printf("MovieFrames %d\n",MovieFrames);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			//Scanning in SN information
			fscanf(CheckIn,"%lf %ld %d %d",&doublevar,&intvar0, &intvar, &intvar2);
			*t = (long double) doublevar;
			*n = intvar0;
			*nc = intvar;
			*nca = intvar2;
			printf("nc %d nca %d\n",*nc, *nca);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			//Scan in total number of sites
			fscanf(CheckIn,"%d",&intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			i = 0;
			Max = intvar;

			while(i<Max){


				fscanf(CheckIn,"%d\t%d\t%d\t%lf\t%d\t%d\t%d\t%lf",&ii,&jj,&kk,&Energy,&Occupancy,&VisFreq,&Vis,&Time);
				printf("%d %d %d %f %d %d %d %g\n",ii,jj,kk,Energy,Occupancy, VisFreq,Vis, Time);
				i++;

				sn = getSN(*snA,ii,jj,kk);	
				setEnergy(sn,Energy);
				setTime(sn,Time);
				setVisFreq(sn,VisFreq);
				setDwelStat(sn,Occupancy);
				setVis(sn,Vis);

			}

			//Calculating Hop Rates for Sites
			electricFieldX = VoltageX/((double)SLength*SiteDistance);
			electricFieldY = VoltageY/((double)SWidth*SiteDistance);
			electricFieldZ = VoltageZ/((double)SHeight*SiteDistance);
			electricEnergyX = SiteDistance*electricFieldX;
			electricEnergyY = SiteDistance*electricFieldY;
			electricEnergyZ = SiteDistance*electricFieldZ;

			printf("electricFieldX %g ElectricEnergyX %g\n",electricFieldX, electricEnergyX);
			KT = kB*TempStart;

			MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrg*kB*300/M_PI,1/2),1/2);
			MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrg*KT),1/2)*exp(-2*gamma*SiteDistance);

			printf("MarcusJ0 %g MarcusCoeff %g\n",MarcusJ0,MarcusCoeff);

			initJumPossibility(electricEnergyX, electricEnergyY, electricEnergyZ,\
					MarcusCoeff, KT,reOrg, *snA,\
					PeriodicX, PeriodicY, PeriodicZ,\
					XElecOn, YElecOn,ZElecOn);

			printf("past init\n");


			//Now we are going to scan in the charges
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%d",&intvar);
			Ntot = intvar;
			printf("Ntot %d\n",Ntot);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			i = 0;

			*chA = newChargeA(Ntot);


			while(i<Ntot){


				fscanf(CheckIn,"%d\t%d\t%d\t%d\t%d\t%lf\t%d",&ii,&jj,&kk,&ChargeID,&Xdist,&dwellTime,&Future);
				printf("%d %d %d %d %d %g %d\n",ii,jj,kk,ChargeID,Xdist, dwellTime,Future);
				ch = getCharge(*chA,i);
				setCx(ch,ii);
				setCy(ch,jj);
				setCz(ch,kk);
				setDwel(ch,dwellTime);
				i++;
				setE(*Sequence,i,1,ChargeID);
				intvar = setE(*FutureSite,ChargeID+1,1,Future);

			}

			printf("Last value of rv %d\n",intvar);
			printf("Load_CheckPt Printing FutureSite");
			//printMatrix(**FutureSite);

			//Read in the number of Charges on an electrode
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elXb = intvar;
			printf("Xb %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elXf = intvar;
			printf("Xf %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elYl = intvar;
			printf("Yl %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elYr = intvar;
			printf("Yr %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elZb = intvar;
			printf("Zb %d\n",intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			*Num_elZa = intvar;
			printf("Za %d\n",intvar);


			fclose(CheckIn);
		}

		//}



}

return 0;
}


int Load_CheckPt_PF( char * FileName, ParameterFrame *PF ){

	if(FileName[0]=='\0'){
		printf("ERROR FileName is \\0\n");
		return -1;
	}
	if(*PF!=NULL){
		printf("ERROR PF was not initialized to NULL\n");
		return -1;
	}
	//int check;
	//unsigned int position;
	int intvar;
	double doublevar;

	struct stat st = {0};
	char bufRead[256];
	char bufTemp[256];
	char bufCheck[256];
	FILE * CheckIn;

	//Check if Directory exists
	if(stat("CHECKPOINT",&st)==-1){
		printf("CHECKPOINT directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufTemp, sizeof bufTemp,"%s%s","CHECKPOINT/",FileName);
		snprintf(bufCheck, sizeof bufCheck,"%s%s","CHECKPOINT/",FileName);

		//Now we know which .ckpt file is the last
		if((CheckIn=fopen(bufCheck,"r"))==NULL){
			printf("ERROR: unable to open chosen .ckpt file, there is a problem!\n");
			printf("File %s\n",bufCheck);
			exit(1);
		}else{

			printf("Creating and loading parameter frame data\n");
			*PF = newParamFrame();

			//Scanning in Parameters
			fgets(bufRead,256,CheckIn);
			//fscanf(CheckIn,"%s",bufRead);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("SLength %d\n",intvar);
			PFset_Len(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("SWidth %d\n",intvar);
			PFset_Wid(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("SHeight %d\n",intvar);
			PFset_Hei(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("PeriodicX %d\n",intvar);
			PFset_Px(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("PeriodicY %d\n",intvar);
			PFset_Py(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("PeriodicZ %d\n",intvar);
			PFset_Pz(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("EndX %d\n",intvar);
			PFset_EndX(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("EndY %d\n",intvar);
			PFset_EndY(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_EndZ(*PF,intvar);
			printf("EndZ %d\n",intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("XElecOn %d\n",intvar);
			PFset_XElecOn(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("YElecOn %d\n",intvar);
			PFset_YElecOn(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			printf("ZElecOn %d\n",intvar);
			PFset_ZElecOn(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphaxb %g\n",doublevar);
			PFset_alphaxb(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphaxf %g\n",doublevar);
			PFset_alphaxf(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphayl %g\n",doublevar);
			PFset_alphayl(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphayr %g\n",doublevar);
			PFset_alphayr(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphazb %g\n",doublevar);
			PFset_alphazb(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("alphaza %g\n",doublevar);
			PFset_alphaza(*PF,doublevar*1E9);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("RelativePerm %g\n",doublevar);
			PFset_RelativePerm(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("vX %g\n",doublevar);
			PFset_vX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("vY %g\n",doublevar);
			PFset_vY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("vZ %g\n",doublevar);
			PFset_vZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("XFermiB %g\n",doublevar);
			PFset_XFermiB(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("XFermiF %g\n",doublevar);
			PFset_XFermiF(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("YFermiL %g\n",doublevar);
			PFset_YFermiL(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("YFermiR %g\n",doublevar);
			PFset_YFermiR(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("ZFermiB %g\n",doublevar);
			PFset_ZFermiB(*PF, doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("ZFermiA %g\n",doublevar);
			PFset_ZFermiA(*PF, doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("VoltageX %g\n",doublevar);
			PFset_VoltageX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("VoltageY %g\n",doublevar);
			PFset_VoltageY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			printf("VoltageZ %g\n",doublevar);
			PFset_VoltageZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("VStepX %d\n",intvar);
			PFset_VStepX(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("VStepY %d\n",intvar);
			PFset_VStepY(*PF,intvar);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			printf("VStepZ %d\n",intvar);
			PFset_VStepZ(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincX(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincY(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);	
			PFset_VincZ(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_SiteDist(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_D(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_TCount(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_NCh(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TStep(*PF,doublevar);

			PFset_Ntot(*PF,PFget_TCount(*PF)*PFget_NCh(*PF));

			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Nstep_av(*PF,intvar);

			PFset_N_av(*PF,(int) PFget_TCount(*PF)/ PFget_Nstep_av(*PF));

			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Time_check(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Rcount(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_CutOff(*PF,doublevar);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_lambda(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_SeedProt(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_Attempts(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_FracSeed(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_E0(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_sigma(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_FracTrap(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_Etrap(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_Tsigma(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_AttemptToHop(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TempStart(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_TempStep(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_TempInc(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_reOrg(*PF,doublevar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fscanf(CheckIn,"%s %lf",bufRead,&doublevar);
			PFset_gamma(*PF,doublevar*1E9);
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_MovieFrames(*PF,intvar);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			fclose(CheckIn);
		}

	}

	return 0;
}

int Save_CheckPt(char * FileName, int * CheckptNum, SNarray snA,\
		ChargeArray chA, matrix Sequence, matrix FutureSite, long double t,\
		ParameterFrame PF,long int n, int nc, int nca, Electrode elXb, Electrode elXf,\
		Electrode elYl, Electrode elYr, Electrode elZb, Electrode elZa){

	if(snA==NULL || FutureSite==NULL || Sequence==NULL || PF==NULL || chA==NULL){
		return -1;
	}

	int i;
	int j;
	int k;
	int element;
	int ID;

	SiteNode sn;
	Charge ch;

	struct stat st = {0};
	char str[20];
	char bufCheck[256];
	FILE * CheckOut;

	//Check if Directory exists
	if(stat("CHECKPOINT",&st)==-1){
		mkdir("CHECKPOINT",0700);
	}

	sprintf(str,"%d",*CheckptNum);
	snprintf(bufCheck, sizeof bufCheck,"%s%s%s%s","CHECKPOINT/",FileName,str,".ckpt");

	*CheckptNum = *CheckptNum+1;

	if((CheckOut=fopen(bufCheck,"w"))==NULL){
		printf("Error! unable to write to .ckpt file!\n");
	}else{

		//Print parameter data first
		fprintf(CheckOut,"//Number of nodes along the x-axis, y-axis and z-axis\n");
		fprintf(CheckOut,"SLength %d\n",getAlen(snA));
		fprintf(CheckOut,"SWidth %d\n",getAwid(snA));
		fprintf(CheckOut,"SHeight %d\n\n",getAhei(snA));

		fprintf(CheckOut,"//Defines whether sides are periodic or not\n");
		fprintf(CheckOut,"// 0 - means non-periodic\n");
		fprintf(CheckOut,"// 1 - means periodic\n");
		fprintf(CheckOut,"PeriodicX %d\n",PFget_Px(PF));
		fprintf(CheckOut,"PeriodicY %d\n",PFget_Py(PF));
		fprintf(CheckOut,"PeriodicZ %d\n\n",PFget_Pz(PF));

		fprintf(CheckOut,"//Defines over how many sample the system is periodic\n");
		fprintf(CheckOut,"// 0 - means it is an infinite sample\n");
		fprintf(CheckOut,"// 1 - means periodic across 1 sample (not-periodic)\n");
		fprintf(CheckOut,"// 2 - periodic across 2 samples then reaches the edge\n");
		fprintf(CheckOut,"//     3, 4, 5, ... etc\n");
		fprintf(CheckOut,"EndX %d\n",PFget_EndX(PF));
		fprintf(CheckOut,"EndY %d\n",PFget_EndY(PF));
		fprintf(CheckOut,"EndZ %d\n\n",PFget_EndZ(PF));

		fprintf(CheckOut,"//Defines whether electrodes exist on the x,y or z axis\n");
		fprintf(CheckOut,"// 0 - no electrodes\n");
		fprintf(CheckOut,"// 1 - yes electrodes\n");
		fprintf(CheckOut,"XElecOn %d\n",PFget_XElecOn(PF));
		fprintf(CheckOut,"YElecOn %d\n",PFget_YElecOn(PF));
		fprintf(CheckOut,"ZElecOn %d\n\n",PFget_ZElecOn(PF));

		fprintf(CheckOut,"//Define the tunneling constant from the electrode [1/nm]\n");
		fprintf(CheckOut,"alphaXb %f\n",PFget_alphaxb(PF)/1E9);
		fprintf(CheckOut,"alphaXf %f\n",PFget_alphaxf(PF)/1E9);
		fprintf(CheckOut,"alphaYl %f\n",PFget_alphayl(PF)/1E9);
		fprintf(CheckOut,"alphaYr %f\n",PFget_alphayr(PF)/1E9);
		fprintf(CheckOut,"alphaZb %f\n",PFget_alphazb(PF)/1E9);
		fprintf(CheckOut,"alphaZa %f\n\n",PFget_alphaza(PF)/1E9);
		/*
			 fprintf(CheckOut,"//Define the work function of the electrodes [eV]\n");
			 fprintf(CheckOut,"workX %f\n",PFget_);
			 fprintf(CheckOut,"workY %f\n",PFget_);
			 fprintf(CheckOut,"workZ %f\n\n",PFget_);
			 fprintf(CheckOut,"//Define the electron affinity of medium [eV]\n");
			 fprintf(CheckOut,"electronAffin %f\n",electronAffin);
			 fprintf(CheckOut,"//Define the ionization energy of medium [eV]\n");
			 fprintf(CheckOut,"ionizationEnergy %f\n\n",ionizationEnergy);
		 */
		fprintf(CheckOut,"//Relative permittivity of medium\n");
		fprintf(CheckOut,"RelativePermittivity %f\n\n",PFget_RelativePerm(PF));

		fprintf(CheckOut,"//Define the attempt to hop rate from the electrode [1/s]\n");
		fprintf(CheckOut,"vX %g\n",PFget_vX(PF));
		fprintf(CheckOut,"vY %g\n",PFget_vY(PF));
		fprintf(CheckOut,"vZ %g\n\n",PFget_vZ(PF));

		fprintf(CheckOut,"//Define Fermi Energy level of the electrodes so we\n");
		fprintf(CheckOut,"//can calculate the rates on and off electrodes [eV]\n");
		fprintf(CheckOut,"//If the Electrodes are not turned on this value will\n");
		fprintf(CheckOut,"//not be used\n");
		fprintf(CheckOut,"XFermiB %f\n",PFget_XFermiB(PF));
		fprintf(CheckOut,"XFermiF %f\n",PFget_XFermiF(PF));
		fprintf(CheckOut,"YFermiL %f\n",PFget_YFermiL(PF));
		fprintf(CheckOut,"YFermiR %f\n",PFget_YFermiR(PF));
		fprintf(CheckOut,"ZFermiB %f\n",PFget_ZFermiB(PF));
		fprintf(CheckOut,"ZFermiA %f\n\n",PFget_ZFermiA(PF));

		fprintf(CheckOut,"//Voltage across a single system in x, y and z [V]\n");
		fprintf(CheckOut,"VoltageX %g\n",PFget_VoltageX(PF));
		fprintf(CheckOut,"VoltageY %g\n",PFget_VoltageY(PF));
		fprintf(CheckOut,"VoltageZ %g\n\n",PFget_VoltageZ(PF));

		fprintf(CheckOut,"//Voltage ramp steps x, y and z\n");\
			fprintf(CheckOut,"VStepX %d\n",PFget_VStepX(PF));
		fprintf(CheckOut,"VStepY %d\n",PFget_VStepY(PF));
		fprintf(CheckOut,"VStepZ %d\n\n",PFget_VStepZ(PF));

		fprintf(CheckOut,"//Voltage ramp increment x, y and z [V]\n");
		fprintf(CheckOut,"VincX %g\n",PFget_VincX(PF));
		fprintf(CheckOut,"VincY %g\n",PFget_VincY(PF));
		fprintf(CheckOut,"VincZ %g\n\n",PFget_VincZ(PF));

		fprintf(CheckOut,"//Distance between nodes in units of [m]\n");
		fprintf(CheckOut,"SiteDistance %g\n\n",PFget_SiteDist(PF));

		fprintf(CheckOut,"//Dimenion to fit real data\n");
		fprintf(CheckOut,"D %g\n\n",PFget_D(PF));

		fprintf(CheckOut,"//Total number of time steps in which charges are injected\n");
		fprintf(CheckOut,"TCount %d\n\n",PFget_TCount(PF));

		fprintf(CheckOut,"//Number of charges injected per time step\n");
		fprintf(CheckOut,"NCh %d\n\n",PFget_NCh(PF));

		fprintf(CheckOut,"//Total number of charges passed through the system\n");
		fprintf(CheckOut,"Ntot TCount*NCh\n\n");

		fprintf(CheckOut,"//Time between injections of charges [s]\n");
		fprintf(CheckOut,"TStep %g\n\n",PFget_TStep(PF));

		fprintf(CheckOut,"//Number of time steps that pass before data is averaged\n");
		fprintf(CheckOut,"Nstep_av %d\n\n",PFget_Nstep_av(PF));

		fprintf(CheckOut,"//Number of time steps before a checkpoint file is created\n");
		fprintf(CheckOut,"Time_check %d\n\n",PFget_Time_check(PF));

		fprintf(CheckOut,"//Number of iterations with different random seeds\n");
		fprintf(CheckOut,"Rcount %d\n\n",PFget_Rcount(PF));

		fprintf(CheckOut,"//Defines the cutoff radius for the correlation function\n");
		fprintf(CheckOut,"//CutOffDistance CutOff*lambda\n");
		fprintf(CheckOut,"//lambda is used in the correlation function\n");
		fprintf(CheckOut,"CutOff %f\n",PFget_CutOff(PF));
		fprintf(CheckOut,"lambda %g\n\n",PFget_lambda(PF));

		fprintf(CheckOut,"//Seed Protocol is used to determine how the energies are\n");
		fprintf(CheckOut,"//spead out between the sites.\n");
		fprintf(CheckOut,"// 0 - means averaged by surrounding seeds\n");
		fprintf(CheckOut,"// 1 - means averaged with closest seeds\n");
		fprintf(CheckOut,"SeedProt %d\n",PFget_SeedProt(PF));
		fprintf(CheckOut,"//How many numerical iterations should be used to approximate\n");
		fprintf(CheckOut,"//cluster behavior (Default to 15%% should be less than 5%% error)\n");
		fprintf(CheckOut,"Attempts %d\n\n",PFget_Attempts(PF));

		fprintf(CheckOut,"//What fraction of sites act as seeds\n");
		fprintf(CheckOut,"fracSeed %f\n",PFget_FracSeed(PF));
		fprintf(CheckOut,"//Average site energy [eV]\n");
		fprintf(CheckOut,"E0 %f\n",PFget_E0(PF));
		fprintf(CheckOut,"//Standard deviation of the site energy\n");
		fprintf(CheckOut,"sigma %f\n\n",PFget_sigma(PF));

		fprintf(CheckOut,"//What fraction of sites act as traps\n");
		fprintf(CheckOut,"fracTrap %f\n",PFget_FracTrap(PF));
		fprintf(CheckOut,"//Average site energy of trap [eV]\n");
		fprintf(CheckOut,"Etrap %f\n",PFget_Etrap(PF));
		fprintf(CheckOut,"//Standard deviation of trap energy\n");
		fprintf(CheckOut,"Tsigma %f\n\n",PFget_Tsigma(PF));

		fprintf(CheckOut,"//Attempt to hop rate [1/s] equates to Markus at 300K\n");
		fprintf(CheckOut,"AttemptToHop %g\n",PFget_AttemptToHop(PF));
		fprintf(CheckOut,"//Temperature [k]\n");
		fprintf(CheckOut,"TempStart %f\n",PFget_TempStart(PF));
		fprintf(CheckOut,"//Temperature steps\n");
		fprintf(CheckOut,"TemperatureStep %d\n",PFget_TempStep(PF));
		fprintf(CheckOut,"//Temperature increment [K]\n");
		fprintf(CheckOut,"TemperatureInc %f\n",PFget_TempInc(PF));
		fprintf(CheckOut,"//Re-organization energy [eV]\n");
		fprintf(CheckOut,"reOrgEnergy %f\n",PFget_reOrg(PF));
		fprintf(CheckOut,"//Tunneling constant [1/nm]\n");
		fprintf(CheckOut,"gamma %f\n",PFget_gamma(PF)/1E9);
		fprintf(CheckOut,"MovieFrames %d\n\n",PFget_MovieFrames(PF));
		fprintf(CheckOut,"Global Time\tTime Increments\tCharges in System\tActive Charges\n%Lg\t%ld\t%d\t%d\n\n",t,n,nc,nca);

		//Printing sites and energies
		fprintf(CheckOut,"%d\n\n",getAtotal(snA));
		for(i=0;i<getAlen(snA);i++){
			for(j=0;j<getAwid(snA);j++){
				for(k=0;k<getAhei(snA);k++){
					sn = getSN(snA,i,j,k);
					fprintf(CheckOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn));

				}
			}
		}

		//Printing charges
		fprintf(CheckOut,"\n%d\n\n",getChargeA_len(chA));

		for(element=0;element<getChargeA_len(chA);element++){
			ID = getE(Sequence,element+1,1);
			ch = getCharge(chA,element);
			//x y z ChargeID Xdist Dwelltime FutureSite
			fprintf(CheckOut,"%d\t%d\t%d\t%d\t%d\t%g\t%d\n",getCx(ch),getCy(ch),getCz(ch),\
					(int)ID,getXdist(ch),getDwel(ch),(int)getE(FutureSite,(int)ID+1,1));
		}

		//Now we need to print the Electrode information 
		fprintf(CheckOut,"\nElectrodes\n");
		if(PFget_XElecOn(PF)==1){
			fprintf(CheckOut,"elXb %d\n",getElectrode_Charges(elXb));
			fprintf(CheckOut,"elXf %d\n",getElectrode_Charges(elXf));
		}else{
			fprintf(CheckOut,"elXb 0\n");
			fprintf(CheckOut,"elXf 0\n");
		}
		if(PFget_YElecOn(PF)==1){
			fprintf(CheckOut,"elYl %d\n",getElectrode_Charges(elYl));
			fprintf(CheckOut,"elYr %d\n",getElectrode_Charges(elYr));
		}else{
			fprintf(CheckOut,"elYl 0\n");
			fprintf(CheckOut,"elYr 0\n");
		}
		if(PFget_ZElecOn(PF)==1){
			fprintf(CheckOut,"elZb %d\n",getElectrode_Charges(elZb));
			fprintf(CheckOut,"elZa %d\n",getElectrode_Charges(elZa));
		}else{
			fprintf(CheckOut,"elZb 0\n");
			fprintf(CheckOut,"elZa 0\n");
		}

	}

	fclose(CheckOut);

	return 0;
}


int printMovie(int * Movie, char * FileName, SNarray snA){

	int i;
	int j;
	int k;
	SiteNode sn;
	struct stat st = {0};
	char str[20];
	char bufM[256];
	FILE * Mout;

	if(stat("Movies",&st)==-1){
		mkdir("Movies",0700);
	}

	sprintf(str,"%d",*Movie);
	snprintf(bufM, sizeof bufM,"%s%s%s%s","Movies/",FileName,str,".xyz");

	*Movie = *Movie+1;

	if((Mout=fopen(bufM,"w"))==NULL){
		printf("Error! unable to write to .xyz movie file!\n");
	}else{

		fprintf(Mout,"%d\n\n",getAtotal(snA));

		for(i=0;i<getAlen(snA);i++){
			for(j=0;j<getAwid(snA);j++){
				for(k=0;k<getAhei(snA);k++){
					sn=getSN(snA,i,j,k);
					fprintf(Mout,"C \t %f \t %f \t %f \t %f \t %g\n",(double)i,(double)j,(double)k,(double)getDwelStat(sn),getTime(sn));
				}
			}
		}

	}

	fclose(Mout);

	return 0;
}


int printTransportData( matrix System, matrix timeArray, matrix Xcurrent, matrix Ycurrent, matrix Zcurrent,\
		matrix Xelec_Drain, matrix Yelec_Drain, matrix Zelec_Drain,\
		matrix Xelec_Source, matrix Yelec_Source, matrix Zelec_Source,\
		matrix Xvelocity, matrix Yvelocity, matrix Zvelocity,\
		int XElecOn, int YElecOn, int ZElecOn, char * FileName,\
		double ElectricFieldX, double ElectricFieldY, double ElectricFieldZ){

	if(timeArray==NULL || Xcurrent==NULL || Ycurrent==NULL || Zcurrent==NULL ||\
			Xelec_Drain==NULL || Yelec_Drain==NULL || Zelec_Drain==NULL ||\
			Xelec_Source==NULL || Yelec_Source==NULL || Zelec_Source==NULL){
		return -1;
	}
	int i;

	char bufx[256];
	char bufy[256];
	char bufz[256];

	snprintf(bufx, sizeof bufx,"%s%s",FileName,"X.txt");
	snprintf(bufy, sizeof bufy,"%s%s",FileName,"Y.txt");
	snprintf(bufz, sizeof bufz,"%s%s",FileName,"Z.txt");

	FILE * Xout;
	FILE * Yout;
	FILE * Zout;

	if((Xout=fopen(bufx,"a"))==NULL){
		printf("Error! unable to open X.txt\n");
	}else{
		if(XElecOn == 1){
			printf("Number of Rows of timeArray %d\n",getRows(timeArray));
			printf("Number of Rows of Xcurrent %d\n",getRows(Xcurrent));
			printf("Number of Rows of Xvelocity %d\n",getRows(Xvelocity));
			printf("Number of Rows of System %d\n",getRows(System));
			printf("Number of Rows of Xelec_Drain %d\n",getRows(Xelec_Drain));
			printf("Number of Rows of Xelec_Source %d\n",getRows(Xelec_Source));

			for(i=1;i<=getRows(timeArray);i++){
				if(getE(timeArray,i,1)!=0){

					if(ElectricFieldX!=0){
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Xout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),\
								getE(Xelec_Source,i,1),getE(Xelec_Drain,i,1),getE(System,i,1),getE(Xvelocity,i,1),getE(Xvelocity,i,1)/ElectricFieldX*10000);
					}else{
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] 
						fprintf(Xout,"%g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),\
								getE(Xelec_Source,i,1),getE(Xelec_Drain,i,1),getE(System,i,1),getE(Xvelocity,i,1));
					}
				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0){
					if(ElectricFieldX!=0){
						//Time [s] Xcurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Xout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),getE(Xvelocity,i,1),getE(Xvelocity,i,1)/ElectricFieldX*1000);
					}else{
						//Time [s] Xcurrent [Amps]  DriftVelocity [m/s]
						fprintf(Xout,"%g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),getE(Xvelocity,i,1));
					}
				}
			}
		}
		fclose(Xout);
	}

	if((Yout=fopen(bufy,"a"))==NULL){
		printf("Error! unable to open Y.txt\n");
	}else{
		if(YElecOn == 1){
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0){
					if(ElectricFieldY!=0){
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Yout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),\
								getE(Yelec_Source,i,1),getE(Yelec_Drain,i,1),getE(System,i,1),getE(Yvelocity,i,1),getE(Yvelocity,i,1)/ElectricFieldY*1000);
					}else{
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s]
						fprintf(Yout,"%g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),\
								getE(Yelec_Source,i,1),getE(Yelec_Drain,i,1),getE(System,i,1),getE(Yvelocity,i,1));
					}

				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0){
					if(ElectricFieldY!=0){
						//Time [s] Ycurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Yout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),getE(Yvelocity,i,1),getE(Yvelocity,i,1)/ElectricFieldY*1000);
					}else{
						//Time [s] Ycurrent [Amps]  DriftVelocity [m/s]
						fprintf(Yout,"%g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),getE(Yvelocity,i,1));
					}
				}
			}
		}
		fclose(Yout);
	}

	if((Zout=fopen(bufz,"a"))==NULL){
		printf("Error! unable to open Z.txt\n");
	}else{
		if(ZElecOn == 1){
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0){
					if(ElectricFieldZ!=0){
						//Time [s] Zcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Zout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),\
								getE(Zelec_Source,i,1),getE(Zelec_Drain,i,1),getE(System,i,1),getE(Zvelocity,i,1),getE(Zvelocity,i,1)/ElectricFieldZ*1000);
					}else{
						//Time [s] Zcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] 
						fprintf(Zout,"%g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),\
								getE(Zelec_Source,i,1),getE(Zelec_Drain,i,1),getE(System,i,1),getE(Zvelocity,i,1));
					}
				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0){
					if(ElectricFieldZ!=0){
						//Time [s] Zcurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Zout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),getE(Zvelocity,i,1),getE(Zvelocity,i,1)/ElectricFieldZ*1000);
					}else{
						//Time [s] Zcurrent [Amps]  DriftVelocity [m/s]
						fprintf(Zout,"%g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),getE(Zvelocity,i,1));
					}
				}
			}
		}
		printf("Closing Z file identifier\n");
		fclose(Zout);
	}

	printf("Returning with a value of 0\n");
	return 0;
}

int ChargeClosure(Electrode el, Charge * one, matrix * Sequence, int * nc, int * nca, int * TotalCollected, int Ntot){

	if(el ==NULL || one==NULL || Sequence==NULL || nc==NULL || TotalCollected==NULL || Ntot<0){
		return -1;
	}
	//Reduce number of charges in the sample

	*nc= *nc-1;
	*nca = *nca-1;
	Electrode_addCharge(el);
	printf("Charge reached electrode value of nc now %d\n",*nc);
	//Total number of collected charges increases
	*TotalCollected = *TotalCollected+1;
	//Set Dwell time to finished state
	printf("Set Dwel exit\n");
	setDwel(*one,1E6);

	int i;
	int ID;
	int ID2;
	//Sequence contains the ids of the charges which go
	//from 0 to Ntot-1;
	//Placing the charge that just reached the electrode to 
	//the back of the list
	ID2 = getE(*Sequence,1,1);

	for(i=0; i<Ntot-1;i++){
		//Updating the waiting queue of charges
		ID = getE(*Sequence,i+2,1);
		setE(*Sequence,i+1,1,ID);
	}

	setE(*Sequence,Ntot,1,ID2);

	return 0;
}

int ChargeSystem( int * nc, Electrode el, const int nca, matrix Sequence, ChargeArray chA){
	//Reduce number of charges in the electrode and increase the number
	//of charges back in the system Jumpint back into the system

	*nc = *nc+1;
	Electrode_minusCharge(el);
	//double tim = 1/getElectrode_Sum(el);
	//printf("Set Dwel having jumped back into th system %lg\n",tim);
	//setDwel(*one, -log(rand()/RAND_MAX)*tim);
	insertDwelltimePos(nca, chA, &Sequence);

	return 0;
}

int ChargeElectrode(Electrode el, Charge * one, matrix * Sequence, ChargeArray chA, int * nc, const int nca, int flag){
	//Hopping from electrode

	if(one==NULL || Sequence==NULL || nc==NULL || chA==NULL){
		return -1;
	}

	if(flag==0){
		//Reduce number of charges in the sample
		*nc = *nc-1;
		//Increase number of charges on electrode
		Electrode_addCharge(el);
		printf("Charge reached electrode value of nc now %d\n",*nc);
	}
	//Set Dwell time to that of the electrode
	double tim = 1/getElectrode_Sum(el); 
	printf("set dwel Electrode tim %g\n",tim);
	setDwel(*one, -log((double)rand()/RAND_MAX)*tim);

	if(tim>1){
		printf("Electrode set Dwel huge\n");
		exit(1);
	}
	//Sequence contains the ids of the charges which go
	//from 0 to Ntot-1;
	insertDwelltimePos(nca, chA, Sequence);

	return 0;

}

int insertDwelltimePos(const int nca, ChargeArray chA, matrix * Sequence){

	if( chA== NULL || Sequence==NULL || nca<0){
		return -1;
	}

	int i;
	int low;
	int high;
	int middle;
	int ID;
	int newID;
	Charge chi;
	Charge chj;

	//Charge id of first charge is the sequence
	ID = getE(*Sequence,1,1);

	low = 2;
	high = nca;

	while(low<=high){
		middle = (low+high)/2;
		//Chi is a probe into the sequence array
		chi = getCharge(chA,getE(*Sequence,middle,1));
		chj = getCharge(chA,ID);

		if( getDwel(chi) > getDwel(chj)){
			//If the probe is higher than the 
			//actual value reduce high to the middle-1
			high = middle-1;
		}else if(getDwel(chi) == getDwel(chj)){
			high = middle;
			low = middle+1;
		}else{
			//If the probe is less than actual 
			//value than increae low to the middle+1
			low = middle+1;
		}
	}
	//After the while loop is completed the value 
	//stored in middle should refer to a dwelltime
	//that is just less than that of charge chj
	for(i=1; i<high; i++){
		newID = getE( *Sequence,i+1,1);
		setE( *Sequence,i,1,newID);
	}

	setE(*Sequence, high,1,ID);

	return 0;
}

int SaveDataPoint(int * CurrentInc, int * NumAvgVel, int nc, int XElecOn, int YElecOn, int ZElecOn,\
		matrix * Xcurrent, matrix * Ycurrent, matrix * Zcurrent, matrix * timeArray,\
		matrix *Xvelocity, matrix * Yvelocity, matrix * Zvelocity,matrix * System,\
		matrix *Xelec_Drain, matrix * Yelec_Drain, matrix * Zelec_Drain,\
		matrix *Xelec_Source, matrix * Yelec_Source, matrix * Zelec_Source,\
		int * TotalX, int * TotalY, int * TotalZ, int * TrackX, int * TrackY,\
		int * TrackZ,const long double t, long double * TimeTrack1, double TStep, double SiteDistance,\
		long double * TotalVelX,long double * TotalVelY,long double * TotalVelZ,\
		int * Xdrain, int * Ydrain, int * Zdrain, int * Xsource, int * Ysource, int * Zsource,\
		const double ElectricFieldX, const double ElectricFieldY, const double ElectricFieldZ,\
		char * FileName){

	int rv;
	double q = 1.602E-19;
	double val;
	double vel;

	//Might need to resize the matrix if it becomes to large
	if((*CurrentInc)>getRows(*Xcurrent)){
		if (*CurrentInc >= 512 ){

			//Matrices have gotten to a size where they should be printed
			//to a file and emptied to free memory

			printTransportData(*System, *timeArray, *Xcurrent, *Ycurrent, *Zcurrent,\
					*Xelec_Drain, *Yelec_Drain, *Zelec_Drain,\
					*Xelec_Source, *Yelec_Source, *Zelec_Source,\
					*Xvelocity, *Yvelocity, *Zvelocity,\
					XElecOn, YElecOn, ZElecOn, FileName,\
					ElectricFieldX, ElectricFieldY, ElectricFieldZ);


			deleteMatrix(*System);
			deleteMatrix(*timeArray);

			deleteMatrix(*Xcurrent);
			deleteMatrix(*Ycurrent);
			deleteMatrix(*Zcurrent);
			if(XElecOn==1){
				deleteMatrix(*Xelec_Drain);
				deleteMatrix(*Xelec_Source);
			}
			if(YElecOn==1){
				deleteMatrix(*Yelec_Drain);
				deleteMatrix(*Yelec_Source);
			}
			if(ZElecOn==1){
				deleteMatrix(*Zelec_Drain);
				deleteMatrix(*Zelec_Source);
			}

			deleteMatrix(*Xvelocity);
			deleteMatrix(*Yvelocity);
			deleteMatrix(*Zvelocity);
			*System = newMatrix(8,1);
			//Keeps track of the time when each of the matrices are incremented
			*timeArray = newMatrix(8,1);

			//Calculate Current
			*Xcurrent = newMatrix(8,1);
			*Ycurrent = newMatrix(8,1);
			*Zcurrent = newMatrix(8,1);

			//Create source and drain for all electrodes that are turned on
			if(XElecOn==1){
				*Xelec_Drain = newMatrix(8,1);
				*Xelec_Source = newMatrix(8,1);
			}
			if(YElecOn==1){
				*Yelec_Drain = newMatrix(8,1);
				*Yelec_Source = newMatrix(8,1);
			}
			if(ZElecOn==1){
				*Zelec_Drain = newMatrix(8,1);
				*Zelec_Source = newMatrix(8,1);
			}

			//Drift Velocities should be in units of [m/s]
			*Xvelocity = newMatrix(8,1);
			*Yvelocity = newMatrix(8,1);
			*Zvelocity = newMatrix(8,1);

			*CurrentInc = 1;
		}else {

			rv = resizeRow(timeArray,getRows(*timeArray)*2);
			rv = resizeRow(Xcurrent,getRows(*Xcurrent)*2);
			rv = resizeRow(Ycurrent,getRows(*Ycurrent)*2);
			rv = resizeRow(Zcurrent,getRows(*Zcurrent)*2);
			rv = resizeRow(System,getRows(*System)*2);

			rv = resizeRow(Xvelocity,getRows(*Xvelocity)*2);
			rv = resizeRow(Yvelocity,getRows(*Yvelocity)*2);
			rv = resizeRow(Zvelocity,getRows(*Zvelocity)*2);

			if(XElecOn==1){
				rv = resizeRow(Xelec_Drain,getRows(*Xelec_Drain)*2);
				rv = resizeRow(Xelec_Source,getRows(*Xelec_Source)*2);
			}
			if(YElecOn==1){
				rv = resizeRow(Yelec_Drain,getRows(*Yelec_Drain)*2);
				rv = resizeRow(Yelec_Source,getRows(*Yelec_Source)*2);
			}
			if(ZElecOn==1){
				rv = resizeRow(Zelec_Drain,getRows(*Zelec_Drain)*2);
				rv = resizeRow(Zelec_Source,getRows(*Zelec_Source)*2);
			}

		}
	}

	printf("Saving Data \n\n");

	rv = setE(*timeArray,(*CurrentInc),1,(double)t);
	assert(rv==0);
	val = ((double)(*TotalX)*q)/(TStep);// * (rN/NCh);
	rv = setE(*Xcurrent,(*CurrentInc),1,val);
	assert(rv==0);
	val = ((double)(*TotalY)*q)/(TStep);// * (rN/NCh);
	rv = setE(*Ycurrent,(*CurrentInc),1,val);
	assert(rv==0);
	val = ((double)(*TotalZ)*q)/(TStep);// * (rN/NCh);
	rv = setE(*Zcurrent,(*CurrentInc),1,val);
	assert(rv==0);

	//
	if(nc!=0){
		vel =(double) ((*TotalVelX)+(long double)(*TrackX)*(long double)SiteDistance/((t-(*TimeTrack1))*(long double)nc))/((long double)(*NumAvgVel)+1);
		//printf("Xvel %g TotalVelX %Lg TrackX %d t %Lg TimeTrack1 %Lg nc %d\n",vel,*TotalVelX, *TrackX, t,*TimeTrack1, nc);
		rv = setE(*Xvelocity,(*CurrentInc),1,(double)vel);
		vel = (double) ((*TotalVelY)+(long double)(*TrackY)*(long double)SiteDistance/((t-(*TimeTrack1))*(long double)nc))/((long double)(*NumAvgVel)+1);
		rv = setE(*Yvelocity,(*CurrentInc),1,(double)vel);
		vel = (double) ((*TotalVelZ)+(long double)(*TrackZ)*(long double)SiteDistance/((t-(*TimeTrack1))*(long double)nc))/((long double)(*NumAvgVel)+1);
		rv = setE(*Zvelocity,(*CurrentInc),1,(double)vel);
	}else{
		vel=0;
		rv = setE(*Xvelocity,(*CurrentInc),1,(double)vel);
		rv = setE(*Yvelocity,(*CurrentInc),1,(double)vel);
		rv = setE(*Zvelocity,(*CurrentInc),1,(double)vel);
	}
	//Number of charges still in the system
	rv = setE(*System,(*CurrentInc),1,(double)nc);
	printf("System %d \t",nc);

	if(XElecOn==1){
		printf("Drain %d Source %d \n",*Xdrain,*Xsource);
		rv = setE(*Xelec_Drain,(*CurrentInc),1,(double)(*Xdrain));
		rv = setE(*Xelec_Source,(*CurrentInc),1,(double)(*Xsource));
	}

	if(YElecOn==1){
		setE(*Yelec_Drain,(*CurrentInc),1,(double)(*Ydrain));
		setE(*Yelec_Source,(*CurrentInc),1,(double)(*Ysource));
	}
	if(ZElecOn==1){
		setE(*Zelec_Drain,(*CurrentInc),1,(double)(*Zdrain));
		setE(*Zelec_Source,(*CurrentInc),1,(double)(*Zsource));
	}
	*TotalX = 0;
	*TotalY = 0;
	*TotalZ = 0;

	//Variables for tracking velocity
	*TimeTrack1 = t;
	*NumAvgVel= 0;
	*TrackX = 0;
	*TrackY = 0;
	*TrackZ = 0;
	*TotalVelX = 0;
	*TotalVelY = 0;
	*TotalVelZ = 0;

	*Xdrain = 0;
	*Xsource = 0;
	*Ydrain = 0;
	*Ysource = 0;
	*Zdrain = 0;
	*Zsource = 0;

	*CurrentInc = *CurrentInc+1;

	return 0;
}

int CheckIfElecHop(int xx, int yy, int zz,\
		int * Xsource, int * Ysource, int * Zsource,\
		int * Xdrain, int * Ydrain, int * Zdrain,\
		int * TrackX, int * TrackY, int * TrackZ,\
		long double * TotalVelX, long double * TotalVelY, long double * TotalVelZ,\
		const int SLength, const int SWidth, const int SHeight,\
		const int EndX, const int EndY, const int EndZ,\
		int * NumAvgVel, long double * TimeTrack1, long double t, double SiteDistance,\
		int * nc, int * ElecExit,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ){


	long double TrackX2 = (long double) (*TrackX);
	long double TrackY2 = (long double) (*TrackY);
	long double TrackZ2 = (long double) (*TrackZ);
	long double nc2 = (long double)(*nc);
	long double SiteDistance2 = (long double)SiteDistance;
	if(xx == -1 && PeriodicX!=1){
		//Arrived at back electrode
		(*Xsource)=*Xsource+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(xx == SLength*EndX && EndX != 0 ){
		//Arrived at front electrode
		(*Xdrain)=*Xdrain+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;
	}

	//Check if a charge has arrived at the Left or Right electrode
	if(yy == -1 && PeriodicY!=1){
		//Arrived at Left electrode
		(*Ysource)=*Ysource+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(yy == SWidth*EndY && EndY != 0 ){
		//Arrived at front electrode
		(*Ydrain)=*Ydrain+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;
	}

	//Check if a charge has arrived at the top or bottom electrode
	if(zz == -1 && PeriodicZ!=1){
		//Arrived at bottom electrode
		(*Zsource)=*Zsource+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(zz == SHeight*EndZ && EndZ != 0 ){
		//Arrived at top electrode
		(*Zdrain)=*Zdrain+1;
		(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;
	}

	return 0;

}

//Check if next hop will use the cluster hopping algorithm 
int ClusterHopCheck(const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
		const int SLength, const int SWidth, const int SHeight,\
		const int EndX, const int EndY, const int EndZ,\
		const int x, const int y, const int z,SiteNode site,\
		const int Cx, const int Cy, const int Cz,\
		int * CheckX, int * CheckY, int * CheckZ){

	int rv; 
	*CheckX = 0;
	*CheckY = 0;
	*CheckZ = 0;

	if(PeriodicX==1){
		rv = getCluster_elecXid((ClusterLL) getClusterList(site));
		if(x/SLength==0){
			//Check to see if next to the back electrode
			if(rv==2 || rv==0){
				(*CheckX)=1;
			}else{
				(*CheckX)=0;
			}
		}else if(Cx/SLength==(EndX-1) || Cx/SLength==EndX){
			//Check to see if next to the front electrode
			if(rv==2 || rv==1){
				(*CheckX)=1;
			}else{
				(*CheckX)=0;
			}
		}
	}
	if(PeriodicY==1){
		rv = getCluster_elecYid((ClusterLL) getClusterList(site));
		if(y/SWidth==0){
			//Check if next to left electrode
			if(rv==2 || rv==0){
				(*CheckY)=1;
			}else{
				(*CheckY)=0;
			}
		}else if(Cy/SWidth==(EndY-1) || Cy/SWidth==EndY){
			//Check if next to right electrode
			if(rv==2 || rv==1){
				(*CheckY)=1;
			}else{
				(*CheckY)=0;
			}
		}
	}
	if(PeriodicZ==1){
		rv = getCluster_elecZid((ClusterLL) getClusterList(site));

		if(z/SHeight==0){
			if(rv==2 || rv==1){
				(*CheckZ)=1;
			}else{
				(*CheckZ)=0;
			}
		}else if(Cz/SHeight==(EndZ-1) || Cz/SHeight==EndZ){
			if(rv==2 || rv==1){
				(*CheckZ)=1;
			}else{
				(*CheckZ)=0;
			}
		}
	}

	return 0;
}

int HoppingToSurroundingSites(SiteNode site, int codeX, int codeY, int codeZ){

	if(site==NULL || codeX<0 || codeY<0 || codeZ<0 ||\
			codeX>2 || codeY>2 || codeZ>2){
		printf("ERROR incorrect input parameters detected in HoppingToSurroundingSites\n");
		return -1;
	}
	//codeX values:
	//0 - can hop to +x and -X
	//1 - cannot hop to -x
	//2 - cannot hop to +x
	//codeY values:
	//0 - can hop to +y and -y
	//1 - cannot hop to -y
	//2 - cannot hop to +y
	//codeZ values:
	//0 - can hop to +z and -z
	//1 - cannot hop to -z
	//2 - cannot hop to +z
	double value;
	double position;
	int l;

	if(codeX==0 && codeY==0 && codeZ==0){
		position = (double)rand() / RAND_MAX;
		l=0;
		//Cycle through the pvals until the pval
		//is equal or above the random number

		while( getSN_p(site,l)<position){
			l++;
		}

	}else if(codeX>0 && codeY==0 && codeZ==0){
		position = (double) rand() /RAND_MAX;
		//Need to redefine the pvals
		matrix mtxPval = newMatrix(5,1);
		if(codeX==2){
			//Ignore the hop forward
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			if(l==1){
				l=0;
			}
			deleteMatrix(mtxPval);

		}else{
			position = (double) rand() /RAND_MAX;
			//Ignore the hop backward
			value = getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			deleteMatrix(mtxPval);

		}
	}else if(codeX==0 && codeY>0 && codeZ==0){
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(5,1);
		if(codeY==1){
			//Not allowed to hop to the left or -y
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			if(l<=2){
				l--;
			}

			deleteMatrix(mtxPval);

		}else{
			//Not allowed to hop to the right or +y
			position = (double) rand() /RAND_MAX;
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			if(l<=3){
				l--;
			}

			deleteMatrix(mtxPval);

		}
	}else if(codeX==0 && codeY==0 && codeZ>0){
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(5,1);
		if(codeZ==1){
			//Cannot hop Below
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			if(l<=4){
				l--;
			}

			deleteMatrix(mtxPval);

		}else{
			//Cannot hop Above
			position = (double) rand() /RAND_MAX;
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
			value += getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,2,1,value);
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,3,1,value);
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,4,1,value);
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,5,1,value);

			DivideEachElement(&mtxPval,value);

			l=1;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getE(mtxPval,l,1)<position){
				l++;
			}

			if(l<=5){
				l--;
			}

			deleteMatrix(mtxPval);

		}
	}else if(codeX>0 && codeY>0 && codeZ==0){
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(4,1);
		if(codeX==2){
			//Not allowed to hop forward
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}else{
			//Not allowed to hop backward
			value = getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}

		if(codeY==1){
			//Not allowed to hop to the left
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,2,1,value);
		}else{
			//Not allowed to hop to the right
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,2,1,value);
		}
		value += getSN_p(site,4)-getSN_p(site,3);
		setE(mtxPval,3,1,value);
		value += getSN_p(site,5)-getSN_p(site,4);
		setE(mtxPval,4,1,value);

		DivideEachElement(&mtxPval,value);

		l=1;
		//Cycle through the pvals until the pval
		//is equal or above the random number
		while( getE(mtxPval,l,1)<position){
			l++;
		}

		if(l>=3){
			l++;
		}else	if(l==2){
			if(codeY==1){
				l=3;	
			}else{
				l=2;
			}
		}else	if(l==1){
			if(codeX==2){
				l=0;
			}else{
				l=1;
			}
		}

		deleteMatrix(mtxPval);


	}else if(codeX==0 && codeY>0 && codeZ>0){
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(4,1);

		value = getSN_p(site,0);
		setE(mtxPval,1,1,value);
		value += getSN_p(site,1)-getSN_p(site,0);
		setE(mtxPval,2,1,value);

		if(codeY==1){
			//Not allowed to hop to the left
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,3,1,value);
		}else{
			//Not allowed to hop to the right
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,3,1,value);
		}
		if(codeZ==1){
			//Not allowed to hop below
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,4,1,value);
		}else{
			//Not allowed to hop above
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
		}

		DivideEachElement(&mtxPval,value);

		l=1;
		//Cycle through the pvals until the pval
		//is equal or above the random number

		while( getE(mtxPval,l,1)<position){
			l++;
		}

		if(l<3){
			l--;
		}else if(l==3){
			if(codeY==1){
				l=3;
			}else{
				l=2;
			}
		}else{
			if(codeZ==1){
				l=5;
			}else {
				l=4;
			}
		}

		deleteMatrix(mtxPval);

	}else if(codeX>0 && codeY==0 && codeZ>0){
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(4,1);

		if(codeX==2){
			//Not allowed to hop forward
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}else{
			//Not allowed to hop backward
			value = getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}

		value += getSN_p(site,2)-getSN_p(site,1);
		setE(mtxPval,2,1,value);
		value += getSN_p(site,3)-getSN_p(site,2);
		setE(mtxPval,3,1,value);

		if(codeZ==1){
			//Not allowed to hop below
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,4,1,value);
		}else{
			//Not allowed to hop backward
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,4,1,value);
		}

		DivideEachElement(&mtxPval,value);

		l=1;
		//Cycle through the pvals until the pval
		//is equal or above the random number
		while( getE(mtxPval,l,1)<position){
			l++;
		}

		if(l==1){
			if(codeX==2){
				l=0;
			}else{
				l=1;
			}
		}else if(l==4){
			if(codeZ==1){
				l=5;
			}else{
				l=4;
			}
		}else{
			l=l;
		}

		deleteMatrix(mtxPval);

	}else{
		position = (double) rand() /RAND_MAX;
		matrix mtxPval = newMatrix(3,1);

		if(codeX==2){
			//Not allowed to hop forward
			value = getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}else{
			//Not allowed to hop backward
			value = getSN_p(site,1)-getSN_p(site,0);
			setE(mtxPval,1,1,value);
		}
		if(codeY==1){
			//Not allowed to hop to the left
			value += getSN_p(site,3)-getSN_p(site,2);
			setE(mtxPval,2,1,value);
		}else{
			//Not allowed to hop to the right
			value += getSN_p(site,2)-getSN_p(site,1);
			setE(mtxPval,2,1,value);
		}
		if(codeZ==1){
			//Not allowed to hop below
			value += getSN_p(site,5)-getSN_p(site,4);
			setE(mtxPval,3,1,value);
		}else{
			//Not allowed to hop backward
			value += getSN_p(site,4)-getSN_p(site,3);
			setE(mtxPval,3,1,value);
		}

		DivideEachElement(&mtxPval,value);

		l=1;
		//Cycle through the pvals until the pval
		//is equal or above the random number
		while( getE(mtxPval,l,1)<position){
			l++;
		}

		if(l==1){
			if(codeX==2){
				l=0;
			}else{
				l=1;
			}
		}else if(l==2){
			if(codeY==1){
				l=3;
			}else{
				l=2;
			}
		}else{
			if(codeZ==1){
				l=5;
			}else{
				l=4;
			}
		}

		deleteMatrix(mtxPval);

	}

	return l;
}

int printFileEnergy(const_SNarray snA, char * FileName,\
		double ElectricEnergyX, double ElectricEnergyY, double ElectricEnergyZ,\
		ParameterFrame PF) {

	if(snA==NULL || PF==NULL){
		return -1;
	}

	//Constants
  //Units of Coulombs
  const double q = 1.601E-19;
  //Units of Farads/m = (C^2/J) * (1/m)
  const double epsilon0 = 8.8541878176E-12;

	//Declaring variables from Parameter Frame
	int SLength;
	int SWidth;
	int SHeight;
	int EndX;
	int EndY;
	int EndZ;
	double SiteDistance;
	double RelativePerm;

	//Declaring local variables
	double EnergyImageForce;

	double midX;
	double midY;
	double midZ;

	double EX;
	double EY;
	double EZ;

	double val;
	
	//Initializing Parameter Frame variables
	SLength = PFget_Len(PF);
	SWidth = PFget_Wid(PF);
	SHeight = PFget_Hei(PF);
	EndX = PFget_EndX(PF);
	EndY = PFget_EndY(PF);
	EndZ = PFget_EndZ(PF);
	SiteDistance = PFget_SiteDist(PF);
	RelativePerm = PFget_RelativePerm(PF);

	//Initializing local variables
	EnergyImageForce = pow(q,2)/(16*M_PI*epsilon0*RelativePerm*SiteDistance)*1/q;
	//Determine if length width and height are even
	midX =((double) SLength)/2;
	midY =((double) SWidth)/2;
	midZ =((double) SHeight)/2;

	char buf[256];
	snprintf(buf, sizeof buf,"%s%s",FileName,"Energy.xyz");

	int i, j, k;
	SiteNode sn;
	FILE * EnergyOut;
	if((EnergyOut = fopen(buf,"w")) == NULL) {
		printf("Error! unable to open Energy.xyz\n");

	}else{
		fprintf(EnergyOut,"%d\n\n",(getAtotal(snA)));
		double id=0;
		double jd=0;
		double kd=0;

		for (i=0;i<SLength;i++){
			jd=0;

			val = (double)i;
			EX = (midX-val)*ElectricEnergyX;

			for(j=0;j<SWidth;j++){
				kd=0;
				val = (double)j;
				EY = (midY-val)*ElectricEnergyY;

				for(k=0;k<SHeight;k++){

					val = (double)k;
					EZ = (midZ-val)*ElectricEnergyZ;

					sn=getSN(snA,i,j,k);
					if( ((i==0 || i==(SLength-1)) && EndX==1) || \
					 		((j==0 || j==(SWidth-1)) && EndY==1) || \
					 		((k==0 || k==(SHeight-1)) && EndZ==1)){
						fprintf(EnergyOut,"C \t %f \t %f \t %f \t %f \t %f\n",\
							id,jd,kd,getEnergy(sn),getEnergy(sn)+EX+EY+EZ-EnergyImageForce);
					}	else{
						fprintf(EnergyOut,"C \t %f \t %f \t %f \t %f \t %f\n",\
							id,jd,kd,getEnergy(sn),getEnergy(sn)+EX+EY+EZ);
					}
					kd++;
				}
				jd++;
			}
			id++;
		}
		fclose(EnergyOut);
	}

	return 0;
}

