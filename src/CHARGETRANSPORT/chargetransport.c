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
#include <stdbool.h>

#include "chargetransport.h"
#include "../PARAMETERS/read.h"
#include "../CHARGE/charge.h"
#include "../CHARGESITENODE/chargesitenode.h"
#include "../SITENODE/sitenode.h"
#include "../CLUSTER/cluster.h"
#include "../CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../MATRIX/matrix.h"
#include "../MONTECARLO/montecarlo.h"
#include "../CLUSTERSITENODE/clustersitenode.h"
#include "../FUNCTIONS/functions.h"
#include "../ELECTRODE/electrode.h"
#include "../IO/io.h"
#include "../ERROR/error.h"

// Why can we not hop past the edges? What if it is periodic?
// I think the reason for this is in the case that there
// are also electrodes defined in the y and z axis if that
// is the case then we would not be hopping to a site but an
// electrode. This needs to be fixed
int HopToElecX(SNarray snA, Electrode elXb, Charge * one, int * future, int EndY, int EndZ){
  // Warning if you have more than one pair of electrodes
  // you will need to make adjustments to this code

  //Elecid defines which electrode could be hopping too
  //where:
  //	getAtotal(snA) + 0 - (-x)
  //	getAtotal(snA) + 1 - (+x)
  //	getAtotal(snA) + 2 - (-y)
  // 	getAtotal(snA) + 3 - (+y)
  //	getAtotal(snA) + 4 - (-z)
  //	getAtotal(snA) + 5 - (+z)

  #ifdef _ERROR_CHECKING_ON_
  if(elXb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR elXb is NULL in HopToElecX\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  
  if( snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in HopToElecX\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  
  if(*one == NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Charge *one is NULL in HopToElecX\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }	
  
  if( EndY<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR EndY is less than 0 in HopToElecX\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  
  if(EndZ<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR EndZ is less than 0 in HopToElecX\n");
    #endif 
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int SLength = getAlen(snA);
  int SWidth  = getAwid(snA);
  int SHeight = getAhei(snA);
  //First Check is to see if all the neighboring sites are occupied. 
  //If the site is not located at a boundary
  //If they are all occupied set the flag to 1
  int xx = getCx(*one);
  int yy = getCy(*one);
  int zz = getCz(*one);

  int l;
  int i;
  int j;
  int k;

  ProjectChargePositionOntoSiteNodeReferenceFrame(&i,&j,&k,xx,yy,zz,SLength,SWidth,SHeight);

  //   codeY = 1 exclude left hop
  //				 = 2 exclude right hop
  //   codeZ = 1 exclude bottom hop
  //				 = 2 exclude top hop

  int codeX = 0;
  int codeY = 0;
  int codeZ = 0;

  if(xx==0){
    //Posibility to jump to (-x) electrode

    //Have to ensure that charge is not at 
    //the edge of the system in the y and z directon
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
  SNarray snAelec       = (SNarray) getElectrode_AdjacentSites(elXb);
  SiteNode site = getSN(snAelec,0,j,k);
  l             = HoppingToSurroundingSites(site,codeX,codeY,codeZ);

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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR elYl is NULL in HopToElecY\n");
    #endif
		return -1;
	}

	int SLength = getAlen(snA);
	int SWidth  = getAwid(snA);
	int SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	int xx = getCx(*one);
	int yy = getCy(*one);
	int zz = getCz(*one);

	int l;
	int i;
	int j;
	int k;

  ProjectChargePositionOntoSiteNodeReferenceFrame(&i,&j,&k,xx,yy,zz,SLength,SWidth,SHeight);

	//   codeX = 1 exclude beh hop
	//				 = 2 exclude front hop
	//   codeZ = 1 exclude bottom hop
	//				 = 2 exclude top hop

	int codeX = 0;
	int codeY = 0;
	int codeZ = 0;

	//This is the random number used to determine which direction a charge hops

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
	SiteNode site   = getSN(snAelec,i,0,k);
	l               = HoppingToSurroundingSites(site, codeX,codeY,codeZ);

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
    #ifdef _ERROR_
    fprintf(stderr,"ERROR elZb is NULL in HopToElecZ\n");
    #endif
		return -1;
	}

  int SLength = getAlen(snA);
  int SWidth  = getAwid(snA);
	int SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	int xx = getCx(*one);
	int yy = getCy(*one);
	int zz = getCz(*one);

	int l;
	int i;
	int j;
	int k;	

  ProjectChargePositionOntoSiteNodeReferenceFrame(&i,&j,&k,xx,yy,zz,SLength,SWidth,SHeight);

	//   codeX = 1 exclude beh hop
	//				 = 2 exclude front hop
	//   codeY = 1 exclude left hop
	//				 = 2 exclude right hop

	int codeX = 0;
	int codeY = 0;
	int codeZ = 0;

	//This is the random number used to determine which direction a charge hops

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
	SiteNode site   = getSN(snAelec,i,j,0);
	l               = HoppingToSurroundingSites(site,codeX,codeY,codeZ);

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

	int    row;
	int    col;
	double position;
  int    flag;

	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
	position   = (double)rand() / RAND_MAX;
	matrix mtx = (matrix) getElectrode_HopRates(el);
	
	//Cycle through the pvals until the pval
	//is equal or above the random number

	if(mtx==NULL){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR matrix mtx does not exist in ElecHopOffX\n");
		#endif
    exit(1);
	}
  
  flag = 0;

	for(row=1;row<=getRows(mtx);row++){
		for(col=1;col<=getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				flag=1;
        break;
			}
		}
    if(flag==1){
      break;
    }
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to future
	*future = getIndex(snA,0,row-1,col-1); 

	if(*future==-1){
    #ifdef _ERROR_
		fprintf(stderr,"Value or row %d Value of col %d future %d position %g in ElecHopOffX\n",row,col,*future,position);
		fprintf(stderr,"Number of columns mtx has %d Number of rows mtx has %d\n",getCols(mtx),getRows(mtx));
    fprintf(stderr,"ERROR Hopping off Electrode and out of system\n");
		#endif
    exit(1);
	}

	return 0;
}

int ElecHopOffY(Electrode el, int * future, SNarray snA){

	int    row;
	int    col;
	double position;
  int    flag;
	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;
	matrix mtx = (matrix) getElectrode_HopRates(el);
	//Cycle through the pvals until the pval
	//is equal or above the random number
  
  flag = 0;

	for(row=1;row<getRows(mtx);row++){
		for(col=1;col<getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				flag = 1;
        break;
			}
		}
    if(flag==1){
      break;
    }
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to futre

	*future = getIndex(snA,row-1,0,col-1); 

	if(*future==-1){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Hopping off Electrode and out of system\n");
		#endif
    exit(1);
	}

	return 0;
}

int ElecHopOffZ(Electrode el, int * future, SNarray snA){

	int    row;
	int    col;
	double position;
  int    flag;
	//Cycle through the hops off the electrode

	//This is the random number used to determine which direction a charge hops
  position   = (double)rand() / RAND_MAX;
  matrix mtx = (matrix) getElectrode_HopRates(el);
	//Cycle through the pvals until the pval
	//is equal or above the random number

  flag = 0;

	for(row=1;row<getRows(mtx);row++){
		for(col=1;col<getCols(mtx);col++){
			if( getE(mtx,row,col)>=position){
				flag = 1;
        break;
			}
		}
    if(flag==1){
      break;
    }
	}

	//The site the charge has chosen to hop to is located on rows (row) and
	//cols (col)
	//The next step is to grab the id of the site the charge has just hopped to
	//and assign it to futre

	*future = getIndex(snA,row-1,col-1,0); 
	if(*future==-1){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Hopping off Electrode and out of system\n");
		#endif
    exit(1);
	}
	return 0;
}

int SiteHop(SNarray snA        , Charge * one     , SiteNode site      , int * future       ,\
            const int EndX     , const int EndY   , const int EndZ     , const int XElecOn  ,\
            const int YElecOn  , const int ZElecOn, const int PeriodicX, const int PeriodicY,\
            const int PeriodicZ){

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
	double position;
	int    l;
	//Below values define location of the charge in the
	//sample when the periodicity is removed
	int x, xx;
	int y, yy;
	int z, zz;
	int SLength;
	int SWidth;
	int SHeight;

	SLength = getAlen(snA);
	SWidth  = getAwid(snA);
	SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	xx = getCx(*one);
	yy = getCy(*one);
	zz = getCz(*one);

  ProjectChargePositionOntoSiteNodeReferenceFrame(&x,&y,&z,xx,yy,zz,SLength,SWidth,SHeight);

	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;

	if(x!=0 && x!=(SLength-1) && y!=0 && y!=(SWidth-1) && z!=0 && z!=(SHeight-1)){
		OccAllNei(snA, x, y, z);

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


  if(*future==-1){
    printf("Future set to -1 exiting this is in SiteHop\n");
    exit(1);
  }

	return 0;

}

int ClusterHop(SNarray snA, Charge * ch,  double * tim , int *newID){

  #ifdef _ERROR_CHECKING_ON_
	if (snA==NULL){
    fprintf(stderr,"ERROR snA is NULL in ClusterHop\n");
    return return_error_val();
  }
  if(ch==NULL){
    fprintf(stderr,"ERROR ch is NULL in ClusterHop\n");
		return return_error_val();
	}
  #endif

	//Flags determine what options a charge has if it is in a cluster
	int flag1;
	int flag2;

	int rv;
	int x;
	int y;
	int z;
	int ID;

	double position;
	double position2;

	int timeflag;

	flag1 = 0;
	flag2 = 0;

  ProjectChargePositionOntoSiteNodeReferenceFrame(&x,&y,&z,getCx(*ch),getCy(*ch),getCz(*ch),\
                                                  getAlen(snA),getAwid(snA),getAhei(snA));

	//newID does not change unless
	//hop to a different position
	ID     = getIndex(snA, x,y,z);
	*newID = ID;

	position  = (double)((double)rand()/(double)RAND_MAX);
	position2 = (double)((double)rand()/(double)RAND_MAX);

	*tim = 0;
  printf("Setting tim to 0 in ClusterHop\n");
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

  printf("Cluster hop flag1 %d flag2 %d\n",flag1,flag2);

	if(flag1==0 && flag2==0){
		//Hopping within cluster and to neighbors is allowed

		//Value of rv
		//1 - hopped off cluster
		//0 - stayed within cluster
		rv = HopOnOffCluster(snA,ID,position);

		if(rv==1){
			timeflag = 2;

			//Hopped off cluster determining which neighboring 
			//site hopped to as well as the time it takes
			HopOffCluster(snA, ID, position2, newID, tim);
      printf("1 Calling HopOffCluster tim %g\n",*tim);
		}else if(rv==0){
			timeflag = 0;
			//Stayed within cluster now will determine which
			//site the charge moves to within the cluster
			HopWithinCluster(snA, ID, position2, newID);
		}

	}else if(flag1==0 && flag2==1){
		//Hopping within cluster not allowed
		//Hopping to neighbors is permitted
		//SiteNode sn = getSN(snA, x,y,z);
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
      printf("2 Calling HopOffCluster tim %g\n",*tim);
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


    printf("1 HoppingWithinCluster tim %g\n",*tim);
	
  }else if(flag1==1 && flag2==1){
		//Hopping to neighbors or within cluster
		//is not allowed
		timeflag = 0;
		//printf("Hopping is not allowed\n");
		//Should use tcluster to increment the time

    printf("2 HoppingWithinCluster tim %g\n",*tim);
	}else{
		//malformed input or not part of a cluster
    printf("Malformed input  tim %g\n",*tim);
		return -1;
	}


  printf("timeflag %d\n",timeflag);

  printf("ClusterHop timeflag %d\n",timeflag);
	if(timeflag==0){
		SiteNode sn    = getSNwithInd(snA,ID);
		ClusterLL ClLL = (ClusterLL) getClusterList(sn);
		*tim           = getCluster_time(ClLL);
    if(*tim <= 0){
      printf("Cluster Hopping time is less than or equal to 0\n");
      exit(1);
    }
	}else if(timeflag==1){
		SiteNode sn    = getSNwithInd(snA,ID);
		ClusterLL ClLL = (ClusterLL) getClusterList(sn);
		*tim           = getCluster_time(ClLL)*20;
    if(*tim <= 0){
      printf("Cluster Hopping time is less than or equal to 0\n");
      exit(1);
    }
	}

	return 0;

}


//Can use MakeHop when jumping off electrode but not when jumping on because (newID?)
int MakeHop(SNarray snA, int newID             , Charge *ch                    , int * totalX       ,\
    int * totalY       , int * totalZ          , ParameterFrame PF             , double KT          ,\
    double electricEnergyX, double electricEnergyY, double electricEnergyZ,double MarcusCoef        ,\
    const long double t){

 // printf("Make Hop time %Lg\n",t);
  const int PeriodicX = PFget_Px(PF);
  const int PeriodicY = PFget_Py(PF);
  const int PeriodicZ = PFget_Pz(PF);

  const int XElecOn = PFget_XElecOn(PF);
  const int YElecOn = PFget_YElecOn(PF);
  const int ZElecOn = PFget_ZElecOn(PF);

  const int    DecayOn           = PFget_DecayOn(PF);
  const double DecayTime         = PFget_DecayTime(PF);
  const double DecayProb         = PFget_DecayProb(PF);
  const double DecayDisplacement = PFget_DecayDisplacement(PF);
	//Now we need to correctly move the charge to
	//the new position
	int x    , y    , z    ;
	int x1   , y1   , z1   ;
	int XDiff, YDiff, ZDiff;
	int xD1  , yD1  , zD1  ;
	int xD2  , yD2  , zD2  ;

  ProjectChargePositionOntoSiteNodeReferenceFrame(&x,&y,&z,getCx(*ch),getCy(*ch),getCz(*ch),\
                                                  getAlen(snA), getAwid(snA), getAhei(snA));

  // Get the old id of the site the charge was on, is hopping from
  int oldID = getIndex(snA,x,y,z);

  //This is the protocol for refering to the electrodes
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
				xD1         = (x1+SLength)-x;
				xD2         = x1-(x+SLength);
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
      //printf("Checking Decay\n");
      //fflush(stdout);
      checkDecay(PF,
                 snA,
                 oldID,
                 KT,
                 ch,
                 electricEnergyX,
                 electricEnergyY,
                 electricEnergyZ,
                 MarcusCoef,
                 t );

      //printf("Function decay success\n");
      //fflush(stdout);
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

  //This means the charge sucessfully hopped
  //to a new location
	return 0;

}

int CheckPt_Test_TOF(int * CheckPtNum, int CheckFileExist      , char * FileNameCheckPtVersion,\
                     int FileNameSize, const double Vx         , const double Vy              ,\
                     const double Vz , const double Temperature){


	if(CheckFileExist==0){
		*CheckPtNum = CheckPt_Latest_TOF(FileNameCheckPtVersion,FileNameSize,Vx,Vy,Vz,Temperature);
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

int CheckPt_Test_CELIV(int * CheckPtNum,int CheckFileExist, char * FileNameCheckPtVersion, int FileNameSize,\
											 const double Temperature){

	if(CheckFileExist==0){
		*CheckPtNum = CheckPt_Latest_CELIV(FileNameCheckPtVersion,FileNameSize,Temperature);
		//Just because a .chpt file exist it does not mean it exists for these parameters

		if(*CheckPtNum>0){
			*CheckPtNum = *CheckPtNum+1;
			//Confirmed that a checkpt file exists
			printf("CELIV CheckPoint file found for Temperature %g\n",Temperature);
			printf("File Name %s\n",FileNameCheckPtVersion);
			return 1;

		} else if(*CheckPtNum==0){
			//Confirmed that a checkpt file did not exist for this version
			printf("CELIV CheckPoint file does not exist for Temperature %g\n",Temperature);
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




int Pre_randomWalk(const int CheckPtStatus  , char * FileNameCheckPtVersion, char * FileName       ,\
		               long double * t          , matrix * Sequence            , ChargeArray * chA     ,\
                   matrix * FutureSite      , ArbArray * ClArLL            , SNarray * snA         ,\
                   ParameterFrame PF        ,	double electricEnergyX       , double electricEnergyY,\
                   double electricEnergyZ   , int r                        , double Vx             ,\
                   double Vy                , double Vz                    , double Temperature    ,\
                   long int * n             , int * nc                     , int * nca             ,\
		               Electrode * elXb         , Electrode * elXf             , Electrode * elYl      ,\
		               Electrode * elYr         , Electrode * elZb             , Electrode * elZa      ){

	//Declaring constants
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	//Declaring Variables from parameter frame
	int    method;
	int    SLength;
	int    SWidth;
	int    SHeight;
	int    EndX;
	int    EndY;
	int    EndZ;
	int    XElecOn;
	int    YElecOn;
	int    ZElecOn;
	int    Ntot;
	int    NCh;
  int    ClusterAlg;
	double AttemptToHop;
	double reOrgEnergy;
	double gamma;
	double SiteDistance;

	//Declaring local variables
	int    loop;
	int    clusterfileExist;
	int    Num_elXf;
	int    Num_elXb;
	int    Num_elYl;
	int    Num_elYr;
	int    Num_elZb;
	int    Num_elZa;
	double MarcusJ0;
	double MarcusCoeff;
	double KT;
	int    OrderL;

	//Initializing Variables from parameter frame
	method       = PFget_method(PF);
	SLength      = PFget_Len(PF);
	SWidth       = PFget_Wid(PF);
	SHeight      = PFget_Hei(PF);
	EndX         = PFget_EndX(PF);
	EndY         = PFget_EndY(PF);
	EndZ         = PFget_EndZ(PF);
	XElecOn      = PFget_XElecOn(PF);
	YElecOn      = PFget_YElecOn(PF);
	ZElecOn      = PFget_ZElecOn(PF);
	Ntot         = PFget_Ntot(PF);
	NCh          = PFget_NCh(PF);
  ClusterAlg   = PFget_ClusterAlg(PF);
	AttemptToHop = PFget_AttemptToHop(PF);
	reOrgEnergy  = PFget_reOrg(PF);
	gamma        = PFget_gamma(PF);
	SiteDistance = PFget_SiteDist(PF);

	//Initializing Local variables
	KT   = kB*Temperature;
	*snA = newSNarray(SLength, SWidth, SHeight);

	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0    = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	//Calculating full Marcus Coefficient;
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(-2*gamma*SiteDistance);

	//This means we are starting from scratch
	if(method==0){
		
		*Sequence     = newMatrix(Ntot,1);
		(*FutureSite) = newMatrix(Ntot,1);
		//printf("Value of Ntot %d\n",Ntot);
		//if(rv==-1){
		//	printf("ERROR FutureSite problem!\n");
		//	exit(1);
		//}
		
		if(FutureSite==NULL){
			printf("ERROR FutureSite NULL\n");
			exit(1);
		}

		if(CheckPtStatus==0){

			//Lets first make sure that a .cluster file does not exist before we
			//create new energies
			clusterfileExist = CheckPt_Cluster_TOF(Vx, Vy, Vz, Temperature, r);

			if(clusterfileExist==0 && ClusterAlg!=0){
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
			*chA = initCharget0( *Sequence, *snA,  Ntot, NCh,\
					XElecOn, YElecOn, ZElecOn,EndX, EndY, EndZ);
      
      //If the cluster algorithme is turned on we will allow 
      //charges to contain a link list that keeps up with 
      //how often a charge hops to a site, not applicable for
      //CELIV method because rates change too much
      //this is only applicaple for ClusterAlg 2
      if(ClusterAlg==2 && method==0){
        initChargeArrayPath(*chA, PFget_ClusterAlgRec(PF));  
      }
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

// This code is obsolete Ben showed that it was not optimal to 
// determine clusters before running the simulations we will
// keep it turned on though for comparitive purposes

    if(clusterfileExist==-1 && PFget_ClusterAlg(PF)==1){
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

    }
      
        printFileEnergy((*snA), &FileName[0], electricEnergyX,\
						electricEnergyY, electricEnergyZ,PF);
		


		}else if(CheckPtStatus==1) {
			//This means we are starting from a checkpt file that already exists
      printf("ERROR Loading Checkpt files has been disabled\n");

			//Lets first make sure that a .cluster file does not exist before we
			//recalculate everything
			clusterfileExist = CheckPt_Cluster_TOF(Vx, Vy, Vz, Temperature, r);

			//Charge Array chA is created in here
			printf("Loading Charge and Site information from .ckpt\n");
			/*Load_CheckPt_Data_TOF( t, snA, chA, Sequence,\
					FutureSite, FileNameCheckPtVersion, n,nc, nca,\
					&Num_elXb, &Num_elXf, &Num_elYl, &Num_elYr,\
					&Num_elZb, &Num_elZa,Vx,Vy,Vz);
      */
			if(FutureSite==NULL || Sequence==NULL || chA==NULL){
				printf("ERROR in the load_checkPt_Data_TOF function a datastructure\n");
				printf("Has been found to be NULL\n");
				exit(1);
			}

			////////////////////////////////////////////////////////////////////////
			//Initialize Electrodes
			printf("Initializing Electrodes\n");
			initElec(electricEnergyX, electricEnergyY, electricEnergyZ, MarcusCoeff,\
					KT, *snA,elXb, elXf, elYl, elYr, elZb, elZa, PF);


			//printf("Updating number of charges on electrodes based on .ckpt\n");
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

			if(PFget_ClusterAlg(PF)!=0){

				if(clusterfileExist==0 ){
					//The .cluster file exists we have already loaded the site information
					//We only want to load cluster information
					LoadCluster_Only( &FileName[0], &OrderL, snA, electricEnergyX,\
							electricEnergyY, electricEnergyZ, ClArLL, KT);
				}

				ConnectClusterElec( ClArLL,\
						(*elXb), (*elXf), (*elYl), (*elYr),\
						(*elZb), (*elZa) );
			}
		}else{
			printf("ERROR found in Pre_randomWalk value of checkptstatus %d\n",CheckPtStatus);
			exit(1);
		}

	}else if(method==1){

		if(CheckPtStatus==0){

			//It makes no sense to check for or create cluster files when using the 
			//CELIV method this is because we would have to recalculate where the clusters
			//were for each time step

			//Initialize Site Energies
			printf("Initializing Site Energies\n");
			initSite(electricEnergyX, electricEnergyY, electricEnergyZ,\
					KT, *snA, PF);
			
			//Initialize all charges in the Matrix and creates chargearray
			printf("Initializing Charges\n");
			*chA = initCharget0_Thermal( *snA, PF, Temperature,\
					XElecOn, YElecOn, ZElecOn);

			//Create Sequence matrix to store the charges and the order 
			//they should be moved which is based on their dwelltime
			//Charge ids start at 0 and go to Ntot-1
			printf("Initializing Sequence of Charges\n");
			Ntot = PFget_Ntot(PF);
			*Sequence = newMatrix(Ntot,1);
			for(loop = 0; loop<Ntot; loop++){
				setE( *Sequence,loop+1,1,loop);
			}
			//printf("Here is the problem\n");
			//exit(1);
			(*FutureSite) = newMatrix(Ntot,1);
			quickSort(0, Ntot-1, *Sequence, *chA);
			
			//t - global time initially 0 when starting from scratch
			//n - number of steps that charges have been injected starts at 1
			//nc - Number of charges initially in the system equal to the Number
			//		initially injected
			//nca - Number of active charges in the system initiallyequal to the number injected
			*t = 0;
			*n = 1;
			NCh = PFget_NCh(PF);
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

			printFileEnergy((*snA), &FileName[0], electricEnergyX,\
					electricEnergyY, electricEnergyZ,PF);

		}else if(CheckPtStatus!=1) {
      printf("CELIV method has not been setup to run from a chkpt file\n");
      printf("or anything other than the parameter.txt file as of now.\n");
      exit(1);
		}


	}

//printf("Printing Future site matrix at end of Pre_randomWalk\n");
//	if(FutureSite==NULL || rv==-1){
//		printf("FutureSite is NULL\n");
//		exit(1);
//	}

	if(FutureSite==NULL){
		printf("ERROR FutureSite is NULL\n");
		exit(1);
	}
	return 0;

}

int Post_randomWalk(ArbArray ClArLL, SNarray snA, Electrode elXb, Electrode elXf,\
		Electrode elYl, Electrode elYr, Electrode elZb, Electrode elZa,\
    double electricEnergyX, double electricEnergyY, double electricEnergyZ, ParameterFrame PF){


  int ClusterAlg = PFget_ClusterAlg(PF);

	if(snA==NULL){
		printf("ERROR snA found to be NULL\n");
		return -1;
	}

  
  printFileEnergy(snA, "Post", electricEnergyX,\
                  electricEnergyY, electricEnergyZ,PF);
	SNarray snAmini;
	matrix mtxmini;

	printf("Deleting Electrodes\n");
	if (elXb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXb);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXb);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elXb);
		elXb = NULL;
	}
	if (elXf!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXf);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXf);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elXf);
		elXf = NULL;
	}

	if (elYl!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYl);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYl);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elYl);
		elYl = NULL;
	}
	if (elYr!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYr);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYr);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elYr);
		elYr = NULL;
	}
	if (elZb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZb);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZb);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elZb);
		elZb = NULL;
	}

	if (elZa!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZa);
		deleteSNarray(&snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZa);
		deleteMatrix(&mtxmini);
		deleteElectrode(&elZa);
		elZa = NULL;
	}

	if(ClusterAlg>=1 && ClArLL!=NULL){
  //If no clusters are formed ClusterAlg could be
  //NULL
		printf("Deleteting Cluster Arbitrary Array Link List\n");
		deleteArbArray(&ClArLL);
	}
	printf("Deleting SiteNode array\n");
	deleteSNarray(&snA);

	return 0;
}

int initFutureSite( SNarray * snA, matrix * FutureSite,ChargeArray * chA, ParameterFrame PF,\
		Electrode elXb, Electrode elYl, Electrode elZb ){

	//printf("Initializing Future sites\n");
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

	//Declaring local variables
	int loop;
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

    ProjectChargePositionOntoSiteNodeReferenceFrame(&x,&y,&z,getCx(one),getCy(one),getCz(one),\
                                                    SLength, SWidth, SHeight);

		site = getSN(*snA,x,y,z);
		if( x==0  && XElecOn==1 ){
			HopToElecX(*snA, elXb, &one, &future, EndY, EndZ);
		}else if( y==0 && YElecOn==1){
			HopToElecY(*snA, elYl, &one, &future, EndX,EndZ);
		}else if( z==0 && ZElecOn==1){
			HopToElecZ(*snA, elZb, &one, &future, EndX, EndY);
		}else{
			SiteHop(*snA, &one, site, &future, EndX, EndY,EndZ,\
					XElecOn,YElecOn,ZElecOn, PeriodicX, PeriodicY, PeriodicZ);
		}
			
		setE(*FutureSite,loop+1,1,future);
	}

	return 0;
}

// You should know that ch is an active charge, as in it has not reached the 
// escape electrode. 
int checkDecay(ParameterFrame PF,\
               SNarray snA,\
               int SN_ID,\
               double KT,\
               Charge * ch,\
               double electricEnergyX,\
               double electricEnergyY,\
               double electricEnergyZ,\
               double MarcusCoef,
               const long double t){

  //printf("Check decay time %Lg\n",t);
  // Let's make sure that the site_ID is not an electrode and is between
  // 0 the length of the actual snA
  if(SN_ID<0 || SN_ID>=getAtotal(snA)) return 0;

  SiteNode sn2 = getSNwithInd(snA,SN_ID);

  // If site decay is allowed
  if(PFget_DecayOn(PF)==1){
    // Determine if the site the charge hopped to has decayed thus its energy
    // has changed 
    double prob_decay = (double)rand() / RAND_MAX;
    if(prob_decay<PFget_DecayProb(PF)){
      // This means the site potentially will decay
      int updateRate = Decay(sn2,PFget_DecayDisplacement(PF)); 

      // Update the rates around the site that decayed
      if(updateRate){
        updateNeigh_JumPossibility(electricEnergyX,\
                                   electricEnergyY,\
                                   electricEnergyZ,\
                                   PF,\
                                   MarcusCoef,\
                                   KT,\
                                   snA,\
                                   SN_ID);

        // Print the line
        int x, y, z;
        getLoc(&x,&y,&z,SN_ID,snA); 
        printFileDecay(x,y,z,t);
        return 1;
      }
    }
  }else if(PFget_DecayOn(PF)==2){
    // Determine if the site the charge hopped to has decayed thus its energy
    // has changed 
    double dw = getDwel(*ch);
    // Double determine how many how likely it is that the charge did not decay
    // in that time period
    double Time_iter = dw/PFget_DecayTime(PF);
    double NonDecayProb = pow((1-PFget_DecayProb(PF)),Time_iter);
    double DecayProbForDwell = 1-NonDecayProb;
    double prob_decay = (double)rand() / RAND_MAX;
    if(prob_decay<DecayProbForDwell){
      // This means the site potentially will decay
      // unless the site has already decayed
      int updateRate = Decay(sn2,PFget_DecayDisplacement(PF)); 

      // Update the rates around the site that decayed
      if(updateRate){
        updateNeigh_JumPossibility(electricEnergyX, electricEnergyY, electricEnergyZ, PF,MarcusCoef, KT,snA,SN_ID);
        int x, y, z;
        getLoc(&x,&y,&z,SN_ID,snA); 
        //printf("time decay %Lg\n",t);
        printFileDecay(x,y,z,t);
        //printf("Site Charge was on has decayed\n");
        //fflush(stdout);
        return 2;
      }
      
    }
  }
  //printf("No decay\n");
  //fflush(stdout); 
  return 0;
}

int randomWalk( SNarray snA,int CheckptNum,\
		char * FileName, double ElectricFieldX,\
		double ElectricFieldY, double ElectricFieldZ,\
		Electrode elXb, Electrode elXf, Electrode elYl,\
		Electrode elYr, Electrode elZb, Electrode elZa,\
		ParameterFrame PF,long double t,matrix Sequence,\
		matrix FutureSite,ChargeArray * chA,\
		long int n,int nc,int nca, double Temperature,\
    ArbArray * ClArLL, int DecayOn, double DecayProb, double DecayDisplacement){

  #ifdef _ERROR_CHECKING_ON_
	if(FutureSite==NULL){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR Future site matrix found to be NULL on entering randomWalk\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif

	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	time_t now;
	time_t later;
	double seconds;
	time(&now);

	//n - Number of time steps that have passed where charges are injected
	//nc - Number of charges in the system
	//nca - Number of charges that are currently active
	int    method       = PFget_method(PF);
	double Tcv          = PFget_Tcv(PF);
	double Vcv          = PFget_Vcv(PF);
	double Tlag         = PFget_Tlag(PF);
	double TStep        = PFget_TStep(PF);
	int    TCount       = PFget_TCount(PF);
  double R_neigh      = PFget_R_neigh(PF);
	int    Time_check   = PFget_Time_check(PF);
	int    Nstep_av     = PFget_Nstep_av(PF);
	int    NCh          = PFget_NCh(PF);
	double D            = PFget_D(PF);
	int    XElecOn      = PFget_XElecOn(PF);
	int    YElecOn      = PFget_YElecOn(PF);
	int    ZElecOn      = PFget_ZElecOn(PF);
	int    EndX         = PFget_EndX(PF);
	int    EndY         = PFget_EndY(PF);
	int    EndZ         = PFget_EndZ(PF);
	int    PeriodicX    = PFget_Px(PF);
	int    PeriodicY    = PFget_Py(PF);
	int    PeriodicZ    = PFget_Pz(PF);
	double SiteDistance = PFget_SiteDist(PF);
	int    MovieFrames  = PFget_MovieFrames(PF);
  
  long double CutOffTime         = (long double) PFget_CutOffTime(PF);
  int         ClusterAlg         = PFget_ClusterAlg(PF);
	//FILE * EndPtFile = NULL;
	//FILE * PathFile = NULL;
	//FILE * LogFile = NULL;

//if(PFget_EndPtFile(PF)==1){
//		EndPtFile = openEndPtFile(FileName);
//		printf("Opening End Pt File\n");
//		if(EndPtFile == NULL){
//			printf("End Pt File is NULL\n");
//			exit(1);
//		}
//	}
//	if(PFget_PathFile(PF)==1 && PFget_NumChargesTrack(PF)!=0){
//		PathFile = openPathFile(FileName);
//	}
	
	printf("will crash here if MovieFrames is 0\n");
	assert(MovieFrames!=0);

  #ifdef _ERROR_CHECKING_ON_
  /* Check inputs */
	if(snA==NULL){
    fprintf(stderr,"ERROR snA is NULL in randomWalk\n");
    return return_error_val();
  }
	if(t<0){
    fprintf(stderr,"ERROR t is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicZ>1){
    fprintf(stderr,"ERROR PeriodicZ is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicY>1){
    fprintf(stderr,"ERROR PeriodicY is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicX>1){
    fprintf(stderr,"ERROR PeriodicX is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicZ<0){
    fprintf(stderr,"ERROR PeriodicZ is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicY<0){
    fprintf(stderr,"ERROR PeriodicY is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(PeriodicX<0){
    fprintf(stderr,"ERROR PeriodicX is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(EndZ<0){
    fprintf(stderr,"ERROR EndZ is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(EndY<0){
    fprintf(stderr,"ERROR EndY is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(EndX<0){
    fprintf(stderr,"ERROR EndX is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(ZElecOn>1){
    fprintf(stderr,"ERROR ZElecOn is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(YElecOn>1){
    fprintf(stderr,"ERROR YElecOn is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(XElecOn>1){
    fprintf(stderr,"ERROR XElecOn is greater than 1 in randomWalk\n");
    return return_error_val();
  }
	if(ZElecOn<0){
    fprintf(stderr,"ERROR ZElecOn is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(YElecOn<0){
    fprintf(stderr,"ERROR YElecOn is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(XElecOn<0){
    fprintf(stderr,"ERROR XElecOn is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(D<0){
    fprintf(stderr,"ERROR D is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(NCh<0){
    fprintf(stderr,"ERROR number of charges is less than 0"
                   " in randomWalk\n");
    return return_error_val();
  }
	if(Nstep_av<0){
    fprintf(stderr,"ERROR Nstep_av is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(TStep<0){
    fprintf(stderr,"ERROR TStep is less than 0 in randomWalk\n");
    return return_error_val();
  }
	if(TCount<0){
    fprintf(stderr,"ERROR TCount is less than 0 in randomWalk\n");
    return return_error_val();
  }
  #endif  

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
	//printf("NTot %d NCh %d\n",Ntot,NCh);
	int flag;
  int SaveCount;
	//Movie start point
	int Movie = (int) ((double)t/((double)TStep*(double)Nstep_av));

	int x, xx, x1;
	int y, yy, y1;
	int z, zz, z1;

	int j;
	int k;

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
	//Number of failed hops because site is occupied
	long int FailedHop;
	long int TotalHopAttempt;
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
	double      ran;

	int TrackX;
	int TrackY;
	int TrackZ;

	long double TimeTrack1;
	int         NumAvgVel;
	long double TotalVelX;
	long double TotalVelY;
	long double TotalVelZ;

	double electricEnergyX;
	double electricEnergyY;
	double electricEnergyZ;
	double Vramp;

	double MarcusJ0;
	double MarcusCoeff;
	double KT;
  double Vx;
  double Energy;

	double tim;
	long double CELIV_totalT = (long double)(Tlag + Tcv);

	int ID;
	int ID2;
  int Global_ClusterID;   //Keeps track of the ids of the 
                                      //clusters as they are created
  Charge    one;
	Charge    two;
	Charge    three;
	SiteNode  site;

	NumAvgVel = 0;

	SLength = getAlen(snA);
	SWidth  = getAwid(snA);
	SHeight = getAhei(snA);

	//Calculate Current 
	matrix Xcurrent = newMatrix(8,1);
	matrix Ycurrent = newMatrix(8,1);
	matrix Zcurrent = newMatrix(8,1);

	//Keeps track of current in X, Y and Z direction
	//during time increment
	TotalX     = 0;
	TotalY     = 0;
	TotalZ     = 0;
	TrackX     = 0;
	TrackY     = 0;
	TrackZ     = 0;
	CurrentInc = 1;
	TimeTrack1 = 0;

	//Create source and drain for all electrodes that are turned on
	matrix Xelec_Drain  = newMatrix(8,1);
	matrix Xelec_Source = newMatrix(8,1);
	matrix Yelec_Drain  = newMatrix(8,1);
	matrix Yelec_Source = newMatrix(8,1);
	matrix Zelec_Drain  = newMatrix(8,1);
	matrix Zelec_Source = newMatrix(8,1);

	//Drift Velocities should be in units of [m/s]
	matrix Xvelocity = newMatrix(8,1);
	matrix Yvelocity = newMatrix(8,1);
	matrix Zvelocity = newMatrix(8,1);

	//Matrices holding energies of sites next to electrodes if CELIV is called
	matrix Xb1;
	matrix Xb2;
	matrix Xf1;
	matrix Xf2;

	SiteNode sn;

	matrix System = newMatrix(8,1);

	//Keeps track of the time when each of the matrices are incremented
	matrix timeArray = newMatrix(8,1);

	printf("Initialization Complete!\n");
	//Initilize number of charges reaching source
	//and drian to 0
	Xdrain  = 0;
	Xsource = 0;
	Ydrain  = 0;
	Ysource = 0;
	Zdrain  = 0;
	Zsource = 0;
	
	FailedHop       = 0;
	TotalHopAttempt = 0;

	TotalVelX = 0;
	TotalVelY = 0;
	TotalVelZ = 0;

	//Continue looping while n is less than the TCount
	//or if there are still charges in the system
	SaveCount        = 1;
	SaveTime         = 0;
	tim              = TStep;
  Global_ClusterID = 1;
  
  electricEnergyX = SiteDistance*ElectricFieldX;
  electricEnergyY = SiteDistance*ElectricFieldY;
  electricEnergyZ = SiteDistance*ElectricFieldZ;

  KT = kB*Temperature;
  //is equivalent to the marcus coefficient at 300 K
  MarcusJ0 = pow( PFget_AttemptToHop(PF)*hbar*pow(4*PFget_reOrg(PF)*kB*300/M_PI,1/2),1/2);
  //Calculating full Marcus Coefficient;
  MarcusCoeff = pow(MarcusJ0,2)/hbar* 
    pow(M_PI/(4*PFget_reOrg(PF)*KT),1/2)*
    exp(-2*PFget_gamma(PF)*PFget_SiteDist(PF));

  /* If CELIV method is specified calculate ramp rate */
	if(method==1){

		Vramp = Vcv/Tcv;
		//Need to define matrices containing energies of sites next to electrodes
		Xb1 = newMatrix(SWidth,SHeight);
		Xb2 = newMatrix(SWidth,SHeight);

		for(j=0;j<SWidth;j++){
			for(k=0;k<SHeight;k++){

				sn = getSN(snA,0,j,k);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Xb1,j+1,k+1,Energy);
				sn = getSN(snA,1,j,k);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Xb2,j+1,k+1,Energy);
				//printf("Energy %g\n",Energy);
			}
		}

		//Only saving energies of sites next to electrodes
		Xf1 = newMatrix(SWidth,SHeight);
		Xf2 = newMatrix(SWidth,SHeight);

		for(j=0;j<SWidth;j++){
			for(k=0;k<SHeight;k++){
				sn = getSN(snA,SLength-1,j,k);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Xf1,j+1,k+1,Energy);
				sn = getSN(snA,SLength-2,j,k);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Xf2,j+1,k+1,Energy);
				//printf("Energy %g\n",Energy);
			}
		}

	}

  matrix MasterM;
  //This is only needed if the second cluster Alg
  //is turned on
  if (ClusterAlg==2){
    MasterM = CalculateAllHops(snA, 
                               electricEnergyX, 
                               electricEnergyY, 
                               electricEnergyZ,
                               kB*Temperature,
                               PFget_reOrg(PF), 
                               SiteDistance, 
                               PFget_AttemptToHop(PF), 
                               PFget_gamma(PF),
                               PeriodicX, 
                               PeriodicY, 
                               PeriodicZ,
                               R_neigh);
  }

  //printf("Entering Loop\n");
  /*************************************************
   *                     KEY LOOP
   *************************************************/
  /* Charge transport simulation will run as long as:
   * 1 The number of charges in the system is not 0
   *   and we have not finished inserting charges
   * 2 The time has not reached the cutoff time
   * 3 The total CELIV ramp time has not been reached
   *   or we are using the TOF method
   */
	while( (n<TCount || nca>0)             &&  
         t<CutOffTime                    && 
         ((t<=CELIV_totalT && method==1) || 
         method==0)){

    //getchar();
		//If no charges have been inserted in the system we will
		//simply increment the time
		CheckConservationCharges(elXb, 
                             elYl, 
                             elZb,
                             nca, 
                             XElecOn, 
                             YElecOn,
                             ZElecOn, 
                             nc);
    //printf("global t %Lg\n",t); 

    if (nca==0){

			//Should be the Step time
			tim = TStep;
			//t is incremented after the charge hops
			//the global time is increased
			if(t==(t+(long double)tim)){
				printf("Exceeded precision t+tim==t\n");
				printf("t %Le tim %g\n",t,tim);
//				exit(1);
			}
			
			t += (long double) tim;

			SaveTime += (long double) tim;
			//The time the charge has been in the sample is updated
			Plust(one,(long double) tim);

			if(SaveTime >= (long double)TStep){
				//Here we check to see if we record the data
				//Nstep is used to determine how many timesteps pass
				//before recording
        printf("SaveTime %Lg TStep %g\n",SaveTime,TStep);
				//Here we will check the method if CELIV will update the site hop rates
				if(method==1){

					/* This is for CELIV */
					//Adjust Vx
					if((double)t<Tcv){
						Vx = Vramp*((double)t);
					}else{
						Vx = 0;
					}
					//Electric field from voltage
					ElectricFieldX = Vx / (((double)SLength)*SiteDistance);
					//Electrical energy from voltage between two sites
					electricEnergyX = SiteDistance*ElectricFieldX;
					if(electricEnergyX>10){
						printf("electricEnergyX %g Vx %g Vramp %g t %g\n",electricEnergyX,Vx,Vramp,(double)t);
						exit(1);
					}
					//Update hop rates of the core system
					initJumPossibility(electricEnergyX, 
                             electricEnergyY, 
                             electricEnergyZ,
                             MarcusCoeff, 
                             KT,
                             PFget_reOrg(PF), 
                             snA,
                             PeriodicX, 
                             PeriodicY, 
                             PeriodicZ, 
                             XElecOn, 
                             YElecOn, 
                             ZElecOn,
                             SiteDistance,
			                       R_neigh);
          /* Update hop rates of the electrodes */	
					Update_initJumPossibility_ElecX(electricEnergyX, 
                                          electricEnergyY,
                                          electricEnergyZ, 
                                          MarcusCoeff, 
                                          KT, 
                                          Xb1, 
                                          Xb2, 
                                          elXb, 
                                          0 , 
                                          PF);

					Update_initJumPossibility_ElecX(electricEnergyX, 
                                          electricEnergyY,
                                          electricEnergyZ, 
                                          MarcusCoeff, 
                                          KT, 
                                          Xf1, 
                                          Xf2, 
                                          elXf, 
                                          1 , 
                                          PF);

				}

				SaveCount++;

				if((SaveCount%Nstep_av)==0 ){
					SaveCount = 1;
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
					if(Movie<MovieFrames){
						printMovie(&Movie,t,FileName,snA,PF);
					}
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
			ChargeID = (int) getE(Sequence,1,1);
      //printf("\nNew Charge ID %d\n",ChargeID);
      for(int indx=1;indx<=getRows(Sequence);indx++){
        Charge four = getCharge(*chA,(int)getE(Sequence,indx,1));
        //printf("ChargeID %g DwellTime %g\n",getE(Sequence,indx,1),getDwel(four));

      }
      fflush(stdout);

			one = getCharge(*chA, ChargeID);
			if(getDwel(one)>1){
				printf("Time to large from getDwel %g\n",getDwel(one));
				exit(1);
			}

			if(getCx(one)<-1 || getCx(one)>SLength){
				printf("Charge position less than -1 or greater than SLength %d\n",SLength);
				exit(1);
			}

			// *********************REGARDLESS OF HOP OR NOT*******************************************

			//Attempt to make the charge hop to site 
			//that was previously determined. If site is
			//occupied return -1. If site unoccupied hop and
			//return 0

      //printf("Project Charge Position\n");
      //fflush(stdout);
      ProjectChargePositionOntoSiteNodeReferenceFrame(&x1,&y1,&z1,getCx(one),getCy(one),getCz(one),\
                                                      SLength,SWidth, SHeight);

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
//      printf("\nAfter Hop Charge ChargeID %d\n",ChargeID);
//      for(int indx=1;indx<=getRows(Sequence);indx++){
//        Charge four = getCharge(*chA,(int)getE(Sequence,indx,1));
//        printf("ChargeID %g DwellTime %g\n",getE(Sequence,indx,1),getDwel(four));
//    }
      


    

			//Function accounts for hops within system from electrodes
			//and two electrodes.
      //The function MakeHop is very important this is where the charge
      //acctually moves, the potential site is held in getE(FutureSite,ChargeID+1,1)
      //it will however only move if the future site is not already occupied.
      //flag - 0 if sucessful
      //flag - 1 if site is already occupied
      //printf("Making hop\n");
      //fflush(stdout);
			flag = MakeHop(snA, (int) getE(FutureSite,ChargeID+1,1),\
					&one, &TotalXtemp, &TotalYtemp, &TotalZtemp,\
					PF, KT, electricEnergyX, electricEnergyY, electricEnergyZ,
          MarcusCoeff,t);

			//Get the time it took to make the hop
      tim = getDwel(one);
		  //printf("tim %g from charge %g\n",tim,getE(Sequence,1,1));	
			if(tim>1){
				printf("time to large %g for charge located at (%d,%d,%d)\n",tim,PrevX,PrevY,PrevZ);
				printChargeA(*chA);
				exit(1);
			}else if(tim<=0){
				printf("time equal or less than 0 tim %g\n",tim);
				printf("time equal or less than 0 tim %g\n",getDwel(one));
				printf("Charge ID %d\n",ChargeID);
				printChargeA(*chA);
        if(one==NULL){
          printf("Charge is NULL\n");
        }
        exit(1);
			}
			//Cannot simply increase the time by the dwelstat if not
			//all the charges have yet been inserted

			if (((n+1)*(long int)NCh)<=Ntot && tim>TStep){
				//This means at maximum can only increase the time by the 
				//TStep because more charges need to be inserted. 
				tim = TStep;
        //printf("tim changed to %g\n",tim);
				//Because the site does not hop we just make it progress
				//in time a little we set the flag to -1
				flag = -1;
			}

			TotalHopAttempt++;

			if(flag==0){
				//it did hop
				
				//If tracking charges and charge id is one of the ones that is being 
				//tracked record information in the .path file
				if((ChargeID)<(PFget_NumChargesTrack(PF)) && (PFget_PathFile(PF))!=0){
					printToPathFile( FileName,getCx(one),getCy(one),getCz(one),ChargeID,tim,t);
				}
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
        
        /* This is only applicable for the second
         * cluster algorithm 
         */
        if(ClusterAlg==2 && method==0){
         //At this point we want to ignore the electrodes
          //the electrodes all have id's greater than getAtotal(snA)
			//printf("Entering ClusterChargePath\n");
	        ClusterChargePath(one,
                            ChargeID,
                            FutureSite,
                            snA,
                            MasterM,
                            PF,
                            ClArLL,
                            &Global_ClusterID);
           
        }//End of ClusterAlg==2 segment
			}else{
				//Hop failed because site was occupied
				//printf("Failed Hop %d flag %d\n",FailedHop,flag);
        FailedHop++;
			}
				
			//t is incremented after the charge hops
			//the global time is increased
			if(t==(long double)t+tim){
				fprintf(stderr,"ERROR Exceeded precision t==t+tim\n");
				fprintf(stderr,"t %Le tim %g\n",t,tim);
//				exit(1);
			}
			t += (long double) tim;

			if(tim==0){
				fprintf(stderr,"Error tim is 0\n");
				fprintf(stderr,"TimeTrack1 %Lg t %Lg\n",TimeTrack1,t);
//				exit(1);
			}
			SaveTime += (long double)tim;
			//The time the charge has been in the sample is updated
			Plust(one,(long double) tim);

			// If the global time has reached the next time step 
			// calculate current

			//printf("SaveTime %Lg TStep %g\n",SaveTime, TStep);
			//printf("SaveTime %Lg TStep %g SaveCount %d Nstep_av %d Movie %d MovieFrames %d\n",SaveTime,TStep,SaveCount,Nstep_av,Movie,MovieFrames);


			if(SaveTime >= (long double)TStep){
				//Here we check to see if we record the data
				//Nstep is used to determine how many timesteps pass
				//before recording
				SaveCount++;

				//Here we will check the method if CELIV will update the site hop rates
				if(method==1){

					//Adjust Vx
					if((double)t<Tcv){
						Vx = Vramp*((double)t);
					}else{
						Vx = 0;
					}
					//Electric field from voltage
					ElectricFieldX = Vx / (((double)SLength)*SiteDistance);
					//Electrical energy from voltage between two sites
					electricEnergyX = SiteDistance*ElectricFieldX;
					if(electricEnergyX>10){
						printf("electricEnergyX %g Vx %g Vramp %g t %g\n",electricEnergyX,Vx,Vramp,(double)t);
						exit(1);
					}
					//Update hop rates
					//printf("Vramp %g t %g Vx %g SLength %d SiteDistance %g electricField %g electricEnergyX %g\n",Vramp,(double)t,Vx,SLength,SiteDistance,ElectricFieldX,electricEnergyX);				
					//exit(1);
					
					initJumPossibility(electricEnergyX, electricEnergyY, electricEnergyZ,\
							MarcusCoeff, KT,PFget_reOrg(PF), snA,\
							PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn,\
              SiteDistance, R_neigh);
					
					Update_initJumPossibility_ElecX( electricEnergyX, electricEnergyY,\
							electricEnergyZ, MarcusCoeff, KT, Xb1, Xb2, elXb, 0 , PF);
					Update_initJumPossibility_ElecX( electricEnergyX, electricEnergyY,\
							electricEnergyZ, MarcusCoeff, KT, Xf1, Xf2, elXf, 1 , PF);

				}

        if(PFget_AvgChargeEnergyFile(PF)){
          printFileChargeEnergy(snA, *chA, Sequence, nca, t,PF);
        }

				if((SaveCount%Nstep_av)==0){
					//Saving Data
					//printf("Saving Data\n");
					SaveCount = 1;
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

					if(Movie<MovieFrames){
						//printf("Should have entered the printMovie Routine\n");
						printMovie(&Movie,t,FileName,snA,PF);
					}

				}

				//printf("Updating time\n");
				time(&later);
				seconds = difftime(later,now);

				if(seconds>(double)(Time_check*60)){
					time(&now);
					//printf("Saving Checkpt file\n");
					Save_CheckPt(FileName, &CheckptNum, snA, *chA, Sequence,FutureSite, t, PF,n, nc, nca,\
							elXb, elXf, elYl, elYr, elZb, elZa);
				}

				SaveTime = 0;
			}

			//For all the charges still active 
			//need to decrease their dwelltime
			UpdateOccTime(&snA,&one,tim, PF);
      int decayed_sites = 0;
      matrix decayed_Sites = newMatrix(1,1);
      // Do not minus from site 1 because tim is the old hop time whereas
      // site 1 now has the new hop time 
      //if(tim>getDwel(one)){
      //  printf("Site thats moving %g\n",getE(Sequence,1,1));
      //}
			for( i = 1; i<=nca; i++ ){
				two = getCharge(*chA,(int)getE(Sequence,i,1));
        //printf("Grabbing all charges %d\n",(int)getE(Sequence,i,1));
				/*if(getDwel(two)>1){
					printf("Dwelltime excessive %g\n",getDwel(two));
					exit(1);
				}*/
				//printf("tim to be minused %g Charge id %g\n",tim,getE(Sequence,i,1));
				//printChargeA(*chA);
				MinusDwel(two,tim);
        // For every charge we need to check to see if the charge
        // decayed or not. If it does we also need to recalculate 
        // The dwell time
        int xxx, yyy, zzz;
        if(PFget_Px(PF)){
          xxx = (getCx(two)+getAlen(snA))%getAlen(snA);
          if(xxx<0){
            xxx = getAlen(snA)+xxx;
          }
        }else{
          xxx = getCx(two);
        }
        if(PFget_Py(PF)){
          yyy = (getCy(two)+getAwid(snA))%getAwid(snA);
          if(yyy<0){
            yyy = getAwid(snA)+yyy;
          }
        }else{
          yyy = getCy(two);
        }
        if(PFget_Pz(PF)){
          zzz = (getCz(two)+getAhei(snA))%getAhei(snA);
          if(zzz<0){
            zzz = getAhei(snA)+zzz;
          }
        }else{
          zzz = getCz(two);
        }
        if(xxx>=0 && xxx<getAlen(snA) &&
           yyy>=0 && yyy<getAwid(snA) &&
           zzz>=0 && zzz<getAhei(snA)){
          int decay_val = checkDecay(PF,
                                     snA,
                                     getIndex(snA,xxx,yyy,zzz),
                                     KT,
                                     &two,
                                     electricEnergyX,
                                     electricEnergyY,
                                     electricEnergyZ,
                                     MarcusCoeff,
                                     t);
          //printf("decay_val %d\n",decay_val);
          //fflush(stdout);
          if(decay_val==2){
            // This means the site decayed so we need to update the time
            // the charge will spend on the site. 
          
            double tim_new = 1/getsum(getSN(snA,xxx,yyy,zzz));
            //printf("tim_new %g\n",tim_new);
            //if(tim_new<0){
            //  printf("tim_new is less than 0\n");
            //}
            double new_dwell = -log((double) rand()/RAND_MAX)*tim_new;
            //printf("new_dwell %g\n",new_dwell);
            //fflush(stdout);
            setDwel(two,new_dwell );
            //printf("Setting Dwell\n");
            //fflush(stdout);
            //updateSequence(nca, *chA, &Sequence, i);
            decayed_sites++;
            if(decayed_sites>getRows(decayed_Sites)){
              //printf("Resizing\n");
              //fflush(stdout);
              resizeRow(&decayed_Sites,getRows(decayed_Sites)+1);
            }
            setE(decayed_Sites,decayed_sites,1,getE(Sequence,i,1)); 
            //printf("setting decay id of Charge %d with new dwell time %g\n",(int)getE(Sequence,i,1),new_dwell);
            //fflush(stdout);
          } 
        
        }
        //printf("Updating OccTime\n");
        //fflush(stdout);
				UpdateOccTime(&snA, &two,tim, PF);
			}
      if(decayed_sites>0){
        for(int dec=1;dec<=getRows(decayed_Sites);dec++){
          //printf("Updating placement of site %d\n",(int)getE(decayed_Sites,dec,1));
          printf("Site that decayed %g\n",getE(decayed_Sites,dec,1));
          updateSequence(nca,*chA,&Sequence,dec,decayed_Sites);
        }
        checkSequence(nca,*chA,&Sequence);
      }
      deleteMatrix(&decayed_Sites);
      
			// *********************************************************************

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
						//printf("Hopping back into system current x %d\n",getCx(one));
						ChargeSystem( &nc, elXb, nca, Sequence, *chA);
					}else if(PrevY == -1 && EndY!=0){
						//printf("Hopping back into system current y %d\n",getCy(one));
						ChargeSystem( &nc, elYl, nca, Sequence, *chA);
					}else if(PrevZ == -1 && EndZ!=0){
						//printf("Hopping back into system current z %d\n",getCz(one));
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

			//Check if a charge has arrived at the front or back electrode
			//If it has updates the following parameters:
			//	Velocities
			//	NumAvgVel
			//	Drain
			//	Source
			//	Number of charges on respective electrodes
			//	TimeTrack1 
			//	t

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


       ProjectChargePositionOntoSiteNodeReferenceFrame(&x, &y, &z, xx, yy, zz,\
                                                     SLength, SWidth, SHeight);

			////////////////////////////////////////////////////////////////////////
			//Determining Next hop
			//printf("Determining Next Hop\n");
			//Before future hop has been determined, tim of hop is 0
			tim=0;
			//printf("Middle Value of t %Lg\n",t);
			//Did not hop to electrode
			if(ElecExit==0){
        //printf("No hop\n");
        //Site jumped too
        site = getSN(snA, x,y,z);
        //Set the DwelStat of new site with the id of 
        //the charge that just jumped
        setDwelStat(site,ChargeID);

        //Check to see if new site is part of a cluster
        ClusterYes=getType(site);
        //Determine future hopping site
     /*   if(ClusterYes==1){
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
              HopToElecX(snA, elXb, &one, &future, EndY, EndZ);
              getLoc(&x,&y,&z,future,snA);
              tim = 1/getsum(getSN(snA,x,y,z));
          printf("6 ste\n");
            }else if( y==0 && YElecOn==1){
              HopToElecY(snA, elYl, &one, &future, EndX, EndZ);
              getLoc(&x,&y,&z,future,snA);
              tim = 1/getsum(getSN(snA,x,y,z));
          printf("5 ste\n");
            }else if( z==0 && ZElecOn==1){
              HopToElecZ(snA, elZb, &one, &future, EndX, EndY);
              getLoc(&x,&y,&z,future,snA);
              tim = 1/getsum(getSN(snA,x,y,z));
          printf("4 ste\n");
            }else{
              
          printf("3 ste\n");
              SiteHop(snA, &one, site, &future, EndX, EndY,EndZ,\
                  XElecOn,YElecOn,ZElecOn, PeriodicX,PeriodicY, PeriodicZ);
              getLoc(&x,&y,&z,future,snA);
              tim = 1/getsum(getSN(snA,x,y,z));
            }

          }else{
            //Future Hop To Cluster or Site using cluster algorithm
          printf("2 ste\n");
            ClusterHop(snA, &one, &tim, &future);
            getLoc(&x,&y,&z,future,snA);
          }

        }else{*/
         //Future Hop To Cluster or Site
          SiteHop(snA, &one, site, &future, EndX, EndY,EndZ,\
              XElecOn,YElecOn,ZElecOn, PeriodicX,PeriodicY, PeriodicZ);
          getLoc(&x,&y,&z,future,snA);
          tim = 1/getsum(getSN(snA,x,y,z));
       // }

        //Having chosen future site recording it
        if(future==getAtotal(snA)){
          printf("ERROR future must be less than total sites future %d total %d",future,getAtotal(snA));
          exit(1);
        }

        setE(FutureSite,ChargeID+1,1,future);

        do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
        //The waiting time of the charge is updated based on
        //it's location and a random number
        //printf("Setting dwel site %lg\n",tim);
        if(-log((double)ran/RAND_MAX)*tim<=0){
          printf("ERROR dweltime is less than or equal to 0 and tim is %g\n",tim);
          exit(1);
        }
        setDwel(one, -log((double) ran/RAND_MAX)*tim);
        

        if(tim>1){
          printf("tim %g\n",tim);
          printf("WARNING setDwel huge %g\n",getDwel(one));
        }
      
        //The location of the charge in the sequence is updated
        insertDwelltimePos(nca, *chA, &Sequence);
      }else{
				//printf("Yes Hopped to Electrode\n");
				//If the charge jumped to an electrode it will now have to jump
				//from the electrode back into the system
        if(xx == -1 && EndX!=0){
					//Hopping from back Electrode
					//Updating dwell time of charge
					ElecHopOffX(elXb, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					//printf("Hopped to X Source x position %d flag %d\n",xx, flag);
					ChargeElectrode(elXb, &one, &Sequence, *chA, &nc, nca, flag);
					//printf("ChargeElectrode value of t %Lg\n",t);
					//printf("Grabbing dwel %g\n",getDwel(one));

				}else if(xx == EndX*SLength && EndX!=0){
					//Hopped to front electrode
					//Setting dwell time to large number removing charge from system
					
					//Saving information to EndPtFile for later calculating the charge
					//mobility
					if(PFget_EndPtFile(PF)!=0){
						//printf("Print to end Pt file X\n");
						//if(EndPtFile==NULL){
						//	printf("End pt file is found to be NULL\n");
						//	exit(1);
						//}
						printToEndPtFile( FileName, xx, yy, zz, ChargeID, t);
						//printf("Finished Printing to end pt file\n");
					}
					ChargeClosure(elXf, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);
					
					ID2 = (int)getE(Sequence,1,1);

					for(i=0; i<Ntot-1;i++){
						//Updating the waiting queue of charges
						ID = (int) getE(Sequence,i+2,1);
						setE(Sequence,i+1,1,ID);
					}

					setE(Sequence,Ntot,1,ID2);

				}else if(yy == -1 && EndY!=0){
					//Hopping from left Electrode
					ElecHopOffY(elYl, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					//printf("Hopped to Y Source y position %d\n",yy);
					ChargeElectrode(elYl, &one, &Sequence, *chA, &nc, nca, flag);
				}else if(yy == EndY*SWidth && EndY!=0 ){
					//Hopped to right electrode
					
					//Saving information to EndPtFile for later calculating the charge
					//mobility
					if(PFget_EndPtFile(PF)!=0){
						//printf("Print to end Pt file Y\n");
						printToEndPtFile( FileName, xx, yy, zz, ChargeID, t);
						//printf("Finished Printing to end pt file\n");
					}
					ChargeClosure(elYr, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);
					ID2 = (int)getE(Sequence,1,1);

					for(i=0; i<Ntot-1;i++){
						//Updating the waiting queue of charges
						ID = (int) getE(Sequence,i+2,1);
						setE(Sequence,i+1,1,ID);
					}

					setE(Sequence,Ntot,1,ID2);


				}else if(zz == -1 && EndZ!=0){
					//Hopping from bottom Electrode
					ElecHopOffZ(elZb, &future, snA);
					setE(FutureSite,ChargeID+1,1,future);
					//printf("Hopped to Z Source z position %d\n",zz);
					ChargeElectrode(elZb, &one, &Sequence, *chA, &nc, nca, flag);
				}else if(zz == EndZ*SHeight && EndZ!=0){
					//Hopped to top electrode
					
					//Saving information to EndPtFile for later calculating the charge
					//mobility
					if(PFget_EndPtFile(PF)!=0){
						//printf("Print to end Pt file Z\n");
						printToEndPtFile( FileName, xx, yy, zz, ChargeID, t);
						//printf("Finished Printing to end pt file\n");
					}
					
					ChargeClosure(elZa, &one, &Sequence, &nc, &nca, &TotalCollected, Ntot);
					//Sequence contains the ids of the charges which go
					//from 0 to Ntot-1;
					//Placing the charge that just reached the electrode to 
					//the back of the list
					ID2 = (int)getE(Sequence,1,1);

					for(i=0; i<Ntot-1;i++){
						//Updating the waiting queue of charges
						ID = (int) getE(Sequence,i+2,1);
						setE(Sequence,i+1,1,ID);
					}

					setE(Sequence,Ntot,1,ID2);

				}
			}

		}	


		////////////////////////////////////////////////////////////////////////
		//Regardless of whether a charge hopped or not need to 
		//see if new charges need to be injected
		//printf("Regardless of whether a charge hopped or not\n");

		//printf("Bottom Value of t %Lg\n",t);
		if(t >= ((long double)n)*(long double)TStep && ((n+1)*(long int)NCh)<=Ntot){

			//increment n only up to certain point
			if( (long double)n*TStep<=t || ((n+1)*(long int)NCh)<Ntot){
				n++;
			}
			//printf("n %d NCh %d Ntot %d\n",n,NCh,Ntot);	
			//Initialize new charges that are inserted
			initCharge( nca, n, chA, Sequence, snA,\
					Ntot, NCh,\
					XElecOn, YElecOn, ZElecOn,\
					EndX, EndY, EndZ);

      
			//Add to the number of charges that are 
			//already in the system

			TotalVelX += (long double)TrackX*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			TotalVelY += (long double)TrackY*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			TotalVelZ += (long double)TrackZ*(long double)SiteDistance/((t-TimeTrack1)*(long double)nc);
			if(isnan(TotalVelX)){
				printf("TotalVelX is nan! t %g n %ld TStep %g NCh %d Ntot %d\n",(double)t,n,TStep,NCh,Ntot);
				exit(1);
			}	
			TrackX = 0;
			TrackY = 0;
			TrackZ = 0;
			NumAvgVel++;
			//printf("TimeTrack1 set to t\n");
			TimeTrack1 = t;
			nc = nc+NCh;
			nca = nca+NCh;
		}

	}

	//If simulation ended before all charges left the system we
	//need to record their positions for accurately calculating 
	//the charge mobility
	if(PFget_EndPtFile(PF)!=0){
		int loop;	
		for(loop = 0; loop<NCh; loop++){
			three = getCharge(*chA,loop);
			xx = getCx(three);
			yy = getCy(three);
			zz = getCz(three);

			if((xx < EndX*SLength && EndX!=0) || EndX==0){
				if((yy < EndY*SWidth && EndY!=0) || EndY==0){
					if((zz < EndZ*SHeight && EndZ!=0) || EndZ==0){
						printToEndPtFile( FileName, xx,yy,zz,loop,t);
					}
				}
			}
		}
		//closeEndPtFile(EndPtFile);
	}

	//Close PathFile
	//if(PFget_PathFile(PF)!=0){
	//		closePathFile( PathFile);
	//}

	//if(PFget_LogFile(PF)!=0){
	//	LogFile = openLogFile(FileName);
		LogFile_printHops( FileName, TotalHopAttempt,FailedHop);
	//	closeLogFile( LogFile);
	//}

	Movie = Movie+1;

	printTransportData(System, timeArray, Xcurrent, Ycurrent, Zcurrent,\
			Xelec_Drain, Yelec_Drain, Zelec_Drain,\
			Xelec_Source, Yelec_Source, Zelec_Source,\
			Xvelocity, Yvelocity, Zvelocity,\
			XElecOn, YElecOn, ZElecOn, FileName,\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ);

	//printf("Deleting matrices\n");
	if(method==1){
		deleteMatrix(&Xb1);
		deleteMatrix(&Xb2);
		deleteMatrix(&Xf1);
		deleteMatrix(&Xf2);
	}
  if (ClusterAlg==2){
    deleteMatrix(&MasterM);
  }
	deleteChargeA(*chA);
	deleteMatrix(&timeArray);
	deleteMatrix(&Xcurrent);
	deleteMatrix(&Ycurrent);
	deleteMatrix(&Zcurrent);
	deleteMatrix(&Xelec_Drain);
	deleteMatrix(&Yelec_Drain);
	deleteMatrix(&Zelec_Drain);
	deleteMatrix(&Xelec_Source);
	deleteMatrix(&Yelec_Source);
	deleteMatrix(&Zelec_Source);
	deleteMatrix(&Xvelocity);
	deleteMatrix(&Yvelocity);
	deleteMatrix(&Zvelocity);
  deleteMatrix(&System);
  deleteMatrix(&Sequence);
	deleteMatrix(&FutureSite);
	return 0;
}

// Updates the time that a charge spend on a node
int UpdateOccTime(SNarray * snA,Charge * ch, double time, ParameterFrame PF){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL ){
    fprintf(stderr,"ERROR snA is NULL in UpdateOccTime\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(ch==NULL){
    fprintf(stderr,"ERROR ch is NULL in UpdateOccTime\n");
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	int SLength;
	int SWidth;
	int SHeight;
	int i;
	int j;
	int k;
	int XElecOn;
	int YElecOn;
	int ZElecOn;

	SLength = PFget_Len(PF);
	SWidth  = PFget_Wid(PF);
	SHeight = PFget_Hei(PF);

	i = getCx(*ch);
	j = getCy(*ch);
	k = getCz(*ch);

	XElecOn = PFget_XElecOn(PF);
	YElecOn = PFget_YElecOn(PF);
	ZElecOn = PFget_ZElecOn(PF);

	if(XElecOn==0){
		i=(i+getAlen(*snA))%(getAlen(*snA));
    if(i<0){
      i=getAlen(*snA)+i;
    }	
	}
	if(YElecOn==0){
		j=(j+getAwid(*snA))%(getAwid(*snA));
    if(j<0){
      j=getAwid(*snA)+j;
    }	
	}
	if(ZElecOn==0){
		k=(k+getAhei(*snA))%(getAhei(*snA));
    if(k<0){
      k=getAhei(*snA)+k;
    }	
	}

	if(i>=0 && i<SLength &&  j>=0 && j<SWidth && k>=0 && k<SHeight){

		SiteNode sn = getSN(*snA,i,j,k);
		if(sn==NULL){
			printf("ERROR sitenode is NULL %d i %d j %d k\n",i,j,k);
			exit(1);
		}
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
		printf("ERROR FileName has not been initialized to \\0\n");
		return -1;
	}

	struct stat st;
	char   FileName[256];
	char   FileEx[20];
	char   FileNameKeep[256];
	int    len;
	int    i;

	int CheckFileExist;
	//Determines if there is a checkpoint file
	//What is the highest number of the checkpointfile

	//Initialize file to Null pointer
	File[0]        = '\0';
	CheckFileExist = 0;

	if(stat("CHECKPOINT",&st)==0){
		//CHECKPOINT directory does exist
		DIR    *d;
		struct dirent *dir;
		d = opendir("CHECKPOINT");
		if(d) {

			printf("CHECKPOINT directory exists!\n");
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

int CheckPt_Cluster_TOF(const double Vx,const double Vy,const double Vz,const double T, int r){

	if(T < 0){
		printf("ERROR Temperature negative\n");
		return -1;
	}

	struct stat st;
	char   FileName[256];
	FILE   *file;
	
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


int CheckPt_Latest_TOF(char * FileNameFull, int buffersize, const double Vx,const double Vy,const double Vz,const double T){

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
	FileName[0]     = '\0';
	FileNameFull[0] = '\0';
	CheckFileExist  = 0;

	ChVx = 0.0;
	ChVy = 0.0;
	ChVz = 0.0;
	ChT  = 0.0;

	if(stat("CHECKPOINT",&st)==0){
		//CHECKPOINT directory does exist
		DIR *d;
		struct dirent *dir;
		d = opendir("CHECKPOINT");
		if(d) {

			num  = 0;
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

int CheckPt_Latest_CELIV(char * FileNameFull, int buffersize,const double T){

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
	FileName[0]     = '\0';
	FileNameFull[0] = '\0';
	CheckFileExist  = 0;

	ChVx = 0.0;
	ChVy = 0.0;
	ChVz = 0.0;
	ChT  = 0.0;

	if(stat("CHECKPOINT",&st)==0){
		//CHECKPOINT directory does exist
		DIR *d;
		struct dirent *dir;
		d = opendir("CHECKPOINT");
		if(d) {

			num  = 0;
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


						if(num>keep && ChT==T && ChVx==0 && \
								ChVy==0 && ChVz==0){

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
/*
int Load_CheckPt(long double * t     , SNarray * snA      , ChargeArray * chA     ,\
                 matrix * Sequence   , matrix * FutureSite, char * FileName       ,\
                 ParameterFrame *PF  , long int *n        , int * nc              ,\
                 int *nca            , int * Num_elXb     , int * Num_elXf        ,\
                 int * Num_elYl      , int * Num_elYr     , int * Num_elZb        ,\
                 int * Num_elZa      ){


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
	int      intvar;
	int      intvar2;
	double   doublevar;

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
			PFset_method(*PF,intvar);
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
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_ClusterAlg(*PF,intvar);
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
			
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Vcv(*PF,intvar);
			
			fscanf(CheckIn,"%s %d",bufRead,&intvar);	
			PFset_Tcv(*PF,intvar);
			

			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);
			fgets(bufRead,256,CheckIn);

			//Scanning in SN information
			fscanf(CheckIn,"%lf %ld %d %d",&doublevar,&intvar0, &intvar, &intvar2);
			*t = (long double) doublevar;
			*n = intvar0;
			*nc = intvar;
			*nca = intvar2;
			//printf("nc %d nca %d\n",*nc, *nca);
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
					PFget_XElecOn(*PF), PFget_YElecOn(*PF),PFget_ZElecOn(*PF),
          SiteDistance,R_neigh);

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
        
        if(dwellTime<=0){
          printf("ERROR dweltime is less than or equal to 0\n");
          exit(1);
        }
				setDwel(ch,dwellTime);
				i++;
				setE(*Sequence,i,1,ChargeID);
				intvar = setE(*FutureSite,ChargeID+1,1,Future);
			}

			printf("Last value of rv %d\n",intvar);
			printf("Load_CheckPt Printing FutureSite");

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
*/

/*
int Load_CheckPt_Data_TOF(long double * t, SNarray * snA, ChargeArray * chA, matrix * Sequence,\
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
			//printf("nc %d nca %d\n",*nc, *nca);
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
        if(dwellTime<=0){
          printf("ERROR dwelltime is less than or equal to 0\n");
          exit(1);
        }
				setDwel(ch,dwellTime);
				i++;
				setE(*Sequence,i,1,ChargeID);
				intvar = setE(*FutureSite,ChargeID+1,1,Future);

			}

			printf("Last value of rv %d\n",intvar);
			printf("Load_CheckPt Printing FutureSite");

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
			fscanf(CheckIn,"%s %d",bufRead,&intvar);
			PFset_ClusterAlg(*PF,intvar);
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
*/
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
		fprintf(CheckOut,"Method %d\n", PFget_method(PF));
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
	  
		//	 fprintf(CheckOut,"//Define the work function of the electrodes [eV]\n");
		//	 fprintf(CheckOut,"workX %f\n",PFget_);
		//	 fprintf(CheckOut,"workY %f\n",PFget_);
		//	 fprintf(CheckOut,"workZ %f\n\n",PFget_);
		//	 fprintf(CheckOut,"//Define the electron affinity of medium [eV]\n");
		//	 fprintf(CheckOut,"electronAffin %f\n",electronAffin);
		//	 fprintf(CheckOut,"//Define the ionization energy of medium [eV]\n");
		//	 fprintf(CheckOut,"ionizationEnergy %f\n\n",ionizationEnergy);
		 
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

		fprintf(CheckOut,"//Turn on Cluster Algorithm 2, Leave off 0, whatever is optimal 1\n");
		fprintf(CheckOut,"ClusterAlg %d\n\n",PFget_ClusterAlg(PF));
		
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
		fprintf(CheckOut,"MovieFrames %d\n",PFget_MovieFrames(PF));
		fprintf(CheckOut,"Vcv %g\n",PFget_Vcv(PF));
		fprintf(CheckOut,"Tcv %g\n\n",PFget_Tcv(PF));
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
			ID = (int) getE(Sequence,element+1,1);
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


int printMovie(int * Movie,long double t, char * FileName, SNarray snA, ParameterFrame PF){

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

		fprintf(Mout,"%d\t%Lg\t%d\t%d\t%d\t%g\n\n",getAtotal(snA),t,PFget_Len(PF),PFget_Wid(PF),PFget_Hei(PF),PFget_SiteDist(PF));

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
//			printf("Number of Rows of timeArray %d\n",getRows(timeArray));
//			printf("Number of Rows of Xcurrent %d\n",getRows(Xcurrent));
//			printf("Number of Rows of Xvelocity %d\n",getRows(Xvelocity));
//			printf("Number of Rows of System %d\n",getRows(System));
//			printf("Number of Rows of Xelec_Drain %d\n",getRows(Xelec_Drain));
//			printf("Number of Rows of Xelec_Source %d\n",getRows(Xelec_Source));
//      
//   
			for(i=1;i<=getRows(timeArray);i++){

				if(getE(timeArray,i,1)!=0.0){

					if(ElectricFieldX!=0){
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Xout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),\
								getE(Xelec_Source,i,1),getE(Xelec_Drain,i,1),getE(System,i,1),getE(Xvelocity,i,1),getE(Xvelocity,i,1)/ElectricFieldX*10000);
					}else{
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] 
						fprintf(Xout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),\
								getE(Xelec_Source,i,1),getE(Xelec_Drain,i,1),getE(System,i,1),getE(Xvelocity,i,1),0.0);
					}
				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0.0){
					if(ElectricFieldX!=0){
						//Time [s] Xcurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Xout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),getE(Xvelocity,i,1),getE(Xvelocity,i,1)/ElectricFieldX*1000);
					}else{
						//Time [s] Xcurrent [Amps]  DriftVelocity [m/s]
						fprintf(Xout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Xcurrent,i,1),getE(Xvelocity,i,1),0.0);
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
				if(getE(timeArray,i,1)!=0.0){
					if(ElectricFieldY!=0){
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Yout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),\
								getE(Yelec_Source,i,1),getE(Yelec_Drain,i,1),getE(System,i,1),getE(Yvelocity,i,1),getE(Yvelocity,i,1)/ElectricFieldY*1000);
					}else{
						//Time [s] Xcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s]
						fprintf(Yout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),\
								getE(Yelec_Source,i,1),getE(Yelec_Drain,i,1),getE(System,i,1),getE(Yvelocity,i,1),0.0);
					}

				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0.0){
					if(ElectricFieldY!=0){
						//Time [s] Ycurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Yout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),getE(Yvelocity,i,1),getE(Yvelocity,i,1)/ElectricFieldY*1000);
					}else{
						//Time [s] Ycurrent [Amps]  DriftVelocity [m/s]
						fprintf(Yout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Ycurrent,i,1),getE(Yvelocity,i,1),0.0);
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
				if(getE(timeArray,i,1)!=0.0){
					if(ElectricFieldZ!=0){
						//Time [s] Zcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Zout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),\
								getE(Zelec_Source,i,1),getE(Zelec_Drain,i,1),getE(System,i,1),getE(Zvelocity,i,1),getE(Zvelocity,i,1)/ElectricFieldZ*1000);
					}else{
						//Time [s] Zcurrent [Amps] Source [unitless] Drain [unitless] System [unitless] DriftVelocity [m/s] 
						fprintf(Zout,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),\
								getE(Zelec_Source,i,1),getE(Zelec_Drain,i,1),getE(System,i,1),getE(Zvelocity,i,1),0.0);
					}
				}
			}
		}else{
			for(i=1;i<=getRows(timeArray);i++){	
				if(getE(timeArray,i,1)!=0.0){
					if(ElectricFieldZ!=0){
						//Time [s] Zcurrent [Amps]  DriftVelocity [m/s] Mobility [cm2/Vs]
						fprintf(Zout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),getE(Zvelocity,i,1),getE(Zvelocity,i,1)/ElectricFieldZ*1000);
					}else{
						//Time [s] Zcurrent [Amps]  DriftVelocity [m/s]
						fprintf(Zout,"%g \t %g \t %g \t %g\n",getE(timeArray,i,1),getE(Zcurrent,i,1),getE(Zvelocity,i,1),0.0);
					}
				}
			}
		}
		//printf("Closing Z file identifier\n");
		fclose(Zout);
	}

	//printf("Returning with a value of 0\n");
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
	//printf("Charge reached electrode value of nc now %d\n",*nc);
	//Total number of collected charges increases
	*TotalCollected = *TotalCollected+1;
	//Set Dwell time to finished state
	//printf("Set Dwel exit\n");
	setDwel(*one,1E6);

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
		//printf("Charge reached electrode value of nc now %d\n",*nc);
	}
	//Set Dwell time to that of the electrode
	double tim = 1/getElectrode_Sum(el); 
  //if(-log((double)rand()/RAND_MAX)*tim<=0){
  //  printf("ERROR in ChargeElectrode dwelltime less than or equal to 0\n");
  //  exit(1);
  //}
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

/* update the whole sequence */
int updateSequence(const int nca, ChargeArray chA, matrix * Sequence, int dec, matrix Decayed_Sites){

  
	if( chA== NULL || Sequence==NULL || nca<0 || dec<0){
    printf("ERROR in updateSequence\n");
    exit(1);
		return -1;
	}

  int IDCharge = (int)getE(Decayed_Sites,dec,1);
  //printf("Updating Sequence\n");
  //fflush(stdout);
	int ID;
  Charge chi;
  Charge chj;
  int Indx_charge_in_Sequence;
  // Step one find the index of the charge
  int i;
  for(i=1;i<=nca;i++){
    ID = (int)getE(*Sequence,i,1); 
    if(ID==IDCharge){
      Indx_charge_in_Sequence = i;
      break;
    }
  }
  if(i>nca){
    printf("We have a problem i is large than the number of active charges and we have not found the chargeID\n");
    exit(1);
  }

  // Id of the charge
  //printf("Gettting Charge index %d active charges %d \n",Indx_charge_in_Sequence, nca);
  //fflush(stdout);

	int high = nca+1;
  //printf("Gettting ID %d ChargeID %d\n",ID,IDCharge);
  //fflush(stdout);
  chi = getCharge(chA,ID);
  //printf("nca %d\n",nca);
  for(int i=1;i<=nca;i++){
  //  printf("i %d\n",i);
 // fflush(stdout);
    if(i!=Indx_charge_in_Sequence){
      //printf("Grabbed Charge id %g\n",getE(*Sequence,i,1));
      //Make sure the charge is not one of the ones that has decayed
      bool ignore = false;    
      for(int k=(dec+1);k<=getRows(Decayed_Sites);k++){
        if(getE(*Sequence,i,1)==getE(Decayed_Sites,k,1)){
          ignore = true;
          break;
        }
      }
      chj = getCharge(chA,(int)getE(*Sequence,i,1));
      // Find the first index that is greater than i and then exit
      if(getDwel(chi)<getDwel(chj) && !ignore){
        //printf("high %d\n",high);
        high = i;
        break;
      }
    }
  }

    //printf("high %d Indx %d \n",high,Indx_charge_in_Sequence);
  if(high<Indx_charge_in_Sequence){
    for(int j=Indx_charge_in_Sequence;j>high;j--){
      setE(*Sequence,j,1,getE(*Sequence,j-1,1));
    }
    setE(*Sequence,high,1,(double)IDCharge);
  }else if(high>Indx_charge_in_Sequence){
    //if(high<nca){
    //printf("high+1 %d Value of charge site %g Charge dwell %g\n",high+1,getE(*Sequence,high+1,1),getDwel(getCharge(chA,(int)getE(*Sequence,high+1,1))));
    //}
    //if(high<=nca){
    //printf("high %d Value of charge site %g Charge dwell %g\n",high,getE(*Sequence,high,1),getDwel(getCharge(chA,(int)getE(*Sequence,high,1))));
    //}
    //printf("high-1 %d Value of charge site %g Charge dwell %g\n",high-1,getE(*Sequence,high-1,1),getDwel(getCharge(chA,(int)getE(*Sequence,high-1,1))));
    //printf("Indx %d value of Charge site %g Charge dwell %g\n",Indx_charge_in_Sequence,getE(*Sequence,Indx_charge_in_Sequence,1),getDwel(getCharge(chA,(int)getE(*Sequence,Indx_charge_in_Sequence,1))));

    for(int j=Indx_charge_in_Sequence;j<(high-1);j++){
      setE(*Sequence,j,1,getE(*Sequence,j+1,1));
    }
    //printf("Setting high-1 %d equal to IDCharge %d\n",high-1,IDCharge);
    setE(*Sequence,(high-1),1,(double)IDCharge);
  }
  //printf("Finished updating\n");
  //fflush(stdout);
	return 0;
}

// Simple check to ensure the sequence is correctly ordered
int checkSequence(const int nca, ChargeArray chA, matrix * Sequence){
  
  Charge chi;
  Charge chj;
  for(int i=1;i<nca;i++){
    chi = getCharge(chA,(int)getE(*Sequence,i,1));
    chj = getCharge(chA,(int)getE(*Sequence,i+1,1));
  
    if(getDwel(chi)>getDwel(chj)){
      printf("ERROR Sequence is now messed up. PrintingSequence\n");
      printf("ERROR Occurred with Charges %g and %g\n",getE(*Sequence,i,1),getE(*Sequence,i+1,1));
      printf("      and times of charges  %g and %g\n",getDwel(chi),getDwel(chj));
      for(int j=1;j<=nca;j++){
        printf("%d Charge ID %g Time %g\n",j,getE(*Sequence,j,1),getDwel(getCharge(chA,(int)getE(*Sequence,j,1))));
      }
      fflush(stdout);
      exit(1);
      return -1;
    }
  }
  return 0;
}

/*// I  believe this updates the sequence
int updateSequence(const int nca, ChargeArray chA, matrix * Sequence, int Indx_charge_in_Sequence){

	if( chA== NULL || Sequence==NULL || nca<0 || Indx_charge_in_Sequence<0){
    printf("ERROR in updateSequence\n");
    exit(1);
		return -1;
	}

  printf("Updating Sequence\n");
  fflush(stdout);
	int i;
	int low;
	int high;
	int middle;
	int ID;
  int newID;
  Charge chi;
  Charge chj;

  // Id of the charge
  printf("Gettting Charge index %d active charges %d \n",Indx_charge_in_Sequence, nca);
  fflush(stdout);
	ID = (int) getE(*Sequence,Indx_charge_in_Sequence,1);

	low = 1;
	high = nca;
  printf("Gettting ID %d\n",ID);
  fflush(stdout);
  chi = getCharge(chA,ID);
  printf("nca %d\n",nca);
  for(int i=1;i<=nca;i++){
    printf("i %d\n",i);
  fflush(stdout);
    if(i!=Indx_charge_in_Sequence){
      printf("Grabbed Charge id %g\n",getE(*Sequence,i,1));
      chj = getCharge(chA,(int)getE(*Sequence,i,1));
      if(getDwel(chi)<getDwel(chj)){
        printf("high %d\n",high);
        high = i;
        break;
      }
    }
  }
  printf("high %d Indx %d\n",high,Indx_charge_in_Sequence);
  if(high<Indx_charge_in_Sequence){
    for(int j=Indx_charge_in_Sequence;j>high;j--){
      setE(*Sequence,j,1,getE(*Sequence,j-1,1));
    }
    setE(*Sequence,high,1,Indx_charge_in_Sequence);
  }else{
    for(int j=Indx_charge_in_Sequence;j<high;j++){
      setE(*Sequence,j,1,getE(*Sequence,j+1,1));
    }
    setE(*Sequence,high,1,Indx_charge_in_Sequence);
  }
  printf("Finished updating\n");
  fflush(stdout);
	return 0;
}
*/
// I  believe this resorts the sequence
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
	ID = (int) getE(*Sequence,1,1);

	low = 2;
	high = nca;

	while(low<=high){
		middle = (low+high)/2;
		//Chi is a probe into the sequence array
		chi = getCharge(chA,(int)getE(*Sequence,middle,1));
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
		newID = (int)getE( *Sequence,i+1,1);
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
  //printf("Calling Save Data Point\n");

	//Might need to resize the matrix if it becomes to large
	if((*CurrentInc)>getRows(*Xcurrent)){
		if (*CurrentInc >= 512 ){
      
      //printf("CurrentInc %d\n",*CurrentInc);
			//Matrices have gotten to a size where they should be printed
			//to a file and emptied to free memory
			printTransportData(*System, *timeArray, *Xcurrent, *Ycurrent, *Zcurrent,\
						*Xelec_Drain, *Yelec_Drain, *Zelec_Drain,\
						*Xelec_Source, *Yelec_Source, *Zelec_Source,\
						*Xvelocity, *Yvelocity, *Zvelocity,\
						XElecOn, YElecOn, ZElecOn, FileName,\
						ElectricFieldX, ElectricFieldY, ElectricFieldZ);


		
			deleteMatrix(System);
			deleteMatrix(timeArray);

			deleteMatrix(Xcurrent);
			deleteMatrix(Ycurrent);
			deleteMatrix(Zcurrent);
			if(XElecOn==1){
				deleteMatrix(Xelec_Drain);
				deleteMatrix(Xelec_Source);
			}
			if(YElecOn==1){
				deleteMatrix(Yelec_Drain);
				deleteMatrix(Yelec_Source);
			}
			if(ZElecOn==1){
				deleteMatrix(Zelec_Drain);
				deleteMatrix(Zelec_Source);
			}

			deleteMatrix(Xvelocity);
			deleteMatrix(Yvelocity);
			deleteMatrix(Zvelocity);
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

	//printf("Saving Data \n\n");

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
		if(isnan(vel)){
			printf("Vel is nan t %g TimeTrack1 %g\n",(double)t,(double)*TimeTrack1);
			exit(1);
		}
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
	//printf("System %d \t",nc);

	if(XElecOn==1){
		//printf("Drain %d Source %d \n",*Xdrain,*Xsource);
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
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
		//therefore the velocity will not have changed.
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(xx == SLength*EndX && EndX != 0 ){
		//Arrived at front electrode
		(*Xdrain)=*Xdrain+1;
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
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
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(yy == SWidth*EndY && EndY != 0 ){
		//Arrived at front electrode
		(*Ydrain)=*Ydrain+1;
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
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
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;

	}else if(zz == SHeight*EndZ && EndZ != 0 ){
		//Arrived at top electrode
		(*Zdrain)=*Zdrain+1;
		if( t!=*TimeTrack1){
			(*TotalVelX) = (*TotalVelX) + (TrackX2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelY) = (*TotalVelY) + (TrackY2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
			(*TotalVelZ) = (*TotalVelZ) + (TrackZ2)*(SiteDistance2)/((t-(*TimeTrack1))*(nc2));
		}
		(*TrackX) = 0;
		(*TrackY) = 0;
		(*TrackZ) = 0;
		(*NumAvgVel)=*NumAvgVel+1;
		(*TimeTrack1) = t;
		(*ElecExit)=1;
	}

	if(isnan(*TotalVelX)){
		printf("ERROR *TotalVelX is nan!\n");
		exit(1);
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
  
    //If it is defined that means there is an electrode
    if(checkCluster_elecXdefined((ClusterLL)getClusterList(site))){
      // rv of 0 means not located adjacent to cluster
      // rv of 1 means located next to back electrode
      // rv of 2 means located next to front electrode
      rv = getCluster_elecXid((ClusterLL) getClusterList(site));
      if(((double)x)/((double)SLength)==0.0){
        //Check to see if next to the back electrode
        if(rv==2 || rv==1){
          //Cluster Hop is not allowed
          (*CheckX)=1;
        }else{
          (*CheckX)=0;
        }
      }else if((double)Cx/((double)SLength)==(double)(EndX-1) || 
               (double)Cx/((double)SLength)==(double)EndX){
        //Check to see if next to the front electrode
        if(rv==2 || rv==1){
          //if The charge is on the leading edge we have
          //to check both electrodes to see if cluster
          //hop is allowed if either is true it is not
          //allowed
          (*CheckX)=1;
        }else{
          (*CheckX)=0;
        }
      }
    }
	}
  if(PeriodicY==1){
    if(checkCluster_elecYdefined((ClusterLL)getClusterList(site))){
      rv = getCluster_elecYid((ClusterLL) getClusterList(site));
      if((double)y/((double)SWidth)==0.0){
        //Check if next to left electrode
        if(rv==2 || rv==1){
          (*CheckY)=1;
        }else{
          (*CheckY)=0;
        }
      }else if((double)Cy/((double)SWidth)==(double)(EndY-1) || 
               (double)Cy/((double)SWidth)==(double)EndY){
        //Check if next to right electrode
        if(rv==2 || rv==1){
          (*CheckY)=1;
        }else{
          (*CheckY)=0;
        }
      }
    }
  }
  if(PeriodicZ==1){
    if(checkCluster_elecZdefined((ClusterLL)getClusterList(site))){
      rv = getCluster_elecZid((ClusterLL) getClusterList(site));

      if((double)z/((double)SHeight)==0.0){
        if(rv==2 || rv==1){
          (*CheckZ)=1;
        }else{
          (*CheckZ)=0;
        }
      }else if((double)Cz/((double)SHeight)==(double)(EndZ-1) || 
               (double)Cz/((double)SHeight)==(double)EndZ){
        
        if(rv==2 || rv==1){
          (*CheckZ)=1;
        }else{
          (*CheckZ)=0;
        }
      }
    }
  }

	return 0;
}

int HoppingToSurroundingSites(SiteNode site, int codeX, int codeY, int codeZ){

  #ifdef _ERROR_CHECKING_ON_
	if(site==NULL || codeX<0 || codeY<0 || codeZ<0 ||\
			codeX>2 || codeY>2 || codeZ>2){
    #ifdef _ERROR_
		printf("ERROR incorrect input parameters detected in HoppingToSurroundingSites\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif
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

  // This means we are next to one of the 
  // x electrodes
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
			deleteMatrix(&mtxPval);

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

			deleteMatrix(&mtxPval);

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

			deleteMatrix(&mtxPval);

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

			deleteMatrix(&mtxPval);

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
		
			deleteMatrix(&mtxPval);

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

			deleteMatrix(&mtxPval);

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

		deleteMatrix(&mtxPval);


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

		deleteMatrix(&mtxPval);

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

		deleteMatrix(&mtxPval);

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

		deleteMatrix(&mtxPval);

	}

	return l;
}

int printFileDecay(int x, int y, int z,const long double time){

	char buf[256];
	snprintf(buf, sizeof buf,"%s","SiteDecay.txt");

  FILE * DecayOut;
  if((DecayOut = fopen(buf,"a")) == NULL){
    printf("ERROR! unable to open SiteDecay.txt file\n");
    return -1;
  }else{
    fprintf(DecayOut,"%Lg %d %d %d\n",time,x,y,z);
    fclose(DecayOut);
  }
  return 0;
}

int printFileChargeEnergy(const_SNarray snA, ChargeArray chA, matrix Sequence,
          int nca, long double time, ParameterFrame PF){

  char buf[256];
  snprintf(buf,sizeof buf,"%s","AverageChargeEnergy");

  FILE * ChargeEnergy;
  if((ChargeEnergy = fopen(buf,"a"))==NULL){
    printf("ERROR! unable to open AverageChargeEnergy\n");
    return -1;
  }else{
    double Energy=0;
    // Cycle through all the charges and add all the energies up
    for(int r=1;r<=nca;r++){
      Charge ch = getCharge(chA,(int)getE(Sequence,r,1));
      int x, y, z;
      x = (getCx(ch)+PFget_Len(PF))%PFget_Len(PF);
      y = (getCy(ch)+PFget_Wid(PF))%PFget_Wid(PF);
      z = (getCz(ch)+PFget_Hei(PF))%PFget_Hei(PF);
      if(x<0){
        x = PFget_Len(PF)+x;
      }
      if(y<0){
        y = PFget_Wid(PF)+y;
      }
      if(z<0){
        z = PFget_Hei(PF)+z;
      }
      SiteNode sn = getSN(snA, x,y,z);
      Energy+=getEnergy(sn); 
    }
    Energy = Energy/((double)nca);
    fprintf(ChargeEnergy,"%Lg %g %d\n",time,Energy,nca);
    fclose(ChargeEnergy);
  }
  return 0;
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

int CheckConservationCharges( Electrode elXb,Electrode elYl,Electrode elZb,\
                             int nca        ,int XElecOn   ,int YElecOn   ,int ZElecOn   ,\
                             int nc){

  int ChargeCheck = nc;
  if(XElecOn==1){
    ChargeCheck = ChargeCheck+getElectrode_Charges(elXb);
  }
  if(YElecOn==1){
    ChargeCheck = ChargeCheck+getElectrode_Charges(elYl);
  }
  if(ZElecOn==1){
    ChargeCheck = ChargeCheck+getElectrode_Charges(elZb);
  }
  if( ChargeCheck != nca){
    #ifdef _ERROR_
    fprintf(stderr,"nc %d elXb %d ",nc,getElectrode_Charges(elXb));
    fprintf(stderr,"elYl %d elZb %d\n",getElectrode_Charges(elYl),getElectrode_Charges(elZb));
    fprintf(stderr,"Charges lost somehow ChargeCheck %d nca %d!\n",ChargeCheck,nca);
    #endif
    exit(1);
  }
  return 0;
}

/* This function is responsible for projecting
 * a position (xx,yy,zz) that may be outside the boundaries
 * of (SLength,SWidth,SHeight) so that it is within those
 * defined boundaries. The new position is described
 * with (x,y,z) pointers
 */
int ProjectChargePositionOntoSiteNodeReferenceFrame(int *x     , int *y    , int *z     ,\
                                                    int xx     , int yy    , int zz     ,\
                                                    int SLength, int SWidth, int SHeight){
  int factorX = 0;
  int factorY = 0;
  int factorZ = 0;

  //Finding position of charge on sites
  if(xx<0){
    factorX = -xx/SLength+1;
    *x = (xx+SLength*factorX)%SLength;
  }else{
    *x = (xx+SLength)%SLength;
  }
  if(yy<0){
    factorY = -yy/SWidth+1;
    *y = (yy+SWidth*factorY)%SWidth;
  }else{
    *y = (yy+SWidth)%SWidth;
  }
  if(zz<0){
    factorZ = -zz/SHeight+1;
    *z = (zz+SHeight*factorZ)%SHeight;
  }else{
    *z = (zz+SHeight)%SHeight;
  }

  return 0;
}
