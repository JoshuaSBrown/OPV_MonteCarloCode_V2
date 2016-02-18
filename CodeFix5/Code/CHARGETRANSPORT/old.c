//FET CONFIGURATION
//INITIALISATION OF SITES
//RANDOM WALK
//	INITIALISES INDICES+OUTPUT FILES
//	FOR EACH TIME STEP 
//		INJECTS CHARGES 
//		UPDATES POSITION INDEX 
//		SORTS CHARGES IN A WAITING LIST 
//		FINDS FIRST CHARGE AND ITS WAIT TIME
//		DETERMINES HOP DIRECTION FOR THE FIRST CHARGE (SEMI-RANDOM: HOP DIRECTION IS WEIGHTED BY HOP PROBABILITY FROM ENERGETICS)
//		MOVES THE CHARGE, UPDATE POSITION, CURRENT AND 'RUN DISTANCE' INDICES... 

#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
//#include <math.h>
//#include <assert.h>

//#include "../parameter.h"
#include "../CHARGE/charge.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
//#include "../MATRIX/matrix.h"
//#include "../MONTECARLO/montecarlo.h"
//#include "../FUNCTIONS/functions.h"
//#include "../MEM/mem.h"

/*
void randomWalk(SNarray snA, double rN, int[][2*SLength+1],\
								double[TCount+1], double[TCount+1], double[TCount+1]);

//insert the jumping charge sequence[0] into sequence. the subscript is from 0 to high -1
void insertDwelltimePos(int high, ChargeArray chA, int * sequence);  
*/

/*
int ClusterHop(SNarray snA, Charge one, int * totalX);

char *outname = "t-I.txt";
char buffer[30], buffer2[30], buffer3[30];

int main(){

	mem_init();      
	atexit(mem_term);

	//real number of charges in the real materials, area = 0.11 (cm*cm)
 	double rN;
  double electricField;
	double electricEnergy;
	int r, l, j, n;
  int AV_Xrundist_counter[N_av+1][2*SLength+1];
	int Xrundist_counter[N_av+1][2*SLength+1];
  double Vindex[2*SLength+1];
	double currentR[TCount+1];
	int currentR_av[TCount+1];
	int rv;
	double ncR[TCount+1];
	double ncR_av[TCount+1];
	double nelec_drainR[TCount+1];
	double nelec_drainR_av[TCount+1];

  electricField = voltage / (SLength * SiteDistance);
	//if SiteDistance = 1E-7 and electricField = 6E5; then electricEnergy = 0.06 
	electricEnergy = SiteDistance * electricField; 

  rN = 1.272E10 + 2.0789E10 * pow(voltage, 0.77388);
  printf("rN is %g\n", rN);

  for (j = 0 ; j < 2 * SLength + 1 ; j++) {
		//initialize velocity array
    Vindex[j]= 0;   
		//assign velocity values
		Vindex[j]= ((j-SLength) * SiteDistance) / (Nstep_av*TStep); 
    
		for (l = 1; l < N_av+1; l++) {
      AV_Xrundist_counter[l][j]=0;
    }
  }
  
  for(n=0;n<TCount;n++){
    currentR[n]=0;
    ncR[n]=0;
    nelec_drainR[n]=0;
  }

  for (r=1;r<Rcount+1;r++){

    printf("Seed iteration # = %d\n", r);
		//set seed of random function
    srand(rdtsc());
		//set seed of random function
    srandom(rdtsc()); 

		SNarray snA = newSNarray(SLength,SWidth,SHeight);
    initSite(electricEnergy, snA);
		//Determine Clusters
		ArbArray ClArLL;
		ArbArray NeighClArLL;
		
		//ArbArrayCheck(ClArLL);
		FindCluster( snA, electricEnergy, &ClArLL, &NeighClArLL);

		rv=ArbArrayCheck(ClArLL);
		assert(rv!=-1);
		rv=ArbArrayCheck(NeighClArLL);
		assert(rv!=-1);

    printf("Past initSite\n");
    randomWalk( snA, rN, Xrundist_counter, currentR, ncR, nelec_drainR);
    printf("Past randomWalk\n");
		rv=ArbArrayCheck(ClArLL);
		assert(rv!=-1);
		rv=ArbArrayCheck(NeighClArLL);
		assert(rv!=-1);


    for (l = 1; l < N_av+1; l++) {
      for ( j = 0 ; j < 2 * SLength + 1 ; j++){
        AV_Xrundist_counter[l][j] += Xrundist_counter[l][j];
      }
    }

    for(n=0;n<TCount;n++){
      currentR_av[n]+=currentR[n];
      ncR_av[n]+=ncR[n];
      nelec_drainR_av[n]+=nelec_drainR[n];
    }

    if(r==Rcount){
      
			printf("Rcount is reached, computing average\n");
      for(n=0;n<TCount;n++){
        currentR_av[n]=currentR_av[n]/Rcount;
        ncR_av[n]=ncR_av[n]/Rcount;
        nelec_drainR_av[n]=nelec_drainR_av[n]/Rcount;
      }
      FILE *AV_tI;
      int coincoin = sprintf(buffer3, "AV_tI.txt");
      char *outname4 = buffer3;

      if((AV_tI = fopen(outname4, "w")) == NULL)
        printf("Could not open output file outname4.\n");   

      fprintf(AV_tI, "Length of sample in cm: %g \n",SLength*SiteDistance);
			fprintf(AV_tI, "Time interval in s: %g \n",TStep);
			fprintf(AV_tI, "Number of charges: %d \n",NCh);
			fprintf(AV_tI, "Applied voltage in V : %g \n", voltage);

      for(n=0;n<TCount;n++){
        fprintf(AV_tI, "%g \t %d \t %g \t %g \n", \
						n*TStep, currentR_av[n], ncR_av[n], nelec_drainR_av[n]);
      }
      fclose(AV_tI);


      for (l = 1; l < N_av+1; l++) {
        FILE *AV_velocity_file;
        int coincoin = sprintf(buffer2, "velocity_dist_AV_%d.txt", l);
        char *outname3 = buffer2;
        if((AV_velocity_file = fopen(outname3, "w")) == NULL)
          printf("Could not open output file outname3.\n");   

        fprintf(AV_velocity_file, "Length of sample in cm: %g \n",SLength*SiteDistance);
				fprintf(AV_velocity_file,"Time interval in s: %g \n",TStep);
				fprintf(AV_velocity_file,"Number of charges: %d \n",NCh);
				fprintf(AV_velocity_file,"Applied voltage in V : %g \n",voltage);
				fprintf(AV_velocity_file,"Time over which velocity is averaged in s: %g \n",Nstep_av*TStep);
				fprintf(AV_velocity_file,"Total number of averaging steps: %d \n",N_av);
				fprintf(AV_velocity_file,"Averaging step #: %d/%d \n",l, N_av);
        fprintf(AV_velocity_file, "l%d\tDistribution counter_av_%d \t",l,l);
				fprintf(AV_velocity_file,"Velocity in cm/s_av%d \t",l);
				fprintf(AV_velocity_file,"Number of counts_av%d \n", l); 

        for ( j = 0 ; j < 2 * SLength + 1 ; j++){
          fprintf(AV_velocity_file, "%d \t %d \t %g \t %d \n", l, j, Vindex[j], AV_Xrundist_counter[l][j]);
        }
        fclose(AV_velocity_file);
      }
    }
		printf("Here2\n");
		assert(ClArLL!=NULL);
	//Delete ClArLL and Neigh ClArLL
		printf("Here2.5\n");
		rv=deleteArbArray(ClArLL);
		assert(rv!=-1);
		printf("Here2.6\n");
		rv=deleteArbArray(NeighClArLL);
		assert(rv!=-1);
		printf("Here3\n");
		deleteSNarray(snA);
		printf("Here3.5\n");

  }

  int t;
  for ( t=0;t< SLength;t++) {
    printf("t is %d\t",t);
  }

			printf("Here4\n");

  return 0;
}
*/
int SiteHop(SNarray snA, Charge one, SiteNode site, int * totalX, int * totalY, int * totalZ, const int EndX, const int EndY, const int EndZ, const int XElecOn, const int YElecOn, const int ZElecOn, const int PeriodicX, const int PeriodicY, const int PeriodicZ){

	//totalX is used to measure the total current in the X direction as the charges move around
	//totalY is used to measure the total current in the Y direction as the charges move around
	//totalZ is used to measure the total current in the Z direction as the charges move around
	if((PeriodicX==1 && EndX==-1)	|| (PeriodicX==0 && EndX>0) ||\
     (PeriodicY==1 && EndY==-1)	|| (PeriodicY==0 && EndY>0) ||\
		 (PeriodicZ==1 && EndZ==-1)	|| (PeriodicZ==0 && EndZ>0)){
		return -1;
	}

	if((EndX==0 && XElecOn==1) ||\
		 (EndY==0 && YElecOn==1) ||\
		 (EndZ==0 && ZElecOn==1)){
		return -1;
	}

	if(snA==NULL || one == NULL || site ==NULL){
		return -1;
	}

	//Initially, the flag is set to 0 
	//means the charge is able to hop
	int flag = 0; 
  double position;
  int middle;
	int rv;
	int l;
	//Below values define location of the charge in the
	//sample when the periodicity is removed
	int x;
	int y;
	int z;
	int SLength;
	int SWidth;
	int SHeight;

	SLength = getAlen(snA);
	SWidth = getAwid(snA);
	SHeight = getAhei(snA);
	//First Check is to see if all the neighboring sites are occupied. 
	//If the site is not located at a boundary
	//If they are all occupied set the flag to 1
	x = getCx(one)%SLength;
	y = getCy(one)%SWidth;
	z = getCz(one)%SHeight;
	
	//This is the random number used to determine which direction a charge hops
	position = (double)rand() / RAND_MAX;
	
	if(x!=0 && x!=(SLength-1) && y!=0 && y!=(SWidth-1) && z!=0 && z!=(SHeight-1)){
		rv = OccAllNei(snA, x, y, z);
		if (rv==1){
			//All neighboring sites are occupied
			flag=1;
		}

		if(flag==0){

			l=0;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getSN_p(site,l)<position){
				l++;
			}
			middle = l;
		


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
					codeX==0;
				}else{
					if(x==0){
					//Need to determine if charge is at the edge of the boundaries.
						if(getCx(one)/SLength==0){
						//Means it is at the boundary
							codeX=1;
						}else{
						//Means not at the boundary
							codeX=0;
						}
					}else if(x==(SLength-1)){
						
						if(getCx(one)/SLength==EndX-1){
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
					codeX=1;
				}else if(x==SLength-1){
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
	
				if(YElecOn==1){
					codeY==0;
				}else{
					if(y==0){
						if(getCy(one)/SWidth==0){
							codeY=1;
						}else{
							codeY=0;
						}
					}else if(y==(SWidth-1)){
						if(getCy(one)/SWidth==EndY-1){
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
					codeY=1;
				}else if(y==SWidth-1){
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
				codeZ==0;

			}else{
			//Charge might hit boundary depending on how many samples
			//charge has traveled through
				if(ZElecOn==1){
					codeZ==0;
				}else{
					if(z==0){
						if(getCz(one)/SHeight==0){
							codeZ=1;
						}else{
							codeZ=0;
						}
					}else if(z==(SHeight-1)){
						if(getCz(one)/SWidth==EndZ-1){
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

		if(codeX==0 && codeY==0 && codeZ==0){
			position = (double)rand() / RAND_MAX;
			l=0;
			//Cycle through the pvals until the pval
			//is equal or above the random number
			while( getSN_p(site,l)<position){
				l++;
			}
			middle = l;
		}else if(codeX>0 && codeY==0 && codeZ==0){
			//Need to redefine the pvals
			matrix mtxPval = newMatrix(5,1);
			if(codeX==1){
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
				middle = l;
				deleteMatrix(mtxPval);

			}else{
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
				
				middle = l;
				deleteMatrix(mtxPval);

			}
		}else if(codeX==0 && codeY>0 && codeZ==0){
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

				middle = l;
				deleteMatrix(mtxPval);

			}else{
			//Not allowed to hop to the right or +y
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

				middle = l;
				deleteMatrix(mtxPval);

			}
		}else if(codeX==0 && codeY==0 && codeZ>0){
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

				middle = l;
				deleteMatrix(mtxPval);

			}else{
			//Cannot hop Above
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

				middle = l;
				deleteMatrix(mtxPval);

			}
		}else if(codeX>0 && codeY>0 && codeZ==0){
			matrix mtxPval = newMatrix(4,1);
			if(codeX==1){
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
				if(codeX==1){
					l=0;
				}else{
					l=1;
				}
			}

			middle = l;
			deleteMatrix(mtxPval);


		}else if(codeX==0 && codeY>0 && codeZ>0){
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
			
			middle = l;
			deleteMatrix(mtxPval);
			
		}else if(codeX>0 && codeY==0 && codeZ>0){
			matrix mtxPval = newMatrix(4,1);
			
			if(codeX==1){
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
				if(codeX==1){
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
				
			middle = l;
			deleteMatrix(mtxPval);
		
		}else{
			matrix mtxPval = newMatrix(3,1);
			
			if(codeX==1){
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
				if(codeX==1){
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

			middle = l;
			deleteMatrix(mtxPval);
		}

	}

	//l value
	//0 - backward 	x-
	//1 - forward 	x+
	//2 - left 			y-
	//3 - right 		y+
	//4 - below 		z-
	//5 - above	 		z+

	if(flag!=1){		
		switch(middle){
		//if the charge is on the illuminated electrode (x=0), it recombines there,
		//ie negative current
		case 0: if( OccXneg(snA, x, y, z)==1) {
							//there is no possibility to hopp
							flag = 1;   
							break;
						}
						//Move the charge backward
						CxMinus(one);
						//current decreases in x
						*totalX=*totalX-1;  
						break;

						//if the charge is adjacent to the opposite electrode (x=Length), 
						//it recombines there, ie positive current
		case 1: if( OccXpos(snA, x,y ,z)== 1){  
							//there is no possibility to hop
							flag = 1;   
							break;
						}
						//move the charge forward
						CxPlus(one);   
						//current increases in x
						*totalX=*totalX+1;   
						break;

		//if the site on the left of the charge position is occupied,
		case 2: if( OccYnegPeriodic(snA, x,y,z)==1){ 
							//there is no possibility to hopp
							flag = 1;
							break;
						}
						//move the charge to the left
						CyMinus(one);
						//current decreases in y
						*totalY=*totalY-1;
						break;

		//if the site on the right of the charge position is occupied,
		case 3: if( OccYposPeriodic(snA, x,y,z)==1){   
							//there is no possibility to hopp
							flag = 1;   
							break;
						}
						//move the charge to the right
						CyPlus(one);
						//current increases in y
						*totalY=*totalY+1;
						break;

		//if the charge is at the lower interface or the site in the decided direction 
		//(below the charge position) is occupied,
		case 4: if( OccZneg(snA, x,y,z)==1){  
							//then there is no possibility to hop
							flag = 1;   
							break;
						}
						//move the charge down
						CzMinus(one);
						//current decreases in z
						*totalZ=*totalZ-1;
						break;
		
		//if the charge is at the upper interface or the site in the decided direction 
		//(above the charge position) is occupied,
		case 5: if( OccZpos(snA, x,y,z)==1){
							//there is no possibility to hopp
							flag = 1;
							break;
						}
						//move the charge up
						CzPlus(one);
						//current increases in z
						*totalZ=*totalZ+1;
						break;
		}
	}
	return flag;

}

int ClusterHop(SNarray snA, Charge one, int * totalX){
	
	printf("Hopping to cluster\n");

	//flag is 0 if charge moves and 1 if charge stays put
	int flag = 0;
	//initially, the flag is set as 0 - can hop, 1 is can not hop
	//flag is used to determine if charges can hop to sites outside
	//of the cluster. 
	int flag1 = 0;
	//flag2 accounts for availability within the cluster. 
	int flag2 = 0;
	int rv;
	int ID, NewID;
	long double sum=0.0;
	double pvalHigh;
	double pval;
	int i, j, k;
	int i1, j1, k1;
	double position2;
  double position;
	double time;
	long double sum2=0.0;
 	SiteNode sn;

	ID=getIndex(snA, getCx(one), getCy(one), getCz(one));
	NewID=ID;
	//Check to see if there are sites adjacent to the cluster that are 
	//open
	rv = OccAllNeiCluster(snA, getCx(one), getCy(one), getCz(one));
	if (rv==1){
		//there is no possibility for hopping outside the cluster as all the neighbouring sites 
		//around the first charge position are all occuppied and the cluster
		//does not sit next to an electrode
		printf("Hopping outside Cluster Not permitted\n");
		flag1 = 1;
	}

	rv = OccAllCluster( snA, getCx(one), getCy(one), getCz(one));
	if(rv==1) {
		//This means hopping within the cluster is not permitted
		printf("Hopping Within Cluster Not permitted\n");
		flag2 = 1;
		//But hopping outside may be allowed;
		if(flag1!=1){

			getNeighClusterPvalHigh( snA, getCx(one), getCy(one), getCz(one), &pvalHigh, &sum);
			//pval=pvalHigh*2;
			sn = getSN(snA, getCx(one), getCy(one), getCz(one));
			setp(sn, 0, 0.0, 4);
			//Don't add to sum because hopping within cluster is not allowed
			//Can't stay where it is. hopping within the cluster is not allowed
			//so pval for the cluster is 0.0
		}
			//what happens if can't hop on or off?

	}else{
	
		double pvalLow;
		//Hopping within the cluster is allowed must account for this in the sum
		getClusterPvalLow( snA, getCx(one), getCy(one), getCz(one), &pvalLow, &sum2);
		assert(pvalLow>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
		
		if (flag1!=1){
			//Grabbing the fastest hop rate off the cluster
			getNeighClusterPvalHigh( snA, getCx(one), getCy(one), getCz(one), &pvalHigh, &sum);
			assert(pvalHigh>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
			
			sn = getSN(snA, getCx(one), getCy(one), getCz(one));
			time=getTimeCluster( (ClusterLL) getList(sn));
			assert(time<1/pvalHigh);

			//if(pvalLow<pvalHigh*2){
			//	pval=pvalLow;
				setp(sn, 0, , 4);
			//}else{
			//	pval=pvalHigh*2;
			//	setp(sn, 0, pval, 4);
			//}
		}else {
			pval=pvalLow;
			sn = getSN(snA, getCx(one), getCy(one), getCz(one));
			setp(sn, 0, pval, 4);
		}
		printf("Value of sum Before %Lg\n",sum);
		sum+=(long double) pval;
		printf("Value of sum After %Lg\n",sum);
		printf("Value of pval %g\n",pval);
		printf("pvalLow %g pvalHigh %g\n",pvalLow, pvalHigh);
		assert(pval>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
	}

	if(flag1==0 && flag2==0){
	//There is a possibility for hopping to at least one of the neighbouring 
	//site off the cluster. There is also the possibility that a site within
	//the cluster is occupied. This loop determines in which direction the charge 
	//will jump to by randomly choosing a number between 0 and 1 
	//(random/max_random)
	//This number is compared with the pvalues within the cluster and those 
	//the pvalues of the sites neighboring the cluster

	//Determine if site will hop on or off cluster
		position = (double)rand() / RAND_MAX;

		//rv will be:
		//1 if hop to new neighbor outside of cluster
		//0 if stay where it is because hop site is occupied
		//2 if hopped to Right electrode
		//3 if hopped to Left electrode
		//-1 if Neither of these ~ error
		printf("\nposition %g \n",position);
		rv = HopOnOffCluster( snA, ID, position, sum, &NewID);
		assert(rv!=-1);
	
		if (rv==0){
			//This means could not hop off cluster but might
			//be able to move around in the cluster
			assert(NewID==ID);
			position2 = (double)rand() / RAND_MAX;
			rv = HopWithinCluster(snA, ID, position2, sum2, &NewID);
			assert(rv!=-1);
			if(rv==1){
				//Moved to new site within the cluster
				//Find the distance in positions between the site hopped
				//to and the site hopped from
				i=getCx(one);
				j=getCy(one);
				k=getCz(one);
				getLoc(&i1,&j1, &k1, NewID, snA);
				*totalX=*totalX+(i1-i);
				setXdist(one, getXdist(one)+(i1-i));
				setCx(one, i1);
				setCy(one, j1);
				setCz(one, k1);

			}else{
				//Charge did not move
				flag=1;
			}
		}else if(rv==1){
			assert(NewID!=ID);

			//Find the distance in positions between the site hopped
			//to and the site hopped from
			i=getCx(one);
			j=getCy(one);
			k=getCz(one);

			getLoc(&i1,&j1, &k1, NewID, snA);
			*totalX=*totalX+(i1-i);
			setXdist(one, getXdist(one)+(i1-i));
			setCx(one, i1);
			setCy(one, j1);
			setCz(one, k1);
		}else if(rv==2){
			//hopped to right electrode x=len
			i=getCx(one);
			*totalX=*totalX+(getAlen(snA)-i);
			setXdist(one, getXdist(one)+(getAlen(snA)-i));
			setCx(one,getAlen(snA));

		}else if(rv==3){
			//hopped to left electrode x=-1
			i=getCx(one);
			*totalX=*totalX-i;
			setXdist(one, getXdist(one)-i);
			setCx(one,0);
		}
	}else if(flag1==0 && flag2==1){
		//Hopping within the cluster is not allowed
		position = (double)rand() / RAND_MAX;

		//rv will be:
		//1 if hop to new neighbor outside of cluster
		//0 if stay where it is because hop site is occupied
		//2 if hopped to Right electrode
		//3 if hopped to Left electrode
		//-1 if Neither of these ~ error
		rv = HopOnOffCluster( snA, ID, position, sum, &NewID);
		assert(rv!=-1);
		if(rv==0){
			flag=1;
		}else if(rv==1){
			assert(NewID!=ID);

			//Find the distance in positions between the site hopped
			//to and the site hopped from
			i=getCx(one);
			j=getCy(one);
			k=getCz(one);

			getLoc(&i1,&j1, &k1, NewID, snA);
			*totalX=*totalX+(i1-i);
			setXdist(one, getXdist(one)+(i1-i));
			setCx(one, i1);
			setCy(one, j1);
			setCz(one, k1);
		}else if(rv==2){
			//hopped to right electrode x=len
			i=getCx(one);
			*totalX=*totalX+(getAlen(snA)-i);
			setXdist(one, getXdist(one)+(getAlen(snA)-i));
			setCx(one,getAlen(snA));

		}else if(rv==3){
			//hopped to left electrode x=-1
			i=getCx(one);
			*totalX=*totalX-i;
			setXdist(one, getXdist(one)-i);
			setCx(one,0);

		}

	}else if(flag1==1 && flag2==0){
		//Hopping off the cluster is not allowed
		assert(NewID==ID);
		position2 = (double)rand() / RAND_MAX;
		printf("position2 %g\n",position2);
		rv = HopWithinCluster(snA, ID, position2, sum2, &NewID);
		assert(rv!=-1);
		if(rv==1){
			//Moved to new site within the cluster
			//Find the distance in positions between the site hopped
			//to and the site hopped from
			i=getCx(one);
			j=getCy(one);
			k=getCz(one);
			getLoc(&i1,&j1, &k1, NewID, snA);
			*totalX=*totalX+(i1-i);
			setXdist(one, getXdist(one)+(i1-i));
			setCx(one, i1);
			setCy(one, j1);
			setCz(one, k1);

		}else{
			//Charge did not move
			flag=1;
		}

	}else{
		//No where for charge to hop within or out of the Cluster
		flag=1;
	}
	return flag;
}

void  randomWalk( SNarray snA,double rN, int Xrundist_counter[][2*SLength+1], double currentR[TCount+1],\
									double ncR[TCount+1], double nelec_drainR[TCount+1]){	
  
	printf("Starting Random Walk\n");
	//if flag = 1, it means that the node that the charge 
	//will jump to is occupied; if flag = 0, it means that 
	//that node is empty
	int flag;  
  int i, j, l;
	//number of time step
	int n;
	//nc total number of charge remaining in the material
	//nelec_drain total number of charges to reach drain
  int nc, nelec, nelec_drain; 
	//the total distance moved in the x direction by all the 
	//charge
  int totalX; 
  int X;
	int sequence[Ntot];
	//time (as a clock starting at the beginning of the steps)
	long double t;  
	//receive the random number
  double ran;    
  double current;
  FILE * fp;
	Charge one;
	Charge chtemp;
	SiteNode site;
	SiteNode site2;
  //initialization of time and charge counters	
	//the value of n is from 1 to TCount
	n = 1;
  X = 0 ;
	//initial value of remaining charges is N
  nc = 0;
  nelec=0;
  nelec_drain=0;
	//initial total distance covered by charges in the x direction is 0
  totalX = 0; 
	//initial time is 0
  t = 0.0; 

	//initialize run dist counter
  for (j = 0 ; j < 2 * SLength + 1 ; j++) {
    for (l = 0; l < N_av+1 ; l++) {
      Xrundist_counter[l][j] = 0;   
    }
  }
  
	for(n=1;n<TCount+1;n++){
    currentR[n]=0;
    ncR[n]=0;
    nelec_drainR[n]=0;
  }
  n = 1;

  if((fp = fopen(outname, "w")) == NULL){
    printf("Could not open output file outname.\n");   
  }

  fprintf(fp, "Temps \t current \t nc \n");

  for(i = 0; i < Ntot; i++){
    sequence[i] = i;
  }
  
	printf("Before initCharget() in random walk \n");
  ChargeArray chA = initCharget0( sequence, snA);
  nc = NCh;
  printf("Past initCharget0() in random walk \n");

	//for each time step as long as there are remaining charges in the material
  for(; n <= TCount && nc > 0;){   

		//address of the first charge in the waiting queue
    one = getCharge(chA, sequence[0]);  
		//address of the node where the first charge in the waiting queue is located
		site = getSN(snA, getCx(one), getCy(one), getCz(one));
		printf("\nnv %d \t",nc);
		//printf("Charge number %d \t",sequence[0]);
		//printf("Site number %d \t",getIndex(snA,getCx(one),getCy(one),getCz(one)));
		//printf("Time %0.17Lg\t",t);

		//if time of the first charge in the queue is equal to global time it means
		//the charge has not been properly initialized, there is probably a problem
		//in initCharget0 or initCharge
    if(getDwel(one) <= 0){
      //printf("Charge number %d\n", (sequence[0]));
      //printf("dwelltime %g\n", getDwel(one));
      //printf("1 charge# %d \t dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n",\
							sequence[0], getDwel(one), getCx(one), getCy(one), getCz(one));
      exit(1);  
    }

		//Check if site is on the cluster or off
		if(getType(site)==0){
			//Charge is not located in a cluster
			flag=SiteHop(snA, one, site, &totalX);
			//printf("NotC\t");
		}else{
			assert(getType(site)==1);
			//Charge is located in a cluster
		  flag = ClusterHop( snA,  one,  &totalX);
			//printf("OnC \t");
		}

		//printf("\tcx %d cy %d cz %d \t",getCx(one), getCy(one), getCz(one));
		//printf("Flag %d",flag);
		//hopping procedure: determine hopping direction and hop if destination site is free,
		//the position of one is updated, not its other properties like dwell time and charge 
		//time

		//if flag is still 0 after the hopping procedure, it means the first charge in the 
		//queue has hopped to an adjacent site
		
		if(flag == 0){
			//then set site as unoccupied ie -1
			setDwelStat(site,-1);

			//if time of the first charge in the queue is equal to global time
      if(t < 0){
        printf("t is negative beginning of loop %Lg \t exit(1)\n", t);
        exit(1);  
      }

			//if time of the first charge in the queue is equal to global time
      if(getDwel(one) <= 0){
        printf("2 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n",\
						    getDwel(one), getCx(one), getCy(one),getCz(one));
        exit(1);  
      }

      t += (long double) getDwel(one);  

			//if time of the first charge in the queue is equal to global time
      if(t < 0){
        printf("t is negative middle of loop %Lg \t exit(1)\n", t);
        exit(1);  
      }

			// if time of the first charge in the queue has reached the current time step, this 
			// loop calculates the current
      if(t >= (long double)(n * TStep)){
        for (l = 1; l < N_av+1; l++) {
					//if averaging step is reached,				
          if (n == l * Nstep_av){

						//for all remaining charges
            for ( i = 0; i < nc ; i++ ) {
              Xrundist_counter[l][getXdist(getCharge(chA,i)) + SLength]++;
							//reset counter to 0 at the end of the averaging step
							setXdist(getCharge(chA,i),0);
            }
					}
        }
				//if the averaging step is reached, calculate the velocity distribution

        current = ((q * totalX) / ((TStep) * SLength)) * (rN/NCh);
				currentR[n]=current;
        ncR[n]=nc;
        nelec_drainR[n]=nelec_drain;

        fprintf(fp, "%Lg\t%g\t%d\t%d\n", t, current, nc, nelec_drain);
				//reset the distance run by all charges to zero for next time step calculation
        totalX = 0;
        nelec_drain=0;
      }
			//if time of the first charge in the queue has reached the current time step, this 
			//loop calculates the current and velocity distribution and prints them in output files

			//if the charge has arrived to an electrode,  
			if(getCx(one) == -1 || getCx(one) == SLength){
        //printf("1st charge arrived to an electrode.\n");
				// nb of remaining charges decreases
				nc--;
        printf("Number charges left %d\n",nc);
				// nb of collected charges increases
        nelec++;
        if(getCx(one)==SLength){
          nelec_drain++;
        }

        setDwel(one, 1E6);

        for(i = 0; i < Ntot-1; i++){
					//the waiting queue of charges is updated
          sequence[i] = sequence[i+1];
        }
        sequence[Ntot-1] = Ntot-1;

				//if time of the first charge in the queue is equal to global time
        if(getDwel(one) <= 0){
          printf("3 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n", \
									getDwel(one), getCx(one),getCy(one),getCz(one));
          exit(1);  //
        }


        //the moving charge (=first in the queue) has arrived to an electrode => update nb 
				//of charges and waiting list
			} else  { 
				site2=getSN(snA, getCx(one),getCy(one),getCz(one));
				//if the first charge in the queue has hopped to an adjacent site
        if(flag == 0) {
					//the dwell status of the site onto which the charge has hopped is set to occupied
          setDwelStat(site2,sequence[0]);
        }
        
				incVisFreq(site2);
        setVis(site2,1);
        for(i = 1; i < nc; i++) {
					//the waiting times of following charges in the queue are decreased by the (passed) 
					//waiting time of the first charge in the queue
          //Should not minus dwelltime for charges that have already crossed to the other electrode.
					chtemp = getCharge(chA, sequence[i]);
					MinusDwel(chtemp, getDwel(one));
        }

        do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
				//the waiting time of ONE is randomly updated
				setDwel(one, -6.0 * log(ran/RAND_MAX) / (D * getsum(site)));
				//printf("\t NewDwel %g",getDwel(one));
				//the time of the first charge in the queue is increased by the (passed) waiting time 
				//of the first charge in the queue
        Plust(one,(long double) getDwel(one)); 

				//if time of the first charge in the queue is equal to global time
        if(getDwel(one) <= 0){
          printf("4 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n", getDwel(one), getCx(one), getCy(one), getCz(one));
          exit(1);  
        }

				// ONE is inserted back into the queue at its new position
        insertDwelltimePos(nc, chA, sequence);
				//if time of the first charge in the queue is equal to global time
        if(t < 0){
          printf("t is negative end of loop %Lg \t exit(1)\n", t);
          exit(1); 
        } 
        //charge has moved inside the material; waiting times are updated and the queue is sorted
        
			}

		}else{

			///////////////////////////////////////////////////////////////////////////////////////////
			//In the case flag == 1 means charge did not move
			//if time of the first charge in the queue is equal to global time
      if(t < 0){
        printf("t is negative beginning of loop %Lg \t exit(1)\n", t);
        exit(1);  
      }

			//if time of the first charge in the queue is equal to global time
      if(getDwel(one) <= 0){
        printf("2 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n",\
						    getDwel(one), getCx(one), getCy(one),getCz(one));
        exit(1);  
      }
      t += (long double) getDwel(one);  

			//if time of the first charge in the queue is equal to global time
      if(t < 0){
        printf("t is negative middle of loop %Lg \t exit(1)\n", t);
        exit(1);  
      }
			// if time of the first charge in the queue has reached the current time step, this 
			// loop calculates the current
      if(t >= (long double)(n * TStep)){
        for (l = 1; l < N_av+1; l++) {
					//if averaging step is reached,				
          if (n == l * Nstep_av){

						//for all remaining charges
            for ( i = 0; i < nc ; i++ ) {
              Xrundist_counter[l][getXdist(getCharge(chA,i)) + SLength]++;
							//reset counter to 0 at the end of the averaging step
							setXdist(getCharge(chA,i),0);
            }
					}
        }
				//if the averaging step is reached, calculate the velocity distribution

        current = ((q * totalX) / ((TStep) * SLength)) * (rN/NCh);
				currentR[n]=current;
        ncR[n]=nc;
        nelec_drainR[n]=nelec_drain;

        fprintf(fp, "%Lg\t%g\t%d\t%d\n", t, current, nc, nelec_drain);
				//reset the distance run by all charges to zero for next time step calculation
        totalX = 0;
        nelec_drain=0;
      }
			for(i = 1; i < nc; i++) {
				//the waiting times of following charges in the queue are decreased by the (passed) 
				//waiting time of the first charge in the queue
				//Should not minus dwelltime for charges that have already crossed to the other electrode.
				chtemp = getCharge(chA, sequence[i]);
				MinusDwel(chtemp, getDwel(one));
			}
			do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
			//the waiting time of ONE is randomly updated
			setDwel(one, -6.0 * log(ran/RAND_MAX) / (D * getsum(site)));
			//printf("\t NewDwel %g",getDwel(one));
			//the time of the first charge in the queue is increased by the (passed) waiting time 
			//of the first charge in the queue
			Plust(one,(long double) getDwel(one)); 

			//if time of the first charge in the queue is equal to global time
			if(getDwel(one) <= 0){
				printf("4 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n", getDwel(one), getCx(one), getCy(one), getCz(one));
				exit(1);  
			}		
			// ONE is inserted back into the queue at its new position
			insertDwelltimePos(nc, chA, sequence);
			//if time of the first charge in the queue is equal to global time
			if(t < 0){
				printf("t is negative end of loop %Lg \t exit(1)\n", t);
				exit(1); 
			} 
			//charge has moved inside the material; waiting times are updated and the queue is sorted
		}
			//CHARGE INJECTION if time has reached a time step
		if(t >= n * TStep && ((n+1)*NCh)<=Ntot ){
			// n incremented (main loop counter)
			n++; 
			printf("initCharge(nc,n)"); 
			//places N charges on the lattice at x=0, at each time step n
			initCharge( nc, n, chA, sequence, snA);
			nc = nc + NCh;
			printf("Number charges left %d\n",nc);
		}
		// if time of the first charge in the queue has reached the current time step, this 
		// loop calculates the current and velocity distribution and prints them in output files

		//if time of the first charge in the queue is equal to global time
		if(getDwel(one) <= 0){
			printf("5 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n",\
							getDwel(one), getCx(one),getCy(one), getCz(one));
			exit(1);  
		}
		
		//the first charge in the queue has hopped to an adjacent site => then set site as 
		//unoccupied; update step # and global time and compute current
		

		//if time of the first charge in the queue is equal to global time
      
		if(getDwel(one) <= 0){
        printf("6 dwelltime %g \t x %d \t y %d \t z %d \t exit(1)\n", \
						    getDwel(one), getCx(one), getCy(one), getCz(one));
        exit(1);  
		}
		
	}
    
	fclose(fp);

	printf("Finished Run\n");
	printVisitFreq(snA);
		
	deleteChargeA(chA);
}

void insertDwelltimePos(int nc, ChargeArray chA, int * sequence){
    int i;
    int low, high, middle;
    int num;
		Charge chi;
		Charge chj;

		//number of the charge that is first in the sequence 
    num = sequence[0];

    low = 1;
    high = nc -1;
    while(low <= high){
      middle = (low + high) / 2;
			chi = getCharge(chA,sequence[middle]);
			chj = getCharge(chA, num);
      if( getDwel(chi) > getDwel(chj))
        high = middle - 1;
      else if(getDwel(chi) == getDwel(chj)){
        high = middle;
        break;
      }
      else
        low = middle + 1;
    }

    for(i = 0; i < high; i++) {
      sequence[i] = sequence[i+1];
    }
    sequence[high] = num;
}
