#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../SITENODE/sitenode.h"
#include "montecarlo.h"
#include "../MATRIX/matrix.h"

//Functions
double Rfind( int i1, int j1, int k1, int i2, int j2, int k2, const double SiteDistance) {
	double i1d= (double)i1;
	double i2d= (double)i2;
	double j1d= (double)j1;
	double j2d= (double)j2;
	double k1d= (double)k1;
	double k2d= (double)k2;
	return pow(pow((i1d-i2d)*SiteDistance,2.0)+pow((j1d-j2d)*SiteDistance,2.0)+pow((k1d-k2d)*SiteDistance,2.0),0.5);
}

int CorrCal(const_matrix A,const int i,const int j,const int k,const double Rad, const double SiteEnergy, \
		        const double SiteDistance, const_SNarray snA, double * SumCor, double * SumEcor,double * distanceij, \
						const int SeedProt,const int PeriodicX,const int PeriodicY,const int PeriodicZ,const double lambda){

	//matrix A contains the x y and z positions of all the traps
	//snA contains all the sites in the material

	if(snA==NULL || A==NULL || SumCor==NULL || SumEcor==NULL){
		return -1;
	}
	if(i<0 || i>=getAlen(snA) || j<0 || j>=getAwid(snA) || k<0 || k>=getAhei(snA)){
		return -1;
	}
	if(Rad<0 || SiteDistance<=0 || lambda<=0){
		return -1;
	}
	if(SeedProt>2 || SeedProt<0){
		return -1;
	}
	if(PeriodicX<0 || PeriodicX>1 ||\
		 PeriodicY<0 || PeriodicY>1 ||\
		 PeriodicZ<0 || PeriodicZ>1){
		return -1;
	}

	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	int m;
	int ii, jj, kk;
	double hx, lx, hy, ly, hz, lz;
	double m_x, m_y, m_z;
	double temp;
	double Ediff;
	double CorrFactor;
	double SeedEnergy;

	//Only used if SeedProt is 1 
	//M is used to store all the seeds that are the same 
	//distance from the site
	int inc;
	int chosen;
	double distancetemp;
	matrix M;
	//Initially matrix is just 2 elements
	M = newMatrixSet(2,1,-1.0);
	*distanceij=-1;
	inc=0;
	
	//Cycling through all the seeds
	for (m=0; m<getRows(A);m++){
	
		hx=0, lx=0, hy=0, ly=0, hz=0, lz=0;

		if(PeriodicX==1){
			//Grabbing x position of trap m
			m_x=getE(A,m+1,1)*SiteDistance;
			//Determine boundaries Low x means trap site will be to the left of the box
			temp=(m_x+Rad)/((double)SLength*SiteDistance);
			if ( (temp>lx) !=0) {
				lx=-abs(floor(temp));
			}
			
			//Determine boundaries High x means trap site will be to the right of box
			temp=(m_x-Rad)/((double)SLength*SiteDistance);
			if ( (temp<hx) !=0) {
				hx=abs(floor(temp));
			}
		}

		if(PeriodicY==1) {
			//Grabbing y position of trap m
			m_y=getE(A,m+1,2)*SiteDistance;
			//Determine boundaries low y
			temp=(m_y+Rad)/((double)SWidth*SiteDistance);
			
			if ( (temp>ly) !=0) {
				ly=-abs(floor(temp));
			}

			//Determine boundaries high y
			temp=(m_y-Rad)/((double)SWidth*SiteDistance);
			if ( (temp<hy) !=0) {
				hy=abs(floor(temp));
			}
		}

		if(PeriodicZ==1){
			//Grabbing z position of trap m
			m_z=getE(A,m+1,3)*SiteDistance;
			//Determine boundaries low z
			temp=(m_z+Rad)/((double)SHeight*SiteDistance);
			if ( (temp>lz) !=0) {
				lz=-abs(floor(temp));
			}

			//Determine boundaries high z
			temp=(m_z-Rad)/((double)SHeight*SiteDistance);
			if ( (temp<hz) !=0) {
				hz=abs(floor(temp));
			}
		}

		//This is in the case that the energy of site (i,j,k) is adjusted
		//based on the average of the seeds within Rad
		if( SeedProt==0) {
			for(ii=lx; ii<=hx; ii++){
				for(jj=ly; jj<=hy;jj++){
					for(kk=lz; kk<=hz; kk++){

						// ( i, j, k) is the site we are calculating the 
						// correlation at.

						//When using getE must remember the matrix
						//starts at row 1 col 1 NOT row 0 col 0
						*distanceij = Rfind( i, j, k, \
								(int)getE(A,m+1,1)+ii*SLength, \
								(int)getE(A,m+1,2)+jj*SWidth, \
								(int)getE(A,m+1,3)+kk*SHeight, \
								SiteDistance);

						if(*distanceij<Rad){
			
							//Calculating difference of energy
							//between the seed and the (i,j,k)
							
							SeedEnergy = getEnergy(getSN(snA,\
											(int)getE(A,m+1,1),\
											(int)getE(A,m+1,2),\
											(int)getE(A,m+1,3)));
							
							Ediff = SeedEnergy-SiteEnergy;

							//Multiplying the Difference by the correlation
							CorrFactor = corr(*distanceij, lambda);
						
							*SumEcor += Ediff*CorrFactor;
							*SumCor += CorrFactor;
						}
					}
				}	
			}
		}
	
		//This is in the case that the energy of site (i,j,k) is adjusted
	  //based on the closest seed within Rad
		if((int)SeedProt==1 || (int)SeedProt==2){
		
			for(ii=lx; ii<=hx; ii++){
				for(jj=ly; jj<=hy;jj++){
					for(kk=lz; kk<=hz; kk++){

						//(Calculating the distance between site
						//(i,j,k) and the seed
						distancetemp = Rfind( i, j, k, \
								(int)getE(A,m+1,1)+ii*SLength, \
								(int)getE(A,m+1,2)+jj*SWidth, \
								(int)getE(A,m+1,3)+kk*SHeight, \
								SiteDistance);
					
						//distancetemp has not yet been correctly initialized
						if(*distanceij==-1){
							//Storing the id of the trap site in matrix M
							inc++;
							setE(M,inc,1,(double) m);
							*distanceij=distancetemp;
						}else if(distancetemp<*distanceij) {
							setAll(M,-1.0);
							inc=1;
							setE(M,inc,1,(double) m);
							*distanceij=distancetemp;
						}else if(*distanceij==distancetemp) {
							inc++;
							if(inc>getRows(M)){
								//If there are more than 2 sites 
								//the closest distance apart 
								//make the matrix larger by 2 elements
								resizeRow(&M,getRows(M)+2);
							}
							setE(M,inc,1,(double) m);
						}
						
					}	
				}
			}
		}

	}
	if ((int)SeedProt==1 || (int)SeedProt==2){
	
		if(inc>1){

			//Randomly choosing from the sites
			//that are equal distance away
			chosen=rand()%inc +1;		
			
			//We are simply resuing the temp value
			temp = getE(M,chosen,1);
			m=(int)temp;

							
			SeedEnergy = getEnergy(getSN(snA,\
											(int)getE(A,m+1,1),\
											(int)getE(A,m+1,2),\
											(int)getE(A,m+1,3)));
			Ediff = SeedEnergy-SiteEnergy;
			
			CorrFactor = corr(*distanceij, lambda);
							
			*SumEcor = Ediff*CorrFactor;
			//*SumCor = CorrFactor;
			*SumCor = 1;
		}
		
		if(inc==1) {
		
			
			temp = getE(M,1,1);
			m=(int)temp;
		
							
			SeedEnergy = getEnergy(getSN(snA,\
											(int)getE(A,m+1,1),\
											(int)getE(A,m+1,2),\
											(int)getE(A,m+1,3)));
			Ediff = SeedEnergy-SiteEnergy;

			CorrFactor = corr(*distanceij, lambda);
			
			*SumEcor=Ediff*CorrFactor;
			//*SumCor=CorrFactor;
			*SumCor=1;
		}
	}
	
	
	
	deleteMatrix(&M);

	return 0;
}


///////////////////////////
double corr(const double distance,const double lambda ) {
	return exp(-distance/lambda);
}

int SiteNodeSeedCompatabilityTest(const_SNarray snA, const_matrix Seed){

	if( snA==NULL || Seed==NULL){
		return -1;
	}

	if(getRows(Seed)>getAwid(snA)*getAlen(snA)*getAhei(snA)){
		return -1;
	}

	int m;
	int i;
	int j;
	int k;

	for(m=0;m<getRows(Seed);m++){

		i =	(int)getE(Seed,m+1,1);
		j = (int)getE(Seed,m+1,2);
		k = (int)getE(Seed,m+1,3);

		if(	i<1 || j<1 ||	k<1 || i>getAlen(snA) || j>getAwid(snA) || k>getAhei(snA)){
			return -1;
		}
	}										
	return 0;
}

//////////////////////////
double grn(double m, double s){
	static int iset = 0;
	static double fac, v1, v2;
	double sum;

	if (iset == 0){
		do {
			v1 = 2.0 * (double)rand() / RAND_MAX - 1.0;
			v2 = 2.0 * (double)rand() / RAND_MAX - 1.0;
			sum = v1 * v1 + v2 * v2;
		} while (sum >= 1.0 || sum == 0.0);

		fac = sqrt(-2.0 * log(sum) / sum);
		iset = 1;
		return m+v1*fac*s;

	} else {
		iset=0;
		return m+v2*fac*s;
	}
}

///////////////////////////////
double hoppingRate(double gapEnergy, double KT, double reOrgEnergy){
	return (double) exp(-(gapEnergy + reOrgEnergy) * (gapEnergy + reOrgEnergy) / (4 * KT * reOrgEnergy));  //here using equation "delat = kt" to replace the 'kt'
}

double hoppingRateMillerAbraham(double Energy, double SiteDistance, double alpha, double KT){
	if(Energy<0){
		return exp(-SiteDistance*alpha*2);
	}else{
		return exp(-SiteDistance*alpha*2)*exp(-Energy/KT);;
	}
}

///////////////////////////////
SiteNode getRandomSite(SNarray snA){

	if(snA==NULL){
		return NULL;
	}

	int i, j, k;
	int num;

	srand(rdtsc());
	srandom(rdtsc());
	do{
		num = random() % (getAtotal(snA));
		getLoc(&i,&j,&k,num, snA);
	}while( getInitE( getSN(snA ,i,j,k)) != 0);

	return getSN(snA, i, j, k);
}

////////////////////////////////////
SiteNode getRandomSitePos(int *i, int *j, int *k, SNarray snA){
	int num;
	if (!i || !j || !k) return NULL;

	srand(rdtsc());
	srandom(rdtsc());

	do{
		num = random() % (getAtotal(snA));
		getLoc(&(*i),&(*j),&(*k),num, snA);
	}while( getInitE(getSN(snA,*i,*j,*k)) != 0);

	return getSN(snA,*i,*j,*k);
}

//int rdtsc() {
//	__asm__ __volatile__("rdtsc");	
//	return 0;
//}

uint64_t rdtsc(){
	uint32_t hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ((( uint64_t)lo)|((uint64_t)hi)<<32);
	
}

