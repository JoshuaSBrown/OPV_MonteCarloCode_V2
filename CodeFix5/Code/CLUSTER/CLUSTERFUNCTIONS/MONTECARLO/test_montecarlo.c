#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"
#include "montecarlo.h"

int main() {

	int len = 11;
	int wid = 11;
	int hei = 11;
	double Rad = 5E-9;

	int SeedProt = 0;
	int PeriodicX = 0;
	int PeriodicY = 0;
	int PeriodicZ = 0;

	double SumEcor = 0;
	double SumCor = 0;
	double distanceij = 0;

	double SiteEnergy = 0;
	double SiteDistance = 1E-9;
	double lambda = 1E-9;

	//Creating Cubic site node array
	SNarray snA = newSNarray(len, wid, hei);

	//Placing Seed in the middle of sitenodes
	matrix seed = newMatrix(1,3);
	//The Sites are labeled from 0 to length-1, width-1, height-1	
	setE(seed,1,1,6);
	setE(seed,1,2,6);
	setE(seed,1,3,6);
	printf("Seed Matrix\n");
	printMatrix(seed);

	//Set Energy of the seed to 1.0 eV
	int rv = setEnergy(getSN(snA,6,6,6),1.0);
	assert(rv==0);

	//Position of site to correlate
	//Site is just to the left of the
	//seed
	int i = 5;
	int j = 6;
	int k = 6;
	//Notice that the seed is at position 6,6,6 in the
	//matrix notation this corresponds to 6,6,6 in the 
	//regular array notation and corresponds to the position
	//in the lattice

	printf("Testing:CorrCal\n");
	rv = CorrCal(NULL, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, -1, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, len, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, -1, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, wid, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, -1, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);
	
	rv = CorrCal(seed, i, j, hei, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, -1, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 0, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, NULL, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, NULL, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, NULL, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, -1,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij,  2,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 -1, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 2, PeriodicY, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, -1, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij,SeedProt,\
					 PeriodicX, 2, PeriodicZ, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, -1, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, 2, lambda);
	assert(rv==-1);

	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, 0);
	assert(rv==-1);

	//Testing if site is one site away from seed
	printf("Testing if site i,j,k is one site away from seed\n");
	printMatrix(seed);
	rv = CorrCal(seed, i, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-SiteDistance/lambda));
	
	SumCor = 0;
	SumEcor = 0;

	//Testing if site is two sites away from seed
	printf("Testing if site i,j,k is two sites away from seed\n");
	rv = CorrCal(seed, i-1, j, k, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-(2*SiteDistance)/lambda));
	
	SumCor = 0;
	SumEcor = 0;

	//Testing if site is outside radius (Rad=5);
	rv = CorrCal(seed, 0, 0, 0, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, SeedProt,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);
	
	//Placing Two Seeds removed from the center
	//node by two sites
	matrix seed2 = newMatrix(2,3);
	setE(seed2,1,1,8);
	setE(seed2,1,2,6);
	setE(seed2,1,3,6);

	setE(seed2,2,1,4);
	setE(seed2,2,2,6);
	setE(seed2,2,3,6);
	//Resetting Energy of sites
	setEnergy(getSN(snA,6,6,6),0.0);
	setEnergy(getSN(snA,4,6,6),1.0);
	setEnergy(getSN(snA,8,6,6),1.0);
	
	rv = CorrCal(seed2, 6, 6, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij,  0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);

	//Values of SumEcor and SumCor should be a combination of 
	//both seed sites for SeedProt = 0
	
	assert(SumCor==SumEcor);
	assert(SumEcor==2*exp(-SiteDistance*2/lambda));

	SumCor = 0;
	SumEcor = 0;

	rv = CorrCal(seed2, 6, 6, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 1 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);

	//Values of SumEcor and SumCor should not be a combination
	//but should be based on one or the other seeds
	//This is for SeedProt = 1
	assert(SumCor==1);
	assert(SumEcor==exp(-SiteDistance*2/lambda));
	
	//Check that the correlation is based on the closest seed
	setEnergy(getSN(snA,4,6,6),0.0);
	setEnergy(getSN(snA,3,6,6),2.0);
	//One of the seeds is three sites away with energy 2 the 
	//other site is 2 sites away with energy 1

	SumCor = 0;
	SumEcor = 0;
	
	matrix seed3 = newMatrix(2,3);
	setE(seed3,1,1,8);
	setE(seed3,1,2,6);
	setE(seed3,1,3,6);

	setE(seed3,2,1,3);
	setE(seed3,2,2,6);
	setE(seed3,2,3,6);
	
	rv = CorrCal(seed3, 6, 6, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, 1 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);

	assert(SumCor==1);
	assert(SumEcor==exp(-SiteDistance*2/lambda));

	//Check that correct information when made periodic in the x
	printf("Testing Periodic in x\n\n");

	SumCor = 0;
	SumEcor = 0;

	matrix seed4 = newMatrix(1,3);
	setE(seed4,1,1,0);
	setE(seed4,1,2,6);
	setE(seed4,1,3,6);

	setEnergy(getSN(snA,3,6,6),0.0);
	setEnergy(getSN(snA,8,6,6),0.0);
	setEnergy(getSN(snA,0,6,6),1);


	rv = CorrCal(seed4, 10, 6, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);

	rv = CorrCal(seed4, 10, 6, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 1, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-SiteDistance/lambda));

	//Check that correct information when made periodic in the y	
	printf("Testing Periodic in y\n\n");

	SumCor = 0;
	SumEcor = 0;
	
	setE(seed4,1,1,6);
	setE(seed4,1,2,0);
	setE(seed4,1,3,6);

	setEnergy(getSN(snA,0,6,6),0);
	setEnergy(getSN(snA,6,0,6),1);

	rv = CorrCal(seed4, 6,10, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);

	rv = CorrCal(seed4, 6, 10, 6, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, 1, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-SiteDistance/lambda));

	//Check that correct information when made periodic in the z
	printf("Testing Periodic in z\n\n");

	SumCor = 0;
	SumEcor = 0;
	
	setE(seed4,1,1,6);
	setE(seed4,1,2,6);
	setE(seed4,1,3,0);

	setEnergy(getSN(snA,6,0,6),0);
	setEnergy(getSN(snA,6,6,0),1);

	rv = CorrCal(seed4, 6, 6, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij, 0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);

	rv = CorrCal(seed4, 6, 6, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, PeriodicY, 1, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-SiteDistance/lambda));

	//Check that correct information when made periodic in xyz
	printf("Testing Periodic in xyz\n\n");

	SumCor = 0;
	SumEcor = 0;
	
	setE(seed4,1,1,6);
	setE(seed4,1,2,6);
	setE(seed4,1,3,0);

	setEnergy(getSN(snA,6,0,6),0);
	setEnergy(getSN(snA,6,6,0),1);

	printf("Test1\n");
	rv = CorrCal(seed4, 6, 6, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);

	printf("Test2\n");
	rv = CorrCal(seed4, 6, 6, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 PeriodicX, PeriodicY, 1, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==exp(-SiteDistance/lambda));

	//Check that correct information when made periodic in xyz
	
	SumCor = 0;
	SumEcor = 0;

	matrix seed5 = newMatrix(3,3);
	setE(seed5,1,1,0);
	setE(seed5,1,2,10);
	setE(seed5,1,3,10);
	setE(seed5,2,1,10);
	setE(seed5,2,2,0);
	setE(seed5,2,3,10);
	setE(seed5,3,1,10);
	setE(seed5,3,2,10);
	setE(seed5,3,3,0);

	setEnergy(getSN(snA,6,6,0),0);
	//Will exist same xy position will be above the box
	rv = setEnergy(getSN(snA,10,10,0),1);
	assert(rv==0);
	//Will exist same xz position will be above the box
	setEnergy(getSN(snA,10,0,10),1);
	//Will exist same yz position will be above the box
	setEnergy(getSN(snA,0,10,10),1);

	rv = SiteNodeSeedCompatabilityTest( snA, seed5);
	assert(rv==-1);

	printf("Test3\n");
	rv = CorrCal(seed5, 10, 10, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor, &distanceij,0 ,\
					 PeriodicX, PeriodicY, PeriodicZ, lambda);
	assert(rv==0);
	assert(SumCor==0);
	assert(SumEcor==0);

	printf("Test4\n");
	rv = CorrCal(seed5, 10, 10, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 0 ,\
					 1, 1, 1, lambda);
	assert(rv==0);
	assert(SumCor==SumEcor);
	assert(SumEcor==3*exp(-SiteDistance/lambda));
	
	printf("Test5\n");
	rv = CorrCal(seed5, 10, 10, 10, Rad, SiteEnergy, \
					 SiteDistance, snA, &SumCor, &SumEcor,&distanceij, 1 ,\
					 1, 1, 1, lambda);
	assert(rv==0);
	assert(SumCor==1);
	assert(SumEcor==exp(-SiteDistance/lambda));

	printf("Testing:grn\n");

	SNarray snA2 = newSNarray(50, 50, 50);
	double TempE=0;

	char * FileName = "Testgrn.xyz";
	FILE * Out;
	if(( Out = fopen( FileName,"w"))!=NULL){

		fprintf(Out,"%d\n\n",getAtotal(snA2));
		for (i=0;i<getAlen(snA2);i++){
			for (j=0;j<getAwid(snA2);j++){
				for (k=0;k<getAhei(snA2);k++){
					TempE = grn(0.0,0.07);		
					setEnergy(getSN(snA2,i,j,k),TempE);
					fprintf(Out,"C\t%f\t%f\t%f\t%f\n",(double)i,(double)j,(double)k,getEnergy(getSN(snA2,i,j,k)));
				}
			}
		}
		fclose(Out);
	}

	printf("Testing:hoppingRate\n");

	double rvd = hoppingRate( 1, 1, 1);

	assert(rvd==exp(-(1+1)*(1+1)/(4*1*1)));

	printf("Testing:getRandomSite\n");
	SNarray snA3 = newSNarray(3, 3, 3);
	
	SiteNode sn1 = getRandomSite(NULL);
	assert(sn1==NULL);

	//Randomly sampling a 3x3x3 sitenode array
	//Ensuring that the values are all correct
	//For a 100 random samples
	for(i=0;i<100;i++){
		sn1 = getRandomSite(snA3);
		assert(getInitE(sn1)==0);
		assert(getDwelStat(sn1)==-1);
		assert(getVisFreq(sn1)==0);
		assert(getVis(sn1)==0);
		assert(getEnergy(sn1)==0);
		assert(getType(sn1)==0);
		assert(getPoint(sn1)!=NULL);
		assert(getClusterList(sn1)==NULL);
	}

	printf("Testing:getRandomSitePos\n");

	int m;

	for(m=0;m<100;m++){
		sn1 = getRandomSitePos(&i,&j,&k,snA3);
		assert(getInitE(sn1)==0);
		assert(getDwelStat(sn1)==-1);
		assert(getVisFreq(sn1)==0);
		assert(getVis(sn1)==0);
		assert(getEnergy(sn1)==0);
		assert(getType(sn1)==0);
		assert(getPoint(sn1)!=NULL);
		assert(getClusterList(sn1)==NULL);
		assert(i>=0);
		assert(i<getAlen(snA3));
		assert(j>=0);
		assert(j<getAwid(snA3));
		assert(k>=0);
		assert(k<getAhei(snA3));
	}


	deleteMatrix(seed);
	deleteMatrix(seed2);
	deleteMatrix(seed3);
	deleteMatrix(seed4);
	deleteMatrix(seed5);
	deleteSNarray(snA);	
	deleteSNarray(snA2);	
	deleteSNarray(snA3);	
	return 0;
}
