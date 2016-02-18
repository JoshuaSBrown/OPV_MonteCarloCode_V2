#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include <stdint.h>

#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"

/* Rfind determines the the magnitude of the distance given
   the two sites. Site 1 (i1, j1, k1) and Site 2 is 
   ( i2, j2, k2). 
*/
double Rfind( int i1, int j1, int k1, int i2, int j2, int k2, const double SiteDistance);

/*Calculates the correlation between site i,j,k with all the sites 
	in the matrix A. Uses the Radius to determine which sites will 
	contribute to the correlation. The SNarray is used to determine
	the actual energetic contribution from the sites listed in the 
	matrix. The site energy is the energy of site i,j,k before
	correlation has been factored in.
	*/
int CorrCal(const_matrix A,const int i,const int j,const int k,const double Rad,const double SiteEnergy,\
		 const double SiteDistance, const_SNarray snA, double * SumCor, double * SumEcor,double * distanceij, \
		 const int SeedProt, const int PeriodicX, const int PeriodicY, const int PeriodicZ, const double lambda); 

/* Returns the correlation value given the distance
   between two sites. Calculates it using an exponential 
   of the form e^(distance/lambda).
*/
double corr( const double distance, const double lambda );

/* Checks to ensure the seed matrix has the appropriate data
	 and is compatable with the SNarray 
	 returns a 0 if compatable and a -1 if not
*/
int SiteNodeSeedCompatabilityTest(const_SNarray snA, const_matrix Seed);

/* Here grn is the gaussian rand number generator which 
   generates numbers based on a gaussian distribution it is
   used to determine the site energies of the seeds
   m is the energy about which the distribution is focused
   s is sigma
*/
double grn(double m, double s);

/* The Marcus equation is used to determine the hopping rate
*/
double hoppingRate( double gapEnergy, double KT, double reOrgEnergy);

/* The Miller & Abrahams expression for hopping rate
*/
double hoppingRateMillerAbraham( double Energy, double SiteDistance, double alpha, double KT);

/* This function returns the pointer to a random node
*/
SiteNode getRandomSite(SNarray snA);

/* This function returns the pointer to a random node as 
   well as the (i, j, k) location of the site in the 
   lattice
*/
SiteNode getRandomSitePos( int *i, int *j, int *k, SNarray snA);

uint64_t rdtsc(void);

#endif
