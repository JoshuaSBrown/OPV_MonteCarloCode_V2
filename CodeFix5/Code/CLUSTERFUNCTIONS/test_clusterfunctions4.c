#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "clusterfunctions.h"
#include "../CLUSTER/cluster.h"
#include "../SITENODE/sitenode.h"
#include "../MATRIX/matrix.h"
#include "../PARAMETERS/read.h"

int main() {
	
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;
	//printf("Value of fracSeed %e.\n",fracSeed);

	int attempts = 15;
	int rv;

	double SiteDistance = 1;     	//units of [nm]
	double AttemptToHop = 1E-13;
	double gamma = 2;							//Units of [1/nm]
	int XElecOn = 1;
	int YElecOn = 0;
	int ZElecOn = 0;

	double reOrgEnergy = 1;
	double KT = 1;

	double MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	//Calculating full Marcus Coefficient;
	double MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(2*gamma*SiteDistance);


  printf("COMPLETED test_cluster4 success\n");
	return 0;
}
