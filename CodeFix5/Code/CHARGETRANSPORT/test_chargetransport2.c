#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "../CHARGE/charge.h"
#include "../PARAMETERS/read.h"
//#include "../MEM/mem.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERSITENODE/clustersitenode.h"
#include "../FUNCTIONS/functions.h"
#include "chargetransport.h"

int main(void){


	//mem_init();
	
	int rv;
	long int n;
	int nc;
	int nca;

	int Num_elXb;
	int Num_elXf;
	int Num_elYl;
	int Num_elYr;
	int Num_elZb;
	int Num_elZa;

/*	SiteNode sn;
	int future;

	int SLength;
	int SWidth;
	int SHeight;

	double electricEnergyX;
	double electricEnergyY;
	double electricEnergyZ;

	double ElectricFieldX;
	double ElectricFieldY;
	double ElectricFieldZ;

	int XElecOn;
	int YElecOn;
	int ZElecOn;

	Electrode elXb = NULL;
	Electrode elXf = NULL;
	Electrode elYl = NULL;
	Electrode elYr = NULL;
	Electrode elZb = NULL;
	Electrode elZa = NULL;

	int EndX;
	int EndY;
	int EndZ;

	int PeriodicX;
	int PeriodicY;
	int PeriodicZ;

	int OrderL;

	//Number of time increments charges are injected
	int Tcount;
	//Size of time increments
	double TStep;
	//How many steps pass before we record the data
	int Nstep_av;
	//How many charges are injected at a time
	int NCh;
	//The distance between sites
	double D;
	//rN some normalization value
	int rN;
	//Attempt to Hop
	double AttemptToHop = 1E13;
	//gamma tunneling constant
	double gamma = 2;
	//SiteDistance is the distance between sites [m]
	double SiteDistance = 1E-9;
	double SiteDistanceNM;
	SiteDistanceNM = SiteDistance*1E9;
	//Testing ChargeClosure

	//Relative Permittivity
	double RelativePerm = 3.9;

	//Hopping Rates from Electrodes
	double vX = 1E13;
	double vY = 1E13;
	double vZ = 1E13;

	//Fermi level of Electrode
	double FermiExb = 0;
	double FermiEyl = 0;
	double FermiEzb = 0;
	double FermiExf = 0;
	double FermiEyr = 0;
	double FermiEza = 0;
	double alphaxb = 2;
	double alphayl = 2;
	double alphazb = 2;
	double alphaxf = 2;
	double alphayr = 2;
	double alphaza = 2;

	//ReOrgEnergy [eV]
	double reOrgEnergy = 1;
	//Ntot is the total number of charges
	int Ntot = 20;
	//nc is the number of charges in the system
	int nc1 = Ntot;
	int nca1 = nc1;
	//Total Charges collected at an electrode
	int TotalCollected = 0;

	double MarcusCoeff;
	double MarcusJ0;
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;

	double KT;

	KT = 1;

	//Calculating Marcus J0 coefficient assuming the Attempt to hop Rate
	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	//Calculating full Marcus Coefficient;
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(-2*gamma*SiteDistanceNM);
*/
	long double t;
	int Total;
	SNarray snA;
	matrix Sequence;
	matrix FutureSite;
	ChargeArray chA;
	ParameterFrame PF;
	char FileName[50] = "DataT300Vx0.5Vy0Vz0R11.ckpt";

	printf("Testing: Load_CheckPt\n");
	rv = Load_CheckPt(&t, &snA, &chA, &Sequence, \
	&FutureSite, FileName, &PF,&n, &nc, &nca,\
							 &Num_elXb, &Num_elXf, &Num_elYl, &Num_elYr, &Num_elZb, &Num_elZa);

	assert(rv==0);
	//printSNarray_Detailed(snA);
	//printMatrix(Sequence);
	//printMatrix(FutureSite);
	//printChargeA(chA);

	deleteMatrix(&Sequence);
	deleteMatrix(&FutureSite);
	deleteChargeA(chA);
	deleteParamFrame(&PF);
	deleteSNarray(snA);

	printf("Testing: CheckPt_exist\n");
	char FileName2[256];
	FileName2[0]='\0';
	rv = CheckPt_exist(FileName2,sizeof(FileName2));
	printf("This is the file the function found %s\n",FileName2);
	assert(rv==0);

	char FileName3[256];
	FileName3[0]='\0';
	rv = CheckPt_Latest(FileName3,sizeof(FileName3), 0.5,0,0,300);
	printf("Value of rv %d\n",rv);
	assert(rv>0);
	printf("Lastest version found %s\n",FileName3);
	
	rv = CheckPt_Latest(FileName3,sizeof(FileName3), 1, 0,0,300);
	assert(rv==0);
	printf(".ckpt file exist but it is not the one you want \n");
	
	/*printf("rv %d\n",rv);
	printf("FileName %s\n",FileName2);
	
	
	Load_CheckPt(&t, &Total, &snA, &chA, &Sequence, &FutureSite, FileName, &PF, &nc, &nca,\
							 &Num_elXb, &Num_elXf, &Num_elYl, &Num_elYr, &Num_elZb, &Num_elZa);
	printf("Sucess!\n");
	printChargeA(chA);

	deleteMatrix(Sequence);
	deleteMatrix(FutureSite);
	deleteChargeA(chA);
	deleteParamFrame(PF);
	deleteSNarray(snA);
	*/
	//atexit(mem_term);

	return 0;

}
