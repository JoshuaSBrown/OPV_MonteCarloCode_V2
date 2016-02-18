#include <stdio.h>
#include <stdlib.h>

#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../PARAMETERS/read.h"
//#include "../MEM/mem.h"
#include "functions.h"

int main(void){

	//mem_init();
	//Parameter List
	//
	ParameterFrame PF = newParamFrame();
	int SLength = 20;
	PFset_Len(PF,SLength);
	int SWidth = 20;
	PFset_Wid(PF,SWidth);
	int SHeight = 20;
	PFset_Hei(PF,SHeight);
	
	double RelativePerm = 3.9;
	PFset_RelativePerm(PF,RelativePerm);
	double electricEnergyX = 0;
	double electricEnergyY = 0;
	double electricEnergyZ = 0;
	double KT = 1;
	double reOrgEnergy = 1;
	PFset_reOrg(PF,reOrgEnergy);
	double SiteDistance = 1E-9;
	PFset_SiteDist(PF,SiteDistance);
	int CutOff = 5;
	PFset_CutOff(PF,CutOff);
	double lambda = 1E-9;
	PFset_lambda(PF,lambda);
	int SeedProt = 1;
	PFset_SeedProt(PF,SeedProt);
	int PeriodicX = 0;
	PFset_Px(PF,PeriodicX);
	int PeriodicY = 0;
	PFset_Py(PF,PeriodicY);
	int PeriodicZ = 0;
	PFset_Pz(PF,PeriodicZ);
	
	int XElecOn = 1;
	PFset_XElecOn(PF,XElecOn);
	int YElecOn = 0;
	PFset_YElecOn(PF,YElecOn);
	int ZElecOn = 0;
	PFset_ZElecOn(PF,ZElecOn);
	
	int EndX = 1;
	PFset_EndX(PF,EndX);
	int EndY = 1;
	PFset_EndY(PF,EndY);
	int EndZ = 1;
	PFset_EndZ(PF,EndZ);
	
	double E0 = -5.2;
	PFset_E0(PF,E0);
	double Etrap = -4.9;
	PFset_Etrap(PF,Etrap);
	//double q = 1.602E-19;
	double fracSeed = 0.01;
	PFset_FracSeed(PF,fracSeed);
	double fraction = 0.005;
	PFset_FracTrap(PF,fraction);
	double Tsigma = 0.07;
	PFset_Tsigma(PF,Tsigma);
	double sigma = 0.07;
	PFset_sigma(PF,sigma);
	double gamma = 2;
	PFset_gamma(PF,gamma);
	double AttemptToHop = 10E13;
	PFset_AttemptToHop(PF,AttemptToHop);
	
	double vX = 10E13;
	PFset_vX(PF,vX);
	double vY = 10E13;
	PFset_vY(PF,vY);
	double vZ = 10E13;
	PFset_vZ(PF,vZ);
	double MarcusCoeff = 10E13;
	
	double FermiExb = -2;
	PFset_XFermiB(PF,FermiExb);
	double FermiEyl = -2;
	PFset_YFermiL(PF,FermiEyl);
	double FermiEzb = -2;
	PFset_ZFermiB(PF,FermiEzb);
	double FermiExf = -2;
	PFset_XFermiF(PF,FermiExf);
	double FermiEyr = -2;
	PFset_YFermiR(PF,FermiEyr);
	double FermiEza = -2;
	PFset_ZFermiA(PF,FermiEza);

	double alphaxb = 2;
	PFset_alphaxb(PF,alphaxb);
	double alphayl = 2;
	PFset_alphayl(PF,alphayl);
	double alphazb = 2;
	PFset_alphazb(PF,alphazb);
	double alphaxf = 2;
	PFset_alphaxf(PF,alphaxf);
	double alphayr = 2;
	PFset_alphayr(PF,alphayr);
	double alphaza = 2;
	PFset_alphaza(PF,alphaza);

	Electrode elXb = NULL;
	Electrode elXf = NULL;
	Electrode elYl = NULL;
	Electrode elYr = NULL;
	Electrode elZb = NULL;
	Electrode elZa = NULL;

	const double D = 20E14; 
	PFset_D(PF,D);

	SNarray snA = newSNarray(SLength,SWidth,SHeight);

/*	initJumPossibility(electricEnergy, KT, reOrgEnergy, snA,\
				PeriodicX, PeriodicY, PeriodicZ,\
				XElecOn, YElecOn, ZElecOn);
*/

	//Non - periodic 
	int rv = initSite(electricEnergyX,electricEnergyY,\
			  electricEnergyZ, KT, snA, PF);
	

	printSNarray_Detailed( snA);
	
	char File1[] = "Test1\0";

	printVisitFreq(snA, &File1[0]);

	int Ntot = 100;
	PFset_Ntot(PF,Ntot);	
	int NCh = 20;
	PFset_NCh(PF,NCh);
	int i;
	//Create an array the size of all the charges
	//label the charges starting from 0 to Ntot-1
	matrix sequence = newMatrix(Ntot,1);
	for(i = 0; i<Ntot; i++){
		setE(sequence,i+1,1,i);
	}
	ChargeArray chA = initCharget0( sequence, snA, Ntot, NCh,\
																	D, XElecOn, YElecOn, ZElecOn,\
																	EndX, EndY, EndZ);
	printChargeA(chA);

	rv = initCharge( 20, 2, &chA, sequence, snA,\
					Ntot, NCh, D, XElecOn, YElecOn, ZElecOn,\
					EndX, EndY, EndZ);

	printChargeA(chA);


	printf("Testing initElec\n");

	rv = initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
		      MarcusCoeff, KT, snA,\
		      &elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
		      PF);
	
	
	deleteSNarray(snA);
	deleteChargeA(chA);
	deleteMatrix(sequence);
	
	SNarray snAmini;
	matrix mtxmini;

	if (elXb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXb);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXb);
		deleteMatrix(mtxmini);
		deleteElectrode(&elXb);	
	
	}
	if (elXf!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elXf);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elXf);
		deleteMatrix(mtxmini);
		deleteElectrode(&elXf);	
	}

	if (elYl!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYl);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYl);
		deleteMatrix(mtxmini);
		deleteElectrode(&elYl);	
	}
	if (elYr!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elYr);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elYr);
		deleteMatrix(mtxmini);
		deleteElectrode(&elYr);	
	}
	if (elZb!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZb);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZb);
		deleteMatrix(mtxmini);
		deleteElectrode(&elZb);	
	}

	if (elZa!=NULL){
		snAmini = (SNarray) getElectrode_AdjacentSites(elZa);
		deleteSNarray(snAmini);
		mtxmini = (matrix) getElectrode_HopRates(elZa);
		deleteMatrix(mtxmini);
		deleteElectrode(&elZa);	
	}
	//atexit(mem_term);
	
	
	/*
	//Fully Periodic
	SNarray snA2 = newSNarray(SLength,SWidth,SHeight);
	
	rv = initSite(electricEnergy, KT, reOrgEnergy,\
							lambda, CutOff,\
							fracSeed, fraction, SiteDistance,\
							Etrap, Tsigma, E0, sigma,\
							SeedProt, 1, 1, 1,\
							XElecOn, YElecOn, ZElecOn, snA2);
	
	char File2[] = "Test2\0";

	printVisitFreq(snA2, &File2[0]);

	//Periodic Only in the x
	SNarray snA3 = newSNarray(SLength, SWidth, SHeight);

	rv = initSite(electricEnergy, KT, reOrgEnergy,\
							lambda, CutOff,\
							fracSeed, fraction, SiteDistance,\
							Etrap, Tsigma, E0, sigma,\
							SeedProt, 1, 0, 0,\
							XElecOn, YElecOn, ZElecOn, snA3);
	
	
	char File3[] = "Test3\0";

	printVisitFreq(snA3, &File3[0]);
*/
	
	deleteParamFrame(&PF);

	return 0;
}
