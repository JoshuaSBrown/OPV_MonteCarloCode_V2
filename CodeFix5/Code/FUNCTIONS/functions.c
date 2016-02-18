#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

#include "../PARAMETERS/read.h"
//#include "../MEM/mem.h"
#include "../CHARGE/charge.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////


int initSite(const double electricEnergyX, const double electricEnergyY,\
		const double electricEnergyZ,double KT,SNarray snA,ParameterFrame PF){
	
	printf("Entering initsite procedure.\n");
	if(snA==NULL ){
		return -1;
	}

	//Grabbing parameters from parameter frame
	double reOrgEnergy = PFget_reOrg(PF);
	double AttemptToHop = PFget_AttemptToHop(PF);
	double gamma = PFget_gamma(PF);
	double lambda = PFget_lambda(PF);
	double CutOff = PFget_CutOff(PF);
	double fracSeed = PFget_FracSeed(PF);
	double fraction = PFget_FracTrap(PF);
	double SiteDistance = PFget_SiteDist(PF);
	double Etrap = PFget_Etrap(PF);
	double Tsigma = PFget_Tsigma(PF);
	double E0 = PFget_E0(PF);
	double sigma = PFget_sigma(PF);
	int SeedProt = PFget_SeedProt(PF);
	int PeriodicX = PFget_Px(PF);
	int PeriodicY = PFget_Py(PF);
	int PeriodicZ = PFget_Pz(PF);
	int XElecOn = PFget_XElecOn(PF);
	int YElecOn = PFget_YElecOn(PF);
	int ZElecOn = PFget_ZElecOn(PF);

	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	int traps, seeds, sites;
	int i, j, k, m;
	int ii, jj, kk;
	int rv;
	int index;
	double seed_dist;
	double trap_dist;
	double SiteEnergy;
	double SumEcor, SumEcor2;
	double SumCor, SumCor2;
	double MarcusJ0;
	double MarcusCoeff;	
	double CorRad=lambda*CutOff;
	double CorRadtemp;
	double percent;
	double SiteDistanceNM;

	//Converting to [nm]
	SiteDistanceNM = SiteDistance*1E9;
	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;
	//printf("Value of fracSeed %e.\n",fracSeed);
	SiteNode pN;
	sites= SLength * SWidth * SHeight;
	traps = (int)((double)sites * fraction);
	seeds = (int)((double)sites * fracSeed);

	//Calculating Marcus J0 coefficient assuming the Attempt to hop Rate
	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	
	//printf("Value of MarcusJ0 %g AttemptToHop %g kB %g\n",MarcusJ0, AttemptToHop, kB);
	//Calculating full Marcus Coefficient;
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(-2*gamma*SiteDistance);

	//printf("Value of Marcus Coeff at begginning %g hbar %g reOrgEnergy %g KT %g gamma %g SiteDistanceNM %g\n",MarcusCoeff, hbar, reOrgEnergy, KT, gamma, SiteDistanceNM);

	if (CorRad < SiteDistance) {
		printf("WARNING CorRad is less than SiteDistance!!!\n");
		printf("Neighbors are not within Correlation distance\n");
		int temp=1;
		while (temp*lambda< SiteDistance) {
			temp++;
		}
		printf("Adjusted CutOff to %d\n", temp);
		CorRad=temp*lambda;
	}

	printf("Value of Seeds %d.\n",seeds);
	printf("Value of Sites %d.\n",sites);
	printf("Value of Traps %d.\n",traps);

	if (seeds+traps > (sites) ){
		seeds=(sites) - traps;
		printf("Total Sites %d.\n",sites);
		printf("Traps %d.\n",traps);
		printf("Sites to be Seeded %d.\n",seeds);
	}

	matrix AsTr=newMatrix(traps,3);

	printf("Created Array\n");

	if (AsTr==NULL) {
		printf("WARNING AsTr returned NULL\n");
		//assert(AsTr!=NULL);
	}

	setDefaultSNa(snA);


	if(traps!=0){
		printf("Randomly determining which sites are traps & assigning energies.\n");
		percent=0;
		i=0;
		//Assign energies for traps
		while(i < traps){
			pN = getRandomSitePos( &ii, &jj, &kk, snA);
			if (getInitE(pN)==0) {

				if ( ((double)i)/((double)traps)>percent){
					printf("Percent Complete %ld\n",(long int)(percent*100));
					percent+=0.1;
				}
				setE(AsTr,i+1,1,(double)ii);
				setE(AsTr,i+1,2,(double)jj);
				setE(AsTr,i+1,3,(double)kk);

				SiteEnergy = grn(Etrap, Tsigma);
				setEnergy(pN,SiteEnergy);
				setInitE(pN,1);
				i++;
			}
		}
		printf("Percent Complete %ld\n",(long int)(100));
	}

	//In the case that there are more correlated energies than assigned energies
	//
	if ( seeds+traps < sites/2 ){

		matrix As=newMatrix(seeds,3);

		if (As==NULL) {
			printf("WARNING Malloc returned NULL for Assigned Matrix\n");
		}
		i=0;
		printf("Calculating Energies for Seeds\n");
		percent=0;
		while( i < seeds) {
			pN = getRandomSitePos( &ii, &jj, &kk, snA);
			//if Site has not already been assigned energy
			if ( getInitE(pN) == 0 ) {         
				//Store the pointer in an array (Array contains sites that have been assigned energies)
				if (((double)i)/((double)seeds)>percent) {
					printf("Percent Complete %ld\n",(int long)(percent*100));
					percent+=0.1;
				}
				setE(As,i+1,1,(double)ii);
				setE(As,i+1,2,(double)jj);
				setE(As,i+1,3,(double)kk);
				//printf("Assigned locations to array\n");
				if(SeedProt<2){
					SiteEnergy = grn(E0, sigma);
				}else if(SeedProt==2){
					SiteEnergy = E0;
				}

				setEnergy(pN,SiteEnergy);
				setInitE(pN,1);
				i++;
			}
		}
		printf("Percent Complete %ld\n",(int long)(100));
		//printVisitFreq(snA);

		m=0;
		printf("Calculating Energies for correlated sites1\n");
		percent=0;
		for (i= 0; i < SLength; i++){
			for(j = 0; j < SWidth; j++){
				for(k = 0; k < SHeight; k++){
					SumEcor=0;
					SumCor=0;
					CorRadtemp=CorRad;
					//printf("Initial ID %d initE %d\n",getIndex(snA,i,j,k), getInitE(getSN(snA,i,j,k)));
					while (SumEcor==0 && getInitE(getSN(snA,i,j,k))==0){

						SiteEnergy=grn(E0, sigma);
						//assert(As!=NULL);
						//Accounting for correlation from seeds
						rv = CorrCal(As, i,j,k, CorRadtemp, SiteEnergy, SiteDistance,snA, &SumCor, &SumEcor, &seed_dist,\
										SeedProt, PeriodicX, PeriodicY, PeriodicZ, lambda);
							
						//Accounting for correlation from traps
						if (traps!=0){
							if(SeedProt==0){
								printf("WARNING: Effect of both seeds and traps has not been correctly accounted for SeedProt = 0!\n");
							}
								
							CorrCal(AsTr, i,j,k, CorRadtemp, SiteEnergy, SiteDistance, snA, &SumCor2, &SumEcor2, &trap_dist,\
										SeedProt, PeriodicX, PeriodicY, PeriodicZ, lambda);
						
							if(trap_dist<seed_dist && SumEcor2!=0){
								SumCor = SumCor2;
								SumEcor = SumEcor2;
							}
						}
						if (SumEcor==0) {
							CorRadtemp=CorRadtemp*2;
						}else{
							setEnergy(getSN(snA,i,j,k),SiteEnergy+SumEcor/SumCor);
							setInitE(getSN(snA,i,j,k),1);
							if ( ((double)m)/((double)sites-(seeds+traps))>percent){
								printf("Percent Complete %ld\n",(int long)(percent*100));
								percent+=0.1;
							}
							
							m++;
							
						}
						
						if(CorRadtemp>1){
							exit(1);
						}
					
					}
					
				}
				
			}
			
		}
		printf("Percent Complete %ld\n",(int long)(100));

		printf("Deleting Matrix As.\n");
		deleteMatrix(&As);
	}

	if ( seeds+traps >= sites/2 ) {
		//In the case that there are more assigned energies than correlated energies
		//Assign all energies

		matrix UnAs=newMatrix((sites-seeds-traps),3);
		matrix As=newMatrix(seeds,3);

		i=0;

		if (As == NULL ){
			printf("ERROR Malloc returned NULL for As Matrix\n");
		  exit(1);
		}
		if( UnAs == NULL) {
			printf("WARNING Malloc returned NULL for UnAs Matrix\n");
		}
		
		//pick unoccupied sites
		
		if(UnAs!=NULL){
			printf("Randomly determining which sites to correlate.\n");
			percent=0;
			while( i<(sites-(seeds+traps))) {
				pN = getRandomSitePos( &ii, &jj, &kk, snA);
				//printf("Value of pN->initE %d\n",pN->initE);
				if ( getInitE(pN) == 0 ) {         
					if (((double)i)/((double)(sites-(seeds+traps)))>percent) {
						printf("Percent Complete %ld\n",(long int) (percent*100));
						percent+=0.1;
					}
					//Choosing sites that have not already been assigned energy
					//For these sites the correlation function will be used to calculate their 
					//energy. 
					setE(UnAs,i+1,1,(double)ii);
					setE(UnAs,i+1,2,(double)jj);
					setE(UnAs,i+1,3,(double)kk);
					setInitE(pN,1);
					i++;
				}
			}
			printf("Percent Complete %ld\n",(long int)(100));
		}

		printf("Calculating Energies for Seeds\n");
		percent=0;
		m=0;
		for(i = 0; i < SLength; i++){
			for(j = 0; j < SWidth; j++){
				for(k = 0; k < SHeight; k++){
					if( getInitE(getSN(snA,i,j,k))==0){   //if site is not considered a trap

						if (((double)m)/((double)(seeds))>percent){
							printf("Percent Complete %ld\n",(long int)(percent*100));
							percent+=0.1;
						}
						SiteEnergy  = grn(E0, sigma);  
						setEnergy(getSN(snA,i,j,k),SiteEnergy);
						setInitE(getSN(snA,i,j,k),1); //initialize  energy for seeds
						setE(As,m+1,1,(double)i);
						setE(As,m+1,2,(double)j);
						setE(As,m+1,3,(double)k);
						m++;
					}
				}
			}
		}
		printf("Percent Complete %ld\n",(long int)(100));

		//Here we are cycling through the sites that have not been assigned energies
		if(UnAs!=NULL){
			printf("Calculating Energies for correlated sites2\n");
			percent=0;
			for( index=0; index<(sites-(seeds+traps));index++) {

				
				printf("index %d\n",index);
				SiteEnergy=grn(E0, sigma);
				CorRadtemp=CorRad;
				//Calculate Correlation Energies for the UnAssigned positions
				SumEcor=0;
				SumCor=0;

				if (((double)index)/((double)(sites-(seeds+traps)))>percent) {
					printf("Percent Complete %ld\n",(long int) (percent*100));
					percent+=0.1;
				}	

				while (SumEcor==0) {

					i=(int)getE(UnAs,index+1,1);
					j=(int)getE(UnAs,index+1,2);
					k=(int)getE(UnAs,index+1,3);

					//Accounting for correlation from seeds
					CorrCal(As, i, j, k, CorRadtemp, SiteEnergy, SiteDistance, snA, &SumCor, &SumEcor, &seed_dist,\
							SeedProt, PeriodicX, PeriodicY, PeriodicZ, lambda);

					//Accounting for correlation from traps
					if (traps!=0) {
						CorrCal(AsTr, i,j,k , CorRadtemp, SiteEnergy, SiteDistance, snA, &SumCor, &SumEcor, &trap_dist,\
								SeedProt, PeriodicX, PeriodicY, PeriodicZ, lambda);
						if(trap_dist<seed_dist && SumEcor2!=0){
							SumEcor=SumEcor2;
							SumCor=SumCor2;
						}
					}
					if (SumEcor==0) {
						CorRadtemp=CorRadtemp*2;
					}else {
						setEnergy(getSN(snA,i,j,k),SiteEnergy+SumEcor/SumCor);
						setInitE(getSN(snA,i,j,k),1);
					}
				}
			}
			printf("Percent Complete %ld\n",(long int) (100));

		}
		deleteMatrix(&As);
		deleteMatrix(&UnAs);
	}

	rv=deleteMatrix(&AsTr);
	
	printf("Finished initializing.\n");
	//initialize jumping possibility; that is initializing sum and p[6] inside this function
	//For all the sitenodes
	initJumPossibility(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  MarcusCoeff, KT,reOrgEnergy, snA,\
				PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);

	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////
int initElec(const double electricEnergyX, const double electricEnergyY,\
		  const double electricEnergyZ, const double MarcusCoeff,\
		  const double KT, SNarray snA,\
			Electrode * elXb, Electrode * elXf, Electrode * elYl,\
			Electrode * elYr, Electrode * elZb, Electrode * elZa,\
			ParameterFrame PF){

	//Grabbing parameters from parameter frame
	double reOrgEnergy = PFget_reOrg(PF);
	double RelativePerm = PFget_RelativePerm(PF);
	double SiteDistance = PFget_SiteDist(PF);
	int PeriodicX = PFget_Px(PF);
	int PeriodicY = PFget_Py(PF);
	int PeriodicZ = PFget_Pz(PF);
	int XElecOn = PFget_XElecOn(PF);
	int YElecOn = PFget_YElecOn(PF);
	int ZElecOn = PFget_ZElecOn(PF);
	double FermiExb = PFget_XFermiB(PF);
	double FermiExf = PFget_XFermiF(PF);
	double FermiEyl = PFget_YFermiL(PF);
	double FermiEyr = PFget_YFermiR(PF);
	double FermiEzb = PFget_ZFermiB(PF);
	double FermiEza = PFget_ZFermiA(PF);
	double alphaxb = PFget_alphaxb(PF);
	double alphaxf = PFget_alphaxf(PF);
	double alphayl = PFget_alphayl(PF);
	double alphayr = PFget_alphayr(PF);
	double alphazb = PFget_alphazb(PF);
	double alphaza = PFget_alphaza(PF);
	double vX = PFget_vX(PF);
	double vY = PFget_vY(PF);
	double vZ = PFget_vZ(PF);

	if((XElecOn+YElecOn+ZElecOn)>1){
		printf("ERROR: Only equipped to deal with one electrode on at a time!\n");
		exit(1);
	}

	if(snA==NULL || KT<0 || SiteDistance<0 || RelativePerm < 0){
		return -1;
	}

	if(XElecOn>1 || YElecOn>1 || ZElecOn>1 ||\
		 XElecOn<0 || YElecOn<0 || ZElecOn<0){
		printf("ERROR: Electron On parameters incorrec in initElec!\n");
		return -1;
	}

	if(PeriodicX<0 || PeriodicY<0 || PeriodicZ<0 ||\
		 PeriodicX>1 || PeriodicY>1 || PeriodicZ>1){
		printf("ERROR: Periodic Conditions incorrect in initElec!\n");
		return -1;
	}


	int i;
	int j;
	int k;
	double Energy;
	
	int SLength;
	int SWidth;
	int SHeight;

	matrix Xb1;
	matrix Xb2;
	matrix Xf1;
	matrix Xf2;

	matrix Yl1;
	matrix Yl2;
	matrix Yr1;
	matrix Yr2;

	matrix Zb1;
	matrix Zb2;
	matrix Za1;
	matrix Za2;

	SNarray snAelec;
	SiteNode sn;

	SLength = getAlen(snA);
	SWidth = getAwid(snA);
	SHeight = getAhei(snA);

	if(XElecOn==1){
		printf("Initializing X Electrodes\n");
		//Creating the back electrode on the x axis
		*elXb = newElectrode();
		setElectrode_alpha(*elXb, alphaxb);
		setElectrode_FermiEnergy(*elXb, FermiExb);
		//Only saving energies of sites next to electrodes
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

		//printf("j %d k %d\n",j,k);
		*elXf = newElectrode();
		setElectrode_alpha(*elXf, alphaxf);
		setElectrode_FermiEnergy(*elXf, FermiExf);
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
	
		//printf("j %d k %d\n",j,k);
		printf("Initializing jump rates to back electrode\n");
		//Calculating jump rates for back electrode
		initJumPossibility_ElecX(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Xb1, Xb2, *elXb,\
			 	RelativePerm, vX, SWidth, SHeight, PeriodicY, PeriodicZ,\
				YElecOn, ZElecOn, 0);

	snAelec = (SNarray) getElectrode_AdjacentSites(*elXb);
	//printf("Length %d Width %d Height %d\n",getAlen(snAelec),getAwid(snAelec),getAhei(snAelec));
		
	deleteMatrix(&Xb1);
		deleteMatrix(&Xb2);
		
		printf("Initializing jump rates to front electrode\n");
		//Calculating jump rates for forward electrode
		initJumPossibility_ElecX(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Xf1, Xf2, *elXf,\
			 	RelativePerm, vX,SWidth, SHeight, PeriodicY, PeriodicZ,\
				YElecOn, ZElecOn, 1);
		
		deleteMatrix(&Xf1);
		deleteMatrix(&Xf2);


	snAelec = (SNarray) getElectrode_HopRates(*elXb);
	//printf("Length %d Width %d Height %d\n",getAlen(snAelec),getAwid(snAelec),getAhei(snAelec));

	}

	if(YElecOn==1){
		printf("Initializing Y Electrodes\n");
		*elYl = newElectrode();
		setElectrode_alpha(*elYl, alphayl);
		setElectrode_FermiEnergy(*elYl, FermiEyl);
		//Only saving energies of sites next to electrodes
		Yl1 = newMatrix(SLength,SHeight);
		Yl2 = newMatrix(SLength,SHeight);
		
		for(i=0;i<SLength;i++){
			for(k=0;k<SHeight;k++){
				sn = getSN(snA,i,0,k);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Yl1,i+1,k+1,Energy);
				sn = getSN(snA,i,1,k);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Yl2,i+1,k+1,Energy);	
			}
		}

		*elYr = newElectrode();
		setElectrode_alpha(*elYr, alphayr);
		setElectrode_FermiEnergy(*elYr, FermiEyr);
		//Only saving energies of sites next to electrodes
		Yr1 = newMatrix(SLength,SHeight);
		Yr2 = newMatrix(SLength,SHeight);
		
		for(i=0;i<SLength;i++){
			for(k=0;k<SHeight;k++){
				sn = getSN(snA,i,SWidth-1,k);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Yr1,i+1,k+1,Energy);
				sn = getSN(snA,i,SWidth-2,k);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Yr2,i+1,k+1,Energy);	

				}
		}

		printf("Initializing jump rates to left electrode\n");
		//Calculating jump rates for forward electrode
		initJumPossibility_ElecY(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Yl1, Yl2, *elYl,\
			 	RelativePerm, vY, SLength,SHeight, PeriodicX, PeriodicZ,\
				XElecOn, ZElecOn, 0);
		
		deleteMatrix(&Yl1);
		deleteMatrix(&Yl2);
		
		printf("Initializing jump rates to right electrode\n");
		//Calculating jump rates for forward electrode
		initJumPossibility_ElecY(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Yr1, Yr2, *elYr,\
			 	RelativePerm, vY, SLength, SHeight, PeriodicX, PeriodicZ,\
				XElecOn, ZElecOn, 1);
	
		deleteMatrix(&Yr1);
		deleteMatrix(&Yr2);
	}

	if(ZElecOn==1){
		printf("Initializing Z Electrodes\n");
		
		*elZb = newElectrode();
		setElectrode_alpha(*elZb, alphazb);
		setElectrode_FermiEnergy(*elZb, FermiEzb);
		//Only saving energies of sites next to electrodes
		Zb1 = newMatrix(SLength,SWidth);
		Zb2 = newMatrix(SLength,SWidth);
		
		for(i=0;i<SLength;i++){
			for(j=0;j<SWidth;j++){
				sn = getSN(snA,i,j,0);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Zb1,i+1,j+1,Energy);
				sn = getSN(snA,i,j,1);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Zb2,i+1,j+1,Energy);	
			}
		}

		*elZa = newElectrode();
		setElectrode_alpha(*elZa, alphaza);
		setElectrode_FermiEnergy(*elZa, FermiEza);
		//Only saving energies of sites next to electrodes
		Za1 = newMatrix(SLength,SWidth);
		Za2 = newMatrix(SLength,SWidth);
		
		for(i=0;i<SLength;i++){
			for(j=0;j<SWidth;j++){
				sn = getSN(snA,i,j,SHeight-1);
				Energy = getEnergy(sn);
				//Energies of first plane
				setE(Za1,i+1,j+1,Energy);
				sn = getSN(snA,i,j,SHeight-2);
				Energy = getEnergy(sn);
				//Energies of second plane
				//This is needed so we can calculate
				//the hops not just to the electrode
				//but further into the system as well
				setE(Za2,i+1,j+1,Energy);	
			}
		}

		printf("Initializing jump rates to bottom electrode\n");
		//Calculating jump rates for forward electrode
		initJumPossibility_ElecZ(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Zb1, Zb2, *elZb,\
			 	RelativePerm, vZ, SLength, SWidth, PeriodicX, PeriodicY,\
				XElecOn, YElecOn, 0);
		
		deleteMatrix(&Zb1);
		deleteMatrix(&Zb2);

		printf("Initializing jump rates to top electrode\n");
		//Calculating jump rates for forward electrode
		initJumPossibility_ElecZ(electricEnergyX, electricEnergyY, electricEnergyZ,\
			  SiteDistance, MarcusCoeff, KT,reOrgEnergy, Za1, Za2, *elZa,\
			 	RelativePerm, vZ, SLength, SWidth, PeriodicX, PeriodicY,\
				XElecOn, YElecOn, 1);
		
		deleteMatrix(&Za1);
		deleteMatrix(&Za2);
	}

  snAelec = (SNarray) getElectrode_AdjacentSites(*elXb);
	//printf("Length %d Width %d Height %d\n",getAlen(snAelec),getAwid(snAelec),getAhei(snAelec));

	

	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
int initJumPossibility_ElecX( const double electricEnergyX,\
											const double electricEnergyY, const double electricEnergyZ,\
											const double SiteDistance,const double MarcusCoeff,\
		  								const double KT,const double reOrgEnergy, matrix X1, matrix X2,\
											Electrode elX,\
											const double RelativePermittivity,const double vX,\
											const int SWidth, const int SHeight,\
											const int PeriodicY, const int PeriodicZ,\
											const int YElecOn, const int ZElecOn, const int BorF){


	if(elX==NULL){
		printf("ERROR elX is NULL in initJumPossibility_ElecX\n");
		exit(1);
	}
	//This function as it is is not equiped to deal with multiple electrodes when they are 
	//turned on 

	//BorF - Back or front electrode back - 0 front - 1
	//X1 - contains the energies of the sites next to the electrodes
	//X2 - contains the energies further into the system at position 1 or position Slength-2

	//This function is for initializing the electrodes on the x 
	//axis, two sitenode arrays are passed which contains the energies
	//of the electrodes(Fermi energy) and the energies of the sites
	//next to the electrode. A matrix is then created that stores
	//the hop rates to and from the electrode only in the x directio

	//Constants
	//Units of Coulombs
	const double q = 1.601E-19;	
	//Units of Farads/m = (C^2/J) * (1/m)
	const double epsilon0 = 8.8541878176E-12;

	//Used to capture the tunneling constant of the electrode [1/nm]
	double alpha;

	//Image force energy [eV]
	double EnergyImageForce;
	//Fermi Energy of the electrode [eV]
	double ElecFermi;
	//Barrier Height between the electrode and the site [eV]
	double W;
	//Total Energy Difference between site and site hopping too [eV]
	double Energy;
	//Energy of current site [eV]
	double SiteEnergy;
	//Energy of Neighbor site [eV]
	double SiteEnergyNeigh;

	double rate;
	double SumRate;
	double pval;

	int j, k, l;
	//one node have 6 hopping rate for 6 neighbor node
	double v[6];
	double sum;
	SiteNode sn;
	//int divisor;

	//Shockley barrier height [eV] this is the difference in energy between 
	//the electrode and the conduction band (LUMO) of the semiconductor
	//If it is positive electrons must overcome energy to move to the semiconductor

 	EnergyImageForce = pow(q,2)/(16*M_PI*epsilon0*RelativePermittivity*SiteDistance)*1/q;
	ElecFermi = getElectrode_FermiEnergy(elX);
	alpha = getElectrode_alpha(elX);
	//Hops from sites within the system should use marcus formalism
	//Hops from electrode should use miller and Abrahams theory

	//Create a SNarray for all the nodes next to the electrode
	//It will be attached to the electrode Datastructure
	SNarray snAX = newSNarray( 1, SWidth, SHeight);
	//Create a matrix that will store all the hops off the electrode
	//will also be attached to a the Electrode datastructure
	matrix mtxX = newMatrix(SWidth,SHeight);
	SumRate = 0;

	//printf("SWidth %d SHeight %d\n",SWidth,SHeight);
	//printf("According to snAX length %d width %d height %d\n",getAlen(snAX),getAwid(snAX),getAhei(snAX));

	for(j = 0; j < SWidth; j++){
		for(k = 0; k < SHeight; k++){

			//Dealing with hops from front and back of the plane 
			if ( BorF == 0 ){
				//assuming plane is infront of back electrode
				SiteEnergy = getE(X1, j+1, k+1);
				sn = getSN(snAX,0,j,k);
				setEnergy(sn,SiteEnergy);
				//printf("SiteEnergy %g\n",SiteEnergy);
				//Hop from Electrode to site forward hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W - electricEnergyX;
				//printf("Energy %g SiteDistance %g alpha %g KT %g value of vX %g\n",Energy,SiteDistance,alpha, KT,vX);
				rate = vX*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				//printf("value of rate %g\n",rate);
				setE(mtxX,j+1,k+1,rate);

				//Hop from site to neighboring site forward hop
				SiteEnergyNeigh = getE(X2,j+1,k+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W - electricEnergyX;
				v[1] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode backward hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W + electricEnergyX;
				v[0] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}else{
				//assuming plane is behind front electrode
				SiteEnergy = getE(X1, j+1, k+1);
				sn = getSN(snAX,0,j,k);
				setEnergy(sn,SiteEnergy);

				//Hop from Electrode to site backwards hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W + electricEnergyX;
				rate = vX*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				setE(mtxX,j+1,k+1,rate);

				//Hop from site to neighboring site backwards hop
				SiteEnergyNeigh = getE(X2,j+1,k+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W + electricEnergyX;
				v[0] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode forward hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W - electricEnergyX;
				v[1] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}

			//No need to include image force energy all these sites 
			//are along the same plane

			//Grabbing energy of site neighboring (0,j,k) to the left
			SiteEnergyNeigh = getE(X1, ((j-1+SWidth)%SWidth)+1, k+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyY;
			v[2] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (0,j,k) to the right
			SiteEnergyNeigh = getE(X1, ((j+1)%SWidth)+1,k+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyY;
			v[3] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (0,j,k) below
			SiteEnergyNeigh = getE(X1, j+1, ((k-1+SHeight)%SHeight)+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyZ;
			v[4] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (0,j,k) above
			SiteEnergyNeigh = getE(X1, j+1, ((k+1)%SHeight)+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyZ;
			v[5] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			if( j == 0 && PeriodicY==0 && YElecOn==0)
					v[2] = 0;

			if( j == SWidth-1 && PeriodicY==0 && YElecOn==0)
					v[3] = 0;

			if( k == 0 && PeriodicZ==0 && ZElecOn==0)
					v[4] = 0;

			if( k == SHeight-1 && PeriodicZ==0 && ZElecOn==0)
					v[5] = 0;

       sum=0;
			 for(l = 0; l < 6; l++) {
				 sum = sum + v[l];
			 }

			 setsum(sn, sum);
			
			 pval=0;
			 for(l=0; l< 6;l++){
				 setSN_p(sn,l,v[l]/sum+pval);
				 pval = getSN_p(sn,l);
			 }
			SumRate += rate;		 
		}
	}

	setElectrode_Sum(elX, SumRate);
	if(SumRate<0.1){
		//printf("j %d k %d\n",j,k);
		printf("MarcusCoeff %g hoppingRate %g\n",MarcusCoeff,hoppingRate(Energy,KT,reOrgEnergy));
		printf("WARNING SumRate way is very small may need to adjust electrode fermi levels %g\n",SumRate);
	}
	
	pval = 0;
	for(j = 0; j< SWidth;j++){
		for(k=0;k<SHeight;k++){
			pval = getE(mtxX,j+1,k+1)/SumRate + pval;
			setE(mtxX,j+1,k+1, pval);
		}
	}

	if(mtxX==NULL || snAX == NULL){
		printf("ERROR mtx or snA NULL in initJumPossibility_ElecX\n");
		exit(1);
	}
	printf("Attaching mtX and snAX to elX\n");
	printf("Printing Hop Rates matrix\n");
	printMatrix(mtxX);
	setElectrode_HopRates(elX, (void *) mtxX);
	setElectrode_AdjacentSites(elX, (void *) snAX);
	printf("End of Init ElectrodeX\n");
	return 0;
}

int initJumPossibility_ElecY( const double electricEnergyX,\
											const double electricEnergyY, const double electricEnergyZ,\
											const double SiteDistance,const double MarcusCoeff,\
		  								const double KT,const double reOrgEnergy, matrix Y1, matrix Y2,\
											Electrode elY,\
											const double RelativePermittivity, const double vY,\
											const int SLength, const int SHeight,\
											const int PeriodicX, const int PeriodicZ,\
											const int XElecOn, const int ZElecOn, const int LorR){

	//This function as it is is not equiped to deal with multiple electrodes when they are 
	//turned on 

	//LorR - Left or Right electrode left - 0 right - 1
	//Y1 - contains the energies of the sites next to the electrodes
	//Y2 - contains the energies further into the system at position 1 or position SWidth-2

	//This function is for initializing the electrodes on the y 
	//axis, two sitenode arrays are passed which contains the energies
	//of the electrodes(Fermi energy) and the energies of the sites
	//next to the electrode. A matrix is then created that stores
	//the hop rates to and from the electrode only in the y direction

	//Constants
	//Units of Coulombs
	const double q = 1.601E-19;	
	//Units of Farads/m = (C^2/J) * (1/m)
	const double epsilon0 = 8.8541878176E-12;

	//Used to capture the tunneling constant of the electrode [1/nm]
	double alpha;

	//Image force energy [eV]
	double EnergyImageForce;
	//Fermi Energy of the electrode [eV]
	double ElecFermi;
	//Barrier Height between the electrode and the site [eV]
	double W;
	//Total Energy Difference between site and site hopping too [eV]
	double Energy;
	//Energy of current site [eV]
	double SiteEnergy;
	//Energy of Neighbor site [eV]
	double SiteEnergyNeigh;

	double rate;
	double SumRate;
	double pval;

	int i, k, l;
	//one node have 6 hopping rate for 6 neighbor node
	double v[6];
	double sum;
	SiteNode sn;
	//int divisor;

	//Shockley barrier height [eV] this is the difference in energy between 
	//the electrode and the conduction band (LUMO) of the semiconductor
	//If it is positive electrons must overcome energy to move to the semiconductor

 	EnergyImageForce = pow(q,2)/(16*M_PI*epsilon0*RelativePermittivity*SiteDistance)*1/q;
	ElecFermi = getElectrode_FermiEnergy(elY);
	alpha = getElectrode_alpha(elY);
	//Hops from sites within the system should use marcus formalism
	//Hops from electrode should use miller and Abrahams theory

	//Create a SNarray for all the nodes next to the electrode
	//It will be attached to the electrode Datastructure
	SNarray snAY = newSNarray( SLength,1, SHeight);
	//Create a matrix that will store all the hops off the electrode
	//will also be attached to a the Electrode datastructure
	matrix mtxY = newMatrix(SLength,SHeight);
	SumRate = 0;

	for(i = 0; i < SLength; i++){
		for(k = 0; k < SHeight; k++){

			//Dealing with hops from left and right of the plane 
			if ( LorR == 0 ){
				//assuming plane is to the right of the of left electrode
				SiteEnergy = getE(Y1, i+1, k+1);
				sn = getSN(snAY,i,0,k);
				setEnergy(sn,SiteEnergy);

				//Hop from Electrode to site right hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W - electricEnergyY;
				rate = vY*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				setE(mtxY,i+1,k+1,rate);

				//Hop from site to neighboring site right hop
				SiteEnergyNeigh = getE(Y2,i+1,k+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W - electricEnergyY;
				v[3] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode backward hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W + electricEnergyY;
				v[2] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}else{
				//assuming plane is to the left of the right electrode
				SiteEnergy = getE(Y1, i+1, k+1);
				sn = getSN(snAY,i,0,k);
				setEnergy(sn,SiteEnergy);

				//Hop from Electrode to site left hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W + electricEnergyY;
				rate = vY*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				setE(mtxY,i+1,k+1,rate);

				//Hop from site to neighboring site left hop
				SiteEnergyNeigh = getE(Y2,i+1,k+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W + electricEnergyX;
				v[3] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode forward hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W - electricEnergyY;
				v[2] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}

			//No need to include image force energy all these sites 
			//are along the same plane

			//Grabbing energy of site neighboring (i,0,k) behind
			SiteEnergyNeigh = getE(Y1, ((i-1+SLength)%SLength)+1, k+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyX;
			v[0] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,0,k) infront
			SiteEnergyNeigh = getE(Y1, ((i+1)%SLength)+1,k+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyX;
			v[1] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,0,k) below
			SiteEnergyNeigh = getE(Y1, i+1, ((k-1+SHeight)%SHeight)+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyZ;
			v[4] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,0,k) above
			SiteEnergyNeigh = getE(Y1, i+1, ((k+1)%SHeight)+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyZ;
			v[5] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			if( i == 0 && PeriodicX==0 && XElecOn==0)
					v[0] = 0;

			if( i == SLength-1 && PeriodicX==0 && XElecOn==0)
					v[1] = 0;

			if( k == 0 && PeriodicZ==0 && ZElecOn==0)
					v[4] = 0;

			if( k == SHeight-1 && PeriodicZ==0 && ZElecOn==0)
					v[5] = 0;

       sum=0;
			 for(l = 0; l < 6; l++) {
				 sum = sum + v[l];
			 }

			 setsum(sn, sum);
			
			 pval=0;
			 for(l=0; l< 6;l++){
				 setSN_p(sn,l,v[l]/sum+pval);
				 pval = getSN_p(sn,l);
			 }
			 
			SumRate += rate;		 
		}
	}

	setElectrode_Sum(elY, SumRate);
	pval = 0;
	for(i = 0; i< SLength;i++){
		for(k=0;k<SHeight;k++){
			pval = getE(mtxY,i+1,k+1)/SumRate + pval;
			setE(mtxY,i+1,k+1, pval);
		}
	}
	setElectrode_HopRates(elY, (void *) mtxY);
	setElectrode_AdjacentSites(elY, (void *) snAY);
	printf("End of Init ElectrodeY\n");
	return 0;
}

int initJumPossibility_ElecZ( const double electricEnergyX,\
											const double electricEnergyY, const double electricEnergyZ,\
											const double SiteDistance,const double MarcusCoeff,\
		  								const double KT, const double  reOrgEnergy, matrix Z1, matrix Z2,\
											Electrode elZ,\
											const double RelativePermittivity, const double vZ,\
											const int SLength, const int SWidth,\
											const int PeriodicX, const int PeriodicY,\
											const int XElecOn, const int YElecOn, const int AorB){

	//This function as it is is not equiped to deal with multiple electrodes when they are 
	//turned on 

	//AorB - above or below electrode below - 0 above - 1
	//Z1 - contains the energies of the sites next to the electrodes
	//Z2 - contains the energies further into the system at position 1 or position SHeight-2

	//This function is for initializing the electrodes on the z 
	//axis, two sitenode arrays are passed which contains the energies
	//of the electrodes(Fermi energy) and the energies of the sites
	//next to the electrode. A matrix is then created that stores
	//the hop rates to and from the electrode only in the z direction

	//Constants
	//Units of Coulombs
	const double q = 1.601E-19;	
	//Units of Farads/m = (C^2/J) * (1/m)
	const double epsilon0 = 8.8541878176E-12;

	//Used to capture the tunneling constant of the electrode [1/nm]
	double alpha;

	//Image force energy [eV]
	double EnergyImageForce;
	//Fermi Energy of the electrode [eV]
	double ElecFermi;
	//Barrier Height between the electrode and the site [eV]
	double W;
	//Total Energy Difference between site and site hopping too [eV]
	double Energy;
	//Energy of current site [eV]
	double SiteEnergy;
	//Energy of Neighbor site [eV]
	double SiteEnergyNeigh;

	double rate;
	double SumRate;
	double pval;

	int i, j, l;
	//one node have 6 hopping rate for 6 neighbor node
	double v[6];
	double sum;
	SiteNode sn;
	//int divisor;

	//Shockley barrier height [eV] this is the difference in energy between 
	//the electrode and the conduction band (LUMO) of the semiconductor
	//If it is positive electrons must overcome energy to move to the semiconductor

 	EnergyImageForce = pow(q,2)/(16*M_PI*epsilon0*RelativePermittivity*SiteDistance)*1/q;
	ElecFermi = getElectrode_FermiEnergy(elZ);
	alpha = getElectrode_alpha(elZ);
	//Hops from sites within the system should use marcus formalism
	//Hops from electrode should use miller and Abrahams theory

	//Create a SNarray for all the nodes next to the electrode
	//It will be attached to the electrode Datastructure
	SNarray snAZ = newSNarray(SLength, SWidth, 1);
	//Create a matrix that will store all the hops off the electrode
	//will also be attached to a the Electrode datastructure
	matrix mtxZ = newMatrix(SLength,SWidth);
	SumRate = 0;

	for(i = 0; i < SLength; i++){
		for(j = 0; j < SWidth; j++){

			//Dealing with hops from top (above) or bottom of the plane 
			if ( AorB == 0 ){
				//assuming plane is above of bottom electrode
				SiteEnergy = getE(Z1, i+1, j+1);
				sn = getSN(snAZ,i,j,0);
				setEnergy(sn,SiteEnergy);

				//Hop from Electrode to site above hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W - electricEnergyZ;
				rate = vZ*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				setE(mtxZ,i+1,j+1,rate);

				//Hop from site to neighboring site above hop
				SiteEnergyNeigh = getE(Z2,i+1,j+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W - electricEnergyZ;
				v[5] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode bottom hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W + electricEnergyZ;
				v[4] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}else{
				//assuming plane is below top electrode
				SiteEnergy = getE(Z1, i+1, j+1);
				sn = getSN(snAZ,i,j,0);
				setEnergy(sn,SiteEnergy);

				//Hop from Electrode to site below hop
				W = (SiteEnergy-EnergyImageForce) - ElecFermi;
				Energy = W + electricEnergyZ;
				rate = vZ*hoppingRateMillerAbraham( Energy, SiteDistance, alpha, KT);
				setE(mtxZ,i+1,j+1,rate);

				//Hop from site to neighboring site below hop
				SiteEnergyNeigh = getE(Z2,i+1,j+1);
				W = SiteEnergyNeigh - (SiteEnergy-EnergyImageForce);
				Energy = W + electricEnergyZ;
				v[3] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

				//Hop from site to Electrode above hop
				W = ElecFermi - (SiteEnergy - EnergyImageForce);
				Energy = W - electricEnergyZ;
				v[4] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);

			}

			//No need to include image force energy all these sites 
			//are along the same plane

			//Grabbing energy of site neighboring (i,j,0) behind
			SiteEnergyNeigh = getE(Z1, ((i-1+SLength)%SLength)+1, j+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyX;
			v[0] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,j,0) infront
			SiteEnergyNeigh = getE(Z1, ((i+1)%SLength)+1,j+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyX;
			v[1] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,j,0) to the left
			SiteEnergyNeigh = getE(Z1,i+1, ((j-1+SWidth)%SWidth)+1);
			Energy = SiteEnergyNeigh - SiteEnergy + electricEnergyY;
			v[2] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			//Grabbing energy of site neighboring (i,j,0) to the right
			SiteEnergyNeigh = getE(Z1,i+1, ((j+1)%SWidth)+1);
			Energy = SiteEnergyNeigh - SiteEnergy - electricEnergyY;
			v[3] = MarcusCoeff*hoppingRate(Energy, KT, reOrgEnergy);	

			if( i == 0 && PeriodicX==0 && XElecOn==0)
					v[0] = 0;

			if( i == SLength-1 && PeriodicX==0 && XElecOn==0)
					v[1] = 0;

			if( j == 0 && PeriodicY==0 && YElecOn==0)
					v[2] = 0;

			if( j == SWidth-1 && PeriodicY==0 && YElecOn==0)
					v[3] = 0;

       sum=0;
			 for(l = 0; l < 6; l++) {
				 sum = sum + v[l];
			 }

			 setsum(sn, sum);
			
			 pval=0;
			 for(l=0; l< 6;l++){
				 setSN_p(sn,l,v[l]/sum+pval);
				 pval = getSN_p(sn,l);
			 }
			 
			SumRate += rate;		 
		}
	}

	setElectrode_Sum(elZ, SumRate);
	pval = 0;
	for(i = 0; i< SLength;i++){
		for(j=0;j<SWidth;j++){
			pval = getE(mtxZ,i+1,j+1)/SumRate + pval;
			setE(mtxZ,i+1,j+1, pval);
		}
	}
	setElectrode_HopRates(elZ, (void *) mtxZ);
	setElectrode_AdjacentSites(elZ, (void *) snAZ);
	printf("End of Init ElectrodeZ\n");
	return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////

int initJumPossibility(const double electricEnergyX,const double electricEnergyY,\
		  const double electricEnergyZ, const double MarcusCoeff,\
		  const double KT,const double reOrgEnergy, SNarray snA,\
			const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
			const int XElecOn, const int YElecOn, const int ZElecOn){

	if(snA==NULL || PeriodicX<0 || PeriodicX>1 ||\
									PeriodicY<0 || PeriodicY>1 ||\
									PeriodicZ<0 || PeriodicZ>1){
		return -1;
	}

	int SLength = getAlen(snA);
	int SWidth = getAwid(snA);
	int SHeight = getAhei(snA);
	int i, j, k, l;
	double pval;
	//one node have 6 hopping rate for 6 neighbor node
	double v[6];
	double sum;
	SiteNode SNi;
	SiteNode SNj;
	//int divisor;

	//printf("KT %g ReOrgEnergy %g electricEnergyX %g\n",KT,reOrgEnergy,electricEnergyX);
	//printf("SLength %d SWidth %d SHeight %d\n",SLength,SWidth,SHeight);
	for(i = 0; i < getAlen(snA); i++){
		for(j = 0; j < getAwid(snA); j++){
			for(k = 0; k < getAhei(snA); k++){

				//divisor=6;
				SNi = getSN(snA, i, j, k);
				SNj = getSN(snA, (i-1+SLength)%SLength,j,k);																						//Behind
				//if sigma = 0.09, generally the maximum difference Ej and Ei is approximate to 0.54
				v[0] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyX,KT,reOrgEnergy); 

				SNj = getSN(snA, (i+1)%SLength,j,k);																										//Front
				v[1] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyX,KT,reOrgEnergy); 
				SNj = getSN(snA,i, (j-1+SWidth)%SWidth,k);																							//Left
				v[2] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyY,KT,reOrgEnergy);
				SNj = getSN(snA,i, (j+1)%SWidth,k);																											//Right
				v[3] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyY,KT,reOrgEnergy);
				SNj = getSN(snA,i, j, (k-1+SHeight)%SHeight);																						//Below
				v[4] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyZ,KT,reOrgEnergy);
				SNj = getSN(snA,i, j, (k+1)%SHeight);																										//Above
				v[5] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyZ,KT,reOrgEnergy);

				//Will want to ignore hop rates at the edges
				//when non periodic and when the Electrodes are off
				
				if( i == 0 && PeriodicX==0 && XElecOn==0){
					v[0] = 0;
				}
				
				if( i == SLength-1 && PeriodicX==0 && XElecOn==0){
					v[1] = 0;
				}

				if( j == 0 && PeriodicY==0 && YElecOn==0){
					v[2] = 0;
				}

				if( j == SWidth-1 && PeriodicY==0 && YElecOn==0){
					v[3] = 0;
				}

				if( k == 0 && PeriodicZ==0 && ZElecOn==0){
					v[4] = 0;
				}

				if( k == SHeight-1 && PeriodicZ==0 && ZElecOn==0){
					v[5] = 0;
				}

				sum=0;
				for(l = 0; l < 6; l++) {
					sum = sum + v[l];
				}
		
				//Don't want the average but the sum
				setsum(SNi,sum);
				//This is WRONG: use the average probability to place charges in the queue--average 
				//rather than sum to avoid "trapping" charges at the interface
				pval=0;
				for(l = 0; l < 6; l++){
					setSN_p(SNi,l,v[l]/sum+pval);
					pval = getSN_p(SNi,l);
				}

			}

		}
	}
	printf("End of Init\n");
	return 0;
}

//////////////////////////////////////////////////////////////////

ChargeArray initCharget0( matrix Sequence, const_SNarray snA,const int Ntot,const int NCh,\
					const double D,const int XElecOn,const int YElecOn,const int ZElecOn, \
					const int EndX, const int EndY, const int EndZ){

	//D is the dimension to fit the real data
	//NCh is the total number of Charges injected per time step
	//Ntot is the total number of charges that will be injected


	if(snA==NULL || NCh>Ntot || D<0 || Sequence==NULL ||\
		 XElecOn<0 || YElecOn<0 || ZElecOn<0 ||\
		 XElecOn>1 || YElecOn>1 || ZElecOn>1 ){
		return NULL;
	}

	int loop;
	int i, j, k;
	int num1;
	int num2;
	int unOccYZ;
	int unOccXZ;
	int unOccXY;
	double ran;
	SiteNode site;
	ChargeArray chA = newChargeA(Ntot);
	Charge ch;

	printf("In initCharget0\n");

	if((XElecOn+YElecOn+ZElecOn)>1){
		printf("WARNING Only equipped to deal with one electrode turned on at a time!!!\n\n");
		return NULL;
	}

	if(XElecOn==1){
		unOccYZ = getUnOccYZplane(snA,0);	
		printf("Number of Unoccuped Sites on back Elec: %d\n",unOccYZ);

		if (unOccYZ<=NCh) {
			printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
		}
	}

	if(YElecOn==1){
		unOccXZ = getUnOccXZplane(snA,0);	
		printf("Number of Unoccuped Sites on back Elec: %d\n",unOccXZ);

		if (unOccXZ<=NCh) {
			printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
		}
	}
	
	if(ZElecOn==1){
		unOccXY = getUnOccXYplane(snA,0);	
		printf("Number of Unoccuped Sites on back Elec: %d\n",unOccXY);

		if (unOccXY<=NCh) {
			printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
		}
	}
	
	//initialization of charge[Ntot] and sequence[Ntot]
	//Basically initilizing all the charges
	for(loop = 0; loop < Ntot; loop++){
		//initializes the waiting time for each charge with a waiting 
			//time=(random nb from 0 to +infty)*average hopping rate of 
		//printf("Within loop\n");
		//printf("Value of loop %d\n",loop);
		//the site WHY NOT SMALLEST HOPPING RATE?
		setDwel(getCharge(chA,loop),1E6);
		//initializes the waiting queue
		//printf("Past set Dwel\n");
		setE( Sequence,loop+1,1, loop); 
		//printf("At setElement\n");
		ch = getCharge(chA,loop);
		//printf("After getCharge\n");
		if(ch==NULL){
			printf("Charge is null\n");
		}
		setCx(ch,0);
		setCy(ch,0);
		setCz(ch,0);
	}
	
	//printf("Length of Charge array %d\n",getChargeA_len(chA));
	//printf("Rows of Sequence matrix %d\n",getRows( Sequence ));
	//printf("After setting initial x y and z\n");
	//printf("Value of loop %d\n",loop);
	//printf("Value of Ntot %d\n",Ntot);
	//printf("Location of charge %d %d %d\n",getCx(ch), getCy(ch), getCz(ch));
	//}

	printf("Past initilizing Sequence\n");
	//Now setting the times for all the charges in the first timestep
	for(loop = 0; loop < NCh; loop++){
	
		//printf("In next loop\n");
		if(XElecOn==1){
			i = 0;
			do{
				//generate a number from 1 to Width
				num1 = rand() % getAwid(snA); 
				//generate a number from 1 to Height
				num2 = rand() % getAhei(snA);
				//Have to minus 1 because sitenodes are labeled
				//starting at 0 to (length,width,height)-1
				j = num1;
				k = num2;
				
			}while( getDwelStat(getSN(snA,i,j,k))!= -1);
		}else if(YElecOn==1){
			j = 0;
			do{
				//generate a number from 1 to Length
				num1 = rand() % getAlen(snA); 
				//generate a number from 1 to Height
				num2 = rand() % getAhei(snA);
				i = num1;
				k = num2;
			}while( getDwelStat(getSN(snA,i,j,k))!= -1);
		}else if(ZElecOn==1){
			k = 0;
			do{
				//generate a number from 1 to length
				num1 = rand() % getAlen(snA); 
				//generate a number from 1 to width
				num2 = rand() % getAwid(snA);
				i = num1;
				j = num2;
			}while( getDwelStat(getSN(snA,i,j,k))!= -1);
		}
		//Initilize information of site charge is occuping
		site = getSN(snA,i,j,k);
		setDwelStat(site,loop);
		setVisFreq(site,getVisFreq(getSN(snA,i,j,k))+1);
		setVis(site,1);
		//printf("After Vis\n");
		ch=getCharge(chA,loop);
		//Initialize the position of charge[loop]
		if(XElecOn==1){
			//X position will not change will start
			//at the electrode
			setCx(ch,i);

			if(EndY==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCy(ch,j);
			}else{
				num1 = rand()%EndY;
				setCy(ch,j+num1*getAwid(snA));
			}
			if(EndZ==0){
				setCz(ch,k);
			}else{
				num1 = rand()%EndZ;
				setCz(ch,k+num1*getAhei(snA));
			}
		}
		if(YElecOn==1){
			//Y position will not change will start
			//at the electrode
			setCy(ch,j);

			if(EndX==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCx(ch,i);
			}else{
				num1 = rand()%EndX;
				setCx(ch,i+num1*getAlen(snA));
			}
			if(EndZ==0){
				setCz(ch,k);
			}else{
				num1 = rand()%EndZ;
				setCz(ch,k+num1*getAhei(snA));
			}
		}
		if(ZElecOn==1){
			//Y position will not change will start
			//at the electrode
			setCz(ch,k);

			if(EndX==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCx(ch,i);
			}else{
				num1 = rand()%EndX;
				setCx(ch,i+num1*getAlen(snA));
			}
			if(EndY==0){
				setCy(ch,j);
			}else{
				num1 = rand()%EndY;
				setCy(ch,j+num1*getAwid(snA));
			}
		}
		//initializes the waiting time for each charge with a waiting 
		//time=(random nb from 0 to +infty)*average hopping rate of 
		do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
		setDwel(ch,-log(ran/RAND_MAX)/(getsum(site)));

		if isinf(getDwel(ch)) {
			printf("Value of getsum(site) %g\n",getsum(site));
			printf("ERROR Dweltime is infinite!\n");
			printf("You probably failed to call initsite. This means that getsum(site)\n");
			printf("will return a 0. -6.0*log(ran/RAND_MAX)/(D*0) ->inf\n");
		}

	}

	printf("Before quickSort in initcharet0\n");
	printf("Ntot is the total number of charges tracked %d \n",Ntot);

	//use quicksort algorithm to sort the disordered array sequence[]
	quickSort(0, Ntot - 1, Sequence , chA);
	printf("After quickSort in initcharget0\n");

	return chA;
}

////////////////////////////////////////////////////////////////////

int initCharge(int nca,long int n, ChargeArray *chA, matrix Sequence, SNarray  snA, \
			const int Ntot, const int NCh, const double D,\
			const int XElecOn, const int YElecOn, const int ZElecOn,\
			const int EndX, const int EndY, const int EndZ){

	if ((nca+NCh-1)>Ntot || (nca+NCh-1)<0 ) {
		printf("Value of Ntot %d value of (nca+NCh-1) %d value of nca %d value of NCh %d\n",\
						Ntot,(nca+NCh-1), nca, NCh);
		printf("ERROR nca+N must be less than Ntot and greater than -1!\n");
		exit(1);
		return -1;
	}

	if(snA==NULL || Sequence==NULL){
		return -1;
	}

	if((XElecOn+YElecOn+ZElecOn)>1){
		printf("WARNING not accounted for when more than 1 electrode is turned on at a time!!!\n\n");
		return -1;
	}

	//queueloops
	int loop;
	int i, j, k;
	int num1;
	int num2;
	int tem;
	double ran;
	SiteNode site;
	Charge ch;

	printf("Time step %ld\n",n);

	if((XElecOn+YElecOn+ZElecOn)>1){
		printf("WARNING Only equipped to deal with one electrode turned on at a time!!!\n");
	}
	//N is the number of charges injected per time step
	if (XElecOn==1){
		printf("Number of Unoccupied Sites on Left Elec: %d New Charges to insert %d\n",getUnOccYZplane(snA,0),NCh);
	}else if(YElecOn==1){
		printf("Number of Unoccupied Sites on Left Elec: %d New Charges to insert %d\n",getUnOccXZplane(snA,0),NCh);
	}else if(ZElecOn==1){
		printf("Number of Unoccupied Sites on Left Elec: %d New Charges to insert %d\n",getUnOccXYplane(snA,0),NCh);
	}

	//printf("initCharge %d %d \n",(n-1)*NCh,n*NCh);
	//exit(1);
	for(loop = (n-1)*(long int)NCh; loop < n*(long int)NCh; loop++){
		
		if (XElecOn==1){
			i = 0;
			if (getUnOccYZplane(snA, 0)<NCh) {
				printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
			}
			tem=0;
			while (tem!=-1) {
				//generate a number from 0 to Width*Height-1
				num1 = rand() % getAwid(snA);
				num2 = rand() % getAhei(snA);
				j = num1;
				k = num2;
				tem = getDwelStat(getSN(snA,i,j,k));
			}
		}else if(YElecOn==1){
			j = 0;
			if (getUnOccXZplane(snA, 0)<NCh) {
				printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
			}
			tem=0;
			while (tem!=-1) {
				//generate a number from 0 to Width*Height-1
				num1 = rand() % getAlen(snA);
				num2 = rand() % getAhei(snA);
				i = num1;
				k = num2;
				tem = getDwelStat(getSN(snA,i,j,k));
			}
		}else if(ZElecOn==1){
			k = 0;
			if (getUnOccXYplane(snA, 0)<NCh) {
				printf("Warning: if there are not enough sites open will become stuck in endless loop\n");
			}
			tem=0;
			while (tem!=-1) {
				//generate a number from 0 to Width*Height-1
				num1 = rand() % getAlen(snA);
				num2 = rand() % getAwid(snA);
				i = num1;
				j = num2;
				tem = getDwelStat(getSN(snA,i,j,k));
			}
		}

		site=getSN(snA,i,j,k);	
		//if(site==NULL){
		//	printf("Zeros initilization\n");
		//	exit(1);
		//}
		setDwelStat(site,loop);
		setVisFreq(site,getVisFreq(getSN(snA,i,j,k))+1);
		setVis(site,1);

		//Initialize the position of charge[loop]
		ch=getCharge(*chA,loop);
		//Initialize the position of charge[loop]
		if(XElecOn==1){
			//X position will not change will start
			//at the electrode
			setCx(ch,i);

			if(EndY==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCy(ch,j);
			}else{
				num1 = rand()%EndY;
				setCy(ch,j+num1*getAwid(snA));
			}
			if(EndZ==0){
				setCz(ch,k);
			}else{
				num1 = rand()%EndZ;
				setCz(ch,k+num1*getAhei(snA));
			}
		}
		if(YElecOn==1){
			//Y position will not change will start
			//at the electrode
			setCy(ch,j);

			if(EndX==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCx(ch,i);
			}else{
				num1 = rand()%EndX;
				setCx(ch,i+num1*getAlen(snA));
			}
			if(EndZ==0){
				setCz(ch,k);
			}else{
				num1 = rand()%EndZ;
				setCz(ch,k+num1*getAhei(snA));
			}
		}
		if(ZElecOn==1){
			//Y position will not change will start
			//at the electrode
			setCz(ch,k);

			if(EndX==0){
				//If the charge is periodic in y there is
				//no reason to spread charges out across
				//sample because there are no walls
				setCx(ch,i);
			}else{
				num1 = rand()%EndX;
				setCx(ch,i+num1*getAlen(snA));
			}
			if(EndY==0){
				setCy(ch,j);
			}else{
				num1 = rand()%EndY;
				setCy(ch,j+num1*getAwid(snA));
			}
		}

		//Initializes the waiting time for each charge with a waiting 
		//time=(random nb from 0 to +infty)*average hopping rate of 
		do{ran = rand();}while(ran == 0 || ran == RAND_MAX);
		setDwel(ch,-log(ran/RAND_MAX)/(getsum(site)));

		if isinf(getDwel(ch)) {
			printf("Value of setsum %g\n",getsum(site));
			printf("ERROR Dweltime is infinite!\n");
			printf("You probably failed to call initsite. This means that getsum(site)\n");
			printf("will return a 0. -6.0*log(ran/RAND_MAX)/(D*0) ->inf\n");
		}

	}

	printf("About to start quicksort\n");
	//use quicksort algorithm to sort the disordered array sequence[]
	quickSort(0, nca+NCh - 1, Sequence, *chA);
	printf("Past Quicksort\n");
	return 0;
}

///////////////////////////////////////////////////////////////////
void quickSort(int low, int high, matrix  Sequence, ChargeArray chA){
	int pivotloc;
	if(low < high){
		pivotloc = partition(low, high, Sequence, chA);
		quickSort(low, pivotloc -1, Sequence, chA);
		quickSort(pivotloc + 1, high, Sequence, chA);
	}
}

////////////////////////////////////////////////////////////////////
int partition(int low, int high, matrix  Sequence, ChargeArray chA){
	int num; //ArrayS;
	double pivotkey;

	num = getE(Sequence,low+1,1);
	pivotkey = getDwel(getCharge(chA,num));

	while(low < high){
		while(low < high && getDwel(getCharge(chA,getE(Sequence,high+1,1))) >= pivotkey)
			high--;

		setE(Sequence,low+1,1,getE(Sequence,high+1,1));
		while(low < high && getDwel(getCharge(chA,getE(Sequence,low+1,1))) <= pivotkey)
			low++;
		
		setE(Sequence,high+1,1,getE(Sequence,low+1,1));
	}

	setE(Sequence,low+1,1,num);
	return low;
}
