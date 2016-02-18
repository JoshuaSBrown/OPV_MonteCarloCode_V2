#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "../CHARGE/charge.h"
#include "../MEM/mem.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERSITENODE/clustersitenode.h"
#include "../FUNCTIONS/functions.h"
#include "chargetransport.h"

int main(void){

	srand(time(NULL));

	ParameterFrame PF;
	SiteNode sn;
	int future;
	int rv;

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

	//Setting up PF
	PF = newParamFrame();
	PFset_reOrg(PF, reOrgEnergy);
	PFset_alphaxf(PF, alphaxf);
	PFset_alphaxb(PF, alphaxb);
	PFset_alphayr(PF, alphayr);
	PFset_alphayl(PF, alphayl);
	PFset_alphazb(PF, alphazb);
	PFset_alphaza(PF, alphaza);
	PFset_XFermiF(PF, FermiExf);
	PFset_XFermiB(PF, FermiExb);
	PFset_YFermiR(PF, FermiEyr);
	PFset_YFermiL(PF, FermiEyl);
	PFset_ZFermiA(PF, FermiEza);
	PFset_ZFermiB(PF, FermiEzb);
	PFset_gamma(PF, gamma);
	PFset_vX(PF,vX);
	PFset_vY(PF,vY);
	PFset_vZ(PF,vZ);
	PFset_Ntot(PF,Ntot);
	PFset_RelativePerm(PF,RelativePerm);
	PFset_SiteDist(PF,SiteDistance);
	PFset_AttemptToHop(PF,AttemptToHop);
	//Calculating Marcus J0 coefficient assuming the Attempt to hop Rate
	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);
	//Calculating full Marcus Coefficient;
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(-2*gamma*SiteDistanceNM);


	//Creating Charge array
	ChargeArray chA = newChargeA(nc1);
	//Creating Sequence
	matrix Sequence1 = newMatrix(nc1,1);
	Charge ch;
	//Creating Electrode
	elXf = newElectrode();
	setElectrode_alpha(elXf,gamma);
	int i;
	//The Charges in chA have ids
	//starting at 0 and going to nc-1
	for(i=0;i<nc1;i++){
		ch = getCharge(chA,i);
		setDwel(ch,i*2+1);
		setE(Sequence1,i+1,1,i);
	}

	int ChargeID = getE(Sequence1,1,1);	
	ch = getCharge(chA, ChargeID);

	ChargeClosure(elXf, &ch, &Sequence1, &nc1, &nca1, &TotalCollected, Ntot);

	printChargeA(chA);
	printMatrix(Sequence1);
	printf("Charges Left %d Charges Collected %d Total Charges %d\n",nc1,TotalCollected, Ntot);

	//Testing InsertDwelltimePos

	//Creating Sequence
	int nc = Ntot;
	matrix Sequence = newMatrix(nc,1);

	for(i=0;i<nc;i++){
		ch = getCharge(chA,i);
		setDwel(ch,i*2+1);
		setE(Sequence,i+1,1,i);
	}

	//Grabbing first charge in sequence
	ch = getCharge(chA,0);
	//setting dwelltime somwhere in between
	setDwel(ch,22);

	printf("Before insertDwell\n");
	printChargeA(chA);
	printMatrix(Sequence);

	insertDwelltimePos( nc, chA, &Sequence); 

	printf("After insertDwell function with Charge 0 changed\n");
	printChargeA(chA);
	printMatrix(Sequence);

	//Charge 0 should appear in the sequence between 10 and 11;

	//Grab the charge at the beginning of the sequence set it so it
	//should reappear at the end
	ChargeID = getE(Sequence,1,1);
	ch = getCharge(chA,ChargeID);
	setDwel(ch,50);
	insertDwelltimePos( nc, chA, &Sequence); 

	printf("After insertDwell function with Charge 1 changed\n");
	printf("Charge 1 should appear at the end\n");
	printChargeA(chA);
	printMatrix(Sequence);

	//Grabbing first charge in sequence should be charge 2
	ChargeID = getE(Sequence,1,1);
	ch = getCharge(chA, ChargeID);
	setDwel(ch,0.1);
	insertDwelltimePos(nc,chA, &Sequence);

	printf("After insertDwell function with Charge 2 changed\n");
	printf("Charge 2 should remain at the start\n");
	printChargeA(chA);
	printMatrix(Sequence);

	//Final test setting Charge 2 dwelltime equal to 13 which is
	//the same as another charge
	ChargeID = getE(Sequence,1,1);
	ch = getCharge(chA, ChargeID);
	setDwel(ch,13);
	insertDwelltimePos(nc,chA, &Sequence);
	printf("ChargeID %d\n",ChargeID);
	printf("After insertDwell function with Charge 2 changed\n");
	printf("Charge 2 should appear next to charge 6\n");
	//printChargeA(chA);
	printMatrix(Sequence);

	//Testing ClusterHop & MakeHop
	
	deleteMatrix(&Sequence1);
	deleteMatrix(&Sequence);
	deleteElectrode(&elXf);
	deleteChargeA(chA);
	///////////////////////////////////////////////////////////////////
	//Must first create a cluster
	//Parameters
	int PeriodicX2=0;
	int PeriodicY2=0;
	int PeriodicZ2=0;
	int attempts=15;

	int XElecOn2 = 1;
	int YElecOn2 = 0;
	int ZElecOn2 = 0;

	int TotalX2=0;
	int TotalY2=0;
	int TotalZ2=0;

	electricEnergyX = 0;
	electricEnergyY = 0;
	electricEnergyZ = 0;
	
	PFset_Px(PF,PeriodicX2);
	PFset_Py(PF,PeriodicY2);
	PFset_Pz(PF,PeriodicZ2);
	PFset_XElecOn(PF,XElecOn2);
	PFset_YElecOn(PF,YElecOn2);
	PFset_ZElecOn(PF,ZElecOn2);
	//variables
	int OrderL2;
	int OrderH2;
	int MidPtsTotal2;
	int TotalOrders2;
	double time;
	int flag;

	time = 0;
	future = 0;

	SNarray snA2 = newSNarray(5,5,5);
	//Creating trap sites for cluster
	SiteNode sn2 = getSN(snA2,1,1,1);
	setEnergy(sn2, -4);
	sn2 = getSN(snA2,2,1,1);
	setEnergy(sn2, -4);

	//Calculate MasterM
	matrix MasterM2 = CalculateAllHops(snA2, electricEnergyX, electricEnergyY, electricEnergyZ,\
			1, 1, SiteDistance, AttemptToHop, gamma,\
			PeriodicX2, PeriodicY2, PeriodicZ2);
	ArbArray mpA2 = MPsort( &OrderL2, &OrderH2, &MidPtsTotal2, MasterM2, snA2, \
			PeriodicX2, PeriodicY2, PeriodicZ2);
	TotalOrders2 = OrderH2-OrderL2 + 1;

	ArbArray Arb2 = SortOrderMag(TotalOrders2, OrderL2, mpA2);
	printArbArray(Arb2,OrderL2);
	ArbArray ArClLL2 = ClusterSort( TotalOrders2, OrderL2, Arb2);
	//printArbArray(ArClLL2,OrderL2);
	FilterCluster( TotalOrders2, OrderL2, MasterM2, &ArClLL2, snA2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,\
			XElecOn2, YElecOn2, ZElecOn2);
	CalculateNeighNodes(TotalOrders2, OrderL2, &ArClLL2, snA2, PeriodicX2, PeriodicY2, PeriodicZ2);
	CalculateSumAndP(TotalOrders2, snA2, &ArClLL2, MasterM2, attempts, PeriodicX2, PeriodicY2, PeriodicZ2);
	ConnectClusterSN(TotalOrders2, snA2, ArClLL2);

	//printArbArray(ArClLL2,OrderL2);
	//SiteNodes have been connected with the cluster
	//Now need to place a charge on the cluster
	ch = newCharge();

	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);
	//Make sure the site reads as occupied
	sn2 = getSN(snA2,1,1,1);
	setDwelStat(sn2,1);

	//Hopping to neighbor and within cluster is allowed
	ClusterHop(snA2, &ch, &time, &future);
	flag = MakeHop(snA2, future, &ch, &TotalX2, &TotalY2, &TotalZ2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,\
			XElecOn2, YElecOn2, ZElecOn2);

	printf("Value of flag should be 0 which means the charge hopped\n");
	printf("Value of flag: %d\n",flag);
	printf("Location of Charge: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge was %d\n",getDwelStat(sn2));
	/*sn2 = getSN(snA2,getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge moved to %d\n",getDwelStat(sn2));
	printf("Time of hop %g\n",time);
	printf("Total Distance moved %d, %d, %d\n",TotalX2,TotalY2,TotalZ2);
	*/
	printf("Deleting ch\n");
	deleteCharge(ch);
	printf("Deleting MasterM2\n");
	deleteMatrix(&MasterM2);
	printf("Deleting ArClLL2\n");
	deleteArbArray(&ArClLL2);	
	printf("Deleting Arb2\n");
	deleteArbArray(&Arb2);	
	printf("Deleting snA2\n");
	deleteSNarray(snA2);
	printf("Deleting mpA2\n");
	deleteAllMidPointArray(&mpA2);
	
	///////////////////////////////////////////////////////////////////
	//Study 2
	//hopping within cluster is not allowed. 
	snA2 = newSNarray(5,5,5);
	//Creating trap sites for cluster
	sn2 = getSN(snA2,1,1,1);
	setEnergy(sn2, -4);
	sn2 = getSN(snA2,2,1,1);
	setEnergy(sn2, 2);
	sn2 = getSN(snA2,1,1,1);
	setDwelStat(sn2,1);

	//Placing charge back on site 1,1,1
	ch = newCharge();
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	TotalX2=0;
	TotalY2=0;
	TotalZ2=0;
	///Calculate MasterM
	MasterM2 = CalculateAllHops(snA2, electricEnergyX, electricEnergyY, electricEnergyZ,\
			1, 1, SiteDistance, AttemptToHop, gamma,\
			PeriodicX2, PeriodicY2, PeriodicZ2);
	mpA2 = MPsort( &OrderL2, &OrderH2, &MidPtsTotal2, MasterM2, snA2, \
			PeriodicX2, PeriodicY2, PeriodicZ2);
	TotalOrders2 = OrderH2-OrderL2 + 1;

	Arb2 = SortOrderMag(TotalOrders2, OrderL2, mpA2);
	printArbArray(Arb2,OrderL2);
	ArClLL2 = ClusterSort( TotalOrders2, OrderL2, Arb2);
	//printArbArray(ArClLL2,OrderL2);
	FilterCluster( TotalOrders2, OrderL2, MasterM2, &ArClLL2, snA2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,\
			XElecOn2, YElecOn2, ZElecOn2);
	CalculateNeighNodes(TotalOrders2, OrderL2, &ArClLL2, snA2, PeriodicX2, PeriodicY2, PeriodicZ2);
	CalculateSumAndP(TotalOrders2, snA2, &ArClLL2, MasterM2, attempts, PeriodicX2, PeriodicY2, PeriodicZ2);
	ConnectClusterSN(TotalOrders2, snA2, ArClLL2);

	//Other site in the cluster is now occupied
	printf("\nLocation of Charge before hop: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	ClusterHop(snA2, &ch, &time, &future);
	flag = MakeHop(snA2, future, &ch, &TotalX2, &TotalY2, &TotalZ2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,\
			XElecOn2, YElecOn2, ZElecOn2);

	printf("Value of flag should be 0 which means the charge hopped\n");
	printf("Value of flag: %d\n",flag);
	printf("Location of Charge: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge was %d\n",getDwelStat(sn2));
	sn2 = getSN(snA2,getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge moved to %d\n",getDwelStat(sn2));
	printf("Time of hop %g\n",time);
	printf("Total Distance moved %d, %d, %d\n",TotalX2,TotalY2,TotalZ2);

	deleteCharge(ch);
	deleteMatrix(&MasterM2);
	deleteArbArray(&ArClLL2);	
	deleteArbArray(&Arb2);	
	deleteSNarray(snA2);
	deleteAllMidPointArray(&mpA2);

	/*	
	///////////////////////////////////////////////////////////////////
	//Study 3
	//hopping within cluster is allowed.
	//hopping off cluster is not allowed.
	setDwelStat(sn3,-1);

	sn3 = getSN(snA2,3,1,1);
	setDwelStat(sn3,2);
	sn3 = getSN(snA2,2,2,1);
	setDwelStat(sn3,3);
	sn3 = getSN(snA2,2,0,1);
	setDwelStat(sn3,4);
	sn3 = getSN(snA2,2,1,2);
	setDwelStat(sn3,5);
	sn3 = getSN(snA2,2,1,0);
	setDwelStat(sn3,6);

	sn3 = getSN(snA2,0,1,1);
	setDwelStat(sn3,7);
	sn3 = getSN(snA2,1,2,1);
	setDwelStat(sn3,8);
	sn3 = getSN(snA2,1,0,1);
	setDwelStat(sn3,9);
	sn3 = getSN(snA2,1,1,2);
	setDwelStat(sn3,10);
	sn3 = getSN(snA2,1,1,0);
	setDwelStat(sn3,11);

	//Placing charge back on site 1,1,1
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	TotalX2=0;
	TotalY2=0;
	TotalZ2=0;
	//Testing	
	printf("\nLocation of Charge before hop: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	ClusterHop(snA2, &ch, &time, &future);
	flag = MakeHop(snA2, future, &ch, &TotalX2, &TotalY2, &TotalZ2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,\
			XElecOn2, YElecOn2, ZElecOn2);
	printf("Value of flag should be 0 which means the charge hopped\n");
	printf("Value of flag: %d\n",flag);
	printf("Location of Charge: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge was %d\n",getDwelStat(sn2));
	sn2 = getSN(snA2,getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge moved to %d\n",getDwelStat(sn2));
	printf("Time of hop %g\n",time);
	printf("Total Distance moved %d, %d, %d\n",TotalX2,TotalY2,TotalZ2);

	//Study 4
	//Now checking when we can't hop within or off the cluster
	sn3 = getSN(snA2,2,1,1);
	setDwelStat(sn3,12);

	//Placing charge back on site 1,1,1
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	TotalX2=0;
	TotalY2=0;
	TotalZ2=0;
	printf("\nLocation of Charge before hop: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	ClusterHop(snA2, &ch, &time, &future);
	flag = MakeHop(snA2, future, &ch, &TotalX2, &TotalY2, &TotalZ2,\
			PeriodicX2, PeriodicY2, PeriodicZ2,
			XElecOn2, YElecOn2, ZElecOn2);
	printf("Value of flag should be 0 which means the charge hopped\n");
	printf("Value of flag: %d\n",flag);
	printf("Location of Charge: %d, %d, %d \n",getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge was %d\n",getDwelStat(sn2));
	sn2 = getSN(snA2,getCx(ch),getCy(ch),getCz(ch));
	printf("DwelStat of site where charge moved to %d\n",getDwelStat(sn2));
	printf("Time of hop %g\n",time);
	printf("Total Distance moved %d, %d, %d\n",TotalX2,TotalY2,TotalZ2);

	//CASE STUDY 1 Testing when placed in center of Sample and made favorable to jump in positive X direction
	//Parameters
	printf("**************Case Study 1******************\n");
	deleteCharge(ch);
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	int totalX = 0;
	int totalY = 0;
	int totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	//Variables
	//int OrderL;
	//int OrderH;
	//int MidPtsTotal;
	//int TotalOrders;

	//Need to correctly set up SN array
	SNarray snA = newSNarray(SLength, SWidth, SHeight);

	//Lowering the Energy of a site right infront of the
	//charge so that it will have the highest probability of 
	//hopping to this location
	SiteNode site= getSN(snA, 2,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias in x, y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn2, YElecOn2, ZElecOn2);
	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 2,1,1);
	printf("\nMost favorable site infront of charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 2 Testing when placed in center of Sample and made favorable to jump in - X direction
	//Parameters
	printf("\n**************Case Study 2******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site right behind the charge
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 0,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias for x, y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn2, YElecOn2, ZElecOn2);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 0,1,1);
	printf("\nMost favorable site behind charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 3 Testing when placed in center of Sample and made favorable to jump in + Y direction
	//Parameters
	printf("\n**************Case Study 3******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site to the right of the charge
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,2,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias for x, y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn2, YElecOn2, ZElecOn2);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,2,1);
	printf("\nMost favorable site right of charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 4 Testing when placed in center of Sample and made favorable to jump in - Y direction
	//Parameters
	printf("\n**************Case Study 4******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site to the left of the charge
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,0,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias for x, y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn2, YElecOn2, ZElecOn2);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,0,1);
	printf("\nMost favorable site to the left of charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 5 Testing when placed in center of Sample and made favorable to jump in +Z direction
	//Parameters
	printf("\n**************Case Study 5******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site above the charge
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,2);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias for x,y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,1,2);
	printf("\nMost favorable site above charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 6 Testing when placed in center of Sample and made favorable to jump in -Z direction
	//Parameters
	printf("\n**************Case Study 6******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site right below the charge
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,0);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias for x, y and z
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the center 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,1,0);
	printf("\nMost favorable site below charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 7 Testing when charge placed infront of x+ electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 7******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site right infront of the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 0,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the front electrode 
	setCx(ch,3);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 3,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 0,1,1);
	printf("\nMost favorable site infront charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 8 Testing when charge placed infront of x- electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 8******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site right behind the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 3,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the back electrode 
	setCx(ch,0);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 0,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 3,1,1);
	printf("\nMost favorable site behind charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 9 Testing when charge placed infront of y+ electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 9******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site to the right of the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,0,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the right electrode 
	setCx(ch,1);
	setCy(ch,3);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,3,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,0,1);
	printf("\nMost favorable site to the right of the charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 10 Testing when charge placed infront of y- electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 10******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site to the right of the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,3,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the Left electrode 
	setCx(ch,1);
	setCy(ch,0);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,0,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,3,1);
	printf("\nMost favorable site to the left of the charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 11 Testing when charge placed infront of z+ electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 11******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site above the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,0);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the top electrode 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,3);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,3);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,1,0);
	printf("\nMost favorable site to the above the charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 12 Testing when charge placed infront of z- electrode and
	//hop to electrode is favorable
	printf("\n**************Case Study 12******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,3);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the bottom electrode 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,0);

	printf("\nCharge Before\n");
	printCharge(ch);
	printf("\nSite Charge should jump to\n");
	printSN(site);
	site= getSN(snA, 1,1,0);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);
	site= getSN(snA, 1,1,0);
	printf("\nMost favorable site to the above the charge\n");
	printSN(site);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 13 Testing when charge placed next to a side without electrode
	//in the x+ direction
	printf("\n**************Case Study 13******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 0,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the front side where there is no electrode
	setCx(ch,3);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 3,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 14 Testing when charge placed next to a side without electrode
	//in the x- direction
	printf("\n**************Case Study 14******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 3,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the back side where there is no electrode 
	setCx(ch,0);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 0,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 15 Testing when charge placed next to a side without electrode
	//in the y+ direction
	printf("\n**************Case Study 15******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site = getSN(snA, 1,0,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the right side where there is no electrode
	setCx(ch,1);
	setCy(ch,3);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,3,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 16 Testing when charge placed next to a side without electrode
	//in the y- direction
	printf("\n**************Case Study 16******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,3,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the left side where there is no electrode
	setCx(ch,1);
	setCy(ch,0);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,0,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 17 Testing when charge placed next to a side without electrode
	//in the z+ direction
	printf("\n**************Case Study 17******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,0);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the top side where there is no electrode 
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,3);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,1,3);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 18 Testing when charge placed next to a side without electrode
	//in the z- direction
	printf("\n**************Case Study 18******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);
	//Lowering the Energy of a site below the charge
	//Because we have included the Xelectrodes this will be 
	//equivalent to the site on the other side of the sample
	//charge so that it will have the highest probability of 
	//hopping to this location
	site= getSN(snA, 1,1,3);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near the bottom side where there is no electrode
	setCx(ch,1);
	setCy(ch,1);
	setCz(ch,0);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,1,0);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 19 Testing when charge placed next to two sides without electrodes
	//placing on xy corner
	printf("\n**************Case Study 19******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x+y+ corner
	setCx(ch,3);
	setCy(ch,3);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 3,3,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 20 Testing when charge placed next to two sides without electrodes
	//placing on xy corner
	printf("\n**************Case Study 20******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x-y- corner
	setCx(ch,0);
	setCy(ch,0);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 0,0,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 21 Testing when charge placed next to two sides without electrodes
	//placing on xy corner
	printf("\n**************Case Study 21******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x+y- corner
	setCx(ch,3);
	setCy(ch,0);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 3,0,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 22 Testing when charge placed next to two sides without electrodes
	//placing on xy corner
	printf("\n**************Case Study 22******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x+y- corner
	setCx(ch,0);
	setCy(ch,3);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 0,3,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 23 Testing when charge placed next to two sides without electrodes
	//placing on yz corner
	printf("\n**************Case Study 23******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near y+z+ corner
	setCx(ch,1);
	setCy(ch,3);
	setCz(ch,3);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,3,3);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 24 Testing when charge placed next to two sides without electrodes
	//placing on yz corner
	printf("\n**************Case Study 24******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near y-z- corner
	setCx(ch,1);
	setCy(ch,0);
	setCz(ch,0);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,0,0);
	printf("\nSite Charge is located on\n");
	printSN(site);

	flag = SiteHop(snA, &ch, site,&future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 25 Testing when charge placed next to two sides without electrodes
	//placing on yz corner
	printf("\n**************Case Study 25******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near y+z- corner
	setCx(ch,1);
	setCy(ch,3);
	setCz(ch,0);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,3,0);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 26 Testing when charge placed next to two sides without electrodes
	//placing on yz corner
	printf("\n**************Case Study 26******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near y-z+ corner
	setCx(ch,1);
	setCy(ch,0);
	setCz(ch,3);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,0,3);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 27 Testing when charge placed next to two sides without electrodes
	//placing on xz corner
	printf("\n**************Case Study 27******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x+z+ corner
	setCx(ch,3);
	setCy(ch,1);
	setCz(ch,3);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 3,1,3);
	printf("\nSite Charge is located on\n");
	printSN(site);

	flag = SiteHop(snA, &ch, site,&future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 28 Testing when charge placed next to two sides without electrodes
	//placing on xz corner
	printf("\n**************Case Study 28******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 1;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	//Placing charge near x+z+ corner
	setCx(ch,0);
	setCy(ch,1);
	setCz(ch,0);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 0,1,0);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 29 Testing Periodic Conditions in x when sample 
	//is doubled
	printf("\n**************Case Study 29******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 2;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 1;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	site= getSN(snA, 0,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	setCx(ch,7);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 3,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 30 Testing Periodic Conditions in x when sample 
	//is doubled
	printf("\n**************Case Study 30******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 2;
	EndY = 1;
	EndZ = 1;
	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;
	PeriodicX = 1;
	PeriodicY = 0;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	site= getSN(snA, 3,1,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	setCx(ch,4);
	setCy(ch,1);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 0,1,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 31 Testing Periodic Conditions in y when sample 
	//is doubled
	printf("\n**************Case Study 31******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 2;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 1;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	site= getSN(snA, 1,0,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	setCx(ch,1);
	setCy(ch,3);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,3,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	//CASE STUDY 32 Testing Periodic Conditions in y when sample 
	//is doubled
	printf("\n**************Case Study 32******************\n");
	SLength = 4;
	SWidth = 4;
	SHeight = 4;
	totalX = 0;
	totalY = 0;
	totalZ = 0;
	EndX = 1;
	EndY = 2;
	EndZ = 1;
	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;
	PeriodicX = 0;
	PeriodicY = 1;
	PeriodicZ = 0;

	snA = newSNarray(SLength, SWidth, SHeight);

	site= getSN(snA, 1,3,1);
	setEnergy(site, -15);

	//Now we need to calculate pvals for array
	//0 - value of bias
	//1 - value of KT
	//5 - value of reOrgEnergy
	initJumPossibility( 0,0,0, MarcusCoeff,\
			1, 15, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	ch = newCharge();

	setCx(ch,1);
	setCy(ch,4);
	setCz(ch,1);

	printf("\nCharge Before\n");
	printCharge(ch);
	site= getSN(snA, 1,0,1);
	printf("\nSite Charge is located on\n");
	printSN(site);

	SiteHop(snA, &ch, site, &future,\
			EndX, EndY, EndZ,\
			XElecOn, YElecOn, ZElecOn,\
			PeriodicX, PeriodicY, PeriodicZ);
	flag = MakeHop(snA, future, &ch, &totalX, &totalY, &totalZ,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("\nFlag %d\n",flag); 
	printf("\nCharge After Hop\n");
	printCharge(ch);

	deleteCharge(ch);
	deleteSNarray(snA);

	printf("Testing RandomWalk");

	//For non-periodic uniform sample 
	//With electrodes along the x-axis
	//And no electric field
	//Single Charge
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA3 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 50;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 0;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	char fid[] = "RandomWalk1\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 1.0, snA3,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA3);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA3, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA3, CheckPointNum, &fid[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printf("Printing Visit Freq\n");

	printVisitFreq(snA3, &fid[0]);

	printf("Deleting Electrodes\n");

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}

	printf("Deleting SiteNode array\n");
	deleteSNarray(snA3);
		 
	printf("Testing RandomWalk2\n");
	//For non-periodic uniform sample 
	//With electrodes along the x-axis
	//And no electric field
	//Single Charge large electric field
	//in the X
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA4 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 14;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	char fid2[] = "RandomWalk2\0";

	printf("Beginning initjumPossibility\n");
	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
	1.0 , 15.0, snA4,\
	PeriodicX, PeriodicY, PeriodicZ,\
	XElecOn, YElecOn, ZElecOn);

	printSNarray(snA4);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA4, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	printf("Beginning randomWalk\n");
	rv = randomWalk(snA4, CheckPointNum, &fid2[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA4, &fid2[0]);

	deleteSNarray(snA4);
	
	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}

	printf("Testing RandomWalk3");
	//Periodic in X only uniform sample
	//With electrodes along the x-axis
	//Single Charge large electric field
	//in the X
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA5 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 14;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 2;
	EndY = 1;
	EndZ = 1;

	PeriodicX=1;
	PeriodicY=0;
	PeriodicZ=0;

	char fid3[] = "RandomWalk3\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA5,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA5);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA5, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	rv = randomWalk(snA5, CheckPointNumm, &fid3[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA5, &fid3[0]);

	deleteSNarray(snA5);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk4");
	//Periodic in X only uniform sample
	//With electrodes along the x-axis
	//Electric Field in x and y
	//Single Charge large electric field
	//in the X
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA6 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 14;
	electricEnergyY = 7;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 2;
	EndY = 1;
	EndZ = 1;

	PeriodicX=1;
	PeriodicY=0;
	PeriodicZ=0;

	char fid4[] = "RandomWalk4\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA6,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA6);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA6, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	rv = randomWalk(snA6, CheckPointNum, &fid4[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA6, &fid4[0]);

	deleteSNarray(snA6);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk5");
	//Periodic in X only uniform sample
	//With electrodes along the x-axis
	//Electric Field in x and z
	//Single Charge large electric field
	//in the X
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA7 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 14;
	electricEnergyY = 0;
	electricEnergyZ = 6;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 2;
	EndY = 1;
	EndZ = 1;

	PeriodicX=1;
	PeriodicY=0;
	PeriodicZ=0;

	char fid5[30] = "RandomWalk5\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA7,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA7);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA7, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

rv = randomWalk(snA7, CheckPointNum, &fid5[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA7, &fid5[0]);

	deleteSNarray(snA7);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk6");
	//Non periodic sample
	//With electrodes along the y-axis only
	//Electric Field in y
	//Single Charge large electric field
	//in the y
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA8 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 0;
	electricEnergyY = 14;
	electricEnergyZ = 0;

	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	char fid6[30] = "RandomWalk6\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA8,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA8);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA8, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	rv = randomWalk(snA8, CheckPointNum, &fid6[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA8, &fid6[0]);

	deleteSNarray(snA8);
	
	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk7");
	//Periodic sample in y only
	//With electrodes along the y-axis only
	//Electric Field in y and in x
	//Single Charge large electric field
	//in the y and smaller in the x
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA9 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 5.5;
	electricEnergyY = 14;
	electricEnergyZ = 0;

	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;

	EndX = 1;
	EndY = 2;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=1;
	PeriodicZ=0;

	char fid7[30] = "RandomWalk7\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA9,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA9);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA9, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA9, CheckPointNum, &fid7[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA9, &fid7[0]);

	deleteSNarray(snA9);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk8");
	//Periodic sample in y only
	//With electrodes along the y-axis only
	//Electric Field in y and and z
	//Single Charge large electric field
	//in the y and smaller in the z
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA10 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 0;
	electricEnergyY = 14;
	electricEnergyZ = 7;

	XElecOn = 0;
	YElecOn = 1;
	ZElecOn = 0;

	EndX = 1;
	EndY = 2;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=1;
	PeriodicZ=0;

	char fid8[30] = "RandomWalk8\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA10,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA10);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA10, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA10, CheckPointNum, &fid8[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA10, &fid8[0]);

	deleteSNarray(snA10);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	printf("Testing RandomWalk9");
	//non-periodic 
	//With electrodes along the z-axis only
	//Electric Field in z and only
	//Single Charge large electric field
	//in the z
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA11 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 0;
	electricEnergyY = 0;
	electricEnergyZ = 14;

	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	char fid9[30] = "RandomWalk9\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA11,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray_Detailed(snA11);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA11, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA11, CheckPointNum, &fid9[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA11, &fid9[0]);

	deleteSNarray(snA11);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	printf("Testing RandomWalk10");
	//Periodic in z only 
	//With electrodes along the z-axis only
	//Electric Field in z and x only
	//Single Charge large electric field
	//in the z small in x
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA12 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 8;
	electricEnergyY = 0;
	electricEnergyZ = 14;

	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;

	EndX = 1;
	EndY = 1;
	EndZ = 2;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=1;

	char fid10[30] = "RandomWalk10\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA12,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA12);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA12, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA12, CheckPointNum, &fid10[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA12, &fid10[0]);

	deleteSNarray(snA12);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	printf("Testing RandomWalk11");
	//Periodic in z only 
	//With electrodes along the z-axis only
	//Electric Field in z and y only
	//Single Charge large electric field
	//in the z small in y
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA13 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 0;
	electricEnergyY = 7;
	electricEnergyZ = 14;

	XElecOn = 0;
	YElecOn = 0;
	ZElecOn = 1;

	EndX = 1;
	EndY = 1;
	EndZ = 2;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=1;

	char fid11[30] = "RandomWalk11\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA13,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA13);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA13, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	rv = randomWalk(snA13, CheckPointNum, &fid11[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA13, &fid11[0]);

	deleteSNarray(snA13);

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	printf("Testing RandomWalk12");
	//Periodic in x y and z 
	//With electrodes along the x-axis only
	//Electric Field in z and x only
	//Single Charge large electric field
	//in the x and small in the y and z
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA14 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 14;
	electricEnergyY = 8;
	electricEnergyZ = 8;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 2;
	EndY = 2;
	EndZ = 2;

	PeriodicX=1;
	PeriodicY=1;
	PeriodicZ=1;

	char fid12[30] = "RandomWalk12\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA14,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printSNarray(snA14);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA14, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);


	rv = randomWalk(snA14, CheckPointNum, &fid12[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA14, &fid12[0]);

	deleteSNarray(snA14);


	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	printf("Testing RandomWalk13");
	//Non periodic 
	//With electrodes along the x-axis only
	//Electric Field in the x
	//Single Charge medium electric field
	//placing cluster in the middle of the sample
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA15 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 8;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	printf("Creating Cluster");

	sn = getSN(snA15,5,2,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,2);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,2);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,3);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,3);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,4);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,4);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,5);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,5);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,6);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,6);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,7);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,7);
	setEnergy(sn,-100);

	sn = getSN(snA15,5,2,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,3,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,4,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,5,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,6,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,7,8);
	setEnergy(sn,-100);
	sn = getSN(snA15,5,8,8);
	setEnergy(sn,-100);


	char fid13[30] = "RandomWalk13\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA15,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	ArbArray ClArLL13;

	FindCluster(&OrderL, snA15, electricEnergyX, electricEnergyY, electricEnergyZ,\
			&ClArLL13, 1.0, 15.0, 15, AttemptToHop, gamma, SiteDistance, PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	printSNarray(snA15);
	printArbArray(ClArLL13);

	char fid13a[30] = "RandomWalk13cluster\0";
	PrintFile_xyz(OrderL,snA15, &ClArLL13, &fid13a[0]);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA15, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA15, CheckPointNum, &fid13[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA15, &fid13[0]);

	deleteSNarray(snA15);
	deleteArbArray(ClArLL13);	

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	printf("Testing RandomWalk14");
	//Non periodic 
	//With electrodes along the x-axis only
	//Electric Field in the x
	//Single Charge medium electric field
	//placing cluster right next to the back electrode 
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA16 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 1;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 1;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 8;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	printf("Creating Cluster");

	sn = getSN(snA16,0,2,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,2);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,2);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,3);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,3);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,4);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,4);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,5);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,5);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,6);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,6);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,7);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,7);
	setEnergy(sn,-100);

	sn = getSN(snA16,0,2,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,3,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,4,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,5,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,6,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,7,8);
	setEnergy(sn,-100);
	sn = getSN(snA16,0,8,8);
	setEnergy(sn,-100);


	char fid14[30] = "RandomWalk14\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 16.0, snA16,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	ArbArray ClArLL14;

	FindCluster(&OrderL, snA16, electricEnergyX, electricEnergyY, electricEnergyZ,\
			&ClArLL14, 1.0, 16.0, 16, AttemptToHop, gamma, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	printSNarray(snA16);
	printArbArray(ClArLL14);

	char fid14a[30] = "RandomWalk14cluster\0";
	PrintFile_xyz(OrderL,snA16, &ClArLL14, &fid14a[0]);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA16, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA16, CheckPointNum,&fid14[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence, FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA16, &fid14[0]);

	deleteSNarray(snA16);
	deleteArbArray(ClArLL14);	

	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}

	printf("Testing RandomWalk15");
	//Non periodic 
	//With electrodes along the x-axis only
	//Electric Field in the x
	//Single Charge medium electric field
	//placing cluster right next to the back electrode 
	SLength = 10;
	SWidth = 10;
	SHeight = 10;

	SNarray snA17 = newSNarray(SLength,SWidth,SHeight);
	//Number of time increments charges are injected
	Tcount = 2;
	//Size of time increments
	TStep = 1;
	//How many steps pass before we record the data
	Nstep_av = 1;
	//How many charges are injected at a time
	NCh = 30;
	//The distance between sites
	D = 1;
	//rN some normalization value
	rN = 1;

	electricEnergyX = 8;
	electricEnergyY = 0;
	electricEnergyZ = 0;

	XElecOn = 1;
	YElecOn = 0;
	ZElecOn = 0;

	EndX = 1;
	EndY = 1;
	EndZ = 1;

	PeriodicX=0;
	PeriodicY=0;
	PeriodicZ=0;

	printf("Creating Cluster");

	sn = getSN(snA17,5,2,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,2);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,2);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,3);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,3);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,4);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,4);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,5);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,5);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,6);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,6);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,7);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,7);
	setEnergy(sn,-20);

	sn = getSN(snA17,5,2,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,3,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,4,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,5,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,6,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,7,8);
	setEnergy(sn,-20);
	sn = getSN(snA17,5,8,8);
	setEnergy(sn,-20);


	char fid15[30] = "RandomWalk15\0";

	initJumPossibility( electricEnergyX,electricEnergyY,electricEnergyZ, MarcusCoeff,\
			1.0 , 15.0, snA17,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	ArbArray ClArLL15;

	FindCluster(&OrderL, snA17, electricEnergyX, electricEnergyY, electricEnergyZ,\
			&ClArLL15, 1.0, 15.0, 15, AttemptToHop, gamma, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);


	printSNarray(snA17);
	printArbArray(ClArLL15);

	char fid15a[30] = "RandomWalk15cluster\0";
	PrintFile_xyz(OrderL,snA17, &ClArLL15, &fid15a[0]);

	initElec(electricEnergyX, electricEnergyY, electricEnergyZ,\
			MarcusCoeff, KT, reOrgEnergy, snA17, RelativePerm, SiteDistance,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn,\
			FermiExb, FermiEyl, FermiEzb,\
			FermiExf, FermiEyr, FermiEza,\
			alphaxb, alphayl, alphazb,\
			alphaxf, alphayr, alphaza,\
			&elXb, &elXf, &elYl, &elYr,&elZb,&elZa,\
			vX, vY, vZ);

	rv = randomWalk(snA17, CheckPointNum, &fid15[0],\
			ElectricFieldX, ElectricFieldY, ElectricFieldZ,\
			elXb, elXf, elYl, elYr, elZb, elZa,\
			PF, t, Sequence,FutureSite, &chA,\
			n, nc, nca);

	printVisitFreq(snA17, &fid15[0]);

	deleteSNarray(snA17);
	deleteArbArray(ClArLL15);	
		
	if(elXb!=NULL){
		if((matrix) getElectrode_HopRates(elXb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXb));
		}
		if((SNarray) getElectrode_AdjacentSites(elXb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXb));
		}
		deleteElectrode(elXb);
		elXb=NULL;
	}
	if(elXf!=NULL){
		if((matrix) getElectrode_HopRates(elXf)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elXf));
		}
		if((SNarray) getElectrode_AdjacentSites(elXf)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elXf));
		}
		deleteElectrode(elXf);
		elXf=NULL;
	}
	if(elYl!=NULL){
		if((matrix) getElectrode_HopRates(elYl)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYl));
		}
		if((SNarray) getElectrode_AdjacentSites(elYl)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYl));
		}
		deleteElectrode(elYl);
		elYl=NULL;
	}
	if(elYr!=NULL){
		if((matrix) getElectrode_HopRates(elYr)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elYr));
		}	
		if((SNarray) getElectrode_AdjacentSites(elYr)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elYr));
		}
		deleteElectrode(elYr);
		elYr=NULL;
	}
	if(elZb!=NULL){
		if((matrix) getElectrode_HopRates(elZb)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZb));
		}
		if((SNarray) getElectrode_AdjacentSites(elZb)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZb));
		}
		deleteElectrode(elZb);
		elZb=NULL;
	}
	if(elZa!=NULL){
		if((matrix) getElectrode_HopRates(elZa)!=NULL){
			deleteMatrix((matrix) getElectrode_HopRates(elZa));
		}
		if((SNarray) getElectrode_AdjacentSites(elZa)!=NULL){
			deleteSNarray((SNarray)getElectrode_AdjacentSites(elZa));
		}
		deleteElectrode(elZa);
		elZa=NULL;
	}
	
	//atexit(mem_term);
	*/
	deleteParamFrame(&PF);
	return 0;

}
