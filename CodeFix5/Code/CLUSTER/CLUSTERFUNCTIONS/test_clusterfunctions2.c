#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "DATASTRUCT/cluster.h"
#include "SITENODE/sitenode.h"
#include "MATRIX/matrix.h"
#include "clusterfunctions.h"
#include "../../PARAMETERS/read.h"
//#include "../../MEM/mem.h"

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

	printf("Testing:CalculateAllHops\n");

	SNarray snA = newSNarray( 3,3,3);
	assert(snA!=NULL);
	matrix mtx = CalculateAllHops(NULL, 1,0,0, 1, 1,SiteDistance,AttemptToHop,gamma, 0, 0, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for NULL snA returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, -1,1, SiteDistance,AttemptToHop,gamma,0, 0, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for -1 KT returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,-1, 0, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for -1 PeriodicX returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,2, 0, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for 2 PeriodicX returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, -1, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for -1 PeriodicY returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 2, 0);
	assert(mtx==NULL);
	printf("CalculateAllHops for 2 PeriodicY returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, -1);
	assert(mtx==NULL);
	printf("CalculateAllHops for -1 PeriodicZ returns NULL PASS\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 2);
	assert(mtx==NULL);
	printf("CalculateAllHops for 2 PeriodicZ returns NULL PASS\n");

	printf("Testing Non periodic conditions\t");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 0);
	assert(mtx!=NULL);

	int i;
	int j;
	double elem;
	double elem2;
	double val = -0.25;
	double val2 = MarcusCoeff*exp( val );

	//Count the number of non zero elements because there are a total of 3x3x3 sites
	//and the sites on the edges cannot hop (162 hops total) - (3x3*6)= 108

	int count;
	count = 0;
	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);

			if(elem!=0){
				count = count+1;
				assert(elem== val2 );
			}

		}
	}
	if(count==108){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	deleteMatrix(mtx);
	/////////////////////////////////////////////////////////////////////////////////
	printf("Testing full periodic conditions in x, y and z\t");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);

	count = 0;
	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			assert(elem == val2 );
			assert(elem2 == 17);
			if(elem!=0){
				count=count+1;
			}
		}
	}

	if(count==162){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	deleteMatrix(mtx);

	//(162 hops total) - (3x3*4)= 144
	printf("Testing periodic in x\t");
	count = 0;
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 0, 0);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==7 || j==8){
				count=count+1;
				assert(elem == val2 );
				assert(elem2 == 17);
			}else if(elem!=0){
				count++;
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	if(count==126){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	deleteMatrix(mtx);

	printf("Testing periodic in y\t");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 1, 0);
	assert(mtx!=NULL);
	count = 0;

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==9 || j==10){
				assert(elem == val2 );
				assert(elem2 == 17);
				count++;
			}else if(elem!=0){
				count++;
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}
	
	if(count==126){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	deleteMatrix(mtx);

	printf("Testing periodic in z\t");
	count = 0;
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 1);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==11 || j==12){
				assert(elem == val2 );
				assert(elem2 == 17);
				count++;
			}else if(elem!=0){
				count++;
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	if(count==126){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
	}

	deleteMatrix(mtx);
	printf("Testing reverse bias\t");
	count = 0;
	mtx = CalculateAllHops(snA, -1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);

	double val3 = MarcusCoeff*exp(-1);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==7 ){
				assert(elem == MarcusCoeff );
				assert(elem2 == 17);
				count++;
			}else if(j==8){
				count++;
				assert(elem==val3);
				assert(elem2 == 16);
			}else if(elem!=0){
				count++;
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	if(count==162){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	deleteMatrix(mtx);
	
	////////////////////////////////////////////////////////////////////////////////
	printf("Testing forward bias\t");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);
	count = 0;

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==8 ){
				assert(elem == MarcusCoeff );
				assert(elem2 == 17);
				count++;
			}else if(j==7){
				assert(elem==val3);
				assert(elem2 == 16);
				count++;
			}else if(elem!=0){
				count++;
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}
	if(count==162){
		printf("PASS\n");
	}else{
		printf("FAIL\n");
		exit(1);
	}

	printf("Testing:MPsort\n");

	int orderL=0;
	int orderH=0;
	int MidPtsTotal=0;

	printf("MPsort with NULL orderL should return NULL\t");
	ArbArray mpA = MPsort( NULL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with NULL orderH should return NULL\t");
	mpA = MPsort( &orderL, NULL, &MidPtsTotal, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with NULL MidPtsTotal should return NULL\t");
	mpA = MPsort( &orderL, &orderH, NULL, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with NULL mtx should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, NULL, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with NULL snA should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, NULL, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with -1 PeriodicX should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, -1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with -1 PeriodicY should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, -1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with -1 PeriodicZ should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, -1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with 2 PeriodicX should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 2, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with 2 PeriodicY should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 2, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	deleteArbArray(&mpA);
	printf("MPsort with 2 PeriodicZ should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 2);
	assert(mpA==NULL);
	printf("PASS\n");

	SNarray snA2 = newSNarray( 3,4,3);
	deleteArbArray(&mpA);
	printf("MPsort with different site node array that does not match the input matrix should return NULL\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA2, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	matrix mtx2 = newMatrix(2,12);
	deleteArbArray(&mpA);
	printf("MPsort with different size matrix that does not match the site node array\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx2, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");
	matrix mtx3 = newMatrix(getAtotal(snA),11);
	deleteArbArray(&mpA);
	printf("MPsort with different size matrix that does not match the site node array test 2\t");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx3, snA, 1, 1, 1);
	assert(mpA==NULL);
	printf("PASS\n");

	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 81 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==81);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);

	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 1, 1);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 72 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 0, 1);

	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 72 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 0);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 72 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 0, 1);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 72 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 1);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 63 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 1, 0);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 63 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 1);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 63 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 0);
	printf("MPsort testing that correctly calculates:\n");
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 54 %d\n",MidPtsTotal);
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==54);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	printf("PASS\n");

	deleteMatrix(mtx);
	mtx = CalculateAllHops(snA, 2,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);
	//printMatrix(mtx);

	deleteAllMidPointArray(&mpA);
	printf("MPsort testing with new snA\n");
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);
	printf("orderL 16 %d\n",orderL);
	printf("orderH 17 %d\n",orderH);
	printf("MidPtsTotal 81 %d\n",MidPtsTotal);
	
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==81);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	printf("Testing:SortOrderMag\n");

	printf("SortOrderMag when Total Orders is less than 0 should return NULL\t");
	ArbArray Arb2 = SortOrderMag( 0, orderL, mpA);
	assert(Arb2==NULL);
	printf("PASS\n");
	deleteArbArray(&Arb2);
	printf("SortOrderMag when mpA is NULL should return NULL\t");
	Arb2 = SortOrderMag( orderH-orderL+1, orderL, NULL);
	assert(Arb2==NULL);
	printf("PASS\n");
	deleteArbArray(&Arb2);
	printf("SortOrderMag when correct arguments are passed should return non NULL\t");
	Arb2 = SortOrderMag( orderH-orderL+1, orderL, mpA);
	assert(Arb2!=NULL);
	printf("PASS\n");

	printf("SortOrderMag testing that total elements in first Order of magnitude link list are 27\t");
	OrderMagLL OMLL = (OrderMagLL) getArbElement( Arb2, 0);
	assert(getOMLL_size(OMLL)==27);
	printf("PASS\n");
	printf("SortOrderMag testing that starting midpoint of first order of magnitude link list is non NULL\t");
	MidPoint mp3 = getOMLLstartMP( OMLL );
	assert(mp3!=NULL);
	printf("PASS\n");
	printf("SortOrderMag testing that midpoint has correct order of magnitude associated with it should be 16\t");
	assert(getMP_order(mp3)==16);
	printf("PASS\n");

	OMLL = (OrderMagLL) getArbElement( Arb2, 1);
	assert(getOMLL_size(OMLL)==54);
	mp3 = getOMLLstartMP( OMLL );
	assert(mp3!=NULL);
	assert(getMP_order(mp3)==17);

	printf("Testing:ClusterSort\n");
	ArbArray Arb3 = ClusterSort( orderH-orderL+1, orderL, NULL);
	assert(Arb3==NULL);
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( 0, orderL, Arb2);
	assert(Arb3==NULL);
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( 3, orderL, Arb2);
	assert(Arb3==NULL);
	
	int TotalOrders = orderH-orderL+1;
	//Sorts Midpoints into clusters based on their order of
	//magnitude and stores them as node data structures
	//Does not ween out degenerate states
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( orderH-orderL+1, orderL, Arb2);
	assert(Arb3!=NULL);

	ClusterLL clLL = (ClusterLL) getArbElement( Arb3, 0);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(getCluster_numNodes(clLL)==3);
	clLL = getNextClusterLL(clLL);
	assert(clLL==NULL);
	ClusterLL clLL2 = (ClusterLL) getArbElement( Arb3, 1);
	assert(getCluster_numNodes(clLL2)==9);
	clLL2 = getNextClusterLL(clLL2);
	assert(getCluster_numNodes(clLL2)==9);
	clLL2 = getNextClusterLL(clLL2);
	assert(getCluster_numNodes(clLL2)==9);
	clLL2 = getNextClusterLL(clLL2);
	assert(clLL2==NULL);

	printf("Testing:FilterCluster\n");
	//printArbArray(mpA);
	//printSNarray(snA);

	int PeriodicX = 1;
	int PeriodicY = 1;
	int PeriodicZ = 1;

	rv = FilterCluster( 0, orderL, mtx, &Arb3, snA, PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New1\n");
	rv = FilterCluster( TotalOrders, orderL, NULL, &Arb3, snA, PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New2\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, NULL, snA, PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New3\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, NULL, PeriodicX, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New4\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, -1, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New5\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, PeriodicX, -1, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New6\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, PeriodicX, PeriodicY, -1, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New7\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, 2, PeriodicY, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New8\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, PeriodicX, 2, PeriodicZ, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);
	printf("New9\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, PeriodicX, PeriodicY, 2, XElecOn, YElecOn, ZElecOn);
	assert(rv==-1);

	printf("New10\n");
	rv = FilterCluster( TotalOrders, orderL, mtx, &Arb3, snA, 1, 1, 1, XElecOn, YElecOn, ZElecOn);
	assert(rv==0);
	assert(getElementsUsed(Arb3)==0);
	assert((ClusterLL) getArbElement(Arb3,0)==NULL);

	printf("Testing:PrintCheck\n");
	rv = PrintCheck( 0, orderL, Arb3, snA, mtx);
	assert(rv==-1);
	rv = PrintCheck( TotalOrders, orderL, NULL, snA, mtx);
	assert(rv==-1);
	rv = PrintCheck( TotalOrders, orderL, Arb3, NULL, mtx);
	assert(rv==-1);
	rv = PrintCheck( TotalOrders, orderL, Arb3, snA, NULL);
	assert(rv==-1);
	//rv = PrintCheck( TotalOrders, orderL, Arb3, snA, mtx);
	//assert(rv==0);

	//Deleting datastructures
	deleteAllMidPointArray(&mpA);
	deleteMatrix(mtx);
	deleteMatrix(mtx2);
	deleteMatrix(mtx3);
	deleteSNarray(snA);
	deleteSNarray(snA2);
	deleteArbArray(&Arb2);
	deleteArbArray(&Arb3);
	
	return 0;
}
