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

	//mem_init();

	printf("Testing:CalculateAllHops\n");

	SNarray snA = newSNarray( 3,3,3);
	assert(snA!=NULL);
	matrix mtx = CalculateAllHops(NULL, 1,0,0, 1, 1,SiteDistance,AttemptToHop,gamma, 0, 0, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, -1,1, SiteDistance,AttemptToHop,gamma,0, 0, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,-1, 0, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,2, 0, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, -1, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 2, 0);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, -1);
	assert(mtx==NULL);
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 2);
	assert(mtx==NULL);

	printf("Testing Non periodic conditions\n");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 0);
	assert(mtx!=NULL);

	int i;
	int j;
	double elem;
	double elem2;
	double val = -0.25;
	double val2 = MarcusCoeff*exp( val );

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);

			if(elem!=0){
				printf("Index %d elem %g val2 %g\n",i,elem,val2);
				assert(elem== val2 );
			}

		}
	}
	
	deleteMatrix(&mtx);
	
  /////////////////////////////////////////////////////////////////////////////////
	printf("Testing full periodic conditions in x, y and z\n");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			assert(elem == val2 );
			assert(elem2 == 17);
		}
	}

	deleteMatrix(&mtx);

	printf("Testing periodic in x\n");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 0, 0);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==7 || j==8){
				assert(elem == val2 );
				assert(elem2 == 17);
			}else if(elem!=0){
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	deleteMatrix(&mtx);

	printf("Testing periodic in y\n");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 1, 0);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==9 || j==10){
				assert(elem == val2 );
				assert(elem2 == 17);
			}else if(elem!=0){
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	deleteMatrix(&mtx);

	printf("Testing periodic in z\n");
	mtx = CalculateAllHops(snA, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,0, 0, 1);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==11 || j==12){
				assert(elem == val2 );
				assert(elem2 == 17);
			}else if(elem!=0){
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	deleteMatrix(&mtx);
	printf("Testing reverse bias\n");
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
			}else if(j==8){
				assert(elem==val3);
				assert(elem2 == 16);
			}else if(elem!=0){
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	deleteMatrix(&mtx);
	
	////////////////////////////////////////////////////////////////////////////////
	printf("Testing forward bias\n");
	mtx = CalculateAllHops(snA, 1,0,0, 1,1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);

	for(i=1;i<=getRows(mtx);i++){
		for(j=7;j<=getCols(mtx);j++){
			elem = getE(mtx,i,j);
			elem2 = getE(mtx,i,j-6);
			if(j==8 ){
				assert(elem == MarcusCoeff );
				assert(elem2 == 17);
			}else if(j==7){
				assert(elem==val3);
				assert(elem2 == 16);
			}else if(elem!=0){
				assert(elem == val2 );
				assert(elem2 == 17);
			}
		}
	}

	printf("Testing:MPsort\n");

	int orderL=0;
	int orderH=0;
	int MidPtsTotal=0;


	ArbArray mpA = MPsort( NULL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, NULL, &MidPtsTotal, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, NULL, mtx, snA, 1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, NULL, snA, 1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, NULL, 1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, -1, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, -1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, -1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 2, 1, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 2, 1);
	assert(mpA==NULL);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 2);
	assert(mpA==NULL);

	SNarray snA2 = newSNarray( 3,4,3);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA2, 1, 1, 1);
	assert(mpA==NULL);
	matrix mtx2 = newMatrix(2,12);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx2, snA, 1, 1, 1);
	assert(mpA==NULL);
	matrix mtx3 = newMatrix(getAtotal(snA),11);
	deleteArbArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx3, snA, 1, 1, 1);
	assert(mpA==NULL);

	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==81);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));
	
	deleteAllMidPointArray(&mpA);

	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 1, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 0, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 0);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 0, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==72);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 1, 0);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==63);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 0, 0, 0);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==54);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	deleteMatrix(&mtx);
	mtx = CalculateAllHops(snA, 2,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,1, 1, 1);
	assert(mtx!=NULL);
	printMatrix(mtx);

	deleteAllMidPointArray(&mpA);
	mpA = MPsort( &orderL, &orderH, &MidPtsTotal, mtx, snA, 1, 1, 1);

	printf("orderL %d orderH %d MidPtsTotal %d TotalsnA %d\n",orderL,orderH,MidPtsTotal,getAtotal(snA));
	assert(mpA!=NULL);
	assert(orderL==16);
	assert(orderH==17);
	assert(MidPtsTotal==81);
	assert(getElementsReserved(mpA)==getElementsUsed(mpA));

	printf("Testing:SortOrderMag\n");

	ArbArray Arb2 = SortOrderMag( 0, orderL, mpA);
	assert(Arb2==NULL);
	deleteArbArray(&Arb2);
	Arb2 = SortOrderMag( orderH-orderL+1, orderL, NULL);
	assert(Arb2==NULL);
	deleteArbArray(&Arb2);
	Arb2 = SortOrderMag( orderH-orderL, orderL, mpA);
	assert(Arb2==NULL);
	deleteArbArray(&Arb2);
	Arb2 = SortOrderMag( orderH-orderL+1, orderL, mpA);
	assert(Arb2!=NULL);

	OrderMagLL OMLL = (OrderMagLL) getArbElement( Arb2, 0);
	assert(getOMLL_size(OMLL)==27);
	MidPoint mp3 = getOMLLstartMP( OMLL );
	assert(mp3!=NULL);
	assert(getMP_order(mp3)==16);

	OMLL = (OrderMagLL) getArbElement( Arb2, 1);
	assert(getOMLL_size(OMLL)==54);
	mp3 = getOMLLstartMP( OMLL );
	assert(mp3!=NULL);
	assert(getMP_order(mp3)==17);

	printf("Testing:ClusterSort\n");
	ArbArray Arb3 = ClusterSort( orderH-orderL+1, NULL);
	assert(Arb3==NULL);
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( 0, Arb2);
	assert(Arb3==NULL);
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( 3, Arb2);
	assert(Arb3==NULL);
	
	int TotalOrders = orderH-orderL+1;
	//Sorts Midpoints into clusters based on their order of
	//magnitude and stores them as node data structures
	//Does not ween out degenerate states
	deleteArbArray(&Arb3);
	Arb3 = ClusterSort( orderH-orderL+1, Arb2);
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
	printArbArray(mpA);
	printSNarray(snA);

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
	rv = PrintCheck( TotalOrders, orderL, Arb3, snA, mtx);
	assert(rv==0);

	//Deleting datastructures
	deleteAllMidPointArray(&mpA);
	deleteMatrix(&mtx);
	deleteMatrix(&mtx2);
	deleteMatrix(&mtx3);
	deleteSNarray(&snA);
	deleteSNarray(&snA2);
	deleteArbArray(&Arb2);
	deleteArbArray(&Arb3);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//There should be no clusters
	
	//Case Study one clusters in the center Periodic
	int OrderLT;
	int OrderHT;
	int MidPtsTotalT;
	int PeriodicXT = 1;
	int PeriodicYT = 1;
	int PeriodicZT = 1;
	int XElecOnT = 1;
	int YElecOnT = 0;
	int ZElecOnT = 0;
	SNarray snAT = newSNarray( 3,4,3);
	//Create an artificial trap in snAT which in increment length is
	//length 0-2 width 0-3 and height 0-2
	SiteNode snT = getSN(snAT,1,1,1);
	printSNarray(snAT);
	rv = setEnergy(snT, -3);
	assert(rv==0);
	snT = getSN(snAT,1,2,1);
	rv = setEnergy(snT, -3);
	assert(rv==0);
	printSNarray(snAT);
	//Testing Periodic x,y and z 0 biasX 0 biasY 0 biasZ KT=1 and reOrgEnergy=1
	printf("\n\n\nCase Study T Periodic Conditions\n\n\n");
	matrix mtxT = CalculateAllHops(snAT, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,PeriodicXT, PeriodicYT, PeriodicZT);
	printMatrix(mtxT);
	//Sorts all the data into midpoints
	ArbArray mpAT = MPsort( &OrderLT, &OrderHT, &MidPtsTotalT, mtxT, snAT, PeriodicXT, PeriodicYT, PeriodicZT);
	//Sorts midpoints into Nodes in link lists. Each element of mpAT is composed of a
	//link list all of the same order of magnitude
	printf("Order Low %d Order High %d\n",OrderLT, OrderHT);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAT);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbT2 = SortOrderMag( OrderHT-OrderLT+1, OrderLT, mpAT);
	assert(ArbT2!=NULL);
	int TotalOrdersT = OrderHT-OrderLT+1;
	//Sorts nodes into clusters
	printf("********************************************\n");
	printArbArray(ArbT2, OrderLT);
	printf("********************************************\n");
	ArbArray ArbT = ClusterSort( TotalOrdersT, ArbT2);
	printf("TotalOrdersT %d OrderLowT %d Elements used by ArbT2 %d\n",TotalOrdersT, OrderLT, getElementsUsed(ArbT2));
	assert(ArbT!=NULL);
	//Filtering Cluster
	rv = FilterCluster( TotalOrdersT, OrderLT, mtxT, &ArbT, snAT, PeriodicXT, PeriodicYT, PeriodicZT, XElecOnT, YElecOnT, ZElecOnT);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersT, OrderLT, ArbT, snAT, mtxT);
	//Because the cluster is smaller energy than the rest of the nodes the rest
	//of the nodes are not considered to be a cluster

	deleteArbArray(&ArbT2);
	////////////////////////////////////////////////////////////////////////////////////////
	//Case Study no clusters Periodic
	int OrderLU;
	int OrderHU;
	int MidPtsTotalU;
	int PeriodicXU = 1;
	int PeriodicYU = 1;
	int PeriodicZU = 1;
	int XElecOnU = 1;
	int YElecOnU = 0;
	int ZElecOnU = 0;
	SNarray snAU = newSNarray( 3,4,3);
	//Create an artificial trap in snAT which in increment length is
	//length 0-2 width 0-3 and height 0-2
	printSNarray(snAU);
	//Testing Periodic x,y and z 0 biasX 0 biasY 0 biasZ KT=1 and reOrgEnergy=1
	matrix mtxU = CalculateAllHops(snAU, 0,0,0, 1,1, SiteDistance,AttemptToHop,gamma,PeriodicXU, PeriodicYU, PeriodicZU);
	printMatrix(mtxU);
	//Sorts all the data into midpoints
	ArbArray mpAU = MPsort( &OrderLU, &OrderHU, &MidPtsTotalU, mtxU, snAU, PeriodicXU, PeriodicYU, PeriodicZU);
	//Sorts midpoints into Nodes in link lists. Each element of mpAT is composed of a
	//link list all of the same order of magnitude
	printf("Order Low %d Order High %d\n",OrderLU, OrderHU);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAU);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbU2 = SortOrderMag( OrderHU-OrderLU+1, OrderLU, mpAU);
	assert(ArbU2!=NULL);
	int TotalOrdersU = OrderHU-OrderLU+1;
	//Sorts nodes into clusters
	printf("********************************************\n");
	printArbArray(ArbU2, OrderLU);
	printf("********************************************\n");
	ArbArray ArbU = ClusterSort( TotalOrdersU, ArbU2);
	printf("TotalOrdersT %d OrderLowT %d Elements used by ArbT2 %d\n",TotalOrdersU, OrderLU, getElementsUsed(ArbU2));
	printf("****************STARTING FILTER****************************\n");
	assert(ArbU!=NULL);
	//Filtering Cluster
	rv = FilterCluster( TotalOrdersU, OrderLU, mtxU, &ArbU, snAU, PeriodicXU, PeriodicYU, PeriodicZU, XElecOnU, YElecOnU, ZElecOnU);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersU, OrderLU, ArbU, snAU, mtxU);

	deleteArbArray(&ArbU2);	
	deleteArbArray(&ArbU);	
	deleteMatrix(&mtxU);
	deleteSNarray(&snAU);
	deleteAllMidPointArray(&mpAU);

	//Case Study two clusters in middle non periodic
	int OrderLS;
	int OrderHS;
	int MidPtsTotalS;
	int PeriodicXS = 0;
	int PeriodicYS = 0;
	int PeriodicZS = 0;
	int XElecOnS = 1;
	int YElecOnS = 0;
	int ZElecOnS = 0;
	SNarray snAS = newSNarray( 4, 3, 3);
	//Create Trap Sites
	SiteNode snS = getSN(snAS,2,1,1);
	printSNarray(snAS);
	rv = setEnergy(snS, -4);
	assert(rv==0);
	snS = getSN(snAS, 1,1,1);
	rv = setEnergy(snS, -4);
	assert(rv==0);
	printSNarray(snAS);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxS = CalculateAllHops(snAS, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXS, PeriodicYS, PeriodicZS);
	printMatrix(mtxS);
	ArbArray mpAS = MPsort( &OrderLS, &OrderHS, &MidPtsTotalS, mtxS, snAS, PeriodicXS, PeriodicYS, PeriodicZS);
	printf("Order Low %d Order High %d\n",OrderLS, OrderHS);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAS);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbS2 = SortOrderMag( OrderHS-OrderLS+1, OrderLS, mpAS);
	assert(ArbS2!=NULL);
	int TotalOrdersS = OrderHS-OrderLS+1;
	printf("********************************************\n");
	printArbArray(ArbS2, OrderLS);
	printf("********************************************\n");
	ArbArray ArbS = ClusterSort( TotalOrdersS, ArbS2);	
	printf("TotalOrdersS %d OrderLowS %d Elements used by ArbS2 %d\n",TotalOrdersS, OrderLS, getElementsUsed(ArbS2));
	assert(ArbS!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbS, OrderLS);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersS, OrderLS, mtxS, &ArbS, snAS, PeriodicXS, PeriodicYS, PeriodicZS, XElecOnS, YElecOnS, ZElecOnS);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersS, OrderLS, ArbS, snAS, mtxS);

	deleteMatrix(&mtxS);
	deleteSNarray(&snAS);
	deleteArbArray(&ArbS);
	deleteAllMidPointArray(&mpAS);
	deleteArbArray(&ArbS2);

	//Case Study two cluster on oposite edges non periodic
	int OrderLR;
	int OrderHR;
	int MidPtsTotalR;
	int PeriodicXR = 0;
	int PeriodicYR = 0;
	int PeriodicZR = 0;
	int XElecOnR = 1;
	int YElecOnR = 0;
	int ZElecOnR = 0;
	SNarray snAR = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snR = getSN(snAR,0,1,1);
	printSNarray(snAR);
	rv = setEnergy(snR, -4);
	assert(rv==0);
	snS = getSN(snAR, 3,1,1);
	rv = setEnergy(snR, -4);
	assert(rv==0);
	printSNarray(snAR);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxR = CalculateAllHops(snAR, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXR, PeriodicYR, PeriodicZR);
	printMatrix(mtxR);
	ArbArray mpAR = MPsort( &OrderLR, &OrderHR, &MidPtsTotalR, mtxR, snAR, PeriodicXR, PeriodicYR, PeriodicZR);
	printf("Order Low %d Order High %d\n",OrderLR, OrderHR);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAR);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbR2 = SortOrderMag( OrderHR-OrderLR+1, OrderLR, mpAR);
	assert(ArbR2!=NULL);
	int TotalOrdersR = OrderHR-OrderLR+1;
	printf("********************************************\n");
	printArbArray(ArbR2, OrderLR);
	printf("********************************************\n");
	ArbArray ArbR = ClusterSort( TotalOrdersR, ArbR2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersR, OrderLR, getElementsUsed(ArbR2));
	assert(ArbR!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbR, OrderLR);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersR, OrderLR, mtxR, &ArbR, snAR, PeriodicXR, PeriodicYR, PeriodicZR, XElecOnR, YElecOnR, ZElecOnR);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersR, OrderLR, ArbR, snAR, mtxR);

	deleteMatrix(&mtxR);
	deleteSNarray(&snAR);
	deleteAllMidPointArray(&mpAR);
	deleteArbArray(&ArbR2);
	deleteArbArray(&ArbR);



	//Case Study two cluster on oposite edges periodic in x
	int OrderLP;
	int OrderHP;
	int MidPtsTotalP;
	int PeriodicXP = 1;
	int PeriodicYP = 0;
	int PeriodicZP = 0;
	int XElecOnP = 1;
	int YElecOnP = 0;
	int ZElecOnP = 0;
	SNarray snAP = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snP = getSN(snAP,0,1,1);
	printSNarray(snAP);
	rv = setEnergy(snP, -4);
	assert(rv==0);
	snP = getSN(snAP, 3,1,1);
	rv = setEnergy(snP, -4);
	assert(rv==0);
	printSNarray(snAP);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxP = CalculateAllHops(snAP, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXP, PeriodicYP, PeriodicZP);
	printMatrix(mtxP);
	ArbArray mpAP = MPsort( &OrderLP, &OrderHP, &MidPtsTotalP, mtxP, snAP, PeriodicXP, PeriodicYP, PeriodicZP);
	printf("Order Low %d Order High %d\n",OrderLP, OrderHP);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAP);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbP2 = SortOrderMag( OrderHP-OrderLP+1, OrderLP, mpAP);
	assert(ArbP2!=NULL);
	int TotalOrdersP = OrderHP-OrderLP+1;
	printf("********************************************\n");
	printArbArray(ArbP2, OrderLP);
	printf("********************************************\n");
	ArbArray ArbP = ClusterSort( TotalOrdersP, ArbP2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersP, OrderLP, getElementsUsed(ArbP2));
	assert(ArbP!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbP, OrderLP);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersP, OrderLP, mtxP, &ArbP, snAP, PeriodicXP, PeriodicYP, PeriodicZP, XElecOnP, YElecOnP, ZElecOnP);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersP, OrderLP, ArbP, snAP, mtxP);
	deleteArbArray(&ArbP2);

	//Case Study two cluster on oposite edges periodic in y
	int OrderLO;
	int OrderHO;
	int MidPtsTotalO;
	int PeriodicXO = 0;
	int PeriodicYO = 1;
	int PeriodicZO = 0;
	int XElecOnO = 1;
	int YElecOnO = 0;
	int ZElecOnO = 0;
	SNarray snAO = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snO = getSN(snAO,1,0,1);
	printSNarray(snAO);
	rv = setEnergy(snO, -4);
	assert(rv==0);
	snO = getSN(snAO, 1,2,1);
	rv = setEnergy(snO, -4);
	assert(rv==0);
	printSNarray(snAO);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxO = CalculateAllHops(snAO, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXO, PeriodicYO, PeriodicZO);
	printMatrix(mtxO);
	ArbArray mpAO = MPsort( &OrderLO, &OrderHO, &MidPtsTotalO, mtxO, snAO, PeriodicXO, PeriodicYO, PeriodicZO);
	printf("Order Low %d Order High %d\n",OrderLO, OrderHO);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAO);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbO2 = SortOrderMag( OrderHO-OrderLO+1, OrderLO, mpAO);
	assert(ArbO2!=NULL);
	int TotalOrdersO = OrderHO-OrderLO+1;
	printf("********************************************\n");
	printArbArray(ArbO2, OrderLO);
	printf("********************************************\n");
	ArbArray ArbO = ClusterSort( TotalOrdersO, ArbO2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersP, OrderLO, getElementsUsed(ArbO2));
	assert(ArbO!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbO, OrderLO);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersP, OrderLO, mtxO, &ArbO, snAO, PeriodicXO, PeriodicYO, PeriodicZO, XElecOnO, YElecOnO, ZElecOnO);
	assert(rv==0);
	printf("****************CHECK***********************\n");
	rv = PrintCheck( TotalOrdersO, OrderLO, ArbO, snAO, mtxO);

	deleteArbArray(&ArbO2);
	//Case Study two cluster on oposite edges non-periodic in y
	int OrderLM;
	int OrderHM;
	int MidPtsTotalM;
	int PeriodicXM = 0;
	int PeriodicYM = 0;
	int PeriodicZM = 0;
	int XElecOnM = 1;
	int YElecOnM = 0;
	int ZElecOnM = 0;
	SNarray snAM = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snM = getSN(snAM,1,0,1);
	printSNarray(snAM);
	rv = setEnergy(snM, -4);
	assert(rv==0);
	snM = getSN(snAM, 1,2,1);
	rv = setEnergy(snM, -4);
	assert(rv==0);
	printSNarray(snAM);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ, KT=0.024, and reOrgEnergy =1
	matrix mtxM = CalculateAllHops(snAM, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXM, PeriodicYM, PeriodicZM);
	printMatrix(mtxM);
	ArbArray mpAM = MPsort( &OrderLM, &OrderHM, &MidPtsTotalM, mtxM, snAM, PeriodicXM, PeriodicYM, PeriodicZM);
	printf("Order Low %d Order High %d\n",OrderLM, OrderHM);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAM);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbM2 = SortOrderMag( OrderHM-OrderLM+1, OrderLM, mpAM);
	assert(ArbM2!=NULL);
	int TotalOrdersM = OrderHM-OrderLM+1;
	printf("********************************************\n");
	printArbArray(ArbM2, OrderLM);
	printf("********************************************\n");
	ArbArray ArbM = ClusterSort( TotalOrdersM, ArbM2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersM, OrderLM, getElementsUsed(ArbM2));
	assert(ArbM!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbM, OrderLM);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersM, OrderLM, mtxM, &ArbM, snAM, PeriodicXM, PeriodicYM, PeriodicZM, XElecOnM, YElecOnM, ZElecOnM);
	assert(rv==0);
	printf("****************CHECK M***********************\n");
	rv = PrintCheck( TotalOrdersM, OrderLM, ArbM, snAM, mtxM);

	deleteArbArray(&ArbM2);
	//Case Study two cluster on oposite edges periodic in z
	int OrderLN;
	int OrderHN;
	int MidPtsTotalN;
	int PeriodicXN = 0;
	int PeriodicYN = 0;
	int PeriodicZN = 0;
	int XElecOnN = 1;
	int YElecOnN = 0;
	int ZElecOnN = 0;
	SNarray snAN = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snN = getSN(snAN,1,1,0);
	printSNarray(snAN);
	rv = setEnergy(snN, -4);
	assert(rv==0);
	snN = getSN(snAN, 1,1,2);
	rv = setEnergy(snN, -4);
	assert(rv==0);
	printSNarray(snAN);
	//Testing non-periodic x, y and z 0 biasX, 0 biasY, 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxN = CalculateAllHops(snAN, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXN, PeriodicYN, PeriodicZN);
	printMatrix(mtxN);
	ArbArray mpAN = MPsort( &OrderLN, &OrderHN, &MidPtsTotalN, mtxN, snAN, PeriodicXN, PeriodicYN, PeriodicZN);
	printf("Order Low %d Order High %d\n",OrderLN, OrderHN);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAN);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbN2 = SortOrderMag( OrderHN-OrderLN+1, OrderLN, mpAN);
	assert(ArbN2!=NULL);
	int TotalOrdersN = OrderHN-OrderLN+1;
	printf("********************************************\n");
	printArbArray(ArbN2, OrderLN);
	printf("********************************************\n");
	ArbArray ArbN = ClusterSort( TotalOrdersN, ArbN2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersN, OrderLN, getElementsUsed(ArbN2));
	assert(ArbN!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbN, OrderLN);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersN, OrderLN, mtxN, &ArbN, snAN, PeriodicXN, PeriodicYN, PeriodicZN, XElecOnN, YElecOnN, ZElecOnN);
	assert(rv==0);
	printf("****************CHECK N***********************\n");
	rv = PrintCheck( TotalOrdersN, OrderLN, ArbN, snAN, mtxN);

	deleteArbArray(&ArbN2);

	//Case Study two cluster on oposite edges periodic in z
	int OrderLK;
	int OrderHK;
	int MidPtsTotalK;
	int PeriodicXK = 0;
	int PeriodicYK = 0;
	int PeriodicZK = 1;
	int XElecOnK = 1;
	int YElecOnK = 0;
	int ZElecOnK = 0;
	SNarray snAK = newSNarray( 4, 3, 3);
	//Create Trap Sites one at the left electrode
	//One at the right electrode should not consider
	//any clusters in the sample
	SiteNode snK = getSN(snAK,1,1,0);
	printSNarray(snAK);
	rv = setEnergy(snK, -4);
	assert(rv==0);
	snK = getSN(snAK, 1,1,2);
	rv = setEnergy(snK, -4);
	assert(rv==0);
	printSNarray(snAK);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxK = CalculateAllHops(snAK, 0,0,0, 1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXK, PeriodicYK, PeriodicZK);
	printMatrix(mtxK);
	ArbArray mpAK = MPsort( &OrderLK, &OrderHK, &MidPtsTotalK, mtxK, snAK, PeriodicXK, PeriodicYK, PeriodicZK);
	printf("Order Low %d Order High %d\n",OrderLK, OrderHK);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAK);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbK2 = SortOrderMag( OrderHK-OrderLK+1, OrderLK, mpAK);
	assert(ArbK2!=NULL);
	int TotalOrdersK = OrderHK-OrderLK+1;
	printf("********************************************\n");
	printArbArray(ArbK2, OrderLK);
	printf("********************************************\n");
	ArbArray ArbK = ClusterSort( TotalOrdersK, ArbK2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersK, OrderLK, getElementsUsed(ArbK2));
	assert(ArbK!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbK, OrderLK);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersK, OrderLK, mtxK, &ArbK, snAK, PeriodicXK, PeriodicYK, PeriodicZK, XElecOnK, YElecOnK, ZElecOnK);
	assert(rv==0);
	printf("****************CHECK K***********************\n");
	rv = PrintCheck( TotalOrdersK, OrderLK, ArbK, snAK, mtxK);

	deleteArbArray(&ArbK2);
	//Case Study periodic conditions percolation pathway across sample
	int OrderLJ;
	int OrderHJ;
	int MidPtsTotalJ;
	int PeriodicXJ = 1;
	int PeriodicYJ = 1;
	int PeriodicZJ = 1;
	int XElecOnJ = 1;
	int YElecOnJ = 0;
	int ZElecOnJ = 0;
	SNarray snAJ = newSNarray( 4, 3, 3);
	//Create Trap Sites that cross the sample
	SiteNode snJ = getSN(snAJ,0,1,1);
	printSNarray(snAJ);
	rv = setEnergy(snJ, -4);
	assert(rv==0);
	snJ = getSN(snAJ, 1,1,1);
	rv = setEnergy(snJ, -4);
	assert(rv==0);
	snJ = getSN(snAJ, 2,1,1);
	rv = setEnergy(snJ, -4);
	assert(rv==0);
	snJ = getSN(snAJ, 3, 1, 1);
	rv = setEnergy(snJ, -4);
	assert(rv==0);
	printSNarray(snAJ);
	//Testing non-periodic x, y and z 0 biasX 0 biasY 0 biasZ KT=0.024, and reOrgEnergy =1
	matrix mtxJ = CalculateAllHops(snAJ, 0, 0,0,1, 1, SiteDistance,AttemptToHop,gamma,PeriodicXJ, PeriodicYJ, PeriodicZJ);
	printMatrix(mtxJ);
	ArbArray mpAJ = MPsort( &OrderLJ, &OrderHJ, &MidPtsTotalJ, mtxJ, snAJ, PeriodicXJ, PeriodicYJ, PeriodicZJ);
	printf("Order Low %d Order High %d\n",OrderLJ, OrderHJ);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	printArbArray(mpAJ);
	printf("oooooooooooooooooooooooooooooooooooooooooooo\n");
	ArbArray ArbJ2 = SortOrderMag( OrderHJ-OrderLJ+1, OrderLJ, mpAJ);
	assert(ArbJ2!=NULL);
	int TotalOrdersJ = OrderHJ-OrderLJ+1;
	printf("********************************************\n");
	printArbArray(ArbJ2, OrderLJ);
	printf("********************************************\n");
	ArbArray ArbJ = ClusterSort( TotalOrdersJ, ArbJ2);	
	printf("TotalOrdersR %d OrderLowR %d Elements used by ArbR2 %d\n",TotalOrdersJ, OrderLJ, getElementsUsed(ArbJ2));
	assert(ArbJ!=NULL);
	printf("--------------------------------------------\n");
	printArbArray(ArbJ, OrderLJ);
	printf("--------------------------------------------\n");
	rv = FilterCluster( TotalOrdersJ, OrderLJ, mtxJ, &ArbJ, snAJ, PeriodicXJ, PeriodicYJ, PeriodicZJ, XElecOnJ, YElecOnJ, ZElecOnJ);
	assert(rv==0);
	printf("****************CHECJ***********************\n");
	rv = PrintCheck( TotalOrdersJ, OrderLJ, ArbJ, snAJ, mtxJ);


	printf("****************CHECK J NeighNodes*********************\n");
	rv = CalculateNeighNodes(TotalOrdersJ, &ArbJ, snAJ, PeriodicXJ, PeriodicYJ, PeriodicZJ);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersJ, OrderLJ, ArbJ, snAJ, mtxJ); 

	deleteArbArray(&ArbJ2);

	printf("****************CHECK K NeighNodes*********************\n");
	//SNarray snAK = newSNarray( 4, 3, 3);
	//int PeriodicXK = 0;
	//int PeriodicYK = 0;
	//int PeriodicZK = 1;
	//SiteNode snK = getSN(snAK,1,1,0);
	//snK = getSN(snAK, 1,1,2);
	rv = PrintCheck( TotalOrdersK, OrderLK, ArbK, snAK, mtxK); 
	rv = CalculateNeighNodes(TotalOrdersK, &ArbK, snAK, PeriodicXK, PeriodicYK, PeriodicZK);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersK, OrderLK, ArbK, snAK, mtxK); 

	deleteAllMidPointArray(&mpAK);
	deleteMatrix(&mtxK);
	deleteSNarray(&snAK);
	deleteArbArray(&ArbK);

	printf("****************CHECK M NeighNodes*********************\n");
	//int PeriodicXM = 0;
	//int PeriodicYM = 0;
	//int PeriodicZM = 0;
	//SNarray snAM = newSNarray( 4, 3, 3);
	//SiteNode snM = getSN(snAM,1,0,1);
	//snM = getSN(snAM, 1,2,1);
	rv = PrintCheck( TotalOrdersM, OrderLM, ArbM, snAM, mtxM); 
	rv = CalculateNeighNodes(TotalOrdersM, &ArbM, snAM, PeriodicXM, PeriodicYM, PeriodicZM);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersM, OrderLM, ArbM, snAM, mtxM); 

	printf("****************CHECK N NeighNodes*********************\n");
	//int PeriodicXN = 0;
	//int PeriodicYN = 0;
	//int PeriodicZN = 0;
	//SNarray snAN = newSNarray( 4, 3, 3);
	//SiteNode snN = getSN(snAN,1,1,0);
	//snN = getSN(snAN, 1,1,2);
	rv = PrintCheck( TotalOrdersN, OrderLN, ArbN, snAN, mtxN); 
	rv = CalculateNeighNodes(TotalOrdersN, &ArbN, snAN, PeriodicXN, PeriodicYN, PeriodicZN);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersN, OrderLN, ArbN, snAN, mtxN);

	printf("****************CHECK O NeighNodes*********************\n");
	//int PeriodicXO = 0;
	//int PeriodicYO = 1;
	//int PeriodicZO = 0;
	//SNarray snAO = newSNarray( 4, 3, 3);
	//SiteNode snO = getSN(snAO,1,0,1);
	//snO = getSN(snAO, 1,2,1);
	rv = PrintCheck( TotalOrdersO, OrderLO, ArbO, snAO, mtxO); 
	rv = CalculateNeighNodes(TotalOrdersO, &ArbO, snAO, PeriodicXO, PeriodicYO, PeriodicZO);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersO, OrderLO, ArbO, snAO, mtxO);

	printf("****************CHECK P NeighNodes*********************\n");
	//int PeriodicXP = 1;
	//int PeriodicYP = 0;
	//int PeriodicZP = 0;
	//SNarray snAP = newSNarray( 4, 3, 3);
	//SiteNode snP = getSN(snAP,1,1,1);
	//snP = getSN(snAP, 1,2,1);
	rv = PrintCheck( TotalOrdersP, OrderLP, ArbP, snAP, mtxP); 
	rv = CalculateNeighNodes(TotalOrdersP, &ArbP, snAP, PeriodicXP, PeriodicYP, PeriodicZP);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersP, OrderLP, ArbP, snAP, mtxP);

	printf("****************CHECK T NeighNodes*********************\n");
	//int PeriodicXT = 1;
	//int PeriodicYT = 1;
	//int PeriodicZT = 1;
	//SNarray snAP = newSNarray( 3, 4, 3);
	//SiteNode snP = getSN(snAP,0,1,1);
	//snP = getSN(snAP, 3,1,1);
	rv = PrintCheck( TotalOrdersT, OrderLT, ArbT, snAT, mtxT); 
	rv = CalculateNeighNodes(TotalOrdersT, &ArbT, snAT, PeriodicXT, PeriodicYT, PeriodicZT);
	assert(rv==0);
	rv = PrintCheck( TotalOrdersT, OrderLT, ArbT, snAT, mtxT);

	 printf("****************Calculate Pval J*********************\n");
	 rv = PrintCheck( TotalOrdersJ, OrderLJ, ArbJ, snAJ, mtxJ);
	 rv = CalculateSumAndP(TotalOrdersJ, snAJ, &ArbJ, mtxJ, attempts, PeriodicXJ, PeriodicYJ, PeriodicZJ); 
	 printf("\n** After Run **\n");
	 rv = PrintCheck( TotalOrdersJ, OrderLJ, ArbJ, snAJ, mtxJ);
	 
	printf("****************Calculate Pval M*********************\n");
	//int PeriodicXM = 0;
	//int PeriodicYM = 0;
	//int PeriodicZM = 0;
	//SNarray snAM = newSNarray( 4, 3, 3);
	//SiteNode snM = getSN(snAM,1,0,1);
	//snM = getSN(snAM, 1,2,1);
	rv = PrintCheck( TotalOrdersM, OrderLM, ArbM, snAM, mtxM);
	printMatrix(mtxM);
	rv = CalculateSumAndP(TotalOrdersM, snAM, &ArbM, mtxM, attempts, PeriodicXM, PeriodicYM, PeriodicZM); 
	printf("\n** After Run **\n");
	rv = PrintCheck( TotalOrdersM, OrderLM, ArbM, snAM, mtxM);

	deleteAllMidPointArray(&mpAM);
	deleteMatrix(&mtxM);
	deleteSNarray(&snAM);
	deleteArbArray(&ArbM);

	printf("****************Calculate Pval N*********************\n");
	//int PeriodicXN = 0;
	//int PeriodicYN = 0;
	//int PeriodicZN = 0;
	//SNarray snAN = newSNarray( 4, 3, 3);
	//SiteNode snN = getSN(snAN,1,1,0);
	//snN = getSN(snAN, 1,1,2);
	rv = PrintCheck( TotalOrdersN, OrderLN, ArbN, snAN, mtxN);
	printMatrix(mtxN);
	rv = CalculateSumAndP(TotalOrdersN, snAN, &ArbN, mtxN, attempts, PeriodicXN, PeriodicYN, PeriodicZN); 
	printf("\n** After Run **\n");
	rv = PrintCheck( TotalOrdersN, OrderLN, ArbN, snAN, mtxN);

	deleteAllMidPointArray(&mpAN);
	deleteMatrix(&mtxN);
	deleteSNarray(&snAN);
	deleteArbArray(&ArbN);

	printf("****************Calculate Pval O*********************\n");
	//int PeriodicXO = 0;
	//int PeriodicYO = 1;
	//int PeriodicZO = 0;
	//SNarray snAO = newSNarray( 4, 3, 3);
	//SiteNode snO = getSN(snAO,1,0,1);
	//snO = getSN(snAO, 1,2,1);
	rv = PrintCheck( TotalOrdersO, OrderLO, ArbO, snAO, mtxO);
	printMatrix(mtxO);
	rv = CalculateSumAndP(TotalOrdersO, snAO, &ArbO, mtxO, attempts, PeriodicXO, PeriodicYO, PeriodicZO); 
	printf("\n** After Run **\n");
	rv = PrintCheck( TotalOrdersO, OrderLO, ArbO, snAO, mtxO);
	printf("\nPrintArb Array\n");
	printArbArray(ArbO,OrderLO);
			
	printf("****************Calculate Pval P*********************\n");
	//int PeriodicXP = 1;
	//int PeriodicYP = 0;
	//int PeriodicZP = 0;
	//SNarray snAP = newSNarray( 4, 3, 3);
	//SiteNode snP = getSN(snAP,0,1,1);
	//snP = getSN(snAP, 3,1,1);
	//Periodic in the X
	rv = PrintCheck( TotalOrdersP, OrderLP, ArbP, snAP, mtxP);
	printMatrix(mtxP);
	rv = CalculateSumAndP(TotalOrdersP, snAP, &ArbP, mtxP, attempts, PeriodicXP, PeriodicYP, PeriodicZP); 
	printf("\n** After Run **\n");
	rv = PrintCheck( TotalOrdersP, OrderLP, ArbP, snAP, mtxP);
			
	printf("****************Calculate Pval T*********************\n");
	if(ArbT==NULL){
		printf("Arb T is NULL\n");
	}else if(snAT==NULL){
		printf("snAT is NULL\n");
	}else if(mtxT==NULL){
		printf("mtxT is NULL\n");
	}
	printf("Total Orders %d OrderL %d\n",TotalOrdersT, OrderLT);
	rv = PrintCheck( TotalOrdersT, OrderLT, ArbT, snAT, mtxT);
	printf("Printing Matrix T\n");
	printMatrix(mtxT);
	rv = CalculateSumAndP(TotalOrdersT, snAT, &ArbT, mtxT, attempts, PeriodicXT, PeriodicYT, PeriodicZT); 
	printf("\n** After Run **\n");
	rv = PrintCheck( TotalOrdersT, OrderLT, ArbT, snAT, mtxT);
	printf("PrintArb Array\n");
	printArbArray(ArbT,OrderLT);

	deleteAllMidPointArray(&mpAT);
	deleteMatrix(&mtxT);
	deleteSNarray(&snAT);
	deleteArbArray(&ArbT);


	printf("***************Connect Cluster J Check*****************\n");
	//int PeriodicXJ = 1;
	//int PeriodicYJ = 1;
	//int PeriodicZJ = 1;
	//SNarray snAJ = newSNarray( 4, 3, 3);
	//Create Trap Sites that cross the sample
	//SiteNode snJ = getSN(snAJ,0,1,1);
	//snJ = getSN(snAJ, 1,1,1);
	//snJ = getSN(snAJ, 2,1,1);
	//snJ = getSN(snAJ, 3, 1, 1);
	rv = ConnectClusterSN( TotalOrdersJ, snAJ, ArbJ);
	printSNarray(snAJ);
		 
	printf("***************Connect Cluster O Check*****************\n");
	//int PeriodicXO = 0;
	//int PeriodicYO = 1;
	//int PeriodicZO = 0;
	//SNarray snAO = newSNarray( 4, 3, 3);
	//SiteNode snO = getSN(snAO,1,0,1);
	//snO = getSN(snAO, 1,2,1);
	rv = ConnectClusterSN( TotalOrdersO, snAO, ArbO);
	printSNarray(snAO);

	deleteAllMidPointArray(&mpAO);
	deleteMatrix(&mtxO);
	deleteSNarray(&snAO);
	deleteArbArray(&ArbO);

	printf("***************FindCluster J Check*****************\n");
	//int PeriodicXJ = 1;
	//int PeriodicYJ = 1;
	//int PeriodicZJ = 1;
	//SNarray snAJ = newSNarray( 4, 3, 3);
	//Create Trap Sites that cross the sample
	//SiteNode snJ = getSN(snAJ,0,1,1);
	//snJ = getSN(snAJ, 1,1,1);
	//snJ = getSN(snAJ, 2,1,1);
	//snJ = getSN(snAJ, 3, 1, 1);
	ParameterFrame PF = newParamFrame();
	PFset_AttemptToHop(PF,AttemptToHop);
	PFset_gamma(PF,gamma);
	PFset_Px(PF,PeriodicXJ);
	PFset_Py(PF,PeriodicYJ);
	PFset_Pz(PF,PeriodicZJ);
	PFset_XElecOn(PF,XElecOnJ);
	PFset_YElecOn(PF,YElecOnJ);
	PFset_ZElecOn(PF,ZElecOnJ);
	PFset_Attempts(PF,15);
	PFset_reOrg(PF,1);
	PFset_MovieFrames(PF,10);
	PFset_SiteDist(PF,SiteDistance);

	ArbArray ClArLLJJ;
	rv = FindCluster( &OrderLJ, snAJ, 0,0,0, &ClArLLJJ, 1,PF); 
	printf("\nOld\n");
	printArbArray(ArbJ,OrderLJ);
	printf("\nNew\n");
	printArbArray(ClArLLJJ, OrderLJ);

	deleteAllMidPointArray(&mpAJ);
	deleteMatrix(&mtxJ);
	deleteSNarray(&snAJ);
	deleteArbArray(&ArbJ);
	deleteArbArray(&ClArLLJJ);

	printf("***************FindCluster P Check*****************\n");
	//int PeriodicXP = 1;
	//int PeriodicYP = 0;
	//int PeriodicZP = 0;
	//SNarray snAP = newSNarray( 4, 3, 3);
	//SiteNode snP = getSN(snAP,0,1,1);
	PFset_Px(PF,PeriodicXJ);
	PFset_Py(PF,PeriodicYJ);
	PFset_Pz(PF,PeriodicZJ);
	PFset_XElecOn(PF,XElecOnJ);
	PFset_YElecOn(PF,YElecOnJ);
	PFset_ZElecOn(PF,ZElecOnJ);
	//snP = getSN(snAP, 3,1,1);
	
	char FileName[] = "Data";
	ArbArray ClArLLPP;
	rv = FindCluster( &OrderLP, snAP, 0,0,0, &ClArLLPP, 1, PF); 

	printArbArray(ClArLLPP);
	//rv = SaveCluster( FileName,OrderLP, snAP, 0.0,0.0,0.0, ClArLLPP, 1.0, PF);
	printf("\nOld\n");
	printArbArray(ArbP,OrderLP);
	printf("\nNew\n");
	printArbArray(ClArLLPP, OrderLP);

	deleteAllMidPointArray(&mpAP);
	deleteMatrix(&mtxP);
	deleteSNarray(&snAP);
	deleteArbArray(&ArbP);
	deleteArbArray(&ClArLLPP);

	double electricEnergyPx;
	double electricEnergyPy;
	double electricEnergyPz;
	double KbT;

	LoadCluster( FileName, &OrderLP, &snAP,&electricEnergyPx,&electricEnergyPy,&electricEnergyPz,&ClArLLPP,&KbT,&PF);

	printSNarray_Detailed(snAP);
	deleteSNarray(&snAP);
	deleteArbArray(&ClArLLPP);
	deleteParamFrame(&PF);	


	int OrderLB;
	int OrderHB;
	int MidPtsTotalB;
	int PeriodicXB = 1;
	int PeriodicYB = 1;
	int PeriodicZB = 1;
	int XElecOnB = 1;
	int YElecOnB = 0;
	int ZElecOnB = 0;
	SNarray snAB = newSNarray( 3,4,3);

	//Create an artificial trap in snAB which in increment length is
	//length 0-2 width 0-3 and height 0-2
	SiteNode snB = getSN(snAB,1,1,1);
	rv = setEnergy(snB, -3);
	snB = getSN(snAB,1,2,1);
	rv = setEnergy(snB, -3);

	//Testing Periodic x,y and z 0 biasX 0 biasY 0 biasZ KT=1 and reOrgEnergy=1
	matrix mtxB = CalculateAllHops(snAB, 0,0,0, 1,1,SiteDistance, AttemptToHop, gamma, PeriodicXB, PeriodicYB, PeriodicZB);

	//Sorts all the data into midpoints
	ArbArray mpAB = MPsort( &OrderLB, &OrderHB, &MidPtsTotalB, mtxB, snAB, PeriodicXB, PeriodicYB, PeriodicZB);

	//Sorts midpoints into Nodes in link lists. Each element of mpAT is composed of a
	//link list all of the same order of magnitude
	int TotalOrdersB = OrderHB-OrderLB+1;
	ArbArray ArbB2 = SortOrderMag( TotalOrdersB, OrderLB, mpAB);

	//Sorts nodes into clusters
	printf("Sorting Nodes into Clusters\n");
	ArbArray ArbB = ClusterSort( TotalOrdersB, ArbB2);

	//Filtering Cluster
	printf("Filtering Clusters\n");
	rv = FilterCluster( TotalOrdersB, OrderLB, mtxB, &ArbB, snAB,\
			PeriodicXB, PeriodicYB, PeriodicZB,\
			XElecOnB, YElecOnB, ZElecOnB);

	printf("Calculating NeighNodes\n");
	CalculateNeighNodes(TotalOrdersB, &ArbB, snAB, PeriodicXB, PeriodicYB, PeriodicZB);

	printf("Calculating P values\n");
	CalculateSumAndP( TotalOrdersB, snAB, &ArbB, mtxB, 15, PeriodicXB, PeriodicYB, PeriodicZB);

	printf("Connecting nodes to clusters\n");
	ConnectClusterSN( TotalOrdersB, snAB, ArbB);

	//Deleting Every thing
	printf("Deleting Phase\n");
	deleteSNarray(&snAB);
	deleteAllMidPointArray(&mpAB);
	deleteMatrix(&mtxB);
	deleteArbArray(&ArbB2);
	deleteArbArray(&ArbB);

  printf("COMPLETED test_cluster success\n");
	return 0;
}
