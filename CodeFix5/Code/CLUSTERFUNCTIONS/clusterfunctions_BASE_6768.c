#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <assert.h>

#include "MONTECARLO/montecarlo.h"
#include "SITENODE/sitenode.h"
#include "MATRIX/matrix.h"
#include "LINKLIST/linklist.h"
#include "DATASTRUCT/cluster.h"
#include "clusterfunctions.h"
#include "../../PARAMETERS/read.h"

matrix CalculateAllHops(const_SNarray snA,const double electricEnergyX, \
		const double electricEnergyY, const double electricEnergyZ, \
		const double KT,const double reOrgEnergy,const double SiteDistance, \
		const double AttemptToHop, const double gamma,\
		const int PeriodicX,const int PeriodicY,const int PeriodicZ){


	if(snA==NULL || KT<0 || PeriodicX<0 || PeriodicX>1 || PeriodicY<0 || PeriodicY>1 || \
			PeriodicZ<0 || PeriodicZ>1){
		return NULL;
	}

	int i, j, k, l ;
	int i1, j1, k1;
	int rv;

	//Boltzmann constant Units of [eV/K]
	static const double kB = 8.6173324E-5;
	//Planck constant Units of [eV s]
	static const double hbar = 6.58211928E-16;
	//printf("Value of fracSeed %e.\n",fracSeed);

	double MarcusJ0;
	double MarcusCoeff;

	//one node have 6 hopping rate for 6 neighbor node
	double v[6];
	matrix MasterM = newMatrix(getAtotal(snA),12);

	//Calculating Marcus J0 coefficient assuming the Attempt to hop Rate
	//is equivalent to the marcus coefficient at 300 K
	MarcusJ0 = pow( AttemptToHop*hbar*pow(4*reOrgEnergy*kB*300/M_PI,1/2),1/2);

	//Convert SiteDistance to [nm] from [m]
	//double SiteDistanceNM = SiteDistance*1E9;

	//Calculating full Marcus Coefficient;

	//printf("MarcusJ0 %g hbar %g gamma %g SiteDistance %g\n",MarcusJ0,hbar,gamma,SiteDistance);
	MarcusCoeff = pow(MarcusJ0,2)/hbar * pow(M_PI/(4*reOrgEnergy*KT),1/2)*exp(2*gamma*SiteDistance);

	//MasterM values
	//1 - Order of Magnitude of hop behind site i,j,k
	//2 - Order of Magnitude of hop infornt of site i,j,k
	//3 - Order of Magnitude of hop to site left of site i,j,k
	//4 - Order of Magnitude of hop to site right of site i,j,k
	//5 - Order of Magnitude of hop to site below site i,j,k
	//6 - Order of Magnitude of hop to site above site i,j,k
	//7 - Hop rate to site behind i,j,k
	//8 - Hop rate to site infront of i,j,k
	//9 - Hop rate to site left of site i,j,k
	//10 - Hop rate to site right of i,j,k
	//11 - Hop rate to site below i,j,k
	//12 - Hop rate to site above i,j,k

	SiteNode SNi;
	SiteNode SNj;

	for(i = 0; i < getAlen(snA); i++){
		for(j = 0; j < getAwid(snA); j++){
			for(k = 0; k < getAhei(snA); k++){

				SNi = getSN(snA, i, j, k);

				//Calculate the hop rates for all the sites in the 
				//System. If hopping off of SNi
				//Site behind site SNi parallel x axis
				i1=(i-1+getAlen(snA))%getAlen(snA);
				SNj = getSN(snA, i1,j,k);
				v[0] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyX, KT, reOrgEnergy);

				//printf("Value of rage v0 %g MarcusCoeff %g Energy SNj %g Energy SNi %g electricEnergyX %g KT %g reOrgEnergy %g\n",v[0],MarcusCoeff,getEnergy(SNj), getEnergy(SNi), electricEnergyX, KT, reOrgEnergy);

				//Site in front of site SNi
				i1=(i+1)%getAlen(snA);
				SNj = getSN(snA,i1,j,k);
				v[1] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyX, KT, reOrgEnergy);
				//Site left of site SNi parallel y axis
				j1=(j-1+getAwid(snA))%getAwid(snA);
				SNj = getSN(snA,i, j1,k);
				v[2] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyY, KT, reOrgEnergy);
				//Site right of site SNi
				j1=(j+1)%getAwid(snA);
				SNj = getSN(snA,i, j1,k);
				v[3] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyY, KT, reOrgEnergy);
				//Site below site SNi parallel z axis
				k1=(k-1+getAhei(snA))%getAhei(snA);
				SNj = getSN(snA,i, j, k1);
				v[4] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) + electricEnergyZ, KT, reOrgEnergy);
				//Site above site SNi
				k1=(k+1)%getAhei(snA);
				SNj = getSN(snA,i, j, k1);
				v[5] = MarcusCoeff*hoppingRate(getEnergy(SNj) - getEnergy(SNi) - electricEnergyZ, KT, reOrgEnergy);

				if(PeriodicX==0){
					if( i == 0)
						v[0] = 0;
					if( i == getAlen(snA)-1)
						v[1] = 0;
				}

				if(PeriodicY==0){
					if( j == 0)
						v[2] = 0;
					if( j == getAwid(snA)-1)
						v[3] = 0;
				}


				if(PeriodicZ==0){
					if( k == 0 )
						v[4] = 0;
					if( k == getAhei(snA)-1 )
						v[5] = 0;
				}

				for(l=0;l<6;l++){

					//printf("Index %d v[%d] %g log10(v[%d]) %g  round(log10(v[%d])) %g\n",getIndex(snA,i,j,k),l,v[l],l,log10(v[l]),l,round(log10(v[l])));
					rv=setE(MasterM,getIndex(snA,i,j,k)+1,(l+1),round(log10(v[l])));
					//Also record the actual rates
					rv=setE(MasterM,getIndex(snA,i,j,k)+1,6+(l+1), v[l]);
				}
			}
		}
	}

	return MasterM;
}

ArbArray MPsort(int * orderL, int * orderH, int * MidPtsTotal, const_matrix MasterM, const_SNarray snA,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ){

	if( orderL==NULL || orderH==NULL || MidPtsTotal==NULL || MasterM==NULL || snA==NULL){
		return NULL;
	}

	if( PeriodicX>1 || PeriodicX<0 || PeriodicY>1 || PeriodicY<0 || PeriodicZ>1 || PeriodicZ<0){
		return NULL;
	}

	if(getRows(MasterM)!=getAtotal(snA)|| getCols(MasterM)!=12){
		return NULL;
	}
	double x1;
	double x2;
	double y1;
	double y2;
	double z1;
	double z2;
	MidPoint mp;

	*MidPtsTotal=0;
	int nei1;
	int nei2;
	int order;
	int orderLow=0;
	int orderHigh=0;

	int i, j, k;

	int NumE;

	//Number of Midpoints if completely non periodic
	NumE = getAlen(snA)*getAwid(snA)*getAhei(snA)*3; 

	if (PeriodicX==0){
		NumE -= (getAwid(snA))*(getAhei(snA));
	}
	if (PeriodicY==0){
		NumE -= (getAlen(snA))*(getAhei(snA));
	}	
	if (PeriodicZ==0){
		NumE -= (getAlen(snA))*(getAwid(snA));
	}

	//printf("NumE %d\n",NumE);
	//The 2 represents that the mpA array stores MidPoints
	ArbArray mpA = newArbArray(NumE,2);

	int Mid_ID=0;

	//printf("orderLow %d orderHigh %d\n",orderLow,orderHigh);

	int lenX = getAlen(snA);
	int lenY = getAwid(snA);
	int lenZ = getAhei(snA);

	for(i = 0; i < lenX; i++){
		for(j = 0; j < lenY; j++){
			for(k = 0; k < lenZ; k++){
				nei1=getIndex(snA,i,j,k);

				//Check to see if periodic in the X
				if((PeriodicX==0 && i<(lenX-1)) || (PeriodicX==1 )){

					nei2=getIndexPeriodic(snA,i+1,j,k);
					//2 - infront i (nei1)
					x1  = (int)getE(MasterM, nei1+1, 2);
					//1 - behind i+1 (nei2);
					x2 = (int)getE(MasterM, nei2+1, 1);
					//printf("x1 %g x2 %g\n",x1,x2);
					if (x1==x2) {

						order=x1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;

					} else if (x1>x2) {
						//If the hop rates are different
						//The order falls to the lower hop rate
						//Or the slower time
						order=x2;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;

					} else {

						order=x1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;

					}
				}

				//Check to see if periodic in the Y
				if((PeriodicY==0 && j<(lenY-1)) || (PeriodicY==1 )){

					nei2=getIndexPeriodic(snA,i,j+1,k);
					//4 - right j (nei1)
					y1 = (int)getE(MasterM,nei1+1,4);
					//3 - left j+1 (nei2)
					y2 = (int)getE(MasterM,nei2+1,3);
					//printf("y1 %g y2 %g\n",y1,y2);

					if (y1==y2) {

						order=y1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;
					} else if (y1>y2){

						order=y2;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;
					} else {

						order=y1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;
					}
				}

				//Check to see if periodic in the Z
				if((PeriodicZ==0 && k<(lenZ-1)) || (PeriodicZ==1 )) {

					nei2=getIndexPeriodic(snA,i,j,k+1);
					//6 - site above k (nei1)
					z1 = (int)getE(MasterM,nei1+1,6);
					//5 - site below k+1 (nei2)
					z2 = (int)getE(MasterM,nei2+1,5);
					//printf("z1 %g z2 %g\n",z1,z2);
					if (z1==z2) {

						order=z1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;
					} else if (z1>z2) {

						order=z2;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);

						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;

					} else {

						order=z1;
						mp = newMidPoint(order, Mid_ID, nei1, nei2);
						setArbElement( mpA, Mid_ID, (void *)mp);
						if (Mid_ID==0){
							orderLow=order;
							orderHigh=order;
						}else{
							if(orderHigh<order){
								orderHigh=order;
							}else if (orderLow>order){
								orderLow=order;
							}
						}
						Mid_ID++;
					}
				}
			}
		}
	}
	//printf("orderLow %d orderHigh %d\n",orderLow,orderHigh);

	*orderL=orderLow;
	*orderH=orderHigh;
	*MidPtsTotal=Mid_ID;
	return mpA;
}

ArbArray SortOrderMag(const int TotalOrders,const int orderLow, const_ArbArray mpA){

	if(TotalOrders<1 || mpA==NULL){
		return NULL;
	}
	//This array contains all the link lists
	ArbArray ArLL=newArbArray(TotalOrders, 0);
	//Need to set default values for all the OrderMagLL in
	//the ArLL, where all OMLL have values:
	//  OMLL->orderMag=orderMag;
	//  OMLL->size=0;
	//  OMLL->start=NULL;

	//int tempOrder;
	int element, rv, tempOrder, order, Mid_ID;

	//printf("Setting Default values to LL\n");
	for(element=0;element<TotalOrders;element++){
		order=element+orderLow;
		//printf("OrderLow %d Order %d\n",orderLow,order);
		setDefaultArbElem(ArLL,element,order);
	}
	//printf("Sorting into LL\n");
	//cycle through mid points and assign to correct link list
	for(Mid_ID=0;Mid_ID<getElementsUsed(mpA);Mid_ID++){
		tempOrder=getMPOrder(mpA, Mid_ID);
		element=tempOrder-orderLow;
		//Check to see if a LL already exists for that order of Magnitude
		//return value of 0 means it was succesfully added 
		//return value of -1 means it was already added this should not happen

		if (element>=TotalOrders){
			//This is the case that we did not make the 
			//ArbArray big enough
			deleteArbArray(&ArLL);
			return NULL;
		}

		rv = addToOrLL(ArLL, element, getMP(mpA,Mid_ID) );
	}

//	printArbArray(ArLL);

	return ArLL;
}

ArbArray ClusterSort(const int TotalOrders,const int orderLow, const_ArbArray ArLL){

	//The totalOrders corresponds to the difference between the highest and lowest
	//order of magnitude in the sample. It can happen that the there is a jump
	//so that not all the orders is represented between the highest and lowest
	//points. 
	if(TotalOrders<=0 || ArLL==NULL || TotalOrders<getElementsUsed(ArLL)){
		return NULL;
	}

	OrderMagLL TempOMLL;
	int element = 0;
	int Lo;
	int Hi;
	TempOMLL = getOrderLL(ArLL,element);
	Lo = getOMLL_order(TempOMLL); 
	while(TempOMLL !=NULL){
		element++;
		TempOMLL = getOrderLL(ArLL,element);
	}
	TempOMLL = getOrderLL(ArLL,element-1);
	Hi = getOMLL_order(TempOMLL);
	if(TotalOrders!=(Hi-Lo+1)){
		//printf("Value Hi %d Lo %d diff %d\n",Hi,Lo,(Hi-Lo+1));
		return NULL;
	}

	ArbArray ClArLL = newArbArray(TotalOrders, 1);

	if(ClArLL==NULL){
		return NULL;
	}

	//printf("Successfully Created Cluster Array Link List.\n");

	MidPoint tempmp;
	int ID;

	//printf("\nOrganizing Nodes into clusters based on proximity.\n");


	for(element=0;element<TotalOrders;element++){
		//printf("Element %d TotalOrders %d\n",element, TotalOrders);

		ID=1;
		//Need to cycle through the mid points in each Link List in ArLL
		TempOMLL = getOrderLL(ArLL, element);
		tempmp = getOMLLstartMP(TempOMLL);

		if(tempmp!=NULL) {
			//Create ClusterLL id is set to one 
			ClusterLL clLL = newClusterLL(ID);
			setArbElement(ClArLL, element, (void *) clLL);

			while(tempmp!=NULL){
				addNodeToCluster(clLL , tempmp);
				tempmp=getNextMP(tempmp);
			}
		}
	}

	//printf("Printing Cluster Array.\n");
	//printArbArray(ClArLL, orderLow);

	return ClArLL;
}

//MasterM Should span [Number of sites][12] 
//The total number of sites should be stored in snA

int FilterCluster(const int TotalOrders,const int orderLow,const_matrix MasterM,\
		ArbArray * ClArLL, const_SNarray snA, const int PeriodicX, \
		const int PeriodicY, const int PeriodicZ, const int XElecOn,\
		const int YElecOn, const int ZElecOn){

	if (TotalOrders<=0 || MasterM==NULL || ClArLL==NULL || snA == NULL){
		return -1;
	}

	if (PeriodicX<0 || PeriodicY<0 || PeriodicZ<0 || PeriodicX>1 || PeriodicY>1 || PeriodicZ>1 || *ClArLL==NULL){
		return -1;
	}

	//Compare Total Orders with Elements Reserved not used because there may not be
	//nodes in a certain orderof magnitude
	if( getRows(MasterM)!=getAtotal(snA) || getElementsReserved(*ClArLL)!=TotalOrders){
		return -1;
	}

	//This is done by first checking all the hop rates off the cluster
	//Take all hops that are within an order of magnitude
	int Node_ID;
	int Node_IDFro, Node_IDBeh, Node_IDLef, Node_IDRig, Node_IDAbo, Node_IDBel;
	int flagSum;
	int diff;
	int DeleteCluster;
	int element;
	int order;
	int i, j, k;
	ClusterLL TempClLL, holder, brief;
	Node tempNode, tempNode2;
	int count;
	int rv;

	//Verify that hops exist that are within two orders of magnitude 
	for(element=0;element<TotalOrders;element++){

		//Get ClusterLL
		TempClLL=(ClusterLL) getArbElement(*ClArLL,element);
		order=element+orderLow;
		holder=NULL;
		//printf("Order %d\n",order);
		count=0;

		while(TempClLL!=NULL){// && getNode(TempClLL)!=NULL){

			count++;
			//printf("\nCluster id %d cluster Num %d\n",getCluster_id(TempClLL),count);
			//Cycle Clusters
			tempNode = getStartNode(TempClLL);

			DeleteCluster=0;
			while(tempNode!=NULL){
				//Cycle Nodes
				//printf("TempNode Not NULL\n");
				Node_ID=getNode_id(tempNode);
				getLoc( &i, &j, &k, Node_ID, snA);

				//Must be careful when considering periodicity in 
				//the x direction sites may act like a cluster if 
				//periodic and they are clumped at the edges of 
				//the simulation box. As soon as the electrode is
				//attached the cluster no longer acts like a
				//cluster
				Node_IDFro=getIndFroP(snA, i, j, k);
				Node_IDBeh=getIndBehP(snA, i, j, k);
				Node_IDLef=getIndLefP(snA, i, j, k);
				Node_IDRig=getIndRigP(snA, i, j, k);
				Node_IDAbo=getIndAboP(snA, i, j, k);
				Node_IDBel=getIndBelP(snA, i, j, k);
				tempNode2=getNextNode(tempNode);
				//Cycle through Cluster and find which neighboring sites are within the cluster
				while (tempNode2!=NULL){

					//If the charge is hopping to the end electrode we do not want
					//to exclude it. Clusters do not exist if it is easy for the 
					//charge to hop to an electrode at either end, but if it is
					//periodic in the X direction the cluster does need to be included

					//In the code below determing which nodes with respect to x, y and z
					//or tempNode are already in the cluster. 
					//printf("Fro %d Beh %d NodeID %d\n",Node_IDFro,Node_IDBeh,getNode_id(tempNode2));

					if(Node_IDFro==getNode_id(tempNode2) ){
						setFlagFro(tempNode);
						setFlagBeh(tempNode2);
					} else if (Node_IDBeh==getNode_id(tempNode2) ){
						setFlagBeh(tempNode);
						setFlagFro(tempNode2);
					}	else if (Node_IDLef==getNode_id(tempNode2) ){
						setFlagLef(tempNode);
						setFlagRig(tempNode2);
					} else if (Node_IDRig==getNode_id(tempNode2) ){
						setFlagRig(tempNode);
						setFlagLef(tempNode2);
					} else if (Node_IDAbo==getNode_id(tempNode2) ){
						setFlagAbo(tempNode);
						setFlagBel(tempNode2);
					} else if (Node_IDBel==getNode_id(tempNode2) ){
						setFlagBel(tempNode);
						setFlagAbo(tempNode2);
					}

					tempNode2=getNextNode(tempNode2);
				} //Compare nodes within the cluster

				flagSum= getFlagFro(tempNode)+getFlagBeh(tempNode)+\
								 getFlagLef(tempNode)+getFlagRig(tempNode)+\
								 getFlagBel(tempNode)+getFlagAbo(tempNode);
				//This is in the case that all the neighbors are not in the cluster
				//MaterM[ blah blah ][l]
				//l-1 hop behind
				//l-2 hop infront
				//l-3 hop left
				//l-4 hop right
				//l-5 hop below
				//l-6 hop above

				//If flagSum is 6 it means that all the neighboring
				//sites are within the same cluster and nothing is
				//done. If however some are not within the cluster
				//we must check their hoprates.
				if(flagSum!=6){

					//If the site is next to a boundary that is not periodic and
					//is next to an electrode the cluster is terminated
					if (PeriodicX==0 && XElecOn==1 && (i==(getAlen(snA)-1) || i==0)){
						DeleteCluster=1;
						break;
					}else{

						if(getFlagFro(tempNode)==0){
							//Hop infront is not in cluster		
							//If the sample is periodic then it means the cluster
							//is confined and would not be deleted

							if(PeriodicX==1 || i!=(getAlen(snA)-1)){
								//If the sample is either periodic or the charge is

								diff=(order-getE(MasterM,getIndex(snA,i,j,k)+1,2));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									//printf("i %d j %d k %d Part 2\n",i,j,k);
									break;
								}	//not on the boundary then we need to check
							}
						}

						if(getFlagBeh(tempNode)==0){

							if(PeriodicX==1 || i!=0){
								diff=(order-getE(MasterM,getIndex(snA,i,j,k)+1,1));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									break;
								}
							}
						}
					}

					//If the site is next to a boundary that is not periodic and
					//is next to an electrode the cluster is terminated
					if (PeriodicY==0 && YElecOn==1 && (j==(getAwid(snA)-1) || j==0)){
						DeleteCluster=1;
						break;
					}else{
						if(getFlagLef(tempNode)==0){
							if(PeriodicY==1 || j!=0){

								diff=(order-getE(MasterM, getIndex(snA,i,j,k)+1,3));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									break;
								}
							}
						}

						if(getFlagRig(tempNode)==0){
							if(PeriodicY==1 || j!=(getAwid(snA)-1)){
								diff=(order-getE(MasterM,getIndex(snA,i,j,k)+1,4));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									//printf("Order of element to Rig %g 5\n",getE(MasterM,getIndex(snA,i,j,k)+1,4));
									//printf("ID %d\n",getIndexRightPeriodic(snA,getIndex(snA,i,j,k)));
									break;
								}	
							}
						}
					}

					//If the site is next to a boundary that is not periodic and
					//is next to an electrode the cluster is terminated
					if (PeriodicZ==0 && ZElecOn==1 && (k==(getAhei(snA)-1) || k==0)){
						DeleteCluster=1;
						break;
					}else{
						if(getFlagBel(tempNode)==0){
							if(PeriodicZ==1 || k!=0){
								diff=(order-getE(MasterM, getIndex(snA,i,j,k)+1, 5));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									break;
								}
							}
						}

						if(getFlagAbo(tempNode)==0){
							if(PeriodicZ==1 || k!=(getAhei(snA)-1)){
								diff=(order-getE(MasterM,getIndex(snA,i,j,k)+1,6));
								if(diff<1){
									//This would mean that this is not a cluster
									DeleteCluster=1;
									break;
								}
							}
						}
					}
				}
				//If Hit DeleteCluster do not need to check the rest of
				//the nodes in the cluster
				tempNode=getNextNode(tempNode);

			}//Cyle through nodes while tempNode!=NULL

			if (DeleteCluster==0){
				//Check to see if there is a percolation pathway
				//If there is will delete the cluster
				//printf("***************BEGINNING PERCOLATION CHECK***********\n");
				DeleteCluster = FeelPercolation( TempClLL, snA, PeriodicY, PeriodicZ );
			}

			if(DeleteCluster==1){
				//If it is the first clLL in the ArbArray element
				if(holder==NULL){
					//Bypass tempClLL and set first element of the ArbArray to the 
					//next clLL


					if (getNextClusterLL(TempClLL)==NULL){
						NullArbElement(*ClArLL, element);
					}else{
						rv = setArbElement(*ClArLL, element, getNextClusterLL(TempClLL));	
						//Delete TempClLL Cluster and all nodes within the cluster
					}
					//printf("Deleting Cluster\n");
					rv = deleteClusterLLNodes( TempClLL );

					//Start over with TempClLL reinitialized

					TempClLL=(ClusterLL) getArbElement(*ClArLL,element);

				} else {
					//This is if the cluster is not the first 
					//Bypass tempClLL
					brief=(ClusterLL) getNextClusterLL(TempClLL);
					//Delete TempClLL and all nodes within the cluster
					//printf("Deleting Cluster\n");
					deleteClusterLLNodes( TempClLL );
					//Reinitialize TempClLL
					TempClLL=brief;
					//Connect previous CllLL (holder) with new TempClLL
					setNextClusterLL(holder,TempClLL);
				}

				//Reinitialize DeleteCluster

			}else {

				//Reinitialize holder
				holder=TempClLL;
				//Before we can delete the cluster must ensure that
				//the pointers correctly connect

				//Move to next cluster
				TempClLL=getNextClusterLL(TempClLL);

			}
		//While TempClLL!=NULL
		}
		
	} //Next element in ArrayLL

	return 0;
}

int FeelPercolation(ClusterLL ClLL, const_SNarray snA, int PeriodicY, int PeriodicZ ) {

	if(ClLL==NULL || snA==NULL){
		return -1;
	}
	if(PeriodicY<0 || PeriodicY>1 || PeriodicZ<0 || PeriodicZ>1){
		return -1;
	}
	//This function will return a 1 if a percolation pathway exist
	//It will return a 0 if no percolation pathway exists and thus can be 
	//treated as a cluster. 
	//Will return -1 if there is a problem with the inputs

	int i, j, k;

	//This matrix stores the IDs of the nodes that are on the 
	//Right side of the sample
	matrix NearRightElec = newMatrix(1,1);
	int rowR=1;

	//This matrix stores the IDs of the nodes that are on the
	//left side of the sample
	matrix NearLeftElec = newMatrix(1,1);
	int rowL=1;

	//This matrix stores the nodes that are within the cluster
	//in matrix form
	matrix AvailableNodes = newMatrix(getCluster_numNodes(ClLL),1);

	int inc=1;
	Node tempNode = getStartNode(ClLL);
	int N_ID = getNode_id(tempNode);

	while(tempNode!=NULL){

		N_ID = getNode_id(tempNode);
		getLoc( &i, &j, &k,N_ID, snA);
		//Storing all the available nodes in this matrix
		setE(AvailableNodes,inc,1,N_ID);
		inc++;

		//Storing all the IDs of the nodes that are 
		//next to the Left and right electrodes
		if (i==0){
			resizeRow(&NearLeftElec,rowL);
			setE(NearLeftElec, rowL,1, (double) N_ID);
			rowL++;
		}else if(i==(getAlen(snA)-1)){
			resizeRow(&NearRightElec,rowR);
			setE(NearRightElec, rowR,1, (double) N_ID);
			rowR++;
		}

		tempNode=getNextNode(tempNode);

	}

//	printf("Set Matrices\n");
//	printf("Printing LeftELec\n");
//	printMatrix(NearLeftElec);
//	printf("Printing RightElec\n");
//	printMatrix(NearRightElec);
//	printf("Printing Available Nodes Matrix\n");
//	printMatrix(AvailableNodes);

	//Looking at the boundaries determining if cluster
	//spans the length of the sample
	int ii, jj;
	int i1, j1, k1;
	int i2, j2, k2;
	int N_ID1, N_ID2;
	int IsCluster = 1;

	ii=1;
	while (ii<=getRows(NearLeftElec) && IsCluster==1 ){
		N_ID1 = (int)getE(NearLeftElec, ii,1);
		getLoc( &i1, &j1, &k1, N_ID1, snA);

		jj=1;
		while (jj<=getRows(NearRightElec) && IsCluster==1){

			N_ID2 = (int) getE(NearRightElec, jj, 1);
			getLoc( &i2, &j2, &k2, N_ID2, snA);

			//This means we have to do a more detailed
			//analysis there could be a percolation
			//pathway
			if(j2==j1 && k1==k2){
				//We aren't sure it's a cluster
				IsCluster=0;
			}
			jj++;
		}
		ii++;
	}

	//printf("Value of IsCluster %d\n",IsCluster);
	int StartID;
	linklist LL = newBlankLinkList();
	int Lef, Rig, Abo, Bel, Fro, Beh;
	int RElec;
	int StartRow;
	int AvailableID;
	//Using non periodic conditions in the x examining whether
	//the cluster spans the length of the sample 
	if (IsCluster!=0){
		deleteLL(LL);
		deleteMatrix(AvailableNodes);
		deleteMatrix(NearLeftElec);
		deleteMatrix(NearRightElec);
		return 0;
		//Here we are sure it is a cluster
	}else{

		if(PeriodicY==1 && PeriodicZ==1){

			for(StartRow=1;StartRow<=getRows(NearLeftElec);StartRow++){
				//Start at the left Electrode
				StartID = (int)getE(NearLeftElec,StartRow,1);

				//printf("Value of StartID %d\n",StartID);
				//printf("Value of MatchExist %d for ID %d\n",matchExist(AvailableNodes,StartID),StartID);

				if (matchExist(AvailableNodes,StartID)==1){
					//Cross the first ID out of the available ids
					matchReplace(AvailableNodes,StartID,-1);
					//printf("Value of MatchExist after match replace %d\n",matchExist(AvailableNodes,StartID));
					//Grab surrounding Neighbors that are within the cluster
					tempNode = getCluster_Node(ClLL,StartID);
					if(getFlagLef(tempNode)==1){
						//Node is within cluster add to LL
						Lef = getIndexLeftPeriodic(snA,StartID);
						if (matchExist(AvailableNodes, Lef)==1){
							addLLNode(LL,Lef);
						}
					}
					if(getFlagRig(tempNode)==1){
						//Node is within cluster add to LL
						Rig = getIndexRightPeriodic(snA,StartID);
						if (matchExist(AvailableNodes, Rig)==1){
							addLLNode(LL,Rig);
						}
					}
					if( getFlagAbo(tempNode)==1){
						//Node is within cluster add to LL
						Abo = getIndexAbovePeriodic(snA, StartID);
						if (matchExist(AvailableNodes, Abo)==1){
							addLLNode(LL,Abo);
						}
					}
					if(getFlagBel(tempNode)==1){
						//Node is within cluster add to LL
						Bel = getIndexBelowPeriodic(snA, StartID);
						if (matchExist(AvailableNodes, Bel)==1){
							addLLNode(LL,Bel);
						}
					}
					if(getFlagFro(tempNode)==1){
						//Node is within cluster
						Fro = getIndexFront(snA, StartID);
						if (Fro!=-1){
							if (matchExist(AvailableNodes, Fro)==1){
								addLLNode(LL,Fro);
							}
						}
					}
					if(getFlagBeh(tempNode)==1){
						//Node is within cluster
						Beh = getIndexBehind(snA, StartID);
						if (Beh!=-1){
							if (matchExist(AvailableNodes, Beh)==1){
								addLLNode(LL,Beh);
							}
						}
					}
					//printf("Printing linklist\n");
					//printLL(LL);
					//Check to see if the added nodes are attached to the right
					//Electrode if any nodes are the cluster is either
					//to big to use the approximation or it is a percolating
					//pathway
					for(jj=1;jj<=getRows(NearRightElec);jj++){

						RElec = (int)getE(NearRightElec,jj,1);

						if(getFlagLef(tempNode)==1){
							if (RElec==Lef){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagRig(tempNode)==1){
							if (RElec==Rig){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if( getFlagAbo(tempNode)==1){
							if (RElec==Abo){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBel(tempNode)==1){
							if (RElec==Bel){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}

						if(getFlagFro(tempNode)==1){
							if (RElec==Fro){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBeh(tempNode)==1){
							if( RElec==Beh){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}


					}

					//Cycle through the link list

					while (getLLlength(LL)!=0 ){

						//printf("Length of LL %d\n",getLLlength(LL));
						AvailableID = getLLstartID(LL); 
						tempNode = getCluster_Node(ClLL,AvailableID);
						//printf("Value of AvailableID %d\n",AvailableID);
						//Cross the first ID out of the available ids
						//And also cut it out of the linklist
						matchReplace(AvailableNodes,AvailableID,-1);
						removeLLNode(LL, AvailableID);

						//Check the neighbors of this node
						//Grab surrounding Neighbors that are within the cluster
						tempNode = getCluster_Node(ClLL,AvailableID);
						if(getFlagLef(tempNode)==1){
							//Node is within cluster add to LL
							Lef = getIndexLeftPeriodic(snA,AvailableID);
							if (matchExist(AvailableNodes, Lef)==1){
								addLLNode(LL,Lef);
							}
						}
						if(getFlagRig(tempNode)==1){
							//Node is within cluster add to LL
							Rig = getIndexRightPeriodic(snA,AvailableID);
							if (matchExist(AvailableNodes, Rig)==1){
								addLLNode(LL,Rig);
							}
						}
						if( getFlagAbo(tempNode)==1){
							//Node is within cluster add to LL
							Abo = getIndexAbovePeriodic(snA, AvailableID);
							if (matchExist(AvailableNodes, Abo)==1){
								addLLNode(LL,Abo);
							}
						}
						if(getFlagBel(tempNode)==1){
							//Node is within cluster add to LL
							Bel = getIndexBelowPeriodic(snA, AvailableID);
							if (matchExist(AvailableNodes, Bel)==1){
								addLLNode(LL,Bel);
							}
						}
						if(getFlagFro(tempNode)==1){
							//Node is within cluster
							Fro = getIndexFront(snA, AvailableID);
							if (Fro!=-1){
								if (matchExist(AvailableNodes, Fro)==1){
									addLLNode(LL,Fro);
								}
							}
						}
						if(getFlagBeh(tempNode)==1){
							//Node is within cluster
							Beh = getIndexBehind(snA, AvailableID);
							if (Beh!=-1){
								if (matchExist(AvailableNodes, Beh)==1){
									addLLNode(LL,Beh);
								}
							}
						}

						//Check to see if the added nodes are attached to the right
						//Electrode if any nodes are the cluster is either
						//to big to use the approximation or it is a percolating
						//pathway
						for(jj=1;jj<=getRows(NearRightElec);jj++){


							RElec = (int)getE(NearRightElec,jj,1);
							//printf("Value of RElec %d\n",RElec);
							//printf("Value of Lef %d Rig %d Abo %d Bel %d\n",Lef,Rig,Abo,Bel);
							//printLL(LL);

							if(getFlagLef(tempNode)==1){
								if (RElec==Lef){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagRig(tempNode)==1){
								if (RElec==Rig){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if( getFlagAbo(tempNode)==1){
								if (RElec==Abo){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBel(tempNode)==1){
								if (RElec==Bel){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

							if(getFlagFro(tempNode)==1){
								if (RElec==Fro){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBeh(tempNode)==1){
								if( RElec==Beh){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

						}

					}
				}
			}

		}else if(PeriodicY==1 && PeriodicZ==0){

			for(StartRow=1;StartRow<=getRows(NearLeftElec);StartRow++){
				//Start at the left Electrode
				StartID = (int)getE(NearLeftElec,StartRow,1);

				if (matchExist(AvailableNodes,StartID)==1){
					//Cross the first ID out of the available ids
					matchReplace(AvailableNodes,StartID,-1);

					//Grab surrounding Neighbors that are within the cluster
					tempNode = getCluster_Node(ClLL,StartID);
					if(getFlagLef(tempNode)==1){
						//Node is within cluster add to LL
						Lef = getIndexLeftPeriodic(snA,StartID);
						if (matchExist(AvailableNodes, Lef)==1){
							addLLNode(LL,Lef);
						}
					}
					if(getFlagRig(tempNode)==1){
						//Node is within cluster add to LL
						Rig = getIndexRightPeriodic(snA,StartID);
						if (matchExist(AvailableNodes, Rig)==1){
							addLLNode(LL,Rig);
						}
					}
					if( getFlagAbo(tempNode)==1){
						//Node is within cluster add to LL
						Abo = getIndexAbove(snA, StartID);
						if (Abo!=-1){
							if (matchExist(AvailableNodes, Abo)==1){
								addLLNode(LL,Abo);
							}
						}
					}
					if(getFlagBel(tempNode)==1){
						//Node is within cluster add to LL
						Bel = getIndexBelow(snA, StartID);
						if (Bel!=-1){
							if (matchExist(AvailableNodes, Bel)==1){
								addLLNode(LL,Bel);
							}
						}
					}

					if(getFlagFro(tempNode)==1){
						//Node is within cluster
						Fro = getIndexFront(snA, StartID);
						if (Fro!=-1){
							if (matchExist(AvailableNodes, Fro)==1){
								addLLNode(LL,Fro);
							}
						}
					}
					if(getFlagBeh(tempNode)==1){
						//Node is within cluster
						Beh = getIndexBehind(snA, StartID);
						if (Beh!=-1){
							if (matchExist(AvailableNodes, Beh)==1){
								addLLNode(LL,Beh);
							}
						}
					}

					//Check to see if the added nodes are attached to the right
					//Electrode if any nodes are the cluster is either
					//to big to use the approximation or it is a percolating
					//pathway
					for(jj=1;jj<=getRows(NearRightElec);jj++){

						RElec = (int)getE(NearRightElec,jj,1);

						if(getFlagLef(tempNode)==1){
							if (RElec==Lef){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagRig(tempNode)==1){
							if (RElec==Rig){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if( getFlagAbo(tempNode)==1){
							if (RElec==Abo){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBel(tempNode)==1){
							if (RElec==Bel){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}

						if(getFlagFro(tempNode)==1){
							if (RElec==Fro){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBeh(tempNode)==1){
							if( RElec==Beh){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}


					}

					//Cycle through the link list

					while (getLLlength(LL)!=0 ){

						AvailableID = getLLstartID(LL); 
						tempNode = getCluster_Node(ClLL,AvailableID);
						//Cross the first ID out of the available ids
						//And also cut it out of the linklist
						matchReplace(AvailableNodes,AvailableID,-1);
						removeLLNode(LL, AvailableID);

						//Check the neighbors of this node
						//Grab surrounding Neighbors that are within the cluster
						tempNode = getCluster_Node(ClLL,AvailableID);
						if(getFlagLef(tempNode)==1){
							//Node is within cluster add to LL
							Lef = getIndexLeftPeriodic(snA,AvailableID);
							if (matchExist(AvailableNodes, Lef)==1){
								addLLNode(LL,Lef);
							}
						}
						if(getFlagRig(tempNode)==1){
							//Node is within cluster add to LL

							Rig = getIndexRightPeriodic(snA,AvailableID);
							if (matchExist(AvailableNodes, Rig)==1){
								addLLNode(LL,Rig);
							}
						}
						if( getFlagAbo(tempNode)==1){
							//Node is within cluster add to LL
							Abo = getIndexAbove(snA, AvailableID);
							if (Abo!=-1){
								if (matchExist(AvailableNodes, Abo)==1){
									addLLNode(LL,Abo);
								}
							}
						}
						if(getFlagBel(tempNode)==1){
							//Node is within cluster add to LL
							Bel = getIndexBelow(snA, AvailableID);
							if (Bel!=-1){
								if (matchExist(AvailableNodes, Bel)==1){
									addLLNode(LL,Bel);
								}
							}
						}

						if(getFlagFro(tempNode)==1){
							//Node is within cluster
							Fro = getIndexFront(snA, AvailableID);
							if (Fro!=-1){
								if (matchExist(AvailableNodes, Fro)==1){
									addLLNode(LL,Fro);
								}
							}
						}
						if(getFlagBeh(tempNode)==1){
							//Node is within cluster
							Beh = getIndexBehind(snA, AvailableID);
							if (Beh!=-1){
								if (matchExist(AvailableNodes, Beh)==1){
									addLLNode(LL,Beh);
								}
							}
						}

						//Check to see if the added nodes are attached to the right
						//Electrode if any nodes are the cluster is either
						//to big to use the approximation or it is a percolating
						//pathway
						for(jj=1;jj<=getRows(NearRightElec);jj++){

							RElec = (int)getE(NearRightElec,jj,1);

							if(getFlagLef(tempNode)==1){
								if (RElec==Lef){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagRig(tempNode)==1){
								if (RElec==Rig){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if( getFlagAbo(tempNode)==1){
								if (RElec==Abo){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBel(tempNode)==1){
								if (RElec==Bel){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

							if(getFlagFro(tempNode)==1){
								if (RElec==Fro){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBeh(tempNode)==1){
								if( RElec==Beh){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

						}

					}
				}
			}

		}else if(PeriodicY==0 && PeriodicZ==1){

			for(StartRow=1;StartRow<=getRows(NearLeftElec);StartRow++){
				//Start at the left Electrode
				StartID = (int)getE(NearLeftElec,StartRow,1);

				if (matchExist(AvailableNodes,StartID)==1){
					//Cross the first ID out of the available ids
					matchReplace(AvailableNodes,StartID,-1);

					//Grab surrounding Neighbors that are within the cluster
					tempNode = getCluster_Node(ClLL,StartID);
					if(getFlagLef(tempNode)==1){
						//Node is within cluster add to LL
						Lef = getIndexLeft(snA,StartID);
						if (Lef!=-1){
							if (matchExist(AvailableNodes, Lef)==1){
								addLLNode(LL,Lef);
							}
						}
					}
					if(getFlagRig(tempNode)==1){
						//Node is within cluster add to LL
						Rig = getIndexRight(snA,StartID);
						if (Rig!=-1){
							if (matchExist(AvailableNodes, Rig)==1){
								addLLNode(LL,Rig);
							}
						}
					}
					if( getFlagAbo(tempNode)==1){
						//Node is within cluster add to LL
						Abo = getIndexAbovePeriodic(snA, StartID);
						if (matchExist(AvailableNodes, Abo)==1){
							addLLNode(LL,Abo);
						}
					}
					if(getFlagBel(tempNode)==1){
						//Node is within cluster add to LL
						Bel = getIndexBelowPeriodic(snA, StartID);
						if (matchExist(AvailableNodes, Bel)==1){
							addLLNode(LL,Bel);
						}
					}

					if(getFlagFro(tempNode)==1){
						//Node is within cluster
						Fro = getIndexFront(snA, StartID);
						if (Fro!=-1){
							if (matchExist(AvailableNodes, Fro)==1){
								addLLNode(LL,Fro);
							}
						}
					}
					if(getFlagBeh(tempNode)==1){
						//Node is within cluster
						Beh = getIndexBehind(snA, StartID);
						if (Beh!=-1){
							if (matchExist(AvailableNodes, Beh)==1){
								addLLNode(LL,Beh);
							}
						}
					}

					//Check to see if the added nodes are attached to the right
					//Electrode if any nodes are the cluster is either
					//to big to use the approximation or it is a percolating
					//pathway
					for(jj=1;jj<=getRows(NearRightElec);jj++){

						RElec = (int)getE(NearRightElec,jj,1);

						if(getFlagLef(tempNode)==1){
							if (RElec==Lef){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagRig(tempNode)==1){
							if (RElec==Rig){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if( getFlagAbo(tempNode)==1){
							if (RElec==Abo){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBel(tempNode)==1){
							if (RElec==Bel){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}

						if(getFlagFro(tempNode)==1){
							if (RElec==Fro){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBeh(tempNode)==1){
							if( RElec==Beh){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}


					}

					//Cycle through the link list

					while (getLLlength(LL)!=0 ){

						AvailableID = getLLstartID(LL); 
						tempNode = getCluster_Node(ClLL,AvailableID);
						//Cross the first ID out of the available ids
						//And also cut it out of the linklist
						matchReplace(AvailableNodes,AvailableID,-1);
						removeLLNode(LL, AvailableID);

						//Check the neighbors of this node
						//Grab surrounding Neighbors that are within the cluster
						tempNode = getCluster_Node(ClLL,AvailableID);
						if(getFlagLef(tempNode)==1){
							//Node is within cluster add to LL
							Lef = getIndexLeft(snA,AvailableID);
							if (Lef!=-1){
								if (matchExist(AvailableNodes, Lef)==1){
									addLLNode(LL,Lef);
								}
							}
						}
						if(getFlagRig(tempNode)==1){
							//Node is within cluster add to LL
							Rig = getIndexRight(snA,AvailableID);
							if(Rig!=-1){
								if (matchExist(AvailableNodes, Rig)==1){
									addLLNode(LL,Rig);
								}
							}
						}
						if( getFlagAbo(tempNode)==1){
							//Node is within cluster add to LL
							Abo = getIndexAbovePeriodic(snA, AvailableID);
							if (matchExist(AvailableNodes, Abo)==1){
								addLLNode(LL,Abo);
							}
						}
						if(getFlagBel(tempNode)==1){
							//Node is within cluster add to LL
							Bel = getIndexBelowPeriodic(snA, AvailableID);
							if (matchExist(AvailableNodes, Bel)==1){
								addLLNode(LL,Bel);
							}
						}

						if(getFlagFro(tempNode)==1){
							//Node is within cluster
							Fro = getIndexFront(snA, AvailableID);
							if (Fro!=-1){
								if (matchExist(AvailableNodes, Fro)==1){
									addLLNode(LL,Fro);
								}
							}
						}
						if(getFlagBeh(tempNode)==1){
							//Node is within cluster
							Beh = getIndexBehind(snA, AvailableID);
							if (Beh!=-1){
								if (matchExist(AvailableNodes, Beh)==1){
									addLLNode(LL,Beh);
								}
							}
						}

						//Check to see if the added nodes are attached to the right
						//Electrode if any nodes are the cluster is either
						//to big to use the approximation or it is a percolating
						//pathway
						for(jj=1;jj<=getRows(NearRightElec);jj++){

							RElec = (int)getE(NearRightElec,jj,1);

							if(getFlagLef(tempNode)==1){
								if (RElec==Lef){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagRig(tempNode)==1){
								if (RElec==Rig){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if( getFlagAbo(tempNode)==1){
								if (RElec==Abo){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBel(tempNode)==1){
								if (RElec==Bel){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagFro(tempNode)==1){
								if (RElec==Fro){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBeh(tempNode)==1){
								if( RElec==Beh){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

						}

					}
				}
			}

		}else{
			//Non periodic

			for(StartRow=1;StartRow<=getRows(NearLeftElec);StartRow++){
				//Start at the left Electrode
				StartID = (int)getE(NearLeftElec,StartRow,1);

				if (matchExist(AvailableNodes,StartID)==1){
					//Cross the first ID out of the available ids
					matchReplace(AvailableNodes,StartID,-1);

					//Grab surrounding Neighbors that are within the cluster
					tempNode = getCluster_Node(ClLL,StartID);
					if(getFlagLef(tempNode)==1){
						//Node is within cluster add to LL
						Lef = getIndexLeft(snA,StartID);
						if (matchExist(AvailableNodes, Lef)==1){
							addLLNode(LL,Lef);
						}
					}
					if(getFlagRig(tempNode)==1){
						//Node is within cluster add to LL
						Rig = getIndexRight(snA,StartID);
						if (matchExist(AvailableNodes, Rig)==1){
							addLLNode(LL,Rig);
						}
					}
					if( getFlagAbo(tempNode)==1){
						//Node is within cluster add to LL
						Abo = getIndexAbove(snA, StartID);
						if (matchExist(AvailableNodes, Abo)==1){
							addLLNode(LL,Abo);
						}
					}
					if(getFlagBel(tempNode)==1){
						//Node is within cluster add to LL
						Bel = getIndexBelow(snA, StartID);
						if (matchExist(AvailableNodes, Bel)==1){
							addLLNode(LL,Bel);
						}
					}

					if(getFlagFro(tempNode)==1){
						//Node is within cluster
						Fro = getIndexFront(snA, StartID);
						if (Fro!=-1){
							if (matchExist(AvailableNodes, Fro)==1){
								addLLNode(LL,Fro);
							}
						}
					}
					if(getFlagBeh(tempNode)==1){
						//Node is within cluster
						Beh = getIndexBehind(snA, StartID);
						if (Beh!=-1){
							if (matchExist(AvailableNodes, Beh)==1){
								addLLNode(LL,Beh);
							}
						}
					}

					//Check to see if the added nodes are attached to the right
					//Electrode if any nodes are the cluster is either
					//to big to use the approximation or it is a percolating
					//pathway
					for(jj=1;jj<=getRows(NearRightElec);jj++){

						RElec = (int)getE(NearRightElec,jj,1);

						if(getFlagLef(tempNode)==1){
							if (RElec==Lef){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagRig(tempNode)==1){
							if (RElec==Rig){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if( getFlagAbo(tempNode)==1){
							if (RElec==Abo){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBel(tempNode)==1){
							if (RElec==Bel){
								//This means delete the Cluster
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}

						if(getFlagFro(tempNode)==1){
							if (RElec==Fro){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}
						if(getFlagBeh(tempNode)==1){
							if( RElec==Beh){
								deleteLL(LL);
								deleteMatrix(AvailableNodes);
								deleteMatrix(NearLeftElec);
								deleteMatrix(NearRightElec);
								return 1;
							}
						}


					}

					//Cycle through the link list

					while (getLLlength(LL)!=0 ){

						AvailableID = getLLstartID(LL); 
						tempNode = getCluster_Node(ClLL,AvailableID);
						//Cross the first ID out of the available ids
						//And also cut it out of the linklist
						matchReplace(AvailableNodes,AvailableID,-1);
						removeLLNode(LL, AvailableID);

						//Check the neighbors of this node
						//Grab surrounding Neighbors that are within the cluster
						tempNode = getCluster_Node(ClLL,AvailableID);
						if(getFlagLef(tempNode)==1){
							//Node is within cluster add to LL
							Lef = getIndexLeft(snA,AvailableID);
							if (matchExist(AvailableNodes, Lef)==1){
								addLLNode(LL,Lef);
							}
						}
						if(getFlagRig(tempNode)==1){
							//Node is within cluster add to LL
							Rig = getIndexRight(snA,AvailableID);
							if (matchExist(AvailableNodes, Rig)==1){
								addLLNode(LL,Rig);
							}
						}
						if( getFlagAbo(tempNode)==1){
							//Node is within cluster add to LL
							Abo = getIndexAbove(snA, AvailableID);
							if (matchExist(AvailableNodes, Abo)==1){
								addLLNode(LL,Abo);
							}
						}
						if(getFlagBel(tempNode)==1){
							//Node is within cluster add to LL
							Bel = getIndexBelow(snA, AvailableID);
							if (matchExist(AvailableNodes, Bel)==1){
								addLLNode(LL,Bel);
							}
						}

						if(getFlagFro(tempNode)==1){
							//Node is within cluster
							Fro = getIndexFront(snA, AvailableID);
							if (Fro!=-1){
								if (matchExist(AvailableNodes, Fro)==1){
									addLLNode(LL,Fro);
								}
							}
						}
						if(getFlagBeh(tempNode)==1){
							//Node is within cluster
							Beh = getIndexBehind(snA, AvailableID);
							if (Beh!=-1){
								if (matchExist(AvailableNodes, Beh)==1){
									addLLNode(LL,Beh);
								}
							}
						}

						//Check to see if the added nodes are attached to the right
						//Electrode if any nodes are the cluster is either
						//to big to use the approximation or it is a percolating
						//pathway
						for(jj=1;jj<=getRows(NearRightElec);jj++){

							RElec = (int)getE(NearRightElec,jj,1);

							if(getFlagLef(tempNode)==1){
								if (RElec==Lef){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagRig(tempNode)==1){
								if (RElec==Rig){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if( getFlagAbo(tempNode)==1){
								if (RElec==Abo){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBel(tempNode)==1){
								if (RElec==Bel){
									//This means delete the Cluster
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagFro(tempNode)==1){
								if (RElec==Fro){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}
							if(getFlagBeh(tempNode)==1){
								if( RElec==Beh){
									deleteLL(LL);
									deleteMatrix(AvailableNodes);
									deleteMatrix(NearLeftElec);
									deleteMatrix(NearRightElec);
									return 1;
								}
							}

						}

					}
				}
			}

		}
	}

	deleteLL(LL);
	deleteMatrix(AvailableNodes);
	deleteMatrix(NearLeftElec);
	deleteMatrix(NearRightElec);
	return 0;
}

int PrintCheck(int TotalOrders, int orderLow, const_ArbArray ClArLL, const_SNarray snA, const_matrix MasterM){

	if ( MasterM==NULL || ClArLL==NULL || snA==NULL || TotalOrders<=0){
		return -1;
	}

	int element;
	int order;
	int Node_ID;
	int i, j, k;
	ClusterLL TempClLL;
	Node tempNode;
	for(element=0;element<TotalOrders;element++){

		TempClLL=(ClusterLL) getArbElement(ClArLL,element);
		order=element+orderLow;
		printf("Order %d\n",order);
		if(TempClLL==NULL){
			printf("Cluster Empty for current order\n");
		} else { 

			printClusterLL(TempClLL);

			while (TempClLL!=NULL){
				//printNodesClusterLL(TempClLL);		
				tempNode= getStartNode(TempClLL);

				printf("\nCluster Number %d\n",getCluster_id(TempClLL));

				while (tempNode!=NULL){

					Node_ID=getNode_id(tempNode);
					getLoc( &i, &j, &k, Node_ID, snA);
					//l-1 hop behind
					//l-2 hop infront
					//l-3 hop left
					//l-4 hop right
					//l-5 hop below
					//l-6 hop above
					if (i==0){
						printf("x = 0\n");
					}else if(i==getAlen(snA)-1){
						printf("x = %d\n",getAlen(snA));
					} 


					if (j==0){
						printf("y = 0\n");
					} else if(j==getAwid(snA)-1) {
						printf("y = %d\n",getAwid(snA));
					}

					if (k==0){
						printf("z = 0\n");
					} else if(k==getAhei(snA)-1) {
						printf("z = %d\n",getAhei(snA));
					} 

					printf("Hops to Neighbors Node: %d i %d j %d k %d\n",Node_ID,i,j,k);
					printf("Behind \t Front \t Left \t Right \t Below \t Above\n");
					if (i==0){
						printf("\t %d \t ", getIndexFront(snA, Node_ID));
					}else if(i==getAlen(snA)-1){
						printf("%d \t \t ",getIndexBehind(snA, Node_ID));
					} else {
						printf("%d \t %d \t ",getIndexBehind(snA, Node_ID), getIndexFront(snA, Node_ID));
					}


					if (j==0){
						printf("\t %d \t ", getIndexRightPeriodic(snA, Node_ID));
					} else if(j==getAwid(snA)-1) {
						printf("%d \t \t ",getIndexLeftPeriodic(snA, Node_ID));
					}else{

						//Because Periodic in y
						printf("%d \t %d \t ",getIndexLeftPeriodic(snA, Node_ID), getIndexRightPeriodic(snA, Node_ID));
					}

					if (k==0){
						printf("\t %d\n",getIndexAbove(snA, Node_ID));
					} else if(k==getAhei(snA)-1) {
						printf("%d \t \t\n",getIndexBelow(snA, Node_ID));
					} else {
						printf("%d \t %d\n",getIndexBelow(snA, Node_ID), getIndexAbove(snA, Node_ID));
					}
					//Printing the order of the hop in the surrounding directions
					printf("%g \t %g \t %g \t %g \t %g \t %g\n",\
							getE(MasterM,Node_ID+1,1),getE(MasterM, Node_ID+1, 2), getE(MasterM,Node_ID+1,3),\
							getE(MasterM,Node_ID+1,4),getE(MasterM, Node_ID+1, 5), getE(MasterM,Node_ID+1,6));

					tempNode=getNextNode(tempNode);
				}
				TempClLL=getNextClusterLL(TempClLL);
			}

		}

	}

	return 0;

}

int CalculateNeighNodes(int TotalOrders, int orderLow, ArbArray * ClArLL, SNarray snA, int PeriodicX, int PeriodicY, int PeriodicZ ){

	if(TotalOrders<0 || ClArLL==NULL || snA==NULL){
		return -1;
	}

	if(*ClArLL==NULL){
		return -1;
	}

	int element;
	int Node_ID;
	int Node_IDFro, Node_IDBeh, Node_IDLef, Node_IDRig, Node_IDAbo, Node_IDBel;
	int Elec_Fro, Elec_Beh, Elec_Lef, Elec_Rig, Elec_Abo, Elec_Bel;
	int i, j, k;
	int rv;
	int order;
	Node tempNode;
	ClusterLL TempClLL;

	for(element=0;element<TotalOrders;element++){
		TempClLL=(ClusterLL) getArbElement((*ClArLL),element);
		order=element+orderLow;
		//printf("Order %d\n",order);

		if(TempClLL==NULL) {
		//	printf("Cluster Empty CalculateNeighNodes function\n");
		}else{
			if(PeriodicX==1 && PeriodicY==1 && PeriodicZ==1){	
				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFroP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist
							rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBehP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLefP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRigP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA, i, j, k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBelP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA, i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAboP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else if(PeriodicX == 0 && PeriodicY==1 && PeriodicZ==1){

				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFro(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBeh(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLefP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRigP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA, i, j, k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBelP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA, i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAboP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}


			}else if(PeriodicX == 1 && PeriodicY==0 && PeriodicZ==1){
				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFroP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist
							rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBehP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLef(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRig(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA, i, j, k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBelP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA, i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAboP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else if(PeriodicX == 1 && PeriodicY==1 && PeriodicZ==0){

				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFroP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist
							rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBehP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLefP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRigP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA,i,j,k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBel(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA,i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAbo(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}	
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else if(PeriodicX == 0 && PeriodicY==0 && PeriodicZ==1){

				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFro(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBeh(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLef(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRig(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA, i, j, k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBelP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA, i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAboP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else if(PeriodicX == 1 && PeriodicY==0 && PeriodicZ==0){

				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFroP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist
							rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBehP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLef(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRig(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA,i,j,k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBel(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA,i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAbo(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else if(PeriodicX == 0 && PeriodicY==1 && PeriodicZ==0){

				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro=getIndFro(snA, i, j, k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFro(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh=getIndBeh(snA, i, j, k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBeh(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLefP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						} else {
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig=getIndRig(snA, i, j, k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRigP(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel=getIndBel(snA, i, j, k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBel(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo=getIndAbo(snA, i, j, k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAbo(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}else{
				//Non-periodic
				while (TempClLL!=NULL){

					tempNode=(Node) getStartNode(TempClLL);

					while (tempNode!=NULL){

						//0 - means that the Node should be considered a 
						//neighbor
						//1 - means the node is part of the cluster and 
						//is therefore not a neighbor
						Node_ID = getNode_id(tempNode);	
						getLoc( &i, &j, &k, Node_ID, snA);

						Elec_Fro= getIndFro(snA,i,j,k);
						if(getFlagFro(tempNode)==0){
							//This is a neighbor that is infront
							Node_IDFro=getIndFro(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							//Add NeighNode to the linklist

							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDFro);
							}
						}else{
							if(Elec_Fro==-1){
								//Add the Front Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,1);
							}
						}
						Elec_Beh = getIndBeh(snA,i,j,k);
						if(getFlagBeh(tempNode)==0){
							//Neighbor that is behind
							Node_IDBeh=getIndBeh(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries

							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBeh);
							}
						}else{
							if(Elec_Beh==-1){ 
								//Add the left Electrode as a neighbor to the cluster
								setCluster_elecXid(&TempClLL,0);
							}
						}
						Elec_Lef = getIndLef(snA,i,j,k);
						if(getFlagLef(tempNode)==0){
							//Neighbor that is to the left
							Node_IDLef=getIndLef(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDLef);
							}
						}else{
							if(Elec_Lef==-1) {
								setCluster_elecYid(&TempClLL,0);
							}
						}
						Elec_Rig = getIndRig(snA,i,j,k);
						if(getFlagRig(tempNode)==0){
							//Neighbor that is to the right
							Node_IDRig=getIndRig(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDRig);
							}
						}else{
							if(Elec_Rig==-1) {
								setCluster_elecYid(&TempClLL,1);
							}
						}
						Elec_Bel = getIndBel(snA,i,j,k);
						if(getFlagBel(tempNode)==0 && k!=0){
							//Neighbor that is below
							Node_IDBel=getIndBel(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDBel);
							}
						}else{
							if(Elec_Bel==-1){
								setCluster_elecZid(&TempClLL,0);
							}
						}
						Elec_Abo = getIndAbo(snA,i,j,k);
						if(getFlagAbo(tempNode)==0 && k!=getAhei(snA)){
							//Neighbor that is above
							Node_IDAbo=getIndAbo(snA, i, j, k);
							//If ID is -1 means it's outside the boundaries
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}else{
								rv = addNeighNodeToCluster(&TempClLL, Node_IDAbo);
							}
						}else{
							if(Elec_Abo==-1) {
								setCluster_elecZid(&TempClLL,1);
							}
						}
						tempNode=getNextNode(tempNode);
					}

					TempClLL=getNextClusterLL(TempClLL);

				}

			}

		}//end of if and else statement
	}//end of for loop with elements

	return 0; 
}

int CalculateSumAndP(const int TotalOrders, const_SNarray snA, ArbArray * ClArLL,\
		const_matrix MasterM, const int attempts,\
		const int PeriodicX,const int PeriodicY,const int PeriodicZ){

	if(PeriodicX<0 || PeriodicY<0 || PeriodicZ<0 ||\
			PeriodicX>1 || PeriodicY>1 || PeriodicZ>1){
		return -1;
	}

	if(TotalOrders<0 || snA==NULL || ClArLL ==NULL ||\
			MasterM==NULL || attempts<1){
		return -1;
	}
	//Cylce through the nodes in ClArLL and for every node
	//determine the hop rate off of the cluster thet total
	//average of all the hop rates off the cluster will be
	//used to calculate the dwell time for a charge stuck 
	//in the cluster.
	//Individual hoprates off the cluster will determine 
	//the probability a charge will hop to a given site
	int element;
	double count, count2;
	double sum, sum2;
	ClusterLL TempClLL;
	Node tempNode;
	matrix mtxHopOpt;
	matrix mtxProb;
	matrix mtxProbNeigh;
	matrix mtxProbNeighDwell;
	matrix mtxDwellTime;
	int inc;
	double PneighTotal;
	//Counts the number of options to hop off the cluster
	double countNeighOpts;

	for(element=0;element<TotalOrders;element++){
		TempClLL=(ClusterLL) getArbElement(*ClArLL,element);

		//printClusterLL(TempClLL);

		if(TempClLL==NULL) {
			//printf("Cluster Empty\n");
		} else {

			while (TempClLL!=NULL){

				count=0.0;
				count2=0.0;
				sum=0.0;
				sum2=0.0;
				tempNode= getStartNode(TempClLL);

				mtxHopOpt = newMatrix( getCluster_numNodes(TempClLL),3);

				/////////////////////////////////////////////////////////
				//1 - Counting how many options each site has to hop off
				//and within the cluster
				//col 1 - hop within cluster
				//col 2 - hop off the cluster
				//col 3 - ID of the node

				CountOptions( tempNode, &mtxHopOpt, snA);

				//This function calculates the probability a charge
				//will hop to a site within the cluster irrespective of
				//time but based on spacial orientation
				mtxProb = CalculateProb( TempClLL, mtxHopOpt, snA, attempts );

				PneighTotal = 0;
				tempNode= getStartNode(TempClLL);
				inc = 1;
				while (tempNode!=NULL){
					//Multiply the probability on a site by the number of 
					//hops off the site
					PneighTotal+=getE(mtxProb,inc,1)*getE(mtxHopOpt,inc,2);
					inc++;
					tempNode = getNextNode(tempNode);
				}
				mtxProbNeigh = duplicateMatrix(mtxProb);
				DivideEachElementCol( &mtxProbNeigh,1, PneighTotal);

				double rateN = 0.0;
				//This function calculates the amount of time on average
				//a charge will spend on a given site within the cluster
				//mtxDwellTime is normalized so it is actually a ratio and not a time
				countNeighOpts = SumOfCol(mtxHopOpt,2);
				matrix mtxTimes = newMatrix(countNeighOpts,1);
				mtxDwellTime = CalculateDwellTimeAndRateN(&TempClLL,&mtxTimes,mtxProbNeigh, MasterM,\
						snA, &rateN, PeriodicX, PeriodicY, PeriodicZ);

				//The matrix produced from this function calculates the
				//amount of time a charge will take to hop to a given 
				//neighbor site. 
				tempNode = getStartNode(TempClLL);
				mtxProbNeighDwell = CalculateProbNeighDwell(countNeighOpts,mtxDwellTime,mtxProb,\
						tempNode,MasterM,snA, rateN,\
						PeriodicX, PeriodicY, PeriodicZ);

				//Now we calculate the pval for the cluster
				//pval - Stay in Cluster
				//pval - Hop off cluster
				//pval - sites off cluster
				//pval - sites within cluster
				CalculatePvalNeigh(&TempClLL, mtxTimes, mtxProbNeighDwell);

				CalculatePvalNodes(&TempClLL, mtxProb, mtxDwellTime);



				//Delete Matrices
				deleteMatrix(mtxHopOpt);
				deleteMatrix(mtxProb);
				deleteMatrix(mtxProbNeigh);
				deleteMatrix(mtxTimes);
				deleteMatrix(mtxDwellTime);
				deleteMatrix(mtxProbNeighDwell);

				TempClLL = getNextClusterLL(TempClLL);

			}// End of while loop for ClLL


		}//End of if statement


	}//End of For loop
	return 0;
}

int CountOptions( Node tempNode, matrix * mtxHopOpt, const_SNarray snA){
	//Function counts the number of hopping options 
	//a given site has within the cluster if it
	//wants to keep hopping within the cluster and
	//off the cluster
	//matrix mtxHopOpt stores the Id of the site
	//and the corresponding number of options
	//col 1 - hops within cluster
	//col 2 - hops off the cluster
	//col 3 - ID of site
	if ((*mtxHopOpt)==NULL || snA==NULL){
		return -1;
	}

	int countOpt;
	int countOpt2;
	int inc = 1;
	int Node_ID;
	int i, j, k;
	while (tempNode!=NULL){

		countOpt=0;
		countOpt2=0;

		Node_ID = getNode_id(tempNode);	
		getLoc( &i, &j, &k, Node_ID, snA);
		if(getFlagFro(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}
		if(getFlagBeh(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}
		if(getFlagLef(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}
		if(getFlagRig(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}
		if(getFlagBel(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}
		if(getFlagAbo(tempNode)==1){
			countOpt++;
		}else{
			countOpt2++;
		}

		setE((*mtxHopOpt),inc,1,countOpt);
		setE((*mtxHopOpt),inc,2,countOpt2);
		setE((*mtxHopOpt),inc,3,Node_ID);					
		inc++;
		tempNode = getNextNode(tempNode);
	}

	return 0;
}

matrix CalculateProb(const_ClusterLL TempClLL, matrix mtxHopOpt, const_SNarray snA, const int attempts ){


	//This function can be used for periodic or non periodic conditions it does not
	//do anything with neighbors outside of the clusters. 
	int inc;
	int Node_ID;
	int Node_IDFro, Node_IDBeh;
	int Node_IDLef, Node_IDRig;
	int Node_IDBel, Node_IDAbo;
	int i, j, k;
	int Row;
	int Row2;
	double val;

	//The second column of mtxProb will contain the IDs of the respective sites
	//We need to reinitialize the 2nd column to the ID's of the CNTs
	matrix mtxProb = newMatrixSet( getCluster_numNodes(TempClLL),2, (1/((double)getCluster_numNodes(TempClLL))));
	matrix mtxProbNew = newMatrix(6,1);
	Node tempNode;

	//Initilize Node Ids in the mtxProb
	tempNode = getStartNode(TempClLL);

	printf("\n****************Calculating Prob Matrix********************\n");
	inc = 1;
	while(tempNode!=NULL){
		Node_ID = getNode_id(tempNode);
		setE(mtxProb,inc,2,(double) Node_ID);
		tempNode = getNextNode(tempNode);
		inc++;
	}

	for(int attempt=1;attempt<attempts;attempt++){
		tempNode= getStartNode(TempClLL);
		inc = 1;

		while (tempNode!=NULL){

			//1 hop behind     	index - 0 
			//2 hop infront			index - 1
			//3 hop left				index - 2
			//4 hop right				index - 3
			//5 hop below				index - 4 
			//6 hop above				index - 5
			Node_ID = getNode_id(tempNode);

			//printf("Node_id %d inc %d\n",Node_ID,inc);
			getLoc( &i, &j, &k, Node_ID, snA);

			if(getFlagFro(tempNode)==1){
				//Node in front is within the cluster
				Node_IDFro = getIndFroP(snA,i,j,k);
				//printf("Node in front of %d is %d\n",Node_ID, Node_IDFro);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDFro, 3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDFro,2);
				//printf("Row %d Row2 %d\n",Row,Row2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,2,1,val);
			}
			if(getFlagBeh(tempNode)==1){
				//Node behind is within the cluster
				Node_IDBeh = getIndBehP(snA,i,j,k);
				//printf("Node behind %d is %d\n",Node_ID, Node_IDBeh);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDBeh,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDBeh,2);
				//printf("Row %d Row2 %d\n",Row,Row2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,1,1,val);
			}

			if(getFlagLef(tempNode)==1){
				//Node behind is within the cluster
				Node_IDLef = getIndLefP(snA,i,j,k);
				//printf("Node to the left of %d is %d\n",Node_ID, Node_IDLef);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDLef,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDLef,2);
				//printf("Row %d Row2 %d\n",Row,Row2);

				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,3,1,val);
			}
			if(getFlagRig(tempNode)==1){
				//Node behind is within the cluster
				Node_IDRig = getIndRigP(snA,i,j,k);
				//printf("Node to the Right of %d is %d\n",Node_ID, Node_IDRig);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDRig,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDRig,2);
				//printf("Row %d Row2 %d\n",Row,Row2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,4,1,val);
			}
			if(getFlagBel(tempNode)==1){
				//Node behind is within the cluster
				Node_IDBel = getIndBelP(snA,i,j,k);
				//printf("Node below %d is %d\n",Node_ID, Node_IDBel);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDBel,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDBel,2);
				//printf("Row %d Row2 %d\n",Row,Row2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,5,1,val);
			}

			if(getFlagAbo(tempNode)==1){
				//Node behind is within the cluster
				Node_IDAbo = getIndAboP(snA,i,j,k);
				//printf("Node above %d is %d\n",Node_ID, Node_IDAbo);
				Row = FindRowOfMatchInCol(mtxHopOpt, Node_IDAbo,3);
				Row2 = FindRowOfMatchInCol(mtxProb, Node_IDAbo,2);
				//printf("Row %d Row2 %d\n",Row,Row2);
				val = (1/getE(mtxHopOpt,Row,1))*getE(mtxProb,Row2,1);
				setE(mtxProbNew,6,1,val);
			}

			val = (SumOfCol(mtxProbNew,1)+getE(mtxProb,inc,1))/2;

			//printf("Prob of site %d val: %g\n",Node_ID,val);
			setE(mtxProb,inc,1,val);
			inc++;
			tempNode = getNextNode(tempNode);
			setAll(mtxProbNew, 0.0);
		}

		val = SumOfCol(mtxProb,1);
		//Normalize the matrix
		DivideEachElementCol(&mtxProb,1, val);

	}

	deleteMatrix(mtxProbNew);
	return mtxProb;
}

matrix CalculateProbNeighDwell(const int countNeighOpts,\
		matrix mtxDwellTime,const_matrix mtxProb,  Node tempNode,\
		const_matrix MasterM, const_SNarray snA, const double rateN,\
		int PeriodicX, int PeriodicY, int PeriodicZ){

	//The second column stores the id of the NeighSite charge
	//would be hopping too
	matrix mtxProbNeighDwell = newMatrix(countNeighOpts,2);
	//Initially set the values to -1 for the IDs so we know where the periodic
	//and non-periodic conditions kick in.
	setAllRowsInCol(mtxProbNeighDwell,2,-1);
	int inc = 1;
	int inc2 = 1;
	double Val = 0;
	double Val2; 
	int i,j,k;
	double dwelltemp;
	int Node_ID;
	int Check;
	int Row;

	printf("\n**********************Creating mtxProbNeighDwell Matrix **************************\n\n");

	while(tempNode!=NULL){

		Node_ID = getNode_id(tempNode);
		Check = getIndexFront(snA, Node_ID);
		dwelltemp = getE(mtxDwellTime,inc,1);
		getLoc(&i,&j,&k,Node_ID,snA);
		Row = (int) FindRowOfMatchInCol( mtxProb, (double) Node_ID,2);
		//Do if node is inside the boundaries or if Periodic in X
		if(Check!=-1 || PeriodicX!=0){
			if(getFlagFro(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,8)*\
							 (1/rateN)*dwelltemp;

				//printMatrix(mtxProb);
				//printMatrix(MasterM);	
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,8));
				//printf("rateN %g\n",rateN);
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexFrontP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		Check = getIndexBehind(snA, Node_ID);
		if(Check!=-1 || PeriodicX!=0){
			if(getFlagBeh(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,7)*\
							 (1/rateN)*dwelltemp;
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,7));
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexBehindP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		Check = getIndexLeft(snA, Node_ID);
		if(Check!=-1 || PeriodicY!=0){
			if(getFlagLef(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,9)*\
							 (1/rateN)*dwelltemp;
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,9));
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexLeftP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		Check = getIndexRight(snA, Node_ID);
		if(Check!=-1 || PeriodicY!=0){
			if(getFlagRig(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,10)*\
							 (1/rateN)*dwelltemp;
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,10));
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexRightP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		Check = getIndexBelow(snA, Node_ID);
		if(Check!=-1 || PeriodicZ!=0){
			if(getFlagBel(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,11)*\
							 (1/rateN)*dwelltemp;
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,11));
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexBelowP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		Check = getIndexAbove(snA, Node_ID);
		if(Check!=-1 || PeriodicZ!=0){
			if(getFlagAbo(tempNode)!=1){
				Val2 = getE(mtxProb,Row,1)*getE(MasterM, getIndex(snA,i,j,k)+1,12)*\
							 (1/rateN)*dwelltemp;
				//printf("Node ID %d Inc %d\t",Node_ID,inc);
				//printf("Increment %d Val %g dwelltemp %g MasterM %g\t",inc2,Val2,dwelltemp,getE(MasterM, getIndex(snA,i,j,k)+1,12));
				//printf("Prob %g\n",getE(mtxProb,Row,1));
				setE(mtxProbNeighDwell,inc2,1,Val2);
				setE(mtxProbNeighDwell,inc2,2,getIndexAboveP(snA,Node_ID));
				inc2++;
				Val+=Val2;
			}
		}
		tempNode = getNextNode(tempNode);
		inc++;
	}

	//Normalize Matrix
	DivideEachElementCol( &mtxProbNeighDwell,1,Val);
	return mtxProbNeighDwell;
}

matrix CalculateDwellTimeAndRateN(ClusterLL *  TempClLL,matrix * mtxTimes, const_matrix mtxProbNeigh,\
		const_matrix MasterM, const_SNarray snA, double * rateN,\
		const int PeriodicX,\
		const int PeriodicY,\
		const int PeriodicZ){

	//printSNarray(snA);
	//printMatrix(MasterM);
	//printf("printing mtxProbNeigh third time\n");
	//printMatrix(mtxProbNeigh);
	//printClusterLL(*TempClLL);
	//printMatrix(*mtxTimes);
	//printf("PeriodicX %d PeriodicY %d PeriodicZ %d rateN %g\n",PeriodicX,PeriodicY,PeriodicZ,*rateN);

	//Now need to multiply each mtxProbNeigh by the appropriate hoprate
	//Only for the sites that are within the clusters need to cycle 
	//through nodes again. 

	//The matrix mtxTimes only stores the hopping times for
	//hops that are off the cluster.

	//mtxDwellTime is counts the time it takes for a charge
	//to move from a site within a cluster to any other site
	//off or on. Considers all the hop rates. Then the matrix
	//is normalized so it becomes a ratio of probability

	//1 - Order of Magnitude of hop behind site i,j,k
	//2 - Order of Magnitude of hop infornt of site i,j,k
	//3 - Order of Magnitude of hop to site left of site i,j,k
	//4 - Order of Magnitude of hop to site right of site i,j,k
	//5 - Order of Magnitude of hop to site below site i,j,k
	//6 - Order of Magnitude of hop to site above site i,j,k
	//7 - Hop rate to site behind i,j,k
	//8 - Hop rate to site infront of i,j,k
	//9 - Hop rate to site left of site i,j,k
	//10 - Hop rate to site right of i,j,k
	//11 - Hop rate to site below i,j,k
	//12 - Hop rate to site above i,j,k
	int inc = 1;
	int inc2 = 1;
	double tempVal;
	double dwelltemp;
	double total;
	double tcluster;
	double Val;
	int Node_ID;
	int i,j,k;
	matrix mtxRates = newMatrix( getCluster_numNodes((*TempClLL)),1);
	matrix mtxDwellTime = newMatrix( getCluster_numNodes((*TempClLL)),1);
	Node tempNode = getStartNode((*TempClLL));
	total=0;

	Val = -1.0;

	while(tempNode!=NULL){

		dwelltemp=0;
		Node_ID = getNode_id(tempNode);	
		getLoc( &i, &j, &k, Node_ID, snA);

		if(getIndexFront(snA, Node_ID)!=-1 || PeriodicX!=0 ){
			if(getFlagFro(tempNode)!=1){
				//Hopping off the cluster
				tempVal = getE(MasterM, Node_ID+1,8);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 8\n",Node_ID, getE(MasterM, Node_ID+1,8),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			}
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,8);
		}

		if(getIndexBehind(snA, Node_ID)!=-1 || PeriodicX!=0 ){
			if(getFlagBeh(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,7);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 7\n",Node_ID, getE(MasterM, Node_ID+1,7),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			}
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,7);
		}

		if(getIndexLeft(snA, Node_ID)!=-1 || PeriodicY!=0 ){
			if(getFlagLef(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,9);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 9\n",Node_ID, getE(MasterM, Node_ID+1,9),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			}
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,9);
		}

		if(getIndexRight(snA, Node_ID)!=-1 || PeriodicY!=0 ){
			if(getFlagRig(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,10);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 10\n",Node_ID, getE(MasterM, Node_ID+1,10),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			}
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,10);
		}

		if(getIndexBelow(snA, Node_ID)!=-1 || PeriodicZ!=0 ){
			if(getFlagBel(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,11);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 11\n",Node_ID, getE(MasterM, Node_ID+1,11),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			} 
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,11);
		}

		//This is if it is within the boundaries or if it is periodic
		if(getIndexAbove(snA, Node_ID)!=-1 || PeriodicZ!=0 ){

			//If it is not a Neighbor
			if(getFlagAbo(tempNode)!=1){
				tempVal = getE(MasterM, Node_ID+1,12);
				Val+=tempVal*getE(mtxProbNeigh,inc,1);
				setE(mtxRates,inc,1,tempVal+getE(mtxRates,inc,1));
				//printf("Node_ID %d Matrix val %g Time %g 12\n",Node_ID, getE(MasterM, Node_ID+1,12),1/tempVal);
				setE((*mtxTimes),inc2,1,1/tempVal);
				inc2++;
			}		
			dwelltemp += 1/getE(MasterM, getIndex(snA,i,j,k)+1,12);
		}	
		total += dwelltemp;
		tempNode = getNextNode(tempNode);
		setE(mtxDwellTime,inc,1,dwelltemp);
		inc++;
	}

	if(Val<0){
		printf("ERROR Val is less than 0 this will result in tcluster that is negative\n");
		printf("Val %g\n",Val);
		exit(1);
	}

	*rateN = SumOfCol(mtxRates,1);
	printf("Val %g\n",Val);
	tcluster = 1/Val*1/20;
	setCluster_time((*TempClLL),tcluster);
	DivideEachElement(&mtxDwellTime,total);

	deleteMatrix(mtxRates);
	return mtxDwellTime;
}

int CalculatePvalNeigh(ClusterLL * TempClLL,const_matrix mtxTimes, const_matrix mtxProbNeighDwell){

	//Calculate Pval & time for NeighNodes 
	if(TempClLL==NULL || mtxTimes==NULL || mtxProbNeighDwell == NULL){
		return -1;
	}

	if((*TempClLL)==NULL){
		return -1;
	}
	NeighNode NeighNod = getStartNeigh((*TempClLL));
	double val=0;
	double time;
	int inc;
	int inc2;
	int ID;
	int Node_ID;

	inc =1;

	while(NeighNod!=NULL){

		inc2 =1;
		ID = getNeighNode_id(NeighNod);
		for(inc=1;inc<=getRows(mtxProbNeighDwell);inc++){
			Node_ID = getE(mtxProbNeighDwell,inc,2);
			if(ID==Node_ID){
				//printf("rows %d inc2 %d\n",getRows(mtxProbNeighDwell),inc2);
				val += getE(mtxProbNeighDwell,inc,1);
				//printf("\nSetting pval of Node %d with val %g\n",Node_ID, val);
				setNeighNodeNew_p(NeighNod,val);
				time = getE(mtxTimes,inc2,1);
				setNeighNode_t(NeighNod, time, inc2);
				inc2++;
			}
		}
		NeighNod = getNextNeigh(NeighNod);
	}

	return 0;
}

int CalculatePvalNodes(ClusterLL * TempClLL, matrix mtxProb, matrix mtxDwellTime){
	//Calculate Pval & time for Nodes

	//printf("Number of Nodes %d\n",getCluster_numNodes(*TempClLL));
	Node tempNode = getStartNode((*TempClLL));
	int inc=1;
	double val=0.0;
	double total=0;

	while(tempNode!=NULL){

		val = getE(mtxProb,inc,1)*getE(mtxDwellTime,inc,1);
		//printf("Value of val %g\n",val);
		total += val;
		setNode_p(tempNode,val);		
		inc++;
		tempNode = getNextNode(tempNode);
	}

	tempNode=getStartNode((*TempClLL));

	val = 0.0;
	inc = 1;
	//Normalize
	while(tempNode!=NULL){

		val += getNode_p(tempNode)/total;
		setNode_p(tempNode,val);	
		tempNode = getNextNode(tempNode);
	}

	return 0;
}

int ConnectClusterSN( int TotalOrders, SNarray snA, ArbArray ClArLL){

	if(snA==NULL || ClArLL==NULL){
		return -1;
	}
	if( snA==NULL){
		return -1;
	}

	int element;
	int i, j, k;
	int Node_ID;
	ClusterLL TempClLL;
	Node tempNode;
	SiteNode sn;

	for(element=0;element<TotalOrders;element++){
		TempClLL = (ClusterLL) getArbElement(ClArLL,element);

		if(TempClLL==NULL) {
			//printf("Cluster Empty\n");
		} else {

			while (TempClLL!=NULL){

				tempNode=(Node) getStartNode(TempClLL);

				while (tempNode!=NULL){

					Node_ID = getNode_id(tempNode);	
					getLoc( &i, &j, &k, Node_ID, snA);

					// Add pointers to SiteNode
					sn = getSN(snA, i, j, k);
					setDataStruc( sn,1, (void *)TempClLL);

					tempNode=getNextNode(tempNode);
					//tempNeighNode = getNextNeigh( tempNeighNode );
				}

				TempClLL=getNextClusterLL(TempClLL);
			}
		}
	}
	return 0;
}

int FindCluster( int * OrderL, SNarray snA, double electricEnergyX,\
		double electricEnergyY, double electricEnergyZ,\
		ArbArray *ClArLL1, double kT,ParameterFrame PF){


	double reOrgEnergy = PFget_reOrg(PF);
	int Attempts = PFget_Attempts(PF);
	double AttemptToHop = PFget_AttemptToHop(PF);
	double gamma = PFget_gamma(PF);
	double SiteDistance = PFget_SiteDist(PF);
	double PeriodicX = PFget_Px(PF);
	double PeriodicY = PFget_Py(PF);
	double PeriodicZ = PFget_Pz(PF);
	int XElecOn = PFget_XElecOn(PF);
	int YElecOn = PFget_YElecOn(PF);
	int ZElecOn = PFget_ZElecOn(PF);

	matrix MasterM = CalculateAllHops(snA, electricEnergyX, electricEnergyY, electricEnergyZ,\
			kT,reOrgEnergy, SiteDistance, AttemptToHop, gamma,\
			PeriodicX, PeriodicY, PeriodicZ);
	//Compare mid points
	//Some cases exist where there is a hop off the cluster that is on the same
	//order of magnitude to hops with in the cluster. If the hop back to the 
	//cluster from this site is faster than the hop off it might still be 
	//considered a cluster. Must account for this case. The mid points where 
	//this occurs are assigned the order of magnitude corresponding to the 
	//mid points smallest order of magnitude
	printf("\nComparing Mid Points\n");
	int MidPtsTotal;
	int orderLow;
	int orderHigh;
	int rv;
	ArbArray mpA = MPsort(&orderLow, &orderHigh, &MidPtsTotal, MasterM, snA, PeriodicX, PeriodicY, PeriodicZ);

	printf("\nSorting Mid Points into Link Lists\n");
	//First need to create link lists Each list contains
	//all the mid points that are the same order of magnitude
	//magnitude share sites
	//Total Order of magnitude range
	//Will create a link list for each order of magnitude
	int TotalOrders=orderHigh-orderLow+1;
	printf("Total Orders %d\n",TotalOrders);

	ArbArray ArLL = SortOrderMag(TotalOrders, orderLow, mpA );

	//Next we need to cycle through the link list and see if there are midpoints
	//in the same list that share neighbors if there are we need to make them a 
	//cluster.
	//Each cluster will be composed of a list of the sites in the cluster. 

	ArbArray ClArLL = ClusterSort( TotalOrders, orderLow, ArLL);
	//rv = ArbArrayCheck(ClArLL);

	//Some cases exist where there is a hop off the cluster that is on the same
	//order of magnitude to hops with in the cluster. If the hop back to the 
	//cluster from this site is faster than the hop off it might still be 
	//considered a cluster. Must account for this case. 

	printf("Now that the nodes are all organized into their respective clusters\n");
	printf("we will check to see if they are actually lead to charge traps.\n");

	rv = FilterCluster( TotalOrders,orderLow, MasterM, &ClArLL, snA,\
			PeriodicX, PeriodicY, PeriodicZ,\
			XElecOn, YElecOn, ZElecOn);

	printf("Finished removing extraneous clusters from Array!\n");
	//rv = ArbArrayCheck(ClArLL);

	//At this point ClArLL only contains clusters where charges
	//have one orders of magnitude greater probability to stay in
	//the cluster

	//printing out the order of magnitudes of the hops to all the neighbors to ensure this is correct
	//PrintCheck(TotalOrders, orderLow, ClArLL, snA, MasterM);
	rv = ArbArrayCheck(ClArLL);

	printf("Order High %d Order Low %d\n",orderHigh,orderLow);
	//Print out one cluster to determine what is going on

	//At this point we need to figure out which sites are neighboring the Clusters
	//Each Cluster should have a list of neighboring sites
	rv = CalculateNeighNodes(TotalOrders, orderLow, &ClArLL, snA,\
			PeriodicX, PeriodicY, PeriodicZ );
	rv = ArbArrayCheck(ClArLL);

	//printf("Printing Neighbor Array.\n");
	//printArbArray(ClArLL, orderLow);

	printf("Calculating Sum and P of clusters\n");
	CalculateSumAndP( TotalOrders, snA, &ClArLL,  MasterM, Attempts,\
			PeriodicX, PeriodicY, PeriodicZ);
	rv = ArbArrayCheck(ClArLL);

	printf("Connecting Clusters to SiteNodes\n");
	ConnectClusterSN(TotalOrders, snA, ClArLL);

	printf("Deleting MasterM matrix\n");
	deleteMatrix(MasterM);
	rv = ArbArrayCheck(ClArLL);

	printf("Deleting ArLL without deleting mid points\n");
	rv=deleteArbArray(&ArLL);
	rv = ArbArrayCheck(ClArLL);
	printf("\nDeleting mpA and mid points\n");
	rv=deleteAllMidPointArray(&mpA);
	*ClArLL1=ClArLL;
	rv = ArbArrayCheck(*ClArLL1);

	*OrderL = orderLow;
	return 0;
}

int SaveCluster( char * FileName, int OrderL, SNarray snA, double electricEnergyX,\
								double electricEnergyY, double electricEnergyZ,\
								ArbArray ClArLL1, double kT, ParameterFrame PF,\
								Electrode elXb, Electrode elXf, Electrode elYl,\
								Electrode elYr, Electrode elZb, Electrode elZa){

	if(snA==NULL || ClArLL1 == NULL || PF==NULL){
		printf("WARNING could not save Cluster file\n");
		return -1;
	}

	int i,j,k;
	struct stat st = {0};
	char str[20];
	char bufCheck[256];
	FILE * ClusterOut;


	if(stat("CLUSTERFILE",&st)==-1){
		mkdir("CLUSTERFILE",0700);
	}

	snprintf(bufCheck, sizeof bufCheck,"%s%s%s","CLUSTERFILE/",FileName,".cluster");

	if((ClusterOut=fopen(bufCheck,"w"))==NULL){
		printf("ERROR! unable to write to .cluster file!\n");
		exit(1);
	}else{


		printf("Starting to print to file\n");
		//Print electricEnergies, kT and OrderLow first
		fprintf(ClusterOut,"%g\t%g\t%g\t%g\t%d\n",electricEnergyX,electricEnergyY,electricEnergyZ,kT,OrderL);
		//Print parameter data second
		fprintf(ClusterOut,"//Number of nodes along the x-axis, y-axis and z-axis\n");
		fprintf(ClusterOut,"SLength %d\n",getAlen(snA));
		fprintf(ClusterOut,"SWidth %d\n",getAwid(snA));
		fprintf(ClusterOut,"SHeight %d\n\n",getAhei(snA));

		fprintf(ClusterOut,"//Defines whether sides are periodic or not\n");
		fprintf(ClusterOut,"// 0 - means non-periodic\n");
		fprintf(ClusterOut,"// 1 - means periodic\n");
		fprintf(ClusterOut,"PeriodicX %d\n",PFget_Px(PF));
		fprintf(ClusterOut,"PeriodicY %d\n",PFget_Py(PF));
		fprintf(ClusterOut,"PeriodicZ %d\n\n",PFget_Pz(PF));

		fprintf(ClusterOut,"//Defines over how many sample the system is periodic\n");
		fprintf(ClusterOut,"// 0 - means it is an infinite sample\n");
		fprintf(ClusterOut,"// 1 - means periodic across 1 sample (not-periodic)\n");
		fprintf(ClusterOut,"// 2 - periodic across 2 samples then reaches the edge\n");
		fprintf(ClusterOut,"//     3, 4, 5, ... etc\n");
		fprintf(ClusterOut,"EndX %d\n",PFget_EndX(PF));
		fprintf(ClusterOut,"EndY %d\n",PFget_EndY(PF));
		fprintf(ClusterOut,"EndZ %d\n\n",PFget_EndZ(PF));

		fprintf(ClusterOut,"//Defines whether electrodes exist on the x,y or z axis\n");
		fprintf(ClusterOut,"// 0 - no electrodes\n");
		fprintf(ClusterOut,"// 1 - yes electrodes\n");
		fprintf(ClusterOut,"XElecOn %d\n",PFget_XElecOn(PF));
		fprintf(ClusterOut,"YElecOn %d\n",PFget_YElecOn(PF));
		fprintf(ClusterOut,"ZElecOn %d\n\n",PFget_ZElecOn(PF));

		fprintf(ClusterOut,"//Define the tunneling constant from the electrode [1/nm]\n");
		fprintf(ClusterOut,"alphaXb %f\n",PFget_alphaxb(PF)/1E9);
		fprintf(ClusterOut,"alphaXf %f\n",PFget_alphaxf(PF)/1E9);
		fprintf(ClusterOut,"alphaYl %f\n",PFget_alphayl(PF)/1E9);
		fprintf(ClusterOut,"alphaYr %f\n",PFget_alphayr(PF)/1E9);
		fprintf(ClusterOut,"alphaZb %f\n",PFget_alphazb(PF)/1E9);
		fprintf(ClusterOut,"alphaZa %f\n\n",PFget_alphaza(PF)/1E9);

		fprintf(ClusterOut,"//Relative permittivity of medium\n");
		fprintf(ClusterOut,"RelativePermittivity %f\n\n",PFget_RelativePerm(PF));

		fprintf(ClusterOut,"//Define the attempt to hop rate from the electrode [1/s]\n");
		fprintf(ClusterOut,"vX %g\n",PFget_vX(PF));
		fprintf(ClusterOut,"vY %g\n",PFget_vY(PF));
		fprintf(ClusterOut,"vZ %g\n\n",PFget_vZ(PF));

		fprintf(ClusterOut,"//Define Fermi Energy level of the electrodes so we\n");
		fprintf(ClusterOut,"//can calculate the rates on and off electrodes [eV]\n");
		fprintf(ClusterOut,"//If the Electrodes are not turned on this value will\n");
		fprintf(ClusterOut,"//not be used\n");
		fprintf(ClusterOut,"XFermiB %f\n",PFget_XFermiB(PF));
		fprintf(ClusterOut,"XFermiF %f\n",PFget_XFermiF(PF));
		fprintf(ClusterOut,"YFermiL %f\n",PFget_YFermiL(PF));
		fprintf(ClusterOut,"YFermiR %f\n",PFget_YFermiR(PF));
		fprintf(ClusterOut,"ZFermiB %f\n",PFget_ZFermiB(PF));
		fprintf(ClusterOut,"ZFermiA %f\n\n",PFget_ZFermiA(PF));

		fprintf(ClusterOut,"//Voltage across a single system in x, y and z [V]\n");
		fprintf(ClusterOut,"VoltageX %g\n",PFget_VoltageX(PF));
		fprintf(ClusterOut,"VoltageY %g\n",PFget_VoltageY(PF));
		fprintf(ClusterOut,"VoltageZ %g\n\n",PFget_VoltageZ(PF));

		fprintf(ClusterOut,"//Voltage ramp steps x, y and z\n");\
			fprintf(ClusterOut,"VStepX %d\n",PFget_VStepX(PF));
		fprintf(ClusterOut,"VStepY %d\n",PFget_VStepY(PF));
		fprintf(ClusterOut,"VStepZ %d\n\n",PFget_VStepZ(PF));

		fprintf(ClusterOut,"//Voltage ramp increment x, y and z [V]\n");
		fprintf(ClusterOut,"VincX %g\n",PFget_VincX(PF));
		fprintf(ClusterOut,"VincY %g\n",PFget_VincY(PF));
		fprintf(ClusterOut,"VincZ %g\n\n",PFget_VincZ(PF));

		fprintf(ClusterOut,"//Distance between nodes in units of [m]\n");
		fprintf(ClusterOut,"SiteDistance %g\n\n",PFget_SiteDist(PF));

		fprintf(ClusterOut,"//Dimenion to fit real data\n");
		fprintf(ClusterOut,"D %g\n\n",PFget_D(PF));

		fprintf(ClusterOut,"//Total number of time steps in which charges are injected\n");
		fprintf(ClusterOut,"TCount %d\n\n",PFget_TCount(PF));

		fprintf(ClusterOut,"//Number of charges injected per time step\n");
		fprintf(ClusterOut,"NCh %d\n\n",PFget_NCh(PF));

		fprintf(ClusterOut,"//Total number of charges passed through the system\n");
		fprintf(ClusterOut,"Ntot TCount*NCh\n\n");

		fprintf(ClusterOut,"//Time between injections of charges [s]\n");
		fprintf(ClusterOut,"TStep %g\n\n",PFget_TStep(PF));

		fprintf(ClusterOut,"//Number of time steps that pass before data is averaged\n");
		fprintf(ClusterOut,"Nstep_av %d\n\n",PFget_Nstep_av(PF));

		fprintf(ClusterOut,"//Number of time steps before a checkpoint file is created\n");
		fprintf(ClusterOut,"Time_check %d\n\n",PFget_Time_check(PF));

		fprintf(ClusterOut,"//Number of iterations with different random seeds\n");
		fprintf(ClusterOut,"Rcount %d\n\n",PFget_Rcount(PF));

		fprintf(ClusterOut,"//Defines the cutoff radius for the correlation function\n");
		fprintf(ClusterOut,"//CutOffDistance CutOff*lambda\n");
		fprintf(ClusterOut,"//lambda is used in the correlation function\n");
		fprintf(ClusterOut,"CutOff %f\n",PFget_CutOff(PF));
		fprintf(ClusterOut,"lambda %g\n\n",PFget_lambda(PF));

		fprintf(ClusterOut,"//Seed Protocol is used to determine how the energies are\n");
		fprintf(ClusterOut,"//spead out between the sites.\n");
		fprintf(ClusterOut,"// 0 - means averaged by surrounding seeds\n");
		fprintf(ClusterOut,"// 1 - means averaged with closest seeds\n");
		fprintf(ClusterOut,"SeedProt %d\n",PFget_SeedProt(PF));
		fprintf(ClusterOut,"//How many numerical iterations should be used to approximate\n");
		fprintf(ClusterOut,"//cluster behavior (Default to 15%% should be less than 5%% error)\n");
		fprintf(ClusterOut,"Attempts %d\n\n",PFget_Attempts(PF));

		fprintf(ClusterOut,"//What fraction of sites act as seeds\n");
		fprintf(ClusterOut,"fracSeed %f\n",PFget_FracSeed(PF));
		fprintf(ClusterOut,"//Average site energy [eV]\n");
		fprintf(ClusterOut,"E0 %f\n",PFget_E0(PF));
		fprintf(ClusterOut,"//Standard deviation of the site energy\n");
		fprintf(ClusterOut,"sigma %f\n\n",PFget_sigma(PF));

		fprintf(ClusterOut,"//What fraction of sites act as traps\n");
		fprintf(ClusterOut,"fracTrap %f\n",PFget_FracTrap(PF));
		fprintf(ClusterOut,"//Average site energy of trap [eV]\n");
		fprintf(ClusterOut,"Etrap %f\n",PFget_Etrap(PF));
		fprintf(ClusterOut,"//Standard deviation of trap energy\n");
		fprintf(ClusterOut,"Tsigma %f\n\n",PFget_Tsigma(PF));

		fprintf(ClusterOut,"//Attempt to hop rate [1/s] equates to Markus at 300K\n");
		fprintf(ClusterOut,"AttemptToHop %g\n",PFget_AttemptToHop(PF));
		fprintf(ClusterOut,"//Temperature [k]\n");
		fprintf(ClusterOut,"TempStart %f\n",PFget_TempStart(PF));
		fprintf(ClusterOut,"//Temperature steps\n");
		fprintf(ClusterOut,"TemperatureStep %d\n",PFget_TempStep(PF));
		fprintf(ClusterOut,"//Temperature increment [K]\n");
		fprintf(ClusterOut,"TemperatureInc %f\n",PFget_TempInc(PF));
		fprintf(ClusterOut,"//Re-organization energy [eV]\n");
		fprintf(ClusterOut,"reOrgEnergy %f\n",PFget_reOrg(PF));
		fprintf(ClusterOut,"//Tunneling constant [1/nm]\n");
		fprintf(ClusterOut,"gamma %f\n",PFget_gamma(PF)/1E9);
		fprintf(ClusterOut,"MovieFrames %d\n\n",PFget_MovieFrames(PF));

		printf("Finished Parameter Frame to File\n");

		SNarray snAtemp;
		matrix mtx;
		SiteNode sn;
		Point pt;
		ClusterLL ClLL;
		int ClusterID;
		//We know how many electrodes there should be by looking
		//at XElecOn YElecOn and ZElecOn in the parameter frame
		if(PFget_XElecOn(PF)==1){
			if(elXb==NULL || elXf==NULL){
				printf("ERROR electrodes Xb and Xf are NULL when trying to save\n");
				exit(1);
			}
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elXb),getElectrode_alpha(elXb),\
					   getElectrode_Charges(elXb),getElectrode_FermiEnergy(elXb));

			mtx = (matrix) getElectrode_HopRates(elXb);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			snAtemp = (SNarray) getElectrode_AdjacentSites(elXb);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


					}
				}
			}
			
			fprintf(ClusterOut,"\n");
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elXf),getElectrode_alpha(elXf),\
					   getElectrode_Charges(elXf),getElectrode_FermiEnergy(elXf));

			mtx = (matrix) getElectrode_HopRates(elXf);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			snAtemp = (SNarray) getElectrode_AdjacentSites(elXf);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


					}
				}
			}
		}


		if(PFget_YElecOn(PF)==1){
			if(elYl==NULL || elYr==NULL){
				printf("ERROR electrodes Yl and Yr are NULL when trying to save\n");
				exit(1);
			}
			
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elYl),getElectrode_alpha(elYl),\
					   getElectrode_Charges(elYl),getElectrode_FermiEnergy(elYl));

			mtx = (matrix) getElectrode_HopRates(elYl);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			snAtemp = (SNarray) getElectrode_AdjacentSites(elYl);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


					}
				}
			}
			
			fprintf(ClusterOut,"\n");
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elYr),getElectrode_alpha(elYr),\
					   getElectrode_Charges(elYr),getElectrode_FermiEnergy(elYr));

			mtx = (matrix) getElectrode_HopRates(elYr);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			SNarray snAtemp = (SNarray) getElectrode_AdjacentSites(elYr);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


					}
				}
			}
		}
		if(PFget_ZElecOn(PF)==1){
			if(elZb==NULL || elZa==NULL){
				printf("ERROR electrodes Zb and Za are NULL when trying to save\n");
				exit(1);
			}
			
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elZb),getElectrode_alpha(elZb),\
					   getElectrode_Charges(elZb),getElectrode_FermiEnergy(elZb));

			mtx = (matrix) getElectrode_HopRates(elZb);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			snAtemp = (SNarray) getElectrode_AdjacentSites(elZb);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


					}
				}
			}
			
			fprintf(ClusterOut,"\n");
			//Sum alpha Charges FermiEnergy
			fprintf(ClusterOut,"%g\t%g\t%d\t%g\n",getElectrode_Sum(elZa),getElectrode_alpha(elZa),\
					   getElectrode_Charges(elZa),getElectrode_FermiEnergy(elZa));

			mtx = (matrix) getElectrode_HopRates(elZa);
			
			//print rates
			for(i=0;i<getRows(mtx);i++){
				for(j=0;j<getCols(mtx);j++){
					//print i j rate
					fprintf(ClusterOut,"%d %d %g\n",i,j,getE(mtx,i+1,j+1));
				}
			}
			
			fprintf(ClusterOut,"\n");

			SNarray snAtemp = (SNarray) getElectrode_AdjacentSites(elZa);

			for(i=0;i<getAlen(snAtemp);i++){
				for(j=0;j<getAwid(snAtemp);j++){
					for(k=0;k<getAhei(snAtemp);k++){
						
						sn = getSN(snAtemp,i,j,k);

						ClLL = (ClusterLL) getClusterList(sn);
						if(ClLL!=NULL){
							ClusterID = getCluster_id(ClLL);
						}else{
							ClusterID = -1;
						}

						pt = (Point) getPoint(sn);
						if(pt==NULL){
							printf("ERROR point should not be NULL if sn exists\n");
							exit(1);
						}

						fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));

					}
				}
			}
		}

		int totalElem;
		int totalNodes;
		int totalNeighNodes;
		int totalHops;
		NeighNode NeiN;
		Node nd;

		fprintf(ClusterOut,"%d\n\n",getAtotal(snA));
		for(i=0;i<getAlen(snA);i++){
			for(j=0;j<getAwid(snA);j++){
				for(k=0;k<getAhei(snA);k++){

					sn = getSN(snA,i,j,k);

					ClLL = (ClusterLL) getClusterList(sn);
					if(ClLL!=NULL){
						ClusterID = getCluster_id(ClLL);
					}else{
						ClusterID = -1;
					}

					pt = (Point) getPoint(sn);
					if(pt==NULL){
						printf("ERROR point should not be NULL if sn exists\n");
						exit(1);
					}

					fprintf(ClusterOut,"%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i,j,k,getEnergy(sn),getDwelStat(sn),getVisFreq(sn),getVis(sn),getTime(sn),ClusterID,getsum(sn),getSN_p(sn,0),getSN_p(sn,1),getSN_p(sn,2),getSN_p(sn,3),getSN_p(sn,4),getSN_p(sn,5));


				}
			}
		}

		totalElem = getElementsReserved(ClArLL1);
		int totalClusters=0;
		fprintf(ClusterOut,"\n%d\n",totalElem);

		for(i=0;i<totalElem;i++){
			ClLL = (ClusterLL) getArbElement(ClArLL1,i);
			if(ClLL==NULL){
				//For this type of ArbArray none of the elements should be NULL
				printf("ArbArray element is NULL in SaveCluster function\n");
				fprintf(ClusterOut,"0\n");
			} else{

				while(ClLL!=NULL){
					totalClusters++;
					ClLL = getNextClusterLL(ClLL);
				}
				ClLL = (ClusterLL) getArbElement(ClArLL1,i);

				printf("TotalClusters %d\n",totalClusters);
				fprintf(ClusterOut,"%d\n",totalClusters);

				while (ClLL!=NULL){
					totalNodes = getCluster_numNodes(ClLL);	
					totalNeighNodes = getCluster_numNeigh(ClLL);
					//ID of Cluster followed by the cluster p val, time and sum 
					fprintf(ClusterOut,"%d\t%g\t%g\t%g\n",getCluster_id(ClLL),getCluster_time(ClLL),\
							getCluster_Sum(ClLL),getCluster_p(ClLL));
					fprintf(ClusterOut,"%d\n",totalNodes);
					//Print node id, pval of nodes and flags of the nodes

					nd = getStartNode(ClLL);

					while(nd!=NULL){
						fprintf(ClusterOut,"%d\t%g\t%d\t%d\t%d\t%d\t%d\t%d\n",getNode_id(nd),getNode_p(nd),getFlag(nd,0),getFlag(nd,1),getFlag(nd,2),getFlag(nd,3),getFlag(nd,4),getFlag(nd,5));

						nd = getNextNode(nd);
					}
					fprintf(ClusterOut,"%d\n",totalNeighNodes);

					NeiN = getStartNeigh(ClLL);

					while(NeiN!=NULL){
						//Id of node, hoplength of node	

						totalHops = getNeighNode_hoplength(NeiN);
						fprintf(ClusterOut,"%d\t%d\n",getNeighNode_id(NeiN),totalHops);

						for(k=1;k<=totalHops;k++){
							fprintf(ClusterOut,"%g\t%g\n",getNeighNode_p(NeiN,k), getNeighNode_t(NeiN,k));
						}

						NeiN = getNextNeigh(NeiN);
					}

					ClLL = getNextClusterLL(ClLL);
				}
			}
		}
	}

	fclose(ClusterOut);
	return 0;

}



int LoadCluster( char * FileName, int * OrderL, SNarray * snA, double *electricEnergyX,\
		double * electricEnergyY, double * electricEnergyZ,\
		ArbArray *ClArLL1, double * kT,ParameterFrame * PF){

	int totalSites;
	int intvar;
	double doublevar;
	struct stat st = {0};
	char bufCluster[256];
	FILE * ClusterIn;

	//Check if Directory exists
	if(stat("CLUSTERFILE",&st)==-1){
		printf("CLUSTERFILE directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufCluster, sizeof bufCluster, "%s%s%s","CLUSTERFILE/",FileName,".cluster");
		if((ClusterIn=fopen(bufCluster,"r"))==NULL){
			printf("WARNING: unable to open chosen .cluster file!\n");
			return -1;
		}else{
		
				*PF = newParamFrame();
				printf("Loading Parameter Frame from .cluster file\n");
				//Scanning in Parameters
				//Scanning in Parameters
				fscanf(ClusterIn,"%lf %lf %lf %lf %d",electricEnergyX,electricEnergyY,electricEnergyZ,kT,OrderL);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Len(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Wid(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Hei(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Px(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Py(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Pz(*PF,intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_EndX(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_EndY(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_EndZ(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_XElecOn(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_YElecOn(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_ZElecOn(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphaxb(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphaxf(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphayl(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphayr(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphazb(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_alphaza(*PF,doublevar*1E9);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_RelativePerm(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_vX(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_vY(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_vZ(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_XFermiB(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_XFermiF(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_YFermiL(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_YFermiR(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_ZFermiB(*PF, doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_ZFermiA(*PF, doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VoltageX(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VoltageY(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VoltageZ(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_VStepX(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_VStepY(*PF,intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_VStepZ(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VincX(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VincY(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_VincZ(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_SiteDist(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_D(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_TCount(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_NCh(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_TStep(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Nstep_av(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Time_check(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Rcount(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_CutOff(*PF,doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_lambda(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_SeedProt(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_Attempts(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_FracSeed(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_E0(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_sigma(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_FracTrap(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_Etrap(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_Tsigma(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_AttemptToHop(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_TempStart(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_TempStep(*PF,intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_TempInc(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_reOrg(*PF,doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				PFset_gamma(*PF,doublevar*1E9);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				PFset_MovieFrames(*PF,intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);
				totalSites = intvar;
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				printf("Len %d Wid %d Hei %d\n",PFget_Len(*PF),PFget_Wid(*PF),PFget_Hei(*PF));
				*snA = newSNarray(PFget_Len(*PF),PFget_Wid(*PF),PFget_Hei(*PF));
				printf("Creating SiteNode array\n");

				int count;
				int i;
				int j;
				int k;
				double Energy;
				int Dstat;
				int VisF;
				int Vis;
				double Tim;
				int ClID;
				double Sum;
				double P0,P1,P2,P3,P4,P5;
				
				SiteNode sn;
				matrix mtx = newMatrix(totalSites,1);
				printf("Loading Site Node data from .cluster file\n");

				printf("totalSites %d\n",totalSites);
				for(count=0;count<totalSites;count++){
				
					fscanf(ClusterIn,"%d\t%d\t%d\t%lf\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",\
														&i, &j, &k,&Energy,&Dstat,&VisF,&Vis,&Tim,&ClID,&Sum,&P0,&P1,&P2,&P3,&P4,&P5);
					sn = getSN(*snA,i,j,k);
					setInitE(sn,1);
					setEnergy(sn,Energy);
					setDwelStat(sn,Dstat);
					setVisFreq(sn,VisF);
					setVis(sn,Vis);
					setTime(sn,Tim);
					setsum(sn,Sum);
					setSN_p(sn,0,P0);
					setSN_p(sn,1,P1);
					setSN_p(sn,2,P2);
					setSN_p(sn,3,P3);
					setSN_p(sn,4,P4);
					setSN_p(sn,5,P5);
					setE(mtx,getIndex(*snA,i,j,k)+1,1,ClID);
				}


				ClusterLL ClLL;
				Node nd;
				NeighLL NeiLL;
				NeighNode NeighNod;
				NeighNode NeighNod2;
				Hop h;
				int count2;
				int count3;
				int count4;
				int NodeID;
				double NodeP;
				int F0, F1, F2, F3, F4, F5;
				int totalElem;
				int ClusterID;
				int NeighID;
				int totalClusters;
				int totalNodes;
				int totalNeigh;
				int totalHops;
				int clusterNum;
				double ClusterTime;
				double ClusterSum;
				double ClusterP;
				double NeighP;
				double NeighT;

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);	
				totalElem = intvar;
				printf("Creating New ArbArray\n");
				*ClArLL1 = newArbArray(totalElem,1);
				
				for (count=0;count<totalElem;count++){
					fscanf(ClusterIn,"%d",&totalClusters);
					printf("TotalNumber of Clusters %d\n",totalClusters);
					for(clusterNum=0;clusterNum<totalClusters;clusterNum++){
						fscanf(ClusterIn,"%d %lf %lf %lf",&ClusterID,&ClusterTime,&ClusterSum,&ClusterP);
						printf("ClusterID %d\n",ClusterID);
						ClLL = newClusterLL(ClusterID);
						printf("1\n");
						setClusterSum(ClLL,ClusterSum);
						printf("2\n");
						setCluster_time(ClLL,ClusterTime);
						printf("3\n");
						setCluster_p(ClLL,ClusterP);
						printf("TotalNodes %d\n",totalNodes);
						fscanf(ClusterIn,"%d",&totalNodes);
						setArbElement(ClArLL1, count, (void *) ClLL);
						printf("Setting element of ArbArray\n");		
						//Connect with the site nodes
						for(count4=0;count4<getRows(mtx);count4++){
						
							printf("Connecting Sitenodes with Clusters\n");
							if(getE(mtx,count4+1,1)==ClusterID){
								//Need to connect the SiteNode with the cluster
								//The index of the siteNode = count4+1
								sn = getSNwithInd(*snA,count4);
								setDataStruc(sn,1,(void *) (ClLL));
							}
						}

						printf("Loading Nodes\n");
						for(count2=0;count2<totalNodes;count2++){
						
							fscanf(ClusterIn,"%d %lf %d %d %d %d %d %d",&NodeID,&NodeP,&F0,&F1,&F2,&F3,&F4,&F5);
							nd = newNode(NodeID);
							setNode_p(nd,NodeP);
							if(F0==1){
								setFlagBeh(nd);
							}
							if(F1==1){
								setFlagFro(nd);
							}
							if(F2==1){
								setFlagLef(nd);
							}
							if(F3==1){
								setFlagRig(nd);
							}
							if(F4==1){
								setFlagBel(nd);
							}
							if(F5==1){
								setFlagAbo(nd);
							}
							addClusterLLNode(&ClLL,nd);
						}

						printf("Creating NeighLL and loading Neighbor nodes\n");
						fscanf(ClusterIn,"%d",&intvar);
						totalNeigh = intvar;
						NeiLL = newNeighLL();
						setCluster_NeighLL(ClLL, NeiLL);
						for(count2=0;count2<totalNeigh;count2++){
							fscanf(ClusterIn,"%d %d",&NeighID,&totalHops);
							if(count2==0){
								NeighNod = newNeighNode(NeighID);
								setNeighLL_start(NeiLL, NeighNod);
								setNeighLL_numNeigh(NeiLL, totalNeigh);
							}else{
								NeighNod2	= newNeighNode(NeighID);
								setNextNeighNode(&NeighNod,&NeighNod2);
								NeighNod = getNextNeigh(NeighNod);
							}
							setNeighNode_hoplength(NeighNod,totalHops);

							for(count3=0;count3<totalHops;count3++){
								fscanf(ClusterIn,"%lf %lf",&NeighP,&NeighT);
								h = newHop();
								setHop_t(h,NeighT);
								setHop_p(h,NeighP);
								if(count3==0){
									setNeighNode_hopstart(NeighNod,h);
								}
							}
						}
					}
				}


				deleteMatrix(mtx);

		}
	
	}

	return 0;
}


int LoadCluster_Data( char * FileName, int * OrderL, SNarray * snA, double electricEnergyX,\
		double electricEnergyY, double electricEnergyZ,\
		ArbArray *ClArLL1, double kT){

	int totalSites;
	int intvar;
	double doublevar;
	struct stat st = {0};
	char bufCluster[256];
	FILE * ClusterIn;

	int len;
	int wid;
	int hei;

	double Ex, Ey, Ez, T;

	//Check if Directory exists
	if(stat("CLUSTERFILE",&st)==-1){
		printf("CLUSTERFILE directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufCluster, sizeof bufCluster, "%s%s%s","CLUSTERFILE/",FileName,".cluster");
		if((ClusterIn=fopen(bufCluster,"r"))==NULL){
			printf("WARNING: unable to open chosen .cluster file!\n");
			return -1;
		}else{
		
				printf("Loading Parameter Frame from .cluster file\n");
				//Scanning in Parameters
				//Scanning in Parameters
				fscanf(ClusterIn,"%lf %lf %lf %lf %d", &Ex,&Ey,&Ez,&T,OrderL);

				if(Ex!=electricEnergyX || Ey!=electricEnergyY ||\
					 Ez!=electricEnergyZ || kT!=T){

					printf("ERROR the cluster file does not have the same parameters\n");
					printf("      as the current simulation run\n");
					exit(1);

				}

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&len);
				fscanf(ClusterIn,"%s %d",bufCluster,&wid);
				fscanf(ClusterIn,"%s %d",bufCluster,&hei);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);
				totalSites = intvar;
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				printf("Len %d Wid %d Hei %d\n",len,wid,hei);
				*snA = newSNarray(len,wid,hei);
				printf("Creating SiteNode array\n");

				int count;
				int i;
				int j;
				int k;
				double Energy;
				int Dstat;
				int VisF;
				int Vis;
				double Tim;
				int ClID;
				double Sum;
				double P0,P1,P2,P3,P4,P5;
				
				SiteNode sn;
				matrix mtx = newMatrix(totalSites,1);
				printf("Loading Site Node data from .cluster file\n");

				printf("totalSites %d\n",totalSites);
				for(count=0;count<totalSites;count++){
				
					fscanf(ClusterIn,"%d\t%d\t%d\t%lf\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",\
														&i, &j, &k,&Energy,&Dstat,&VisF,&Vis,&Tim,&ClID,&Sum,&P0,&P1,&P2,&P3,&P4,&P5);
					sn = getSN(*snA,i,j,k);
					setInitE(sn,1);
					setEnergy(sn,Energy);
					setDwelStat(sn,Dstat);
					setVisFreq(sn,VisF);
					setVis(sn,Vis);
					setTime(sn,Tim);
					setsum(sn,Sum);
					setSN_p(sn,0,P0);
					setSN_p(sn,1,P1);
					setSN_p(sn,2,P2);
					setSN_p(sn,3,P3);
					setSN_p(sn,4,P4);
					setSN_p(sn,5,P5);
					setE(mtx,getIndex(*snA,i,j,k)+1,1,ClID);
				}


				ClusterLL ClLL;
				Node nd;
				NeighLL NeiLL;
				NeighNode NeighNod;
				NeighNode NeighNod2;
				Hop h;
				int count2;
				int count3;
				int count4;
				int NodeID;
				double NodeP;
				int F0, F1, F2, F3, F4, F5;
				int totalElem;
				int ClusterID;
				int NeighID;
				int totalClusters;
				int totalNodes;
				int totalNeigh;
				int totalHops;
				int clusterNum;
				double ClusterTime;
				double ClusterSum;
				double ClusterP;
				double NeighP;
				double NeighT;

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);	
				totalElem = intvar;
				printf("Creating New ArbArray\n");
				*ClArLL1 = newArbArray(totalElem,1);
				
				for (count=0;count<totalElem;count++){
					fscanf(ClusterIn,"%d",&totalClusters);
					printf("TotalNumber of Clusters %d\n",totalClusters);
					for(clusterNum=0;clusterNum<totalClusters;clusterNum++){
						fscanf(ClusterIn,"%d %lf %lf %lf",&ClusterID,&ClusterTime,&ClusterSum,&ClusterP);
						printf("ClusterID %d\n",ClusterID);
						ClLL = newClusterLL(ClusterID);
						printf("1\n");
						setClusterSum(ClLL,ClusterSum);
						printf("2\n");
						setCluster_time(ClLL,ClusterTime);
						printf("3\n");
						setCluster_p(ClLL,ClusterP);
						printf("TotalNodes %d\n",totalNodes);
						fscanf(ClusterIn,"%d",&totalNodes);
						setArbElement(ClArLL1, count, (void *) ClLL);
						printf("Setting element of ArbArray\n");		
						//Connect with the site nodes
						for(count4=0;count4<getRows(mtx);count4++){
						
							printf("Connecting Sitenodes with Clusters\n");
							if(getE(mtx,count4+1,1)==ClusterID){
								//Need to connect the SiteNode with the cluster
								//The index of the siteNode = count4+1
								sn = getSNwithInd(*snA,count4);
								setDataStruc(sn,1,(void *) (ClLL));
							}
						}

						printf("Loading Nodes\n");
						for(count2=0;count2<totalNodes;count2++){
						
							fscanf(ClusterIn,"%d %lf %d %d %d %d %d %d",&NodeID,&NodeP,&F0,&F1,&F2,&F3,&F4,&F5);
							nd = newNode(NodeID);
							setNode_p(nd,NodeP);
							if(F0==1){
								setFlagBeh(nd);
							}
							if(F1==1){
								setFlagFro(nd);
							}
							if(F2==1){
								setFlagLef(nd);
							}
							if(F3==1){
								setFlagRig(nd);
							}
							if(F4==1){
								setFlagBel(nd);
							}
							if(F5==1){
								setFlagAbo(nd);
							}
							addClusterLLNode(&ClLL,nd);
						}

						printf("Creating NeighLL and loading Neighbor nodes\n");
						fscanf(ClusterIn,"%d",&intvar);
						totalNeigh = intvar;
						NeiLL = newNeighLL();
						setCluster_NeighLL(ClLL, NeiLL);
						for(count2=0;count2<totalNeigh;count2++){
							fscanf(ClusterIn,"%d %d",&NeighID,&totalHops);
							if(count2==0){
								NeighNod = newNeighNode(NeighID);
								setNeighLL_start(NeiLL, NeighNod);
								setNeighLL_numNeigh(NeiLL, totalNeigh);
							}else{
								NeighNod2	= newNeighNode(NeighID);
								setNextNeighNode(&NeighNod,&NeighNod2);
								NeighNod = getNextNeigh(NeighNod);
							}
							setNeighNode_hoplength(NeighNod,totalHops);

							for(count3=0;count3<totalHops;count3++){
								fscanf(ClusterIn,"%lf %lf",&NeighP,&NeighT);
								h = newHop();
								setHop_t(h,NeighT);
								setHop_p(h,NeighP);
								if(count3==0){
									setNeighNode_hopstart(NeighNod,h);
								}
							}
						}
					}
				}


				deleteMatrix(mtx);

		}
	
	}

	return 0;
}

int LoadCluster_Only( char * FileName, int * OrderL, SNarray * snA, double electricEnergyX,\
		double electricEnergyY, double electricEnergyZ,\
		ArbArray *ClArLL1, double kT){

	if(snA==NULL){
		printf("ERROR snA can not be NULL to use the LoadCluster_Only function\n");
		exit(1);
	}

	if(*snA==NULL){
		printf("ERROR *snA can not be NULL to use the LoadCluster_Only function\n");
		printf("			it must have been previously initialized\n");
		exit(1);
	}

	int totalSites;
	int intvar;
	double doublevar;
	struct stat st = {0};
	char bufCluster[256];
	FILE * ClusterIn;

	int len;
	int wid;
	int hei;

	double Ex, Ey, Ez, T;

	//Check if Directory exists
	if(stat("CLUSTERFILE",&st)==-1){
		printf("CLUSTERFILE directory does not exist!\n");
		return -1;
	}else{
		snprintf(bufCluster, sizeof bufCluster, "%s%s%s","CLUSTERFILE/",FileName,".cluster");
		if((ClusterIn=fopen(bufCluster,"r"))==NULL){
			printf("WARNING: unable to open chosen .cluster file!\n");
			return -1;
		}else{
		
				printf("Loading Parameter Frame from .cluster file\n");
				//Scanning in Parameters
				//Scanning in Parameters
				fscanf(ClusterIn,"%lf %lf %lf %lf %d", &Ex,&Ey,&Ez,&T,OrderL);

				if(Ex!=electricEnergyX || Ey!=electricEnergyY ||\
					 Ez!=electricEnergyZ || kT!=T){

					printf("ERROR the cluster file does not have the same parameters\n");
					printf("      as the current simulation run\n");
					exit(1);

				}

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&len);
				fscanf(ClusterIn,"%s %d",bufCluster,&wid);
				fscanf(ClusterIn,"%s %d",bufCluster,&hei);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%s %lf",bufCluster,&doublevar);
				fscanf(ClusterIn,"%s %d",bufCluster,&intvar);

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);
				totalSites = intvar;
				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);

				printf("Creating SiteNode array\n");

				int count;
				int i;
				int j;
				int k;
				double Energy;
				int Dstat;
				int VisF;
				int Vis;
				double Tim;
				int ClID;
				double Sum;
				double P0,P1,P2,P3,P4,P5;
				
				SiteNode sn;
				matrix mtx = newMatrix(totalSites,1);
				printf("Loading Site Node data from .cluster file\n");

				printf("totalSites %d\n",totalSites);
				for(count=0;count<totalSites;count++){
				
					fscanf(ClusterIn,"%d\t%d\t%d\t%lf\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",\
														&i, &j, &k,&Energy,&Dstat,&VisF,&Vis,&Tim,&ClID,&Sum,&P0,&P1,&P2,&P3,&P4,&P5);
					
					if(Energy!=getEnergy(getSN(*snA,i,j,k))){
						printf("ERROR Energies in .cluster file do not line up with energies already in snA\n");
						exit(1);
					}
					setE(mtx,getIndex(*snA,i,j,k)+1,1,ClID);
				}


				ClusterLL ClLL;
				Node nd;
				NeighLL NeiLL;
				NeighNode NeighNod;
				NeighNode NeighNod2;
				Hop h;
				int count2;
				int count3;
				int count4;
				int NodeID;
				double NodeP;
				int F0, F1, F2, F3, F4, F5;
				int totalElem;
				int ClusterID;
				int NeighID;
				int totalClusters;
				int totalNodes;
				int totalNeigh;
				int totalHops;
				int clusterNum;
				double ClusterTime;
				double ClusterSum;
				double ClusterP;
				double NeighP;
				double NeighT;

				fgets(bufCluster,256,ClusterIn);
				fgets(bufCluster,256,ClusterIn);
				fscanf(ClusterIn,"%d",&intvar);	
				totalElem = intvar;
				printf("Creating New ArbArray\n");
				*ClArLL1 = newArbArray(totalElem,1);
				
				for (count=0;count<totalElem;count++){
					fscanf(ClusterIn,"%d",&totalClusters);
					printf("TotalNumber of Clusters %d\n",totalClusters);
					for(clusterNum=0;clusterNum<totalClusters;clusterNum++){
						fscanf(ClusterIn,"%d %lf %lf %lf",&ClusterID,&ClusterTime,&ClusterSum,&ClusterP);
						printf("ClusterID %d\n",ClusterID);
						ClLL = newClusterLL(ClusterID);
						printf("1\n");
						setClusterSum(ClLL,ClusterSum);
						printf("2\n");
						setCluster_time(ClLL,ClusterTime);
						printf("3\n");
						setCluster_p(ClLL,ClusterP);
						printf("TotalNodes %d\n",totalNodes);
						fscanf(ClusterIn,"%d",&totalNodes);
						setArbElement(ClArLL1, count, (void *) ClLL);
						printf("Setting element of ArbArray\n");		
						//Connect with the site nodes
						for(count4=0;count4<getRows(mtx);count4++){
						
							printf("Connecting Sitenodes with Clusters\n");
							if(getE(mtx,count4+1,1)==ClusterID){
								//Need to connect the SiteNode with the cluster
								//The index of the siteNode = count4+1
								sn = getSNwithInd(*snA,count4);
								setDataStruc(sn,1,(void *) (ClLL));
							}
						}

						printf("Loading Nodes\n");
						for(count2=0;count2<totalNodes;count2++){
						
							fscanf(ClusterIn,"%d %lf %d %d %d %d %d %d",&NodeID,&NodeP,&F0,&F1,&F2,&F3,&F4,&F5);
							nd = newNode(NodeID);
							setNode_p(nd,NodeP);
							if(F0==1){
								setFlagBeh(nd);
							}
							if(F1==1){
								setFlagFro(nd);
							}
							if(F2==1){
								setFlagLef(nd);
							}
							if(F3==1){
								setFlagRig(nd);
							}
							if(F4==1){
								setFlagBel(nd);
							}
							if(F5==1){
								setFlagAbo(nd);
							}
							addClusterLLNode(&ClLL,nd);
						}

						printf("Creating NeighLL and loading Neighbor nodes\n");
						fscanf(ClusterIn,"%d",&intvar);
						totalNeigh = intvar;
						NeiLL = newNeighLL();
						setCluster_NeighLL(ClLL, NeiLL);
						for(count2=0;count2<totalNeigh;count2++){
							fscanf(ClusterIn,"%d %d",&NeighID,&totalHops);
							if(count2==0){
								NeighNod = newNeighNode(NeighID);
								setNeighLL_start(NeiLL, NeighNod);
								setNeighLL_numNeigh(NeiLL, totalNeigh);
							}else{
								NeighNod2	= newNeighNode(NeighID);
								setNextNeighNode(&NeighNod,&NeighNod2);
								NeighNod = getNextNeigh(NeighNod);
							}
							setNeighNode_hoplength(NeighNod,totalHops);

							for(count3=0;count3<totalHops;count3++){
								fscanf(ClusterIn,"%lf %lf",&NeighP,&NeighT);
								h = newHop();
								setHop_t(h,NeighT);
								setHop_p(h,NeighP);
								if(count3==0){
									setNeighNode_hopstart(NeighNod,h);
								}
							}
						}
					}
				}


				deleteMatrix(mtx);

		}
	
	}

	return 0;
}
