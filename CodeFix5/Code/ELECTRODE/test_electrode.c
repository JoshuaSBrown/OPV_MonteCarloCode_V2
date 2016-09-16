#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "electrode.h"

int main() {
	int rv;
	double rvd;
	Electrode el;

  //////////////////////////////////////////////////////////////////////

	printf("Testing:newElectrode\n");
	el = newElectrode();
	assert(el!=NULL);

	printf("Testing:deleteElectrode\n");
	rv = deleteElectrode(NULL);
	assert(rv==-1);
	rv = deleteElectrode(&el);
	assert(rv==0);

	printf("Testing:printElectrode\n");
	rv = printElectrode(NULL);
	assert(rv==-1);
	el = newElectrode();
	printf("Should print p of 0.0\n");
	rv = printElectrode(el);
	assert(rv==0);

	printf("Testing:getElectrode_alpha\n");
	assert(getElectrode_alpha(NULL)==-1);
	assert(getElectrode_alpha(el)==0);

	printf("Testing:setElectrode_alpha\n");
	rv = setElectrode_alpha(NULL, 34);
	assert(rv==-1);
	rv = setElectrode_alpha(el, 12);
	assert(rv==0);
	assert(getElectrode_alpha(el)==12);

	printf("Testing:getElectrode_Charges\n");
	rv = getElectrode_Charges(NULL);
	assert(rv==-1);
	rv = getElectrode_Charges(el);
	assert(rv==0);

	printf("Testing:setElectrode_Charges\n");
	rv = setElectrode_Charges(NULL, 1);
	assert(rv==-1);
	rv = setElectrode_Charges(el, -1);
	assert(rv==-1);
	rv = setElectrode_Charges(el, 32);
	assert(rv==0);
	assert(getElectrode_Charges(el)==32);

	printf("Testing:Electrode_addCharge\n");
	rv = Electrode_addCharge(NULL);
	assert(rv==-1);
	rv = Electrode_addCharge(el);
	assert(rv==0);
	assert(getElectrode_Charges(el)==33);

	printf("Testing:Electrode_minusCharge\n");
	rv = Electrode_minusCharge(NULL);
	assert(rv==-1);
	rv = Electrode_minusCharge(el);
	assert(rv==0);
	assert(getElectrode_Charges(el)==32);

	printf("Testing:getElectrode_FermiEnergy\n");
	rvd = getElectrode_FermiEnergy(NULL);
	assert(rvd==-1);
	rvd = getElectrode_FermiEnergy(el);
	assert(rvd==0.0);

	printf("Testing:setElectrode_FermiEnergy\n");
	rv = setElectrode_FermiEnergy(NULL, 2.3);
	assert(rv==-1);
	rv = setElectrode_FermiEnergy(el, -1.1);
	assert(rv==0);
	assert(getElectrode_FermiEnergy(el)==-1.1);

	printf("Testing:getElectrode_Sum\n");
	rv = getElectrode_Sum(NULL);
	assert(rv==-1);
	rv = getElectrode_Sum(el);
	assert(rv==0);

	printf("Testing:setElectrode_Sum\n");
	rv=setElectrode_Sum(NULL,15.4);
	assert(rv==-1);
	rv=setElectrode_Sum(el,-1);
	assert(rv==-1);
	rv=setElectrode_Sum(el,15.4);
	assert(rv==0);
	assert(getElectrode_Sum(el)==15.4);

	printf("Testing:getElectrode_HopRates\n");

	printf("Testing:setElectrode_HopRates\n");
  void * Hops = (void *) 5;
	assert(setElectrode_HopRates(NULL, Hops)==-1);
	assert(setElectrode_HopRates(el, NULL)==-1);
	assert(setElectrode_HopRates(el,Hops ) == 0); 

	printf("Testing:getElectrode_AdjacentSites\n");
	assert(getElectrode_AdjacentSites(NULL)==NULL);
	assert(getElectrode_AdjacentSites(el)==NULL);

	printf("Testing:setElectrode_AdjacentSites\n"); 
  void * AdjacentSites = (void *) 6;
	assert(setElectrode_AdjacentSites(NULL, AdjacentSites)==-1);
	assert(setElectrode_AdjacentSites(el,NULL)==-1);
	assert(setElectrode_AdjacentSites(el, AdjacentSites )==0);
	assert(getElectrode_AdjacentSites(el)!=NULL);
  deleteElectrode(&el);
}
