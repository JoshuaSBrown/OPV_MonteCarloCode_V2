#include <stdio.h>
#include <stdlib.h>

#include "electrode.h"
#include "../ERROR/error.h"

struct _Electrode{
	double sum;
	double alpha;
	int Charges;
	double FermiEnergy;
	void * HopRates;
	void * AdjacentSites;
};

//////////////////////////////////////////////////////////////
Electrode newElectrode(void){

	Electrode el = (Electrode) malloc(sizeof(struct _Electrode));

	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newElectrode\n");
    #endif
  	return NULL;
	}

	el->sum=0.0;
	el->alpha=0.0;
	el->Charges=0.0;
	el->FermiEnergy=0.0;
	el->HopRates = NULL;
	el->AdjacentSites = NULL;
	return el;
}

int deleteElectrode(Electrode * el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in deleteElectrode\n");
    #endif
		return -1;
	}
	if((*el)==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR *el is NULL in deleteElectrode \n");
    #endif
		return -1;
	}

	free(*el);
	*el=NULL;
	return 0;
}

int setElectrode_alpha(Electrode el, double alpha){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode\n");
    #endif
		return -1;
	}
  if(alpha<=0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR alpha is less than 0 in setElectrode\n");
    #endif
		return -1;
  }

	el->alpha = alpha;
	return 0;
}

double getElectrode_alpha(Electrode el){
  if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_alpha\n");
    #endif
    return -1.0;
  }
	return el->alpha;
}

int getElectrode_Charges(Electrode el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_Charges \n");
    #endif
		return -1;
	}

	return el->Charges;
}

int setElectrode_Charges(Electrode el, double NumCharge){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_Charges\n");
    #endif
		return -1;
	}
  if(NumCharge<0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR NumCharge is less than 0 in setElectrode_Charges\n");
    #endif
    return -1;
  }

	el->Charges = NumCharge;
	return 0;
}

int Electrode_addCharge(Electrode el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in Electrode_addCharge\n");
    #endif
		return -1;
	}
	el->Charges++;
	return 0;
}

int Electrode_minusCharge(Electrode el){
	if(el==NULL || el->Charges==0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in Electrode_minusCharge\n");
    #endif
		return -1;
	}
  if(el->Charges==0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR cannot minus charge when electrode has");
    fprintf(stderr," 0 charges in it in function Electrode_minusCharge\n");
    #endif
		return -1;
  }

	el->Charges--;
	return 0;
}


int setElectrode_FermiEnergy(Electrode el, double FermiE){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_FermiEnergy \n");
    #endif
		return -1;
	}

	el->FermiEnergy = FermiE;
	return 0;
}

double getElectrode_FermiEnergy(Electrode el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_FermiEnergy\n");
    #endif
		return -1;
	}
	return el->FermiEnergy;
}

int setElectrode_HopRates(Electrode el, void * HopS){
	if(el==NULL || HopS==NULL ){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_HopRates\n");
    #endif
		return -1;
	}
  if(HopS==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR HopS is NULL in setElectrode_HopRates\n");
    #endif
		return -1;
  }
	el->HopRates=HopS;
	return 0;
}

void * getElectrode_HopRates(Electrode el){ 
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_HopRates\n");
    #endif
		return NULL;
	}
	return el->HopRates;
}

int setElectrode_AdjacentSites(Electrode el, void * AdjacentSites){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_AdjacentSites\n");
    #endif
		return -1;
	}
  if(AdjacentSites==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR AdjacentSites is NULL in setElectrode_AdjacentSites\n");
    #endif
		return -1;
  }
	el->AdjacentSites = AdjacentSites;
	return 0;
}

void * getElectrode_AdjacentSites(Electrode el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_AdjacentSites\n");
    #endif
		return NULL;
	}
	return el->AdjacentSites;
}

int setElectrode_Sum(Electrode el, double sum){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_Sum\n");
    #endif
		return -1;
	}
  if(sum<0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR sum is less than 0 in setElectrode_Sum\n");
    #endif
		return -1;
  }
	el->sum=sum;
	return 0;
}

int setElectrode_AddToSum(Electrode el, double val){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in setElectrode_AddToSum\n");
    #endif
		return -1;
	}
  if((el->sum+val)<0){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR sum is less than 0 after adding val in setElectrode_AddToSum\n");
    #endif
		return -1;
  }
	el->sum+=val;
	return 0;
}

double getElectrode_Sum(Electrode el){
  if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in getElectrode_Sum\n");
    #endif
    return -1.0;
  }
	return el->sum;
}

int printElectrode(const_Electrode el){
	if(el==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR el is NULL in printElectrode\n");
    #endif
		return -1;
	}
	printf("Electrode sum %g\n",el->sum);
	return 0;
}

