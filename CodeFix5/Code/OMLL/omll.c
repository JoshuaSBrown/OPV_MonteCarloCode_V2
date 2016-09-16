#include <stdio.h>
#include <stdlib.h>

#include "omll.h"
#include "../MIDPOINT/midpoint.h"

#include "../ERROR/error.h"

struct _OrderMagLL{
	int orderMag;
	int size;
	MidPoint start;
};

////////////////////////////////////////////////////////////////////////////////////
//Tools for accessing Order of Magnitude Link List
OrderMagLL newOrLL(int orderMag){
	OrderMagLL OMLL= (OrderMagLL) malloc(sizeof(struct _OrderMagLL));
	#ifdef _ERROR_CHECKING_ON_
  if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR newOrLL returned NULL from malloc\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
	}
  #endif
	OMLL->orderMag=orderMag;
	OMLL->size=0;
	OMLL->start=NULL;
	return OMLL;
}

int checkNewOrLL(OrderMagLL OMLL, MidPoint mp){ //use when adding newMidPoint to OMLL
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"OMLL in checkNewOrLL is NULL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1; 
	}
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"mp in checkNewOrLL is NULL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	while(mp != NULL){
		if(OMLL->orderMag != getMP_order(mp) ){
			return -1;
		}
		mp=getMP_next(mp);
	}
	return 0;
}

int deleteOrLL(OrderMagLL * OMLL){ //used only if no midpoints
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in deleteOrLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
	if(*OMLL==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in deleteOrLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
  if((*OMLL)->start!=NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR start of OMLL is not NULL in deleteOrLL");
    fprintf(stderr," this would lead to leaked memory\n");
    #endif
    return -1; 
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  }
  #endif
	free(*OMLL);
  *OMLL=NULL;
	return 0;
}

int deleteOrLLwithOutDeletingMP(OrderMagLL * OMLL){ 
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in deleteOrLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
	if(*OMLL==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in deleteOrLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
  #endif
 	free(*OMLL);
  *OMLL=NULL;
	return 0;
}

int deleteOMLL_startMP(OrderMagLL OMLL){
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
	  #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in deleteOMLL_startMP\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
  	return -1;
	}
  if(OMLL->start!=NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR start of OMLL is not NULL in deleteOMLL_startMP");
    fprintf(stderr,"this would lead to leaked memory\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1; 
  }
  #endif
  
  MidPoint mp = getOMLLstartMP(OMLL);
  if(getMP_next(mp)==NULL){
    deleteMidPoint(mp);
    OMLL->start = NULL;
  }else{
    setOMLL_startMP(OMLL, getMP_next(mp));
    deleteMidPoint(mp);
  }
  return 0;
}

int deleteAllOrLL(OrderMagLL OMLL){
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL in deleteAllOrLL is NULL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//Delete All the MidPoints in the Link List first
	MidPoint tempmp;
	MidPoint tempmp2;
	int rv;

	if (OMLL->start!=NULL) {
		tempmp=OMLL->start;
		if (getMP_next(tempmp)!=NULL){
			tempmp2=getMP_next(tempmp);
			while(tempmp2!=NULL){
				rv = deleteMidPoint(tempmp);
				if(rv==-1){
					return -1;
				}
				tempmp=tempmp2;
				tempmp2=getMP_next(tempmp2);
			}
		}
		rv = deleteMidPoint(tempmp);
		if (rv==-1){
			return -1;
		}
	}
	free(OMLL);
	return 0;
}

int printOrLL(const_OrderMagLL OMLL) {

	#ifdef _ERROR_CHECKING_ON_
	if (OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in printOrLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	MidPoint tempmp;	
	tempmp=OMLL->start;
	printf("\nOrderMag %d Size of LL: %d\n",OMLL->orderMag, OMLL->size);
	if (tempmp) {
		printMP(tempmp);
		while(getMP_next(tempmp)){
			tempmp=getMP_next(tempmp);
			printMP(tempmp);

		}
	}
	return 0;
}

int getOMLL_size(const_OrderMagLL OMLL){
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in getOMLL_size\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	return OMLL->size;
}

int getOMLL_order(const_OrderMagLL OMLL){
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in getOMLL_order\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	return OMLL->orderMag;
}

MidPoint getOMLLstartMP(const_OrderMagLL OMLL){
	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in getOMLLstartMP\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return OMLL->start;
}

int setOMLL_startMP(OrderMagLL OMLL, MidPoint mp){

	#ifdef _ERROR_CHECKING_ON_
	if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in setOMLL_startMP\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
 
	if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in setOMLL_startMP\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
  
  if(OMLL->start==NULL){
    OMLL->start = mp;
  }else{
    MidPoint mp2 = OMLL->start;
    OMLL->start = mp;
    setMP_nextMP(OMLL->start,mp2);
  }
  OMLL->size++;
  return 0;
}

int addMPToOMLL(OrderMagLL OMLL, MidPoint mp){

	#ifdef _ERROR_CHECKING_ON_
  if(OMLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR OMLL is NULL in addMPToOMLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1; 
  }
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in setOMLL_addMPToOMLL\n");
    #endif
		#ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1; 
  }
  #endif
  MidPoint tempmp;
  tempmp=getOMLLstartMP(OMLL);

  if (tempmp==NULL) {
    setOMLL_startMP(OMLL,mp);
    return 0;

  }
  //Cycle through the mid points to find the end of the link list and to
  //ensure that the mid point has not already been added. 
  int rv;
  rv = addMPToEnd(tempmp,mp);
  if(rv!=0){
    return rv;
  }else {
    OMLL->size++;
    return 0;
  }
}
