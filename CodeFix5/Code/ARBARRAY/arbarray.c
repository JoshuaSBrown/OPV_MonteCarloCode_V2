#include <stdio.h>
#include <stdlib.h>

#include "arbarray.h"
#include "../MIDPOINT/midpoint.h"
#include "../OMLL/omll.h"
#include "../CLUSTER/cluster.h"
#include "../ERROR/error.h"

struct _ArbArray{
	int type;
	int reserved;
	int used;
	void * Num[0];
};

/////////////////////////////////////////////////////////////////////////
//Tools for accessing ArbArray
ArbArray newArbArray(int len, int type){

  #ifdef _ERROR_CHECKING_ON_
	if(len<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR len is less than 0 in newArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(type>=3){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR type is greater than 3 in newArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(type<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR type is less than 0 in newArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	ArbArray Arb;	

	//Each element contains an Order of Magnitude Link List
	Arb = (ArbArray) malloc(sizeof(struct _ArbArray)+sizeof(void *)*len);

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	Arb->reserved=len;
	Arb->used=0;
	Arb->type=type;
	int i;
	for(i=0;i<len;i++){
		Arb->Num[i]=NULL;
	}

	return Arb;
}

int ArbArrayCheck(const_ArbArray Arb){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in ArbArrayCheck\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	printf("Arb Type Check %d\n",Arb->type);
	return 0;
}

int deleteAllMidPointArray(ArbArray * Arb){

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in deleteAllMidPointArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if((*Arb)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *Arb is NULL in deleteAllMidPointArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
		
	}
  #endif
	
	if ((*Arb)->type==2){
		
		int i;
		for(i=0;i<(*Arb)->reserved;i++){
			if((*Arb)->Num[i]!=NULL){
				deleteMidPoint((*Arb)->Num[i]);
			}
		}
		free(*Arb);
		*Arb=NULL;
		return 0;
	}
	printf("ERROR not a midpoint array!\n");
	return -1;
}


int deleteArbArray(ArbArray * Arb){

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in deleteArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if((*Arb)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *Arb is NULL in deleteArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//OrderMagLL
  OrderMagLL OMLL;

	if((*Arb)->type==0){
		int i;
		for(i=0;i<(*Arb)->reserved;i++){
			if ((*Arb)->Num[i]!=NULL){
				OMLL = (OrderMagLL) (*Arb)->Num[i];
        deleteOrLLwithOutDeletingMP(&OMLL);
				(*Arb)->Num[i]=NULL;
			}
		}
		free((*Arb));
		*Arb=NULL;
		return 0;
		//ClusterLL
	} else if ((*Arb)->type==1) {
		//Must Cycle through all the elements in the array and get rid of 
		//all the clusters first
		int i;
		ClusterLL clLL;
		for(i=0;i<(*Arb)->reserved;i++){
			if ((*Arb)->Num[i]!=NULL){
				clLL = (*Arb)->Num[i];
				deleteAllClusterLL(&clLL );
				(*Arb)->Num[i]=NULL;
			}
		}
		free(*Arb);
		*Arb = NULL;
		return 0;
	} else if ((*Arb)->type==2){
		
		//WARNING we can not delete the midpoints here
		//because more than one arbitrary array may use them
		//They need to be delted externally
		/*
		int i;
		for(i=0;i<Arb->reserved;i++){
			if(Arb->Num[i]!=NULL){
				deleteMidPoint(Arb->Num[i]);
			}
		}
		*/
		free(*Arb);
		*Arb=NULL;
		return 0;
	}

	return -1;
}

void * getArbElement(const_ArbArray Arb, int element){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is less than 0 in getArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is larger than the elements reserved in Arb in getArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return Arb->Num[element];
}

int NullArbElement(ArbArray Arb, int element){

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in NullArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is less than 0 in NullArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is larger than reserved in Arb in NullArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	if(Arb->Num[element]!=NULL){
		Arb->used--;
	}
	Arb->Num[element]=NULL;

	return 0;

}

int appendArbElement(ArbArray * Arb, void * ptr){

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in appendArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
	if((*Arb)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *Arb is NULL in appendArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
	if( ptr==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ptr is NULL in appendArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif

  int len = (*Arb)->reserved+1;
  //Create new ArbArray with one extra element
  ArbArray ArbNew = newArbArray( len, (*Arb)->type);

  //Copy over elements from old ArbArray to new one  
  int i;
  for(i=0;i<(*Arb)->reserved;i++){
    ArbNew[i]=(*Arb)[i];
  }
  
  //free just the old arbarray data structure without 
  //actually deleting the elements in the arb array
  free(*Arb);
  //Set the pointer equal to the newArbArray
  (*Arb)=ArbNew;
  return 0;
}

int setArbElement(ArbArray Arb, int element, void * ptr){

  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in setArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is less than 0 in setArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is larger than reserved in Arb in setArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(ptr==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ptr is NULL in setArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	if(Arb->Num[element]!=NULL){
		Arb->Num[element]=ptr;
	}else{
		Arb->used++;
		Arb->Num[element]=ptr;
	}
	return 0;
}

int printArbArray(const_ArbArray Arb, ...) {
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in printArbArray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	if(Arb->type==0) {
		int i;
		printf("Printing Array of Link Lists\n");
		printf("Total number of Linked Lists: %d\n",Arb->used);
		for(i=0;i<Arb->used;i++) {
			printOrLL((OrderMagLL) Arb->Num[i]);
		}	
	}else if(Arb->type==1) {
		va_list argument;
		va_start (argument, Arb);
		int orderLow = va_arg(argument, int);
		int element;
		ClusterLL tempclLL;
		for(element=0;element<Arb->reserved;element++){
			printf("\nOrder of Magnitude %d\n",element+orderLow);
			tempclLL = (ClusterLL) getArbElement(Arb, element);
			if (tempclLL==NULL){
				printf("Cluster does not exist.\n");
			} else {
				printClusterLL(tempclLL);
			}
		}
		va_end(argument);
	}else if(Arb->type==2) {
		int i;
		for(i=0;i<Arb->used;i++){
			printMP((MidPoint) Arb->Num[i]);
		}
	}

	return 0;
}

int getElementsUsed(const_ArbArray Arb){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getElementsUsed\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Arb->used;
}

int getElementsReserved(const_ArbArray Arb){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getElementReserved\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Arb->reserved;
}

int removeArbElement(ArbArray Arb, int element){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in removeArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is less than 0 in removeArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is greater than reserved in removeArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	if(Arb->Num[element]!=NULL){
		Arb->used--;
		Arb->Num[element]=NULL;
		return 0;
	}
	return -1;
}

MidPoint getMP(const_ArbArray Arb,int element) {
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getMP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if (Arb->type!=2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb must be of type 2 to use getMP function\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if ( element<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested element is less than 0 in getMP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if (element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is greater than the number of elements reserved in Arb in getMP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	MidPoint mp = (MidPoint) Arb->Num[element];
	return mp;
}

int getMPnei1(const_ArbArray Arb, int Mid_ID){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getMPnei1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( Mid_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Mid_ID is less than 0 in getMPnei1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Mid_ID>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Mid_ID is larger than the number of elements reserved by Arb in getMPnei1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Arb->type!=2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb needs to be of type 2 to call getMPnei1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if( mp!=NULL){
		return getMP_nei1(mp);
	}
	return -1;
}

int getMPnei2(const_ArbArray Arb, int Mid_ID){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getMPnei2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Mid_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Mid_ID is less than 0 in getMPnei2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Mid_ID>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR MID_ID is larger than the number of elements reserved in Arb in getMPnei2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Arb->type!=2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is of the wrong type to call getMPnei2 function\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif


	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if (mp!=NULL){
		return getMP_nei2(mp);
	}
	return -1;
}

int getMPOrder(const_ArbArray Arb, int Mid_ID){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getMPOrder\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Mid_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Mid_ID (Midpoint id) is less than 0 in getMPOrder\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Mid_ID>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested Mid_ID is larger than the elements reserved by Arb in getMPOrder\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Arb->type!=2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is of the wronge type to call function getMPOrder must be of type 2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	MidPoint mp = (MidPoint) Arb->Num[Mid_ID];
	if(mp!=NULL){
		return getMP_order(mp);
	}
	return -1;
}

int addToOrLL(ArbArray Arb, int element , MidPoint * mp){
	//Ensure that input parameters are error free
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested element is less than 0 in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested element is larger than the elements reserved by Arb in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( (*mp)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *mp is NULL in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Arb->type!=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb->type is no 0 in addToOrLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//increment the size of the order of magnitude link list
	OrderMagLL OMLL = (OrderMagLL) Arb->Num[element];
	if(OMLL==NULL){
		Arb->used++;
		int Order = getMP_order(*mp);
		Arb->Num[element]=(void *) newOrLL(Order);
		OMLL = (OrderMagLL) Arb->Num[element];	
	}
  return addMPToOMLL(OMLL, *mp);
}

OrderMagLL getOrderLL(const_ArbArray Arb, int element){
  #ifdef _ERROR_CHECKING_ON_
	if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in getOrderLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested element in getOrderLL is larger than the elements reserved to Arb\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if( element<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR requested element is less than 0 in getOrderLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(Arb->type!=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is of the wrong type when calling getOrderLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return (OrderMagLL) Arb->Num[element];
}

int setDefaultArbElem(ArbArray Arb, int element, int orderMag){
	
  #ifdef _ERROR_CHECKING_ON_
  if(Arb==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Arb is NULL in setDefaultArbElem\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is less than 0 in setDefaultArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(element>=Arb->reserved){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR element is larger than what is reserved in setDefaultArbElement\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Arb->type!=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR type is not equal to 0 in setDefaultArbElem\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

  OrderMagLL temp;
	//Setting default values
	if(Arb->Num[element]==NULL){
		Arb->used++;
	}else{
		//If not NULL need to remove
		//the old version
    temp = (OrderMagLL) Arb->Num[element];
		deleteOrLLwithOutDeletingMP(&temp);
	}
	Arb->Num[element]=newOrLL(orderMag);
	return 0;
}
