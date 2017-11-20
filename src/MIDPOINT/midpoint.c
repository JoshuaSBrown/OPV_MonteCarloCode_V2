
#include <stdio.h>
#include <stdlib.h>
#include "midpoint.h"

#include "../ERROR/error.h"
struct _MidPoint{
  int id;
  int orderMag; //log(HopRate)
  //Neigbors between the mid points
  int nei1;
  int nei2;
  MidPoint next;
};

MidPoint newMidPoint(int order,int Mid_ID, int nei1, int nei2){

  if( nei1<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nei1 is less than 0 in newMidPoint\n");
    #endif
    return NULL;
  } 
  if( nei2<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nei2 is less than 0 in newMidPoint\n");
    #endif
    return NULL;
  }
  if(Mid_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Mid_ID is less than 0 in newMidPoint\n");
    #endif
    return NULL;
  }
  if(nei1==nei2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nei1==nei2 in newMidPoint\n");
    #endif
    return NULL;
  }
  MidPoint mp = (MidPoint) malloc(sizeof(struct _MidPoint));

  if(mp==NULL ){
    return NULL;
  }

  mp->orderMag=order;
  mp->nei1=nei1;
  mp->nei2=nei2;
  mp->next=NULL;
  mp->id=Mid_ID;
  return mp;
} 

int deleteMidPoint(MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in deleteMidPoint\n");
    #endif 
    return -1;
  }
  free(mp);
  return 0;
}

int printMP(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in printMP");
    #endif
    return -1;
  }
  printf("Id: %d \t OrderMag %d \t nei1: %d \t nei2: %d\n",mp->id,mp->orderMag, mp->nei1, mp->nei2);
  return 0;
}

MidPoint getNextMP(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getNext\n");
    #endif
    return NULL;
  }

  return mp->next;
}

int getMP_order(const_MidPoint mp){
  if(mp==NULL){ //added condition
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getMP_order\n");
    #endif
    return -1;
  }
  return mp->orderMag;
}

int getMP_nei1(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getMP_nei1\n");
    #endif
    return -1;
  }
  return mp->nei1;
}

int getMP_nei2(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getMP_nei2\n");
    #endif
    return -1;
  }
  return mp->nei2;
}

int getMP_id(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getMP_id\n");
    #endif
    return -1;
  }
  return mp->id;
}

MidPoint getMP_next(const_MidPoint mp){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in getMP_next\n");
    #endif
    return NULL;
  }

  return mp->next;
}

int setMP_id(MidPoint mp, int ID){
  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in setMP_id\n");
    #endif
    return -1;
  }
  if(ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ID is less than 0 in setMP_id\n");
    #endif
    return -1;
  }
  mp->id=ID;
  return 0;
}

int setMP_nextMP(MidPoint mp, MidPoint nextmp){

  if(mp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp is NULL in setMP_nextMP\n");
    #endif
    return -1;
  }

  if(nextmp==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR nextmp is NULL in setMP_nextMP\n");
    #endif
    return -1;
  }

  mp->next = nextmp;
  return 0;
}

int CompareNeiMidPoint(const_MidPoint mp1, MidPoint mp2) {

  if(mp1==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp1 is NULL in CompareNeiMidPoint\n");
    #endif
    return -1;
  }
  if(mp2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp2 is NULL in CompareNeiMidPoint\n");
    #endif
    return -1;
  }

  if(mp1->nei1==mp2->nei1){
    return 1;
  }
  if(mp1->nei2==mp2->nei2){
    return 1;
  }
  if(mp1->nei1==mp2->nei2){
    return 1;
  }
  if(mp1->nei2==mp2->nei1){
    return 1;
  }

  return 0;
}

int addMPToEnd(MidPoint mp1, MidPoint mp2){

  if(mp1==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp1 is NULL in addMPToEnd\n");
    #endif
    return -1;
  }
  if(mp2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR mp2 is NULL in addMPToEnd\n");
    #endif
    return -1;
  }

  if(getMP_id(mp1)==getMP_id(mp2)){ 
    return -2;
  }

  MidPoint tempmp;
  tempmp = mp1;
  while(tempmp->next){
    if(tempmp->id==mp2->id){
      return -2;
    }
    tempmp=tempmp->next;
  }
  tempmp->next=mp2;
  return 0;
}
