#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "charge.h"
#include "../ERROR/error.h"
#include "../MATRIX_LINKLIST/matrix_linklist.h"
#include "../MATRIX/matrix.h"
#include "../CLUSTERSITENODE/clustersitenode.h"

struct _Charge {
	int x, y, z; 		    //current position
	int tot_x; 		      //the distance run by the charge in the x direction
	double t;  		      //total time the electron uses
	double dwelltime; 	//the time the electron will dwell on current site
  
  //path is used to keep track of where charges have hopped
  //The first column stores the id of the site the second 
	//column stores the number of times the charge has hopped
	//to the site
	matrix_linklist path;		
};				
  
struct _ChargeArray {
	int length;
	Charge C[0];
};

Charge newCharge(void) {
	Charge ch = (Charge) malloc( sizeof(struct _Charge ));
	#ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returns NULL in newCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
  }
  #endif
	//By default set values to 0
	ch->x=0;
	ch->y=0;
	ch->z=0;
	ch->t=0;
	ch->tot_x=0;
	ch->dwelltime=0;
  ch->path=NULL;
	return ch;
}

int initChargePath(Charge * ch, int NumNodes){
	#ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR charge is NULL in setChargePath\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(*ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR *ch is NULL in setChargePath\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(NumNodes<2){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR NumNodes is less than two in setChargePath ");
    fprintf(stderr,"Path must have a length of at least 2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if((*ch)->path!=NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR charge path has already been initiated ");
    fprintf(stderr,"cannot create a new one\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int i;
  //Default setting is for a linklist nodes with three elements
  double data[3];
  //First element stores the id of the site
  data[0] = -1;
  //Second element stores the number of times the charge
  //has visited that particular site. 
  data[1] = 0;
  //Third element stores the id of the cluster if it is part of 
  //one
  data[2] = -1;

  //Creating first node
	matrix_linklist path = newMatrix_LinkList( data, 3);
  //Append the rest of the nodes
  for(i=1;i<NumNodes;i++){
    addLL_MNode( &path, data, 3);
  }
  (*ch)->path = path;

  return 0;
}

int initChargeArrayPath(ChargeArray chA, int NumNodes){
	#ifdef _ERROR_CHECKING_ON_
  if(chA==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR charge array is NULL in initChargeArrayPath\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(NumNodes<2){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR NumNodes is less than two in initChargeArrayPath ");
    fprintf(stderr,"Path must have a length of at least 2\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  int i;
  int len;
  len = chA->length;
  Charge ch; 
  for(i=0;i<len;i++){
    ch = chA->C[i];
    initChargePath(&ch,NumNodes);
  }
  return 0;
}

ChargeArray newChargeA(int len) {
	#ifdef _ERROR_CHECKING_ON_
	if(len < 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR len is less than 0 in newCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
  }
  #endif
	ChargeArray chA=(ChargeArray) malloc( sizeof(struct _ChargeArray)+sizeof(Charge)*(len));
	#ifdef _ERROR_CHECKING_ON_
	if(chA==NULL){	//can't check on test script because created within function
		#ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newChargeA\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
  #endif
	chA->length=len;
	for(int i=0;i<len;i++){ 
		chA->C[i]=newCharge();
		if(chA->C[i]==NULL) //can't check
			return NULL; 
	}
	return chA;
}

int deleteCharge(Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL indeleteCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
  if(ch->path!=NULL){
	  deleteMatrixLL(&(ch->path)); 
  }
	free(ch);
	return 0; 
}

int deleteChargeA(ChargeArray chA) {
	#ifdef _ERROR_CHECKING_ON_
	if(chA==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR chA is NULL in deleteChargeA\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	int i;
	for(i=0;i<chA->length;i++){
		deleteCharge(chA->C[i]);
	}
	free(chA);
	return 0;
}

int getChargeA_len(ChargeArray chA){
	#ifdef _ERROR_CHECKING_ON_
	if(chA==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR chA is NULL in deleteChargeA_len\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return chA->length;
}

Charge getCharge(const_ChargeArray chA, int i){
	#ifdef _ERROR_CHECKING_ON_
	if(chA==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR chA is NULL in getCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
	if(i>= chA->length ){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to chA->length in getCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
	if(i < 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in getCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
  #endif
	return chA->C[i];
}

matrix_linklist getChargePath(Charge ch){
  
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getChargePath\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
	if(ch->path==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch->path is NULL in getChargePath\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL; 
  }
  #endif
  return ch->path;
}

int printCharge(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in printCharge\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	printf("\nx: %d y: %d z: %d\n",ch->x,ch->y,ch->z);
	printf("tot_x: %d\n",ch->tot_x);
	printf("t: %g\n",ch->t);
	printf("dwelltime: %g\n",ch->dwelltime);
  if(ch->dwelltime<=0){
    printf("ERROR dwelltime less than or equal to 0\n");
    exit(1);
  }
  if(ch->path!=NULL){
    printMatrixLL(ch->path);
  }else{
    printf("ch->path is NULL\n");
  }
  printf("Finished with printCharge\n");
	return 0; 
}

int printChargeA(const_ChargeArray chA){
	#ifdef _ERROR_CHECKING_ON_
	if(chA==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR chA is NULL in printChargeA\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
		printf("\nCharge Array\n");
	printf("Length: %d\n",chA->length);
	for(int i=0;i<(chA->length);i++){
		printf("\nCharge %d\n",i);
		printCharge(chA->C[i]);
	}
	return 0; 
}

double getDwel(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getDwel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
	return ch->dwelltime;
}

int getCx(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getCx\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	return ch->x;
}

int getCy(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getCy\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	return ch->y;
}

int getCz(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getCz\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	return ch->z;
}

int setDwel(Charge ch, double dw) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in setDwel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(dw < 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR dw is less than 0 in setDwel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->dwelltime=dw;
	return 0;
}

int MinusDwel(Charge ch, double dw){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in MinusDwel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(dw < 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR dw is less than 0 in MinusDwel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->dwelltime-=dw;
	return 0; 
}

int sett(Charge ch, long double t) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in sett\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(t<= 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR t is less than or equal to 0 in sett\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->t=t;
	return 0; 
}

int Plust(Charge ch, long double t){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in Plust\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(t<=0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR t is less than or equal to 0 in Plust\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->t+=t;
	return 0; 
}

double gett(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in gett\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	return ch->t;
}
int CxPlus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CxPlus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->x++;
	ch->tot_x++;
	return 0; 
}

int CxMinus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CxMinus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->x--;
	ch->tot_x--;
	return 0; 
}

int CyPlus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CyPlus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->y++;
	return 0;
}

int CyMinus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CyMinus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->y--;
	return 0;
}

int CyPlusPeriodic(Charge ch,int w){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CyPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(w <=0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR w is less than or equal to 0 in CyPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->y++;
	ch->y=(ch->y%w);
	if(ch->y>=w) //can't check
		return -2;
	return 0; 
}

int CyMinusPeriodic(Charge ch, int w) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CyMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(w <= 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR w is less than or equal to 0 in CyMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->y--;
	ch->y=(ch->y+w)%w;
	if( ch->y >= w || ch->y<=-1) //can't check
		return -2; 
	return 0;
}

int CzPlusPeriodic(Charge ch,int h){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CzPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(h <= 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR h is less than or equal to 0 in CzPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->z++;
	ch->z=(ch->z%h);
	if(ch->z>=h) //can't check
		return -2;
	return 0; 
}

int CzMinusPeriodic(Charge ch, int h) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CzMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(h <=0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR h is less than or equal to 0 in CzMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->z--;
	ch->z=(ch->z+h)%h;
	if(ch->z>=h || ch->z==-1) //can't check
		return -2; 
	return 0;
}

int CxPlusPeriodic(Charge ch,int l){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CxPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(l <= 0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR l is less than or equal to 0 in CxPlusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->x++;
	ch->tot_x++;
	ch->x=(ch->x%l);
	if(ch->x>=l) //can't check
		return -2;
	return 0; 
}

int CxMinusPeriodic(Charge ch, int l) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CxMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
	if(l <=0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR l is less than or equal to 0 in CxMinusPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->x--;
	ch->tot_x--;
	ch->x=(ch->x+l)%l;
	if(ch->x>=l || ch->x==-1) //can't check
		return -2; 
	return 0;
}

int CzPlus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CzPlus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->z++;
	return 0;
}

int CzMinus(Charge ch){
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in CzMinus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->z--;
	return 0;
}

int setCx(Charge ch, int cx) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL ){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in setCx\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->x=cx;
	ch->tot_x = cx;
	return 0;
}
int setCy(Charge ch, int cy) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL ){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in setCy\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 	
  }
  #endif
	ch->y=cy;
	return 0;
}
int setCz(Charge ch, int cz) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL ){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in setCz\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	ch->z=cz;
	return 0;
}
int getXdist(const_Charge ch) {
	#ifdef _ERROR_CHECKING_ON_
	if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR ch is NULL in getXdist\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1; 
  }
  #endif
	return ch->tot_x;
}

double getChargePathVisits(Charge ch, int seq){
	#ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR can not getChargePathVisits "); 
		fprintf(stderr,"charge is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(seq<1){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR can not getChargePathVisits ");
    fprintf(stderr,"sequence value is less than 1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  double visits = getMatrixLLNodeElem(ch->path,seq,2);
  return visits;
}

int resetChargePathVisit(Charge ch, int visit){
 
	#ifdef _ERROR_CHECKING_ON_
 if(ch==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR can not getChargePathVisitsForSite "); 
		fprintf(stderr,"charge is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif
  
  setMatrixLLValueAtRow(ch->path,2,visit);
  return 0;
}

double getChargePathVisitsForSite(Charge ch, int SiteID){
 
	#ifdef _ERROR_CHECKING_ON_
 if(ch==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR can not getChargePathVisitsForSite "); 
		fprintf(stderr,"charge is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(SiteID<0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR can not getChargePathVisits for  ");
    fprintf(stderr,"SiteID value is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  int seq;
  seq = getMatrixLLLastMatchAtRow(ch->path, SiteID, 1);
  if(seq==-1){
    return -1;
  }else{
    double visits = getMatrixLLNodeElem(ch->path,seq,2);
    return visits;
  }
}

int triggerMatch(Charge ch, double match){
	#ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR can not update charge path "); 
		fprintf(stderr,"charge is NULL.\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif

  return getMatrixLLNumberElemGreaterThanMatchAtRow(ch->path,match,2);
}

int getIDsOfTwoOfMostFrequentlyVisitedSites(Charge ch, int * ID_1, int * ID_2){

	#ifdef _ERROR_CHECKING_ON_
  if(ch==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR in getIDsOfTwoOfMostFrequentlyVisitedSites ");
    fprintf(stderr,"ch is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  int max_visits1;
  int max_visits2;
  int ID1;
  int ID2;
  int seq;
  int flag;

  max_visits1 = 0;
  max_visits2 = 0;
  ID1 = -1;
  ID2 = -1;

  max_visits1 = getMatrixLLNodeElem(ch->path,1,2);
  max_visits1 = getMatrixLLNodeElem(ch->path,1,2);
  ID1         = getMatrixLLNodeElem(ch->path,1,1);
  ID2         = getMatrixLLNodeElem(ch->path,1,1);
 
  if( getMatrixLLNodeElem(ch->path,2,2) >= max_visits1){
    max_visits1 = getMatrixLLNodeElem(ch->path,2,2); 
    ID1         = getMatrixLLNodeElem(ch->path,2,1); 
  }else{
    max_visits2 = getMatrixLLNodeElem(ch->path,2,2); 
    ID2         = getMatrixLLNodeElem(ch->path,2,1); 
  }
  
  flag = 0;

  for(seq=3;seq<=getMatrixLLlength(ch->path);seq++){ 
    if( getMatrixLLNodeElem(ch->path,seq,2)>max_visits1){
      max_visits2 = max_visits1;
      ID2         = ID1;
      max_visits1 = getMatrixLLNodeElem(ch->path,seq,2); 
      ID1         = getMatrixLLNodeElem(ch->path,seq,2);
      flag        = 1; 
    }else if( getMatrixLLNodeElem(ch->path,seq,2)>max_visits2){
      max_visits2 = getMatrixLLNodeElem(ch->path,seq,2); 
      ID2         = getMatrixLLNodeElem(ch->path,seq,2); 
      flag        = 1;
    }
  }
  *ID_1 = ID1;
  *ID_2 = ID2;
  return flag;
}
