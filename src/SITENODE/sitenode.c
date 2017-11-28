#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "sitenode.h"
#include "../ERROR/error.h"

struct _SiteNode{

	//flag to judge whether assign the site energy. If initE = 0,
	//site energy is not initialed. And vice versa for initE = 1;
	int initE; 

	// if 0 <= dwellStatus < N, it means charge[dwellStatus] on 
	//the node, -1 means no electron on this node
	int dwellStatus; 

	// 0 - Normal Node with sum and p value
	// 1 - Cluster 
	int type;
	double visitFreq;
	double visit;

  // Determines if the site has decayed or not
  // 0 - has not decayed
  // 1 - has decayed
  int decay_status;
	//site energy
	double energy;

	//time occupied
	//How much time during the run was the site occupied
	double time;

	// Pointer to the point structure
	void * Point;

	// Cluster link list
	void * Cluster;
};

struct _Point{
	//The total hopping rate
	double sum;
	//the length(<1) of jumping possiblities to the six neighbor sites;
	//0:x-1, 1:x+1, 2:y-1, 3:y+1, 4:z-1, 5:z+1
	double p[6];
};

struct _SNarray{
	int length;
	int width;
	int height;
	int total;

	SiteNode N[0];
};

SiteNode newSN(void) {

	SiteNode sn = (SiteNode) malloc(sizeof(struct _SiteNode));

  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newSN\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	//Initialize attributes
	sn->initE=0;
	sn->dwellStatus=-1;
	sn->visitFreq=0;
	sn->visit=0;
  sn->decay_status = 0;
	sn->time=0;
	sn->energy=0;
	sn->type=0;
	// By default DataStru set to Point
	sn->Point = (void *) newPoint();
	//sn->NeighList=NULL;
	sn->Cluster = NULL;
	return sn;
}

int clearPoint(SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
  if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in clearPoint\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  deletePoint(sn->Point);
  sn->Point = NULL;
  return 0;
}

int printSN(SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in printSN\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	printf("Init Energy     %d\n",sn->initE);
	printf("Dwell Status    %d\n",sn->dwellStatus);
	printf("Type            %d\n",sn->type);
  printf("Decay Status    %d\n",sn->decay_status);
	printf("Visit Frequency %g\n",sn->visitFreq);
	printf("Visit 					%g\n",sn->visit);
	printf("Energy					%g\n",sn->energy);

	if(sn->Point!=NULL){
		printPoint(sn->Point);
	}
	printf("\n");
	return 0;
}

void * getPoint(const_SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getPoint\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return sn->Point;
}

/*
void * getCluster(const_SiteNode sn){
	if(sn==NULL){
		return NULL;
	}
	return sn->Cluster;
}
*/
Point newPoint(void){
	
	Point pt = (Point) malloc(sizeof(struct _Point));
	
  #ifdef _ERROR_CHECKING_ON_
	if(pt==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newPoint\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	pt->sum=0;
	int i;
	for(i=0;i<6;i++){
		pt->p[i]=0;
	}
	return pt;	
}

SNarray newSNarray( int len, int wid, int hei) {

  #ifdef _ERROR_CHECKING_ON_
	if(len<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR len is less than 0 in newSNarray must be positive\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(wid<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR wid is less than 0 in newSNarray must be positive\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
	if(hei<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR hei is less than 0 in newSNarray must be positive\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	printf("Before Creation\n");
	SNarray snA = (SNarray) malloc(sizeof(struct _SNarray)+sizeof(SiteNode)*len*wid*hei);

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newSNarray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	snA->length=len;
	snA->width=wid;
	snA->height=hei;
	snA->total=len*wid*hei;

	printf("After Creation\n");
	int i;
	for(i=0;i<snA->total;i++) {
		snA->N[i]=newSN();	

		//if(i>(snA->total-140000)){
		//printf("Value of i %d total %d\n",i,snA->total);
		//mem_check();
		//}
	}

	printf("After Creation of sites\n");
	return snA;
}

int deletePoint( Point pt){
    #ifdef _ERROR_CHECKING_ON_
	if(pt==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in deletePoint\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
    #endif

	free(pt);
	return 0;
}

int deleteSN(SiteNode sn){
//Does not Delete the Cluster	
    #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in deleteSN\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	if(sn->type==0){
		deletePoint((Point) sn->Point);
	}else if(sn->type==1){
		deletePoint((Point) sn->Point);
	}

	free(sn);
	return 0;
}

int deleteSNarray(SNarray * snA){

	#ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR snA is NULL in deleteSNarray\n");
		#endif
		#ifdef _FORCE_HARD_CRASH_
		exit(1);
		#endif
		return -1;
	}
	if(*snA==NULL){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR *snA is NULL in deleteSNarray\n");
		#endif
		#ifdef _FORCE_HARD_CRASH_
		exit(1);
		#endif
		return -1;
	}
	#endif

	int i;
	for(i=0;i<(*snA)->total;i++){
		deleteSN((*snA)->N[i]);
	}

	free(*snA);
	*snA=NULL;
	return 0;
}

int printPoint(void * vpt){
	
  #ifdef _ERROR_CHECKING_ON_
	if(vpt==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR vpt is NULL in printPoint\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	Point pt = (Point) vpt;
	printf("sum: %g\n",(pt)->sum);
	int i;
	double prev=0;
	for(i=0;i<6;i++){
		printf("p[%d]: %g\t",i,(pt)->p[i]);
	}
	printf("\nDifference in Pval\n");
	for(i=0;i<6;i++){
		printf("%g\t",(pt)->p[i]-prev);
		prev = (pt)->p[i];
	}

	return 0;

}

int printVisitFreq(const_SNarray snA, char * FileName) {

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in printVisitFreq\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	char buf[256];
	snprintf(buf, sizeof buf,"%s%s",FileName,".xyz");
	
	int i, j, k;
	SiteNode sn;
	FILE * FreqOut;
	if((FreqOut = fopen(buf,"w")) == NULL) {
		printf("Error! unable to open VisitFreq.xyz\n");

	}else{
		fprintf(FreqOut,"%d\n\n",(snA->total));
		double id=0;
		double jd=0;
		double kd=0;

		for (i=0;i<snA->length;i++){
			jd=0;
			for(j=0;j<snA->width;j++){
				kd=0;
				for(k=0;k<snA->height;k++){
					sn=getSN(snA,i,j,k);
					fprintf(FreqOut,"C \t %f \t %f \t %f \t %f \t %f \t %g\n",\
							id,jd,kd,getVisFreq(sn),getVis(sn),getTime(sn));
					kd++;
				}
				jd++;
			}
			id++;
		}
		fclose(FreqOut);
	}

	return 0;
}

/*int setNeighList(SiteNode sn, void * ptr){
	
	if(sn==NULL || ptr==NULL || sn->type!=1){
		return -1;
	}
	
	sn->NeighList=ptr;
	return 0;
}
*/
int setDataStruc(SiteNode * sn, int ty, void * ptr){
	
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setDataStruc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(ty<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ty is less than 0 in setDataStruc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(ptr==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ptr is NULL in setDataStruc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
	if(ty>1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ty is greater than 1 in setDataStruc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	//Structure can either have type 0 or 1 
	//if 0 only points to a point if type 1
	//points to both a point and is part of 
	//a cluster
	//Set to point structure
	if(ty==0){
		if((*sn)->type!=1){
			(*sn)->type=0;
		}
		//There is no way to ensure that the variable ptr
		//is of type Point
		//Must use caution
    #ifdef _ERROR_CHECKING_ON_
    if((*sn)->Point!=NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR sn->Point is not NULL in setDataStruc can not set a new point until the old one has been removed from the sitenode\n");
      #endif
      #ifdef _FORCE_HARD_CRASH_
      exit(1);
      #endif
      return -1;
    }
    #endif
		(*sn)->Point=ptr;
		//sn->NeighList=NULL;
		//Set to Cluster
	}else if(ty==1){
		(*sn)->type=1;
		(*sn)->Cluster=ptr;
	}
	return 0;
}

void * getClusterList(const_SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getClusterList\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return sn->Cluster;
}

int Decay(SiteNode sn,double Energy){
  #ifdef _ERROR_CHECKING_ON_
	if( sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setDwelStat\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

  if(sn->decay_status==0){
    sn->energy = sn->energy+Energy;  
    sn->decay_status=1;
  }
  return 0;
}

int setDwelStat(SiteNode sn, int stat) {
	
  #ifdef _ERROR_CHECKING_ON_
	if( sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setDwelStat\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(stat<-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR stat is less than -1 in setDwelStat\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->dwellStatus=stat;
	return 0;
}

int getDwelStat(const_SiteNode sn) {
	return sn->dwellStatus;
}

int setVisFreq(SiteNode sn, int freq) {
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setVisFreq\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(freq<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR freq is less than 0 in setVisFreq\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->visitFreq=freq;
	return 0;
}

void incVisFreq(SiteNode sn) {
	sn->visitFreq++;
}

double getVisFreq(const_SiteNode sn) {
	return sn->visitFreq;
}

int setVis(SiteNode sn, double vis) {
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setVis\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(vis>2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR vis is greater than 2 in setVis\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(vis<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR vis is less than 0 in setVis\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->visit=vis;
	return 0;
}

double getVis(const_SiteNode sn) {
  #ifdef _ERROR_CHECKING_ON_
  if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getVis\n");
    #endif
    exit(1);
  }
  #endif
	return sn->visit;
}

//initE is used to keep track of whether the energy of the 
//site has been initialized yet or not. 
int setInitE(SiteNode sn, int E) {
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setInitE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(E<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR E is less than 0 in setInitE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(E>2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR E is greater than 2 in setInitE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->initE=E;
	return 0;
}

int getInitE(const_SiteNode sn) {
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getInitE\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return sn->initE;
}

int setEnergy(SiteNode sn, double Energy) {
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setEnergy\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(isnan(Energy)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Energy is nan in setEnergy\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->energy=Energy;
	return 0;
}

double getEnergy( const_SiteNode sn) {
	return sn->energy;
}

int setDecayStatus(SiteNode sn, int decay_status){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setDecayStatus\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(!(decay_status==1 || decay_status==0)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR decay_status set to a value that is not 0 or 1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->decay_status=decay_status;
	return 0;
}

int getDecayStatus( const_SiteNode sn) {
	return sn->decay_status;
}

int setTime(SiteNode sn, double time){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(time<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR time is negative in setTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	sn->time=time;
	return 0;
}

int addTime(SiteNode * sn, double time){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in addTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(*sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *sn is NULL in addTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(time==0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR adding time of 0 in addTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
	  return 0;
  }
	if(time<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR added a negative time in addTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	(*sn)->time=(*sn)->time+time;
	return 0;
}

double getTime(SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getTime\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return sn->time;
}

int setsum(SiteNode sn, double s) {

  #ifdef _ERROR_CHECKING_ON_
	if (sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setsum\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if (s<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR s is less than 0 in setsum, s must be positive or 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if (sn->type!=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn->type does not equal 0 in setsum\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(((Point)sn->Point)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in setsum\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Point pt = (Point) sn->Point;
	pt->sum=s;

	return 0;
}

double getsum(SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getsum\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(((Point)sn->Point)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in getsum\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Point pt = (Point) sn->Point;
	return pt->sum;
}

int setSN_p(SiteNode sn, int i, double val){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in setSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in setSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>5){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater than 5 in setSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(val>1.01){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR val is greater than 1.01 in setSN_p supposed to be a probability must have a value between 0 and 1 inclusive\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(((Point)sn->Point)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in setSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Point pt = (Point) sn->Point;
	pt ->p[i]=val;
	return 0;
}

double getSN_p(SiteNode sn, int i){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in getSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>5){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater than 5 in getSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(((Point)sn->Point)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in getSN_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Point pt = (Point) sn->Point;
	return pt ->p[i];

}

double getsumPt( const_Point pt){
  #ifdef _ERROR_CHECKING_ON_
  if(pt==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pt is NULL in getsumPt\n");
    #endif
		exit(1);
  }
  #endif
	return pt->sum;
}

int getType(const_SiteNode sn){
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in getType\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return sn->type;
}


int getLoc(int *i, int *j, int *k, int Index, const_SNarray snA) {

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getLoc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(Index>=(snA->height*snA->width*snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Index is greater or equal to the total number of nodes in snA in getLoc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( Index<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Index is less than 0 in getLoc\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	*i=(Index)/(snA->width*snA->height);
	Index = (Index) % (snA->width*snA->height);
	*j=Index/(snA->height);
	*k=Index % (snA->height);

	return 0;
}

int getIndex( const_SNarray snA, int i, int j, int k) {

  #ifdef _ERROR_CHECKING_ON_
  if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(i>=snA->length){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(j>=snA->width){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(k>=snA->height){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in getIndex\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
		
  return ((snA->width*snA->height)*i+(snA->height)*j+k);
}

SiteNode getSN(const_SNarray snA, int i, int j, int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getSN\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	
	int m=getIndex(snA,i,j,k);
	if (m!=-1){
		return snA->N[m];
	}
	return NULL;
}

SiteNode getSNwithInd(const_SNarray snA, int Ind){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA returned NULL in getSNwithInd\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  if(Ind>=getAtotal(snA)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Ind is greater or equal to getAtotal(snA) in getSNwithInd\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
  }
  if(Ind==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Ind is -1 in getSNwithInd\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
  }
  #endif

	return snA->N[Ind];

}

int getIndexFrontP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	i++;
	return getIndexPeriodic( snA, i, j, k);
}

int getIndexFront( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	i++;
  #ifdef _ERROR_CHECKING_ON_
	if (rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexFront\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if (i>=snA->length){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in getIndexFront\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex( snA, i, j, k);
}

int checkBoundsIndexFront( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	i++;
  #ifdef _ERROR_CHECKING_ON_
	if (rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexFront\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if (i>=snA->length){
    //Out of bounds 
		return 0;
	}else{
    //Within bounds
    return 1;
  }
}
int getIndexBehindP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	i--;
	return getIndexPeriodic( snA, i, j, k);
}

int getIndexBehind( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	i--;
  #ifdef _ERROR_CHECKING_ON_
	if (rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexBehind\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if (i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in getIndexBehind\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex( snA, i, j, k);
}

int checkBoundsIndexBehind( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	i--;
  #ifdef _ERROR_CHECKING_ON_
	if (rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexBehind\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if (i<0){
    //Means it is out of bounds
		return 0;
	}else{
    //It is within bounds
    return 1;
  }
}

int getIndexLeftP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	j--;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexLeft( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	j--;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexLeft\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in getIndexLeft\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int checkBoundsIndexLeft( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	j--;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexLeft\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if(j<0){
    //Out of Bounds
    return 0;
	}else{
    //Within Bounds
    return 1;
  }
	return getIndex(snA,i,j,k);
}

int getIndexRightP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	j++;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexRight( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	j++;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexRight\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=snA->width){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in getIndexRight\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int checkBoundsIndexRight( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	j++;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexRight\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if(j>=snA->width){
  	//Oub of Bounds
    return 0;
	}else{
    //Within Bounds
    return 1;
  }
}

int getIndexBelowP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	k--;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexBelow( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	k--;
  #ifdef _ERROR_CHECKING_ON_
	if(rv ==-1 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexBelow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in getIndexBelow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int checkBoundsIndexBelow( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	k--;
  #ifdef _ERROR_CHECKING_ON_
	if(rv ==-1 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexBelow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if(k<0 ){
    //Out of Bounds
    return 0;
	}else{
    //Within Bounds
    return 1;
  }
}


int getIndexAboveP( const_SNarray snA, int Index){
	int i, j, k;
	getLoc(&i, &j, &k, Index, snA);
	k++;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexAbove( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	k++;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexAbove\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=snA->height){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int checkBoundsIndexAbove( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
	k++;
  #ifdef _ERROR_CHECKING_ON_
	if(rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in checkBoundsIndexAbove\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	if(k>=snA->height){
    //Out of Bounds
		return 0;
	}else{
    //Within Bounds
    return 1;
  }
}


int getIndexLeftPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv = getLoc(&i, &j, &k, Index, snA);
  #ifdef _ERROR_CHECKING_ON_
	if (rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexLeftPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	j--;
	return getIndexPeriodicY(snA,i,j,k);
}

int getIndexRightPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
  #ifdef _ERROR_CHECKING_ON_
	if ( rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexRightPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	j++;
	return getIndexPeriodicY(snA,i,j,k);
}

int getIndexAbovePeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
  #ifdef _ERROR_CHECKING_ON_
	if ( rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexAbovePeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	k++;
	return getIndexPeriodicZ(snA,i,j,k);
}

int getIndexBelowPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
  #ifdef _ERROR_CHECKING_ON_
	if ( rv==-1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getLoc returned -1 in getIndexBelowPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	k--;
	return getIndexPeriodicZ(snA,i,j,k);
}

int getIndexPeriodicX( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL && \
			j<snA->width && k<snA->height && \
			j>-1 && k>-1) {
		i=(i+snA->length)%snA->length;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}
	
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR either snA is NULL or j, k are out of bounds in getIndexPeriodicX\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}
	
int getIndexPeriodicY( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL && \
			i<snA->length && k<snA->height && \
			i>-1 && k>-1) {
		j=(j+snA->width)%snA->width;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}
	
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR either snA is NULL or i, k are out of bounds in getIndexPeriodicY\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getIndexPeriodicZ( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL && \
			i<snA->length && j<snA->width && \
			i>-1 && j> -1){
		
		k=(k+snA->height)%snA->height;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}

  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR either snA is NULL or i,j are out of bounds in getIndexPeriodicZ\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getIndexPeriodic( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL){
		i=(i+snA->length)%snA->length;
		j=(j+snA->width)%snA->width;
		k=(k+snA->height)%snA->height;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}

  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR snA is NULL in getIndexPeriodic\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getIndFro( const_SNarray snA, int i, int j, int k){
	i++;
  #ifdef _ERROR_CHECKING_ON_
	if( i<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than or equal to 0 in getIndFro\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndBeh( const_SNarray snA, int i, int j, int k){
	i--;
  #ifdef _ERROR_CHECKING_ON_
	if( i>=(snA->length-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than or equal to snA->length-1 in getIndBeh\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndLef( const_SNarray snA, int i, int j, int k){
	j--;
  #ifdef _ERROR_CHECKING_ON_
	if( j>=(snA->width-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width-1 in getIndLef\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndRig( const_SNarray snA, int i, int j, int k){
	j++;
  #ifdef _ERROR_CHECKING_ON_
	if( j<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than or equal to 0 in getIndRig\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndBel( const_SNarray snA, int i, int j, int k){
	k--;
  #ifdef _ERROR_CHECKING_ON_
	if(k>=(snA->height-1) ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height-1 in getIndBel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndAbo( const_SNarray snA, int i, int j, int k){
	k++;
  #ifdef _ERROR_CHECKING_ON_
	if(k<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than or equal to 0 in getIndAbo\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndex(snA,i,j,k);
}

int getIndFroP( const_SNarray snA, int i, int j, int k){
	i++;
  #ifdef _ERROR_CHECKING_ON_
	if(i<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than or equal to 0 in getIndFroP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndexPeriodicX( snA, i,j,k);
}

int getIndBehP( const_SNarray snA, int i, int j, int k){
	i--;
  #ifdef _ERROR_CHECKING_ON_
	if(i>=(snA->length-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length-1 in getIndBehP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndexPeriodicX( snA, i,j,k);
}

int getIndLefP( const_SNarray snA, int i, int j, int k){
	j--;
  #ifdef _ERROR_CHECKING_ON_
	if(j>=(snA->width-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width-1 in getIndLefP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndexPeriodicY( snA, i,j,k);
}

int getIndRigP(const_SNarray snA, int i, int j, int k){
	j++;
  #ifdef _ERROR_CHECKING_ON_
	if(j<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less or equal to 0 in getIndRigP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	return getIndexPeriodicY( snA, i,j,k);
}

int getIndBelP( const_SNarray snA, int i, int j, int k){
	k--;
  #ifdef _ERROR_CHECKING_ON_
	if(k>=(snA->height-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height-1 in getIndBelP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndexPeriodicZ( snA, i, j, k);
}

int getIndAboP( const_SNarray snA, int i, int j, int k){
	k++;
  #ifdef _ERROR_CHECKING_ON_
	if(k<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less or equal to 0 in getIndAboP\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return getIndexPeriodicZ( snA, i, j, k);
}

int getAlen(const_SNarray snA) {

	if(snA!=NULL){
		return snA->length;
	}
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR snA is NULL in getAlen\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getAwid(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->width;
	}
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR snA is NULL in getAwid\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getAhei(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->height;
	}
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR snA is NULL in getAhei\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int getAtotal(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->total;
	}
  #ifdef _ERROR_CHECKING_ON_
  #ifdef _ERROR_
  fprintf(stderr,"ERROR snA is NULL in getAtotal\n");
  #endif
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  #endif
	return -1;
}

int checkSNconnectedCluster(const_SiteNode sn){
    
  #ifdef _ERROR_CHECKING_ON_
	if(sn==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR sn is NULL in checkSNconnectedCluster\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
  if(sn->Cluster==NULL){
    //Not connected
    return 0;
  }else{
    //Connected
    return 1;
  }
}

int getUnOccYZplane(const_SNarray snA,const int l){
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getUnOccYZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(l<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is less than 0 in getUnOccYZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(l>=snA->length){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is greater or equal to snA->length in getUnOccYZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	int j;
	int k;
	int count;
	int test;
	count=0;
	SiteNode sn;

	for (j=0;j<snA->width;j++){
		for (k=0;k<snA->height;k++){
	
			test=getIndex(snA,l,j,k);
			if(test!=-1){
				sn = snA->N[getIndex(snA,l,j,k)];
				getSN(snA, l, j, k);		
				if (sn->dwellStatus==-1){
					count++;
				}
			}else{
				return -1;
			}
		}
	}
	return count;
}

int getUnOccXZplane(const_SNarray snA,const int l){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getUnOccXZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(l<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is less than 0 in getUnOccXZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(l>=snA->width){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is greater or equal to snA->width in getUnOccXZplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	int i;
	int k;
	int count;
	int test;
	count=0;
	SiteNode sn;

	for (i=0;i<snA->length;i++){
		for (k=0;k<snA->height;k++){
	
			test=getIndex(snA,i,l,k);
			if(test!=-1){
				sn = snA->N[getIndex(snA,i,l,k)];
				getSN(snA, i, l, k);		
				if (sn->dwellStatus==-1){
					count++;
				}
			}else{
				return -1;
			}
		}
	}
	return count;
}

int getUnOccXYplane(const_SNarray snA,const int l){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in getUnOccXYplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( l<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is less than 0 in getUnOccXYplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(l>=snA->height){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR l is greater or equal to snA->height in getUnOccXYplane\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	int i;
	int j;
	int count;
	int test;
	count=0;
	SiteNode sn;

	for (i=0;i<snA->length;i++){
		for (j=0;j<snA->width;j++){
	
			test=getIndex(snA,i,j,l);
			if(test!=-1){
				sn = snA->N[getIndex(snA,i,j,l)];
				getSN(snA, i, j, l);		
				if (sn->dwellStatus==-1){
					count++;
				}
			}else{
				return -1;
			}
		}
	}
	return count;
}


int OccXpos(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length-1 in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccXpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	SiteNode sn1;
	sn1 = getSN(snA,i+1,j,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccXneg(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<=0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than or greater than 0 in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccXneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	SiteNode sn1;
	sn1 = getSN(snA,i-1,j,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYpos(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width-1 in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccYpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	SiteNode sn1;
	sn1 = getSN(snA,i,j+1,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYneg(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( i<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( j<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 1 in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal snA->width in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccYneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	SiteNode sn1;
	sn1 = getSN(snA,i,j-1,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYposPeriodic(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( i<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccYposPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	
	SiteNode sn1;
	sn1 = getSN(snA,i,(j+1)%snA->width,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYnegPeriodic(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( i<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 0 in OccYnegPeriodic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height in OccYnegPeridic\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	SiteNode sn1;
	sn1 = getSN(snA,i,(j-1+snA->width)%snA->width,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccZpos(const_SNarray snA,int i, int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccZpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccZpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length in OccZpos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccZPos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width OccZPos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 1 in OccZPos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height-1)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height-1 in OccZPos\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	SiteNode sn1;
	sn1 = getSN(snA,i,j,k+1);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccZneg(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( i<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 0 in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(i>=(snA->length)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater than snA->length in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 1 in OccZneg\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(k>=(snA->height)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	SiteNode sn1;
	sn1 = getSN(snA,i,j,k-1);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccAllNei(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if((i-1)<0 ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( i+1>=snA->length ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less or equal to snA->length-1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if((j-1)<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j+1>=snA->width){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width-1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if ((k-1)<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if ( k+1>=snA->height ){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater or equal to snA->height-1 in OccAllNei\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	
	//Need to check six surrounding neiboring sites
	SiteNode sn1;
	SiteNode sn2;
	SiteNode sn3;
	SiteNode sn4;
	SiteNode sn5;
	SiteNode sn6;

	sn1 = getSN(snA,i-1,j,k);
	sn2 = getSN(snA,i+1,j,k);
	sn3 = getSN(snA,i,j-1,k);
	sn4 = getSN(snA,i,j+1,k);
	sn5 = getSN(snA,i,j,k-1);
	sn6 = getSN(snA,i,j,k+1);

	int stat1 = sn1->dwellStatus;
	int stat2 = sn2->dwellStatus;
	int stat3 = sn3->dwellStatus;
	int stat4 = sn4->dwellStatus;
	int stat5 = sn5->dwellStatus;
	int stat6 = sn6->dwellStatus;

	if(stat1!=-1 && stat2!=-1 && stat3!=-1 && stat4!=-1 && stat5!=-1 && stat6!=-1){
		return 1;
	}

	return 0;
}

int OccAllNeiPeriodicY(const_SNarray snA,int i,int j,int k) {
	
  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( (i-1)<0) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is less than 1 in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if( (i+1)>=snA->length) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR i is greater or equal to snA->length-1 in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
		if(	(k-1)<0 ) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is less than 1 in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if((k+1)>=snA->height) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR k is greater than or equal to snA->height-1 in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j<0) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is less than 0 in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
	if(j>=(snA->width)) {
    #ifdef _ERROR_
    fprintf(stderr,"ERROR j is greater or equal to snA->width in OccAllNeiPeriodicY\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//Need to check six surrounding neiboring sites
	SiteNode sn1;
	SiteNode sn2;
	SiteNode sn3;
	SiteNode sn4;
	SiteNode sn5;
	SiteNode sn6;

	sn1 = getSN(snA,i-1,j,k);
	sn2 = getSN(snA,i+1,j,k);
	sn3 = getSN(snA,i,(j-1+snA->width)%(snA->width),k);
	sn4 = getSN(snA,i,(j+1)%(snA->width),k);
	sn5 = getSN(snA,i,j,k-1);
	sn6 = getSN(snA,i,j,k+1);

	int stat1 = sn1->dwellStatus;
	int stat2 = sn2->dwellStatus;
	int stat3 = sn3->dwellStatus;
	int stat4 = sn4->dwellStatus;
	int stat5 = sn5->dwellStatus;
	int stat6 = sn6->dwellStatus;

	if(stat1!=-1 && stat2!=-1 && stat3!=-1 && stat4!=-1 && stat5!=-1 && stat6!=-1){
		return 1;
	}
	return 0;
	
}

int setDefaultSNa(SNarray snA){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in setDefaultSNa\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	int i;
	int j;
	int k;
	SiteNode temp;
	for( i=0;i<snA->length;i++){
		for (j=0;j<snA->width;j++){
			for (k=0;k<snA->height;k++){
				temp=getSN(snA,i,j,k);
				temp->initE=0;
				temp->dwellStatus=-1;
				temp->visitFreq=0;
				temp->visit=0;
			}
		}
	}

	return 0;
}

int printSNarray( const_SNarray snA){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in printSNarray\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif

	printf("\nLength %d Width %d Height %d\n",snA->length,snA->width,snA->height);
	
	for(int i=0;i<snA->total;i++) {
		printf("\nSite %d\n",i);
		printf("initE %d\n",snA->N[i]->initE);
		printf("dwellStatus %d\n",snA->N[i]->dwellStatus);
		printf("type %d\n",snA->N[i]->type);
		printf("visitFreq %g\n",snA->N[i]->visitFreq);
		printf("visit %g\n",snA->N[i]->visit);
		printf("energy %g\n",snA->N[i]->energy);

		if(snA->N[i]->type==1){
			printf("Attached to Cluster\n");
		}else{
			printf("Not attached to Cluster\n");
		}
	}
	return 0;
}

int printSNarray_Detailed( const_SNarray snA){

  #ifdef _ERROR_CHECKING_ON_
	if(snA==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR snA is NULL in printSNarray_Detailed\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	int i1, j1, k1;
	printf("\nLength %d Width %d Height %d\n",snA->length,snA->width,snA->height);
	
	for(int i=0;i<snA->total;i++) {
		printf("\n\nSite %d\n",i);
		printf("initE %d\n",snA->N[i]->initE);
		printf("dwellStatus %d\n",snA->N[i]->dwellStatus);
		printf("type %d\n",snA->N[i]->type);
		printf("visitFreq %g\n",snA->N[i]->visitFreq);
		printf("visit %g\n",snA->N[i]->visit);
		printf("energy %g\n",snA->N[i]->energy);
		getLoc(&i1,&j1,&k1,i,snA);
		printf("Position %d,%d,%d\n",i1,j1,k1);
		if(snA->N[i]->type==1){
			printf("Attached to Cluster\n");
			printPoint(snA->N[i]->Point);
		}else{
			printf("Not attached to Cluster\n");
			printPoint(snA->N[i]->Point);
		}
	}
	return 0;
}
