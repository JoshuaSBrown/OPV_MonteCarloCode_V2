#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "sitenode.h"
//#include "../../../MEM/mem.h"

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

	if(sn==NULL){
		return NULL;
	}

	//Initialize attributes
	sn->initE=0;
	sn->dwellStatus=-1;
	sn->visitFreq=0;
	sn->visit=0;
	sn->time=0;
	sn->energy=0;
	sn->type=0;
	// By default DataStru set to Point
	sn->Point = (void *) newPoint();
	//sn->NeighList=NULL;
	sn->Cluster = NULL;
	return sn;
}

int printSN(SiteNode sn){
	if(sn==NULL){
		return -1;
	}
	printf("Init Energy     %d\n",sn->initE);
	printf("Dwell Status    %d\n",sn->dwellStatus);
	printf("Type            %d\n",sn->type);
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
	if(sn==NULL){
		return NULL;
	}
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
	
	if(pt==NULL){
		return NULL;
	}

	pt->sum=0;
	int i;
	for(i=0;i<6;i++){
		pt->p[i]=0;
	}
	return pt;	
}

SNarray newSNarray( int len, int wid, int hei) {

	if(len<0 || wid<0 || hei<0){
		return NULL;
	}

	printf("Before Creation\n");
	SNarray snA = (SNarray) malloc(sizeof(struct _SNarray)+sizeof(SiteNode)*len*wid*hei);

	if(snA==NULL){
		return NULL;
	}

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
	if(pt==NULL){
		return -1;
	}

	free(pt);
	return 0;
}

int deleteSN(SiteNode sn){
//Does not Delete the Cluster	
	if(sn==NULL){
		return -1;
	}

	if(sn->type==0){
		deletePoint((Point) sn->Point);
	}else if(sn->type==1){
		deletePoint((Point) sn->Point);
	}

	free(sn);
	return 0;
}

int deleteSNarray(SNarray snA){

	if(snA==NULL){
		return -1;
		}

	int i;
	for(i=0;i<snA->total;i++){
		deleteSN(snA->N[i]);
	}

	free(snA);
	return 0;
}

int printPoint(void * vpt){
	
	if(vpt==NULL){
		return -1;
	}

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

	if(snA==NULL){
		return -1;
	}

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
					fprintf(FreqOut,"C \t %f \t %f \t %f \t %f \t %f \t %f\n",\
							id,jd,kd,getVisFreq(sn),getVis(sn),getEnergy(sn));
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
int setDataStruc(SiteNode sn, int ty, void * ptr){
	
	if(sn==NULL || ty<0 || ptr==NULL || ty>1){
		printf("ERROR Unable to set data structure\n");
		return -1;
	}

	//Structure can either have type 0 or 1 
	//if 0 only points to a point if type 1
	//points to both a point and is part of 
	//a cluster
	//Set to point structure
	if(ty==0){
		if(sn->type!=1){
			sn->type=0;
		}
		//There is no way to ensure that the variable ptr
		//is of type Point
		//Must use caution
		sn->Point=ptr;
		//sn->NeighList=NULL;
		//Set to Cluster
	}else if(ty==1){
		sn->type=1;
		sn->Cluster=ptr;
	}
	return 0;
}

void * getClusterList(const_SiteNode sn){
	if(sn==NULL){
		return NULL;
	}
	return sn->Cluster;
}

int setDwelStat(SiteNode sn, int stat) {
	
	if(stat<-1 || sn==NULL){
		return -1;
	}
	sn->dwellStatus=stat;
	return 0;
}

int getDwelStat(const_SiteNode sn) {
	return sn->dwellStatus;
}

int setVisFreq(SiteNode sn, int freq) {
	if(sn==NULL || freq<0){
		return -1;
	}
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
	if(sn==NULL || vis>2 || vis<0){
		return -1;
	}
	sn->visit=vis;
	return 0;
}

double getVis(const_SiteNode sn) {
	return sn->visit;
}

int setInitE(SiteNode sn, int E) {
	if(sn==NULL || E<0 || E>1){
		return -1;
	}
	sn->initE=E;
	return 0;
}

int getInitE(const_SiteNode sn) {
	if(sn==NULL){
		return -1;
	}
	return sn->initE;
}

int setEnergy(SiteNode sn, double Energy) {
	if(sn==NULL || isnan(Energy)){
		return -1;
	}
	sn->energy=Energy;
	return 0;
}

double getEnergy( const_SiteNode sn) {
	return sn->energy;
}

int setTime(SiteNode sn, double time){
	if(sn==NULL || time<0){
		return -1;
	}
	sn->time=time;
	return 0;
}

int addTime(SiteNode * sn, double time){
	if(sn==NULL){
		printf("WARNING sn does not exist\n");
		return -1;
	}
	if(*sn==NULL){
		printf("WARNING sn does not exist\n");
		return -1;
	}
	if(time==0){
		printf("time added to sitenode is 0!\n");
		exit(1);
	}
	if(time<0){
		printf("WARNING trying to set time to negative value\n");
		return -1;
	}
	(*sn)->time=(*sn)->time+time;
	return 0;
}

double getTime(SiteNode sn){
	if(sn==NULL){
		return -1;
	}
	return sn->time;
}

int setsum(SiteNode sn, double s) {

	if (sn==NULL || s<0.0 || sn->type!=0){
		return -1;
	}
	Point pt = (Point) sn->Point;
	if(pt==NULL){
		return -1;
	}
	pt->sum=s;

	return 0;
}

double getsum(SiteNode sn){
	if(sn==NULL){
		return -1;
	}
	Point pt = (Point) sn->Point;
	if(pt==NULL){
		return -1;
	}
	return pt->sum;
}

int setSN_p(SiteNode sn, int i, double val){
	if(sn==NULL || i<0 || i>5 || val>1.01){
		return -1;
	}

	Point pt = (Point) sn->Point;
	if(pt==NULL){
		return -1;
	}
	pt ->p[i]=val;
	return 0;
}

double getSN_p(SiteNode sn, int i){
	if(sn==NULL || i<0 || i>5){
		return -1;
	}

	Point pt = (Point) sn->Point;
	if(pt==NULL){
		return -1;
	}
	return pt ->p[i];

}

double getsumPt( const_Point pt){
	return pt->sum;
}

int getType(const_SiteNode sn){
	if(sn==NULL){
		return -1;
	}
	return sn->type;
}


int getLoc(int *i, int *j, int *k, int Index, const_SNarray snA) {

	if(snA==NULL || Index<0 || Index>=(snA->height*snA->width*snA->length)){
		return -1;
	}

	*i=(Index)/(snA->width*snA->height);
	Index = (Index) % (snA->width*snA->height);
	*j=Index/(snA->height);
	*k=Index % (snA->height);

	return 0;
}

int getIndex( const_SNarray snA, int i, int j, int k) {

	if (snA!=NULL && i<snA->length && j<snA->width && k<snA->height && i>-1 && j>-1 && k>-1){
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}
	return -1;
}

SiteNode getSN(const_SNarray snA, int i, int j, int k) {
	
	if(snA==NULL){
		return NULL;
	}
	
	int m=getIndex(snA,i,j,k);
	if (m!=-1){
		return snA->N[m];
	}
	return NULL;
}

SiteNode getSNwithInd(const_SNarray snA, int Ind){

	if(snA==NULL){
		return NULL;
	}

	if(Ind!=-1 && Ind<getAtotal(snA)){
		return snA->N[Ind];
	}
	return NULL;

}

int getIndexFrontP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	i++;
	return getIndexPeriodic( snA, i, j, k);
}

int getIndexFront( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	i++;
	if (i>=snA->length || rv==-1){
		return -1;
	}
	return getIndex( snA, i, j, k);
}

int getIndexBehindP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	i--;
	return getIndexPeriodic( snA, i, j, k);
}

int getIndexBehind( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	i--;
	if (i<0 || rv==-1){
		return -1;
	}
	return getIndex( snA, i, j, k);
}

int getIndexLeftP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	j--;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexLeft( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	j--;
	if(j<0 || rv==-1){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndexRightP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	j++;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexRight( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	j++;
	if(j>=snA->width || rv==-1){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndexBelowP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	k--;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexBelow( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	k--;
	if(k<0 || rv ==-1 ){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndexAboveP( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	k++;
	return getIndexPeriodic(snA,i,j,k);
}

int getIndexAbove( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	k++;
	if(k>=snA->height || rv==-1){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndexLeftPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	if (rv==-1){
		return -1;
	}
	j--;
	return getIndexPeriodicY(snA,i,j,k);
}

int getIndexRightPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	if ( rv==-1){
		return -1;
	}
	j++;
	return getIndexPeriodicY(snA,i,j,k);
}

int getIndexAbovePeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	if ( rv==-1){
		return -1;
	}
	k++;
	return getIndexPeriodicZ(snA,i,j,k);
}

int getIndexBelowPeriodic( const_SNarray snA, int Index){
	int i, j, k;
	int rv;
	rv = getLoc(&i, &j, &k, Index, snA);
	if ( rv==-1){
		return -1;
	}
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
	
	return -1;
}
	
int getIndexPeriodicY( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL && \
			i<snA->length && k<snA->height && \
			i>-1 && k>-1) {
		j=(j+snA->width)%snA->width;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}
	
	return -1;
}

int getIndexPeriodicZ( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL && \
			i<snA->length && j<snA->width && \
			i>-1 && j> -1){
		
		k=(k+snA->height)%snA->height;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}

	return -1;
}

int getIndexPeriodic( const_SNarray snA, int i, int j, int k){
	
	if (snA!=NULL){
		i=(i+snA->length)%snA->length;
		j=(j+snA->width)%snA->width;
		k=(k+snA->height)%snA->height;
		return ((snA->width*snA->height)*i+(snA->height)*j+k);
	}

	return -1;
}

int getIndFro( const_SNarray snA, int i, int j, int k){
	i++;
	if( i<=0){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndBeh( const_SNarray snA, int i, int j, int k){
	i--;
	if( i>=(snA->length-1)){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndLef( const_SNarray snA, int i, int j, int k){
	j--;
	if( j>=(snA->width-1)){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndRig( const_SNarray snA, int i, int j, int k){
	j++;
	if( j<=0){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndBel( const_SNarray snA, int i, int j, int k){
	k--;
	if(k>=(snA->height-1) ){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndAbo( const_SNarray snA, int i, int j, int k){
	k++;
	if(k<=0){
		return -1;
	}
	return getIndex(snA,i,j,k);
}

int getIndFroP( const_SNarray snA, int i, int j, int k){
	i++;
	if(i<=0){
		return -1;
	}
	return getIndexPeriodicX( snA, i,j,k);
}

int getIndBehP( const_SNarray snA, int i, int j, int k){
	i--;
	if(i>=(snA->length-1)){
		return -1;
	}
	return getIndexPeriodicX( snA, i,j,k);
}

int getIndLefP( const_SNarray snA, int i, int j, int k){
	j--;
	if(j>=(snA->width-1)){
		return -1;
	}
	return getIndexPeriodicY( snA, i,j,k);
}

int getIndRigP(const_SNarray snA, int i, int j, int k){
	j++;
	if(j<=0){
		return -1;
	}

	return getIndexPeriodicY( snA, i,j,k);
}

int getIndBelP( const_SNarray snA, int i, int j, int k){
	k--;
	if(k>=(snA->height-1)){
		return -1;
	}
	return getIndexPeriodicZ( snA, i, j, k);
}

int getIndAboP( const_SNarray snA, int i, int j, int k){
	k++;
	if(k<=0){
		return -1;
	}
	return getIndexPeriodicZ( snA, i, j, k);
}

int getAlen(const_SNarray snA) {

	if(snA!=NULL){
		return snA->length;
	}
	return -1;
}

int getAwid(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->width;
	}
	return -1;
}

int getAhei(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->height;
	}
	return -1;
}

int getAtotal(const_SNarray snA) {
	
	if(snA!=NULL){
		return snA->total;
	}
	return -1;
}

int getUnOccYZplane(const_SNarray snA,const int l){
	
	if(snA==NULL || l<0 || l>=snA->length){
		return -1;
	}
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

	if(snA==NULL || l<0 || l>=snA->width){
		return -1;
	}

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

	if(snA==NULL || l<0 || l>=snA->height){
		return -1;
	}

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
	
	if(snA==NULL || i<0 || i>=(snA->length-1) ||\
		 j<0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i+1,j,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccXneg(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<=0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i-1,j,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYpos(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width-1) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i,j+1,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYneg(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<=0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i,j-1,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYposPeriodic(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}
	
	SiteNode sn1;
	sn1 = getSN(snA,i,(j+1)%snA->width,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccYnegPeriodic(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i,(j-1+snA->width)%snA->width,k);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccZpos(const_SNarray snA,int i, int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width) ||\
		 k<0 || k>=(snA->height-1)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i,j,k+1);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccZneg(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || i<0 || i>=(snA->length) ||\
		 j<0 || j>=(snA->width) ||\
		 k<=0 || k>=(snA->height)){
		return -1;
	}

	SiteNode sn1;
	sn1 = getSN(snA,i,j,k-1);

	int stat1 = sn1->dwellStatus;

	if(stat1!=-1) {
		return 1;
	}

	return 0;
}

int OccAllNei(const_SNarray snA,int i,int j,int k) {
	
	if(snA==NULL || (i-1)<0 || i+1>=snA->length ||\
									(j-1)<0 || j+1>=snA->width  ||\
									(k-1)<0 || k+1>=snA->height ){
		return -1;
	}
	
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
	
	if(snA==NULL || (i-1)<0 || (i+1)>=snA->length || 
									(k-1)<0 || (k+1)>=snA->height || 
									j<0 || j>=(snA->width)) {
		return -1;
	}
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

	if(snA==NULL){
		return -1;
	}
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

	if(snA==NULL){
		return -1;
	}

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

	if(snA==NULL){
		return -1;
	}

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
