#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "charge.h"
//#include "../MEM/mem.h"

struct _Charge {
	int x, y, z; //current position
	int tot_x; //the distance run by the charge in the x direction
	double t;  //total time the electron uses
	double dwelltime; //the time the electron will dwell on current site
};

struct _ChargeArray {
	int length;
	Charge C[0];
};

Charge newCharge(void) {
	Charge ch = (Charge) malloc( sizeof(struct _Charge ));
	if(ch==NULL)	//can't check
		return NULL;
	//By default set values to 0
	ch->x=0;
	ch->y=0;
	ch->z=0;
	ch->t=0;
	ch->tot_x=0;
	ch->dwelltime=0;
	return ch;
}

ChargeArray newChargeA(int len) {
	if(len < 0)
		return NULL;
	ChargeArray chA=(ChargeArray) malloc( sizeof(struct _ChargeArray)+sizeof(Charge)*(len));
	if(chA==NULL)	//can't check on test script because created within function
		return NULL; 
	chA->length=len;
	for(int i=0;i<len;i++){ 
		chA->C[i]=newCharge();
		if(chA->C[i]==NULL) //can't check
			return NULL; 
	}
	return chA;
}

int deleteCharge(Charge ch) {
	if(ch==NULL)
		return -1; 
	free(ch);
	return 0; 
}

int deleteChargeA(ChargeArray chA) {
	if(chA==NULL)
		return -1; 
	int i;
	for(i=0;i<chA->length;i++){
		free(chA->C[i]);
	}
	free(chA);
	return 0;
}

int getChargeA_len(ChargeArray chA){
	if(chA==NULL){
		return -1;
	}
	return chA->length;
}

Charge getCharge(const_ChargeArray chA, int i){
	if(chA==NULL || i>= chA->length || i < 0)
		return NULL; 
	return chA->C[i];
}

int printCharge(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	printf("\nx: %d y: %d z: %d\n",ch->x,ch->y,ch->z);
	printf("tot_x: %d\n",ch->tot_x);
	printf("t: %g\n",ch->t);
	printf("dwelltime: %g\n",ch->dwelltime);
	return 0; 
}

int printChargeA(const_ChargeArray chA){
	if(chA==NULL)
		return -1; 
		printf("\nCharge Array\n");
	printf("Length: %d\n",chA->length);
	for(int i=0;i<(chA->length);i++){
		printf("\nCharge %d\n",i);
		printCharge(chA->C[i]);
	}
	return 0; 
}

double getDwel(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->dwelltime;
}

int getCx(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->x;
}

int getCy(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->y;
}

int getCz(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->z;
}

int setDwel(Charge ch, double dw) {
	if(ch==NULL || dw < 0)
		return -1; 
	ch->dwelltime=dw;
	return 0;
}

int MinusDwel(Charge ch, double dw){
	if(ch==NULL || dw < 0)
		return -1; 
	ch->dwelltime-=dw;
	return 0; 
}

int sett(Charge ch, long double t) {
	if(ch==NULL || t<= 0)
		return -1; 
	ch->t=t;
	return 0; 
}

int Plust(Charge ch, long double t){
	if(ch==NULL || t<=0)
		return -1; 
	ch->t+=t;
	return 0; 
}

double gett(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->t;
}
int CxPlus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->x++;
	ch->tot_x++;
	return 0; 
}

int CxMinus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->x--;
	ch->tot_x--;
	return 0; 
}

int CyPlus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->y++;
	return 0;
}

int CyMinus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->y--;
	return 0;
}

int CyPlusPeriodic(Charge ch,int w){
	if(ch==NULL || w <=0)
		return -1; 
	ch->y++;
	ch->y=(ch->y%w);
	if(ch->y>=w) //can't check
		return -2;
	return 0; 
}

int CyMinusPeriodic(Charge ch, int w) {
	if(ch==NULL || w <= 0)
		return -1; 
	ch->y--;
	ch->y=(ch->y+w)%w;
	if( ch->y >= w || ch->y<=-1) //can't check
		return -2; 
	return 0;
}

int CzPlusPeriodic(Charge ch,int h){
	if(ch==NULL || h <= 0)
		return -1; 
	ch->z++;
	ch->z=(ch->z%h);
	if(ch->z>=h) //can't check
		return -2;
	return 0; 
}

int CzMinusPeriodic(Charge ch, int h) {
	if(ch==NULL || h <=0)
		return -1; 
	ch->z--;
	ch->z=(ch->z+h)%h;
	if(ch->z>=h || ch->z==-1) //can't check
		return -2; 
	return 0;
}

int CxPlusPeriodic(Charge ch,int l){
	if(ch==NULL || l <= 0)
		return -1; 
	ch->x++;
	ch->tot_x++;
	ch->x=(ch->x%l);
	if(ch->x>=l) //can't check
		return -2;
	return 0; 
}

int CxMinusPeriodic(Charge ch, int l) {
	if(ch==NULL || l <=0)
		return -1; 
	ch->x--;
	ch->tot_x--;
	ch->x=(ch->x+l)%l;
	if(ch->x>=l || ch->x==-1) //can't check
		return -2; 
	return 0;
}

int CzPlus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->z++;
	return 0;
}

int CzMinus(Charge ch){
	if(ch==NULL)
		return -1; 
	ch->z--;
	return 0;
}

int setCx(Charge ch, int cx) {
	if(ch==NULL )
		return -1; 
	ch->x=cx;
	ch->tot_x = cx;
	return 0;
}
int setCy(Charge ch, int cy) {
	if(ch==NULL )
		return -1; 	
	ch->y=cy;
	return 0;
}
int setCz(Charge ch, int cz) {
	if(ch==NULL )
		return -1; 
	ch->z=cz;
	return 0;
}
int getXdist(const_Charge ch) {
	if(ch==NULL)
		return -1; 
	return ch->tot_x;
}
