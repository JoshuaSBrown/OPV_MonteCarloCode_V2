#include <stdio.h>
#include <assert.h>

#include "charge.h"
#include "parameter.h"

int main(void){
	Charge chNull = NULL;
	ChargeArray chANull = NULL;
	Charge chrv;
	ChargeArray chArv; 
	int rv; 
	int len = 5; 
	int val = -999; 
	int spot = 3; //spot in array where testing functions
	double drv; 
	//valid:
	//Charge ch
	//ChargeArray chA

	printf("Testing: newCharge\n");
	Charge ch = newCharge();
	assert(ch !=NULL);  
	
	printf("Testing: newChargeA\n");
	ChargeArray chA = newChargeA(len);
	assert(chA !=NULL);  
	chArv = newChargeA(val);
	assert(chArv ==NULL);  

	printf("Testing: setDwel\n");
	rv= setDwel(chNull, 15.75);
	assert(rv==-1);
	rv= setDwel(ch, val);
	assert(rv==-1);
	rv= setDwel(ch, 15.75);
	assert(rv==0);
	
	printf("Testing: MinusDwel\n");
	rv= MinusDwel(chNull, 1.05);
	assert(rv==-1);
	rv= MinusDwel(ch, val);
	assert(rv==-1);
	rv= MinusDwel(ch, 14.25);
	assert(rv==0);
	
	printf("Testing: getDwel\n");
	drv= getDwel(chNull);
	assert(drv==-1);
	drv= getDwel(ch);
	printf("%f\n", drv);
	assert(drv==1.5);

	printf("Testing: printCharge, Expecting: dwelltime = 1.5 \n");
	rv= printCharge(chNull);
	assert(rv==-1); 
	rv= printCharge(ch);
	assert(rv==0); 

	printf("Testing: sett\n"); 
	rv= sett(chNull, 2.451);
	assert(rv==-1);
	rv= sett(ch, val);
	assert(rv==-1);
	rv= sett(ch, 2.451);
	assert(rv==0);
	printf("Testing: Plust\n"); 
	rv= Plust(chNull, 3.5);
	assert(rv==-1);
	rv= Plust(ch, val);
	assert(rv==-1);
	rv= Plust(ch, 3.5);
	assert(rv==0);
	printf("Testing: gett\n");
	drv= gett(chNull);
	assert(drv==-1);
	drv= gett(ch);
	printf("%f\n", drv);
	assert(drv==(2.451+3.5)); 

	printf("Testing: getCharge\n");
	chrv= getCharge(chANull, 4);
	assert(chrv==NULL); 
	chrv= getCharge(chA, -val); // i > len
	assert(chrv==NULL); 
	chrv= getCharge(chA, val);
	assert(chrv==NULL); 
	chrv= getCharge(chA, 4);
	assert(chrv!=NULL); 	

	printf("Testing: setCx\n"); 
	rv= setCx(chNull, 6);
	assert(rv==-1);
	rv= setCx(getCharge(chA, spot), val);
	assert(rv==0);
	rv= setCx(getCharge(chA, spot), 6);
	assert(rv==0);
	printf("Testing: CxPlus\n"); //x=8
	rv= CxPlus(chNull);
	assert(rv==-1);
	rv= CxPlus(getCharge(chA, spot));
	rv= CxPlus(getCharge(chA, spot)); //done twice
	assert(rv==0);
	printf("Testing: CxMinus\n"); //x=7 //tot_x = 7
	rv= CxMinus(chNull);
	assert(rv==-1);
	rv= CxMinus(getCharge(chA, spot)); //done once
	assert(rv==0);
	printf("Testing: CxPlusPeriodic\n"); //x=0  //tot_x = 8
	rv= CxPlusPeriodic(chNull, 2);
	assert(rv==-1);
	rv= CxPlusPeriodic(getCharge(chA, spot), val);
	assert(rv==-1);
	rv= CxPlusPeriodic(getCharge(chA, spot), 2); //add once
	assert(rv==0);
	printf("Testing: CxMinusPeriodic\n"); //x=4 //tot_x = 7
	rv= CxMinusPeriodic(chNull, 5);
	assert(rv==-1);
	rv= CxMinusPeriodic(getCharge(chA, spot), val);
	assert(rv==-1);
	rv= CxMinusPeriodic(getCharge(chA, spot), 5);
	assert(rv==0);
	printf("Testing: getCx\n");
	rv= getCx(chNull);
	assert(rv==-1);
	rv= getCx(getCharge(chA, spot));
	assert(rv==4);
	printf("Testing: getXdist\n");
	rv= getXdist(chNull);
	assert(rv==-1);
	rv= getXdist(getCharge(chA, spot));
	assert(rv==7);
	printf("Testing: printChargeA, Expecting: ch[3]->x = 4\n");
	rv= printChargeA(chANull);
	assert(rv==-1);
	rv= printChargeA(chA);
	assert(rv==0);

	printf("Testing: setCy\n"); 
	rv= setCy(chNull,8);
	assert(rv==-1);
	rv= setCy(ch, val);
	assert(rv==0);
	rv= setCy(ch, 8);
	assert(rv==0);
	printf("Testing: CyPlus\n"); //y=9 (after function successfully called)
	rv= CyPlus(chNull);
	assert(rv==-1);
	rv= CyPlus(ch);
	assert(rv==0); 
	printf("Testing: CyMinus\n");//y=7
	rv= CyMinus(chNull);
	assert(rv==-1);
	rv= CyMinus(ch);
	rv= CyMinus(ch); //subtract twice
	assert(rv==0); 
	printf("Testing: CyPlusPeriodic\n"); //y=0 
	rv= CyPlusPeriodic(chNull, 2);
	assert(rv==-1);
	rv= CyPlusPeriodic(ch, val);
	assert(rv==-1);
	rv= CyPlusPeriodic(ch, 2); //add once
	assert(rv==0);
	printf("Testing: CyMinusPeriodic\n"); //y=4
	rv= CyMinusPeriodic(chNull, 5);
	assert(rv==-1);
	rv= CyMinusPeriodic(ch, val);
	assert(rv==-1);
	rv= CyMinusPeriodic(ch, 5);
	assert(rv==0);
	printf("Testing: getCy\n");
	rv= getCy(chNull);
	assert(rv==-1);
	rv= getCy(ch);
	assert(rv==4);
	printf("Testing: setCz\n"); 
	rv= setCz(chNull,16);
	assert(rv==-1);
	rv= setCz(ch, val);
	assert(rv==0);
	rv= setCz(ch, 16);
	assert(rv==0);
	printf("Testing: CzPlus\n"); //z=18
	rv= CzPlus(chNull);
	assert(rv==-1);
	rv= CzPlus(ch);	
	rv= CzPlus(ch); //add twice
	assert(rv==0);
	printf("Testing: CzMinus\n"); //z=17
	rv= CzMinus(chNull);
	assert(rv==-1);
	rv= CzMinus(ch);
	assert(rv==0);
	printf("Testing: CzPlusPeriodic\n"); //z=0
	rv= CzPlusPeriodic(chNull, 5);
	assert(rv==-1);
	rv= CzPlusPeriodic(ch, val);
	assert(rv==-1);
	rv= CzPlusPeriodic(ch, 6);
	assert(rv==0);
	printf("Testing: CzMinusPeriodic\n"); // z=12
	rv= CzMinusPeriodic(chNull, 13);
	assert(rv==-1); 
	rv= CzMinusPeriodic(ch, val);
	assert(rv==-1); 
	rv= CzMinusPeriodic(ch, 13);
	assert(rv==0);
	printf("Testing: getCz\n");
	rv= getCz(chNull);
	assert(rv==-1);
	rv= getCz(ch);
	assert(rv==12);
	printf("Testing: deleteCharge\n");
	rv= deleteCharge(chNull);
	assert(rv==-1); 
	rv= deleteCharge(ch);
	assert(rv==0); 
	printf("Testing: deleteChargeA\n");
	rv= deleteChargeA(chANull);
	assert(rv==-1); 
	rv= deleteChargeA(chA);
	assert(rv==0);

	printf("Mission Completed\n");
	return 0;
}

