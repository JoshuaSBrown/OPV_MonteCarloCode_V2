#ifndef _CHARGETRANSPORT_H_
#define _CHARGETRANSPORT_H_

#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CHARGE/charge.h"
#include "../PARAMETERS/read.h"

/* Function designed to take care of hops from site neighboring an electrode. Determines
	 if the future site a charge will hop too if the site is still available when the 
	 appropriate amount of time has passed. 
	 IDs for the electrodes are as follows:
	 getAtotal(snA) + 0 - (-x)
	 getAtotal(snA) + 1 - (+x)
	 getAtotal(snA) + 2 - (-y)
	 getAtotal(snA) + 3 - (+y)
	 getAtotal(snA) + 4 - (-z)
	 getAtotal(snA) + 5 - (+z)
 */
int HopToElecX(SNarray snA, Electrode elXb, Charge * one, int * future, int EndY, int EndZ);
int HopToElecY(SNarray snA, Electrode elYl, Charge * one, int * future, int EndX, int EndZ);
int HopToElecZ(SNarray snA, Electrode elZb, Charge * one, int * future, int EndX, int EndY);

/* Function determines which site a charge to when hopping off the electrode
 */
int ElecHopOffX(Electrode el, int * future, SNarray snA);
int ElecHopOffY(Electrode el, int * future, SNarray snA);
int ElecHopOffZ(Electrode el, int * future, SNarray snA);


/* This function is used in the case that the charge is not 
	 attached to a cluster and therefore will only move to 
	 sites neighboring it.

	 EndX, EndY, EndZ defines at what site the sample ends this 
	 is important when using periodic conditions.	
	 A value of:
	 0 - means infinetly periodic can not be used if those electrodes
	 are not turned on because the charges will never be able to escape. 
	 1 - and above describes how many samples are periodic before reaches 
	 the end electrode. (Periodic for 1 sample is the same as non-periodic)
	 2 - and above how many samples periodic for
 */
int SiteHop(SNarray snA, Charge * ch, SiteNode site, int * future, const int EndX, const int EndY, const int EndZ, const int XElecOn, const int YElecOn, const int ZElecOn, const int PeriodicX, const int PeriodicY, const int PeriodicZ);

/* This function is used in the case that a charge is part of a cluster
	 Determines whether chooses a site for the charge to hop to
	 WARNING does not change the DwelStat of the site the charge hopped to
	 or from
	 Will return a value of:
	 -1 - if malformed input
	 0 - if charge hopped to new location
 */
int ClusterHop(SNarray snA, Charge * ch, double * time, int * newID); 

/* This function actually moves the charge will return:
	 newID is the id of the site the charge has chosen to
	 hop to.
	 0 - if sucessful
	 1 - if site is already occupied
 */
int MakeHop(SNarray snA, int newID, Charge * ch, int * totalX, int * totalY, int * totalZ,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
		const int XElecOn, const int YElecOn, const int ZElecOn);

int initFutureSite( SNarray * snA, matrix * FutureSite, ChargeArray *chA, ParameterFrame PF,\
										Electrode elXb, Electrode elYl, Electrode elZb);

/* This function is responsible for correctly generating all the parameters needed to run 
	 the randomWalk function including:
	 t - global time of the simulation
	 Sequence - the order charges are called 
	 chA - data related to the charges and thier locations
	 snA - data related to sites their locations and energies
	 PF - all the parameters input by the user
	 n - the number of steps that have been passed that charges are injected
	 nc - number of charges in the system (excluding the electrodes)
	 nca - total number of active charges 

	 Finally the varias electrodes
*/
int Pre_randomWalk(const int CheckPtStatus, char * CheckPtVersion,char * FileName,\
									 long double * t,\
									 matrix * Sequence, ChargeArray * chA, matrix * FuturSite,\
									 ArbArray * ClArLL, SNarray * snA,\
									 ParameterFrame PF,double electricEnergyX, double ElectricEnergyY,\
									 double ElectricEnergyZ, int r, double Vx, double Vy, double Vz,\
									 double Temperature, long int * n, int * nc, int * nca,\
									 Electrode * elXb, Electrode * elXf, Electrode * elYl,\
									 Electrode * elYr, Electrode * elZb, Electrode * elZa);

/* This is the heart of the program injects charges and allows them to hop
	 through the system till they reach an electrode
 */
int randomWalk( SNarray snA, int CheckPointNum,\
		char * FileName, double ElectricFieldX,\
		double ElectricFieldY, double ElectricFieldZ,\
		Electrode elXb, Electrode elXf, Electrode elYl,\
		Electrode elYr, Electrode elZb, Electrode elZa,\
		ParameterFrame PF, long double t, matrix Sequence,\
		matrix FutureSite, ChargeArray * chA,\
		long int n, int nc, int nca);

int UpdateOccTime(SNarray * snA, Charge * ch,double tim, ParameterFrame PF);

/* This function is responsibly for appropriately removing the ClArLL, the snA and 
	 the electrodes
*/
int Post_randomWalk(ArbArray ClArLL, SNarray snA, Electrode elXb, Electrode elXf,\
									  Electrode elYl, Electrode elYr, Electrode elZb, Electrode elZa,\
										ParameterFrame PF);

/* If a .ckpt file does exist will determine whether it is the correct kind and wheather
	 it is the latest version if it is not it will find the latest version and return it
*/
int CheckPt_Test(int * CheckPtNum, int CheckFileExist, char * FileNameCheckPtVersion,\
								 int FileNameSize, const double Vx, const double Vy, const double Vz, const double Temperature);

/* Check if Checkpoint files exist
	 pass a character array e.g.
	 Passes back the first .ckpt file
	 so we can grab the Parameter frame
	 char File[256]
	 And the size of the character array

	 CheckPt_exist(File,sizeof(File));
*/
int CheckPt_exist(char * File, int buffersize);

/* Check if .cluster file exists in the CLUSTERFILE folder
	 returns 0 if exists
	 returns -1 if does not
*/
int CheckPt_Cluster(const double Vx, const double Vy, const double Vz, const double T, int r);

/* Grabs the latest .ckpt file of the particular version
*/
int CheckPt_Latest(char * File, int buffersize, double Vx, double Vy, double Vz, double T);

/* Load CheckPoint file
	 */
int Load_CheckPt(long double * t, SNarray * snA, ChargeArray * chA, matrix * Sequence,\
								 matrix * FutureSite, char * FileName, ParameterFrame *PF,long int * n, int * nc, int * nca,\
								 int * Num_elXb, int * Num_elXf, int * Num_elYl, int * Num_elYr, int * Num_elZb,\
								 int * Num_elZa);

/* Grab parameter frame from checkpoint file
*/
int Load_CheckPt_PF(char * FileName, ParameterFrame *PF);

/* Grabs only the data from a checkpoint file
*/
int Load_CheckPt_Data(long double * t, SNarray * snA, ChargeArray * chA,matrix * Sequence,\
		      matrix * FutureSite, char * FileName,long int * n,int * nc,int * nca,\
		      int * Num_elXb, int * Num_elXf, int * Num_elYl,\
		      int * Num_elYr, int * Num_elZb, int * Num_elZa,\
		      double Vx, double Vy, double Vz);



/* Checkpoint file
*/
int Save_CheckPt(char * FileName, int * CheckptNum, SNarray snA, ChargeArray chA, matrix Sequence,matrix FutureSite, long double t,\
		 						 ParameterFrame PF, long int n, int nc, int nca, Electrode elXb, Electrode elXf, Electrode elYl, Electrode elYr,\
								 Electrode elZb, Electrode elZa);

/* prints xys movie files
*/
int printMovie(int * Movie, long double t,char * FileName, SNarray snA,ParameterFrame PF);

/* Prints the currents flowing in the x,y and z direction durint each time
	 interval. As well as the number of charges reaching the x,y and z source
	 and drain electrodes (if the electrodes are turned on).
 */
int printTransportData( matrix System, matrix timeArray, matrix Xcurrent, matrix Ycurrent, matrix Zcurrent,\
		matrix Xelec_Drain, matrix Yelec_Drain, matrix Zelec_Drain,\
		matrix Xelec_Source, matrix Yelec_Source, matrix Zelec_Source,\
		matrix Xvelocity, matrix Yvelocity, matrix Zvelocity,\
		int XElecOn, int YElecOn, int ZElecOn, char * FileName,\
		double ElectricFieldX, double ElectricFieldY, double ElectricFieldZ);

/* If charge is hopping to an electrode
*/
int ChargeElectrode(Electrode el, Charge * one, matrix * Sequence, ChargeArray chA, int * nc, const int nca, int flag);

/* If charge is returning from an electrode to the system
 */
int ChargeSystem(int * nc, Electrode el, const int nca, matrix Sequence, ChargeArray chA);

/*	This function keeps track of all the operations that need to be done
		if a charge reaches an electrode.
 */
int ChargeClosure(Electrode el, Charge * one, matrix * Sequence, int * nc,int * nca, int * TotalCollected, int Ntot);

/* After a charge makes an action (tries to hop or does hop)
	 a new dwell time is assigned to it. The order by which charges
	 are moved is dependent on their dwell time and stored in the 
	 Sequence matrix. After the action the charge that was first in
	 the Sequence needs to be placed in a different spot in the 
	 Sequence depending on its new dwell time. That's what this 
	 function does.
WARNING: this function only looks at the first charge in the 
sequence it does not sort the whole sequence
 */
int insertDwelltimePos(int nc, ChargeArray chA, matrix * Sequence);

/* Checks if a charge has hopped to an electrode
 */
int CheckIfElecHop(int xx, int yy, int zz,\
		int * Xsource, int * Ysource, int * Zsource,\
		int * Xdrain, int * Ydrain, int * Zdrain,\
		int * TrackX, int * TrackY, int * TrackZ,\
		long double * TotalVelX, long double * TotalVelY, long double * TotalVelZ,\
		const int SLength, const int SWidth, const int SHeight,\
		const int EndX, const int EndY, const int EndZ,\
		int * NumAvgVel, long double * TimeTrack1, long double t, double SiteDistance,\
		int * nc, int * ElecExit,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ);

/* Check if the next hop should use the cluster hopping algorithm
 */
int ClusterHopCheck(const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
		const int SLength, const int SWidth, const int SHeight,\
		const int EndX, const int EndY, const int EndZ,\
		const int x, const int y, const int z,SiteNode site,\
		const int Cx, const int Cy, const int Cz,\
		int * CheckX, int * CheckY, int * CheckZ);

/* Saves Data point
 */
int SaveDataPoint(int * CurrentInc, int * NumAvgVel, int nc, int XElecOn, int YElecOn, int ZElecOn,\
		matrix * Xcurrent, matrix * Ycurrent, matrix * Zcurrent, matrix *timeArray,\
		matrix * Xvelocity, matrix * Yvelocity, matrix * Zvelocity,matrix * System,\
		matrix * Xelec_Drain, matrix * Yelec_Drain, matrix * Zelec_Drain,\
		matrix * Xelec_Source, matrix * Yelec_Source, matrix * Zelec_Source,\
		int * TotalX, int * TotalY, int * TotalZ, int * TrackX, int * TrackY,\
		int * TrackZ, long double t, long double * TimeTrack1, double TStep, double SiteDistance,\
		long double * TotalVelX,long double * TotalVelY,long double * TotalVelZ,\
		int * Xdrain, int * Ydrain, int * Zdrain, int * Xsource, int * Ysource, int * Zsource,\
	  const double ElectricFieldX, const double ElectricFieldY, const double ElectricFieldZ,\
	  char * FileName);


/* For a SiteNode which has surrouding sites determines which of the 8 
	 surrounding sites a charge is randomly likely to land on
 */
int HoppingToSurroundingSites(SiteNode site, int codeX, int codeY, int codeZ);

/* This function prints out the energies of all the sites showing the effect
	 of the electric field on the energies as well as the image force
*/
int printFileEnergy(const_SNarray snA, char * FileName,\
		double ElectricEnergyX, double ElectricEnergyY, double ElectricEnergyZ,\
		ParameterFrame PF);

#endif
