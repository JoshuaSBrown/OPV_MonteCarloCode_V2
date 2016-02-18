#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.h"
#include "../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.h"
#include "../PARAMETERS/read.h"
#include "../CHARGE/charge.h"

/* Here electricEnergy is the the energy between two sites caused
	 from the electric field. E.g.
	 electricField = voltage / (Length*SiteDistance);
	 electricEnergy = SiteDistance * electricField;

	 This function is responsible for setting up the Site energies of
	 all the sites in SNarray snA. It uses a correlation function and 
	 also takes into account the effect of traps. 
 */
int initSite(const double electricEnergyX, const double electricEnergyY,\
		const double electricEnergyZ, double KT, SNarray snA, ParameterFrame PF);

/* Initializes the Electrode hops onto the electrode and hops off of the 
	 electrode
 */
int initElec(const double electricEnergyX, const double electricEnergyY,\
		const double electricEnergyZ, const double MarcusCoeff,\
		const double KT, SNarray snA,\
		Electrode * elXb, Electrode * elXf, Electrode * elYl,\
		Electrode * elYr, Electrode * elZb, Electrode * elZa,\
		ParameterFrame PF);

/* Initilalizes parameters on the electrodes position on the x axis
 */
int initJumPossibility_ElecX( const double electricEnergyX,\
		const double electricEnergyY, const double electricEnergyZ,\
		const double SiteDistance, const double MarcusCoeff, const double KT, const double ReOrgEnergy,\
		matrix X1, matrix X2, Electrode elX, const double RelativePermittivity,\
		const double vX, const int SWidth, const int SHeight,\
		const int PeriodicY, const int PeriodicZ,\
		const int YElecOn, const int ZElecOn, const int BorF);

/* Initilalizes parameters on the electrodes position on the x axis
 */
int initJumPossibility_ElecY( const double electricEnergyX,\
		const double electricEnergyY, const double electricEnergyZ,\
		const double SiteDistance, const double MarcusCoeff, const double KT, const double ReOrgEnergy,\
		matrix Y1, matrix Y2, Electrode elY, const double RelativePermittivity,\
		const double vY,const int SLength, const int SHeight,\
		const int PeriodicX, const int PeriodicZ,\
		const int XElecOn, const int ZElecOn, const int LorR);

/* Initializes parameters on the electrodes position on the z axis
*/
int initJumPossibility_ElecZ( const double electricEnergyX,\
		const double electricEnergyY, const double electricEnergyZ,\
		const double SiteDistance,const double MarcusCoeff, const double KT, const double ReOrgEnergy,\
		matrix Z1, matrix Z2, Electrode elZ, const double RelativePermittivity,\
		const double vZ, const int SLength, const int SWidth,\
		const int PeriodicX, const int PeriodicY,\
		const int XElecOn, const int YElecOn, const int AorB);


int initJumPossibility(const double electricEnergyX,const double electricEnergyY,\
		const double electricEnergyZ, const double MarcusCoeff,\
		const double KT,const double reOrgEnergy,SNarray snA,\
		const int PeriodicX, const int PeriodicY, const int PeriodicZ,\
		const int XElecOn, const int YElecOn, const int ZElecOn);

/* InitCharget0 initializes the charge data structure
	 It returns the charge data structure
 */
ChargeArray initCharget0( matrix Sequence, const_SNarray snA,  const int Ntot, const int NCh,\
		const double D, const int XElecOn, const int YElecOn, const int ZElecOn,\
		const int EndX, const int EndY, const int EndZ);

/* Init Charge responsible for defining initial conditions of charges
	 as they are added at each timestep. 

	 chA  	contains all the charges and the properties associated with them
	 sequence  	contains the order the charges should be called based on 
	 their dwelltime
	 snA 		contains the properties associated with each of the sites
	 nc 		is the number of charges currently in the styem
	 n			is the time step iterating over. 
 */
int initCharge(int nca, long int n, ChargeArray *chA, matrix Sequence, SNarray snA,\
		const int Ntot, const int NCh, const double D,\
		const int XElecOn, const int YElecOn, const int ZElecOn,\
		const int EndX, const int EndY, const int EndZ);

/* Quicksort sorts the charges based on their dwell times. It
	 calls the partition function to determine the pivot and
	 recursively breaks up the charge data structure in order to sort
	 all the charges quickly. It then stores the correct order in 
	 the sequence array.
	 .
 */
void quickSort( int low, int high, matrix Sequence, ChargeArray chA);

/* The partition function is used to find the pivot of the
	 charge data structure and sort the structure based on
	 the dwelltimes
 */
int partition( int low, int high, matrix Sequence, ChargeArray chA);
#endif
