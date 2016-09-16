#ifndef _ELECTRODE_H_
#define _ELECTRODE_H_

#include <stdarg.h>
#include <stdlib.h>

typedef struct _Electrode * Electrode;
typedef struct _Electrode const * const_Electrode;

//////////////////////////////////////////////////////////
/* Creates a new Electrode
	 ID defined where the electrode is
	 0 - (-(x,y,z))
	 1 - (+(x,y,z))
*/
Electrode newElectrode(void);

/* Deletes the electrode
*/
int deleteElectrode(Electrode *el);

/* Sets the tunneling constant for the electrode
*/
int setElectrode_alpha(Electrode el, double alpha);

/* Grab the tunneling constant for hopping 
	 from the electrode to the medium
*/
double getElectrode_alpha(Electrode el);

/* Sets the fermi energy of the electrode
*/
int setElectrode_FermiEnergy(Electrode el, double FermiE);

/* gets the fermi energy of the electrode
*/
double getElectrode_FermiEnergy(Electrode el);

/* Get number of charges in electrode
*/
int getElectrode_Charges(Electrode el);

/* Sets the number of charges on the electrode
*/
int setElectrode_Charges(Electrode el, double NumCharge);

/* Increases the number of charges on the electrode by 1
*/
int Electrode_addCharge(Electrode el);

/* Decreases the number of charges on the electrode by 1
	 if there are 0 and try to reduce will return -1
*/
int Electrode_minusCharge(Electrode el);

/* connects to a structure that has the pvals of hops to 
	 the site neighboring the electrode.
*/
int setElectrode_HopRates(Electrode el, void * HopS);

/* Grabs the datastructure containing hop rates off
	 the electrode
*/
void * getElectrode_HopRates(Electrode el);

/* Sets the data structure containing the adjacent site
	 information
*/
int setElectrode_AdjacentSites(Electrode el, void * snA);

/* Gets the data structure with the adjacent site information
*/
void * getElectrode_AdjacentSites(Electrode el);

/* Sets the sum of the electrode used when calculating
	 the transition time off the electrode
*/
int setElectrode_Sum(Electrode el, double sum);

/* Adds to the sum value of the electrode el
 */
int setElectrode_AddToSum(Electrode el, double val);

/* Grabs the sum of the electrode
*/
double getElectrode_Sum(Electrode el);

/* Prints contents of electrode id and p value
*/
int printElectrode(const_Electrode el);
#endif
