#ifndef _CHARGE_H_
#define _CHARGE_H_

/* The type declartion of the ADT
 */
typedef struct _Charge * Charge;
typedef struct _Charge const * const_Charge;

typedef struct _ChargeArray * ChargeArray;
typedef struct _ChargeArray const * const_ChargeArray;

/* Create type Charge
 */
Charge newCharge(void);

/* Create charge Array with a total of len elements.
 * Returns NULL if len <0, chA = NULL, or chA->C[i] = NULL. 
 * Returns chA if succesful. 
*/
ChargeArray newChargeA(int len);

/* Delete Charge, Returns 0 if successful, -1 if ch is null
 */
int deleteCharge(Charge ch);

/* Delete Charge Array and all charges within the
   array. Returns 0 if successful, -1 if chA is null.*/
int deleteChargeA(ChargeArray chA);

/* gets the charge arrays length
*/
int getChargeA_len(ChargeArray chA);

/* Returns a Charge i from charge array chA the value of
   i should be 0<=i<chA->length. Returns -1 if i is outside
   the range or chA is NULL. 
 */
Charge getCharge(const_ChargeArray chA, int i);

/* Print properties of Charge ch which includes
   x, y and z position
   run_Xdist
   t
   dwelltime
   Returns 0 if successful and -1 if ch = NULL*/
int printCharge(const_Charge ch);

/* Print all the properties of each charge in Charge 
   Array chA*/
int printChargeA(const_ChargeArray chA);

/* Returns the dwell time of charge ch
	 */
double getDwel(const_Charge ch);

double gett(const_Charge ch);

int getCx(const_Charge ch);

int getCy(const_Charge ch);

int getCz(const_Charge ch);

int setDwel( Charge ch, double dw);

// Minus dw from the current dwelltime
int MinusDwel( Charge ch, double dw);

// Set ch->t to t
int sett(Charge ch, long double t);

// Add time t to current ch->t
int Plust(Charge ch, long double t);

// Move charge forward one increment (one site) in x direction
// also increments x_tot
//Return 0 if succesful, -1 if ch is NULL. 
int CxPlus(Charge ch);

// Move charge back one increment (one site) in x direction.
// also decrements x_tot
//Return 0 if succesful, -1 if ch is NULL. 
int CxMinus(Charge ch);

int CyPlus(Charge ch);
int CyMinus(Charge ch);
int CzPlus(Charge ch);
int CzMinus(Charge ch);

// Move the charge forward one increment (one site) in the y 
// direction with the periodicity defined by w. Returns 0 if
//successful, -1 if ch is NULL or w is negative, and -2 if y = w. 
int CyPlusPeriodic(Charge ch, int w);

// Move the charge back one increment (one site) in the y 
// direction with the periodicity defined by w. Returns 0 if
//successful, -1 if ch is NULL or w is negative, and -2 if y = w 
//or y =-1. 
int CyMinusPeriodic(Charge ch, int w); 

// Move the charge forward one increment (one site) in the z 
// direction with the periodicity defined by w. Returns 0 if
//successful, -1 if ch is NULL or w is negative, and -2 if y = w. 
int CzPlusPeriodic(Charge ch, int w);

// Move the charge back one increment (one site) in the z
// direction with the periodicity defined by w. Returns 0 if
//successful, -1 if ch is NULL or w is negative, and -2 if y = w. 
int CzMinusPeriodic(Charge ch, int w); 

// Move the charge forward one increment (one site) in the x 
// direction and the total Xdist with the periodicity defined by w.
// Returns 0 if successful, -1 if ch is NULL or w is negative, and -2 if
// y = w. 
int CxPlusPeriodic(Charge ch, int w);

// Move the charge back one increment (one site) in the x
// direction and the total Xdist with the periodicity defined by w. 
//Returns 0 if successful, -1 if ch is NULL or w is negative, and -2 if 
//y = w. 
int CxMinusPeriodic(Charge ch, int w); 

//returns -1 if ch is NULL or if cx is negative or zero. 
//If successful, sets x and tot_x to cx, and returns 0. 
int setCx(Charge ch, int cx);

//returns 0 if successful, -1 if ch is NULL or if cx is negative
//or zero
int setCy(Charge ch, int cy);

//returns 0 if successful, -1 if ch is NULL or if cx is negative
//or zero
int setCz(Charge ch, int cz);

//returns the total x distance travelled if successful, -1 if ch is NULL
int getXdist(const_Charge ch);

#endif
