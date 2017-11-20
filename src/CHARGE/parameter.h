#define Length 30  //the number of nodes in x direction, the subscript is from 0 to 99
#define Width 30  //the number of nodes in y direction, the subscript is from 0 to 999
#define Height 20  //the number of nodes in z direction, the subscript is from 0 to 9
#define SiteDistance 1E-9  //unit:cm, the distance between two nodes is 1 nm
#define D 20E14 //dimension, to fit the real data
#define TCount 300  //total number of time steps
#define N 100  //number of charge injected per time step: N/TCount
#define Ntot TCount*N  //total number of charge in this model
#define TStep 1E-11  //the length of time between two data points in seconds NB: the charges move on their own tmes, this value is just the spacing between points in output files and numbers.
#define Nstep_av 20 // number of time steps to average on 
#define N_av TCount/Nstep_av //number of averages to look at
#define Rcount 1 // number of iterations with different random seeds
#define CutOff 4    //Defines the cutoff radius for the Correlation function, The final CutOffDistance is a multiple of CutOff*lambda
#define lambda 5E-9 //Defines lambda which is related is part of the correlation function
#define fracSeed  0.2  //Ratio of sites to be seeded [ 0 - 1 ]. For the rest of them a correlation function is used. Does not count sites allocated with trap sites as seeded
#define  E0 -5.2  //average site energy, HOMO of P3HT is -5.2ev
#define Etrap -5.06  //the meaning energy of trap states
#define q 1.602E-19  //the charge of an electron
#define voltage 0.06  //the voltage between two electrode

#define fraction 0//0.01;//0.01;  //the fraction of localised trap states
#define Tsigma  0.07  //the standard devation of trap states energy
#define sigma 0.07 //the standard deviation of site energy
#define KT 0.01524  //k =8.617E-5 ev/k; T = 293.15 k; KT = k * T 
#define reOrgEnergy 0.1  //reOrgEnergy is reorganisation energy 

