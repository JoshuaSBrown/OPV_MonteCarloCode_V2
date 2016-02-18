#ifndef READ_H_
#define READ_H_

typedef struct _ParameterFrame * ParameterFrame;
typedef struct _ParameterFrame const * const_ParameterFrame;

//Function takes text and finds where the location of pattern 
//in the text.
int match(char text[], char pattern[]);

//Given a string with a word followed by a number which
//is separated by the word with a single space will find
//and convert the number to an integer
int GrabInt(unsigned int position, char * buf);

//Given a string with a word followed by a number which
//is separated by the word with a single space will find
//and convert the number to a double
double GrabDouble(unsigned int position, char * buf);

//Function Reads all the parameters from the parameter 
//text file to read correctly the name of the variable 
//must be placed next on a new line followed by a number
//e.g.
//
//variable 33
int ReadParameter(int * method, int * SLength, int * SWidth, int * SHeight,\
		int * PeriodicX, int * PeriodicY, int * PeriodicZ,\
		int * EndX, int * EndY, int * EndZ,\
		int * XElecOn, int * YElecOn, int *ZElecOn,\
		double * XFermiB, double * XFermiF, double * YFermiL,\
		double * YFermiR, double * ZFermiB, double * ZFermiA,\
		double * alphaxb, double * alphaxf, double * alphayl,\
		double * alphayr, double * alphazb, double * alphaza,\
		double * vX, double * vY, double * vZ,\
		double * VoltageX, double * VoltageY, double * VoltageZ,\
		int * VStepX, int * VStepY, int * VStepZ,\
		double * VincX, double * VincY, double * VincZ,\
		double * SiteDistance, double * D, int * TCount,\
		int * NCh, int * Ntot, double * TStep, int * N_av,\
		int * Nstep_av, int * Time_check, int * Rcount,int * ClusterAlg, double * CutOff,\
		double * lambda, int * SeedProt, int * Attempts,\
		double * fracSeed, double * E0, double * sigma,\
		double * fracTrap, double * Etrap, double * Tsigma,\
		double * TempStart, int * TemperatureStep,\
		double * TemperatureInc, double * reOrgEnergy,\
		double * AttemptToHop, double * gamma,\
		double * RelativePerm, int * MovieFrames,\
		double * Tcv, double * Vcv);

int deleteParamFrame(ParameterFrame * PF);

ParameterFrame newParamFrame_File(void);

ParameterFrame newParamFrame(void);

int PFset_method(ParameterFrame Pf, int method);
int PFset_Len(ParameterFrame PF,int SLength);
int PFset_Wid(ParameterFrame PF,int SWidth);
int PFset_Hei(ParameterFrame PF,int SHeight);
int PFset_Px(ParameterFrame PF,int PeriodicX);
int PFset_Py(ParameterFrame PF,int PeriodicY);
int PFset_Pz(ParameterFrame PF,int PeriodicZ);
int PFset_EndX(ParameterFrame PF,int EndX);
int PFset_EndY(ParameterFrame PF,int EndY);
int PFset_EndZ(ParameterFrame PF,int EndZ);
int PFset_XElecOn(ParameterFrame PF,int XElecOn);
int PFset_YElecOn(ParameterFrame PF,int YElecOn);
int PFset_ZElecOn(ParameterFrame PF,int ZElecOn);
int PFset_XFermiB(ParameterFrame PF,double XFermiB);
int PFset_XFermiF(ParameterFrame PF,double XFermiF);
int PFset_YFermiL(ParameterFrame PF,double YFermiL);
int PFset_YFermiR(ParameterFrame PF,double YFermiF);
int PFset_ZFermiB(ParameterFrame PF,double ZFermiB);
int PFset_ZFermiA(ParameterFrame PF,double ZFermiA);
int PFset_alphaxb(ParameterFrame PF,double alphaxb);
int PFset_alphaxf(ParameterFrame PF,double alphaxf);
int PFset_alphayl(ParameterFrame PF,double alphayl);
int PFset_alphayr(ParameterFrame PF,double alphayr);
int PFset_alphazb(ParameterFrame PF,double alphazb);
int PFset_alphaza(ParameterFrame PF,double alphaza);
int PFset_vX(ParameterFrame PF,double vX);
int PFset_vY(ParameterFrame PF,double vY);
int PFset_vZ(ParameterFrame PF,double vZ);
int PFset_VoltageX(ParameterFrame PF,double VoltageX);
int PFset_VoltageY(ParameterFrame PF,double VoltageY);
int PFset_VoltageZ(ParameterFrame PF,double VoltageZ);
int PFset_VStepX(ParameterFrame PF,int VStepX);
int PFset_VStepY(ParameterFrame PF,int VStepY);
int PFset_VStepZ(ParameterFrame PF,int VStepZ);
int PFset_VincX(ParameterFrame PF,double VincX);
int PFset_VincY(ParameterFrame PF,double VincY);
int PFset_VincZ(ParameterFrame PF,double VincZ);
int PFset_SiteDist(ParameterFrame PF,double SiteDistance);
int PFset_D(ParameterFrame PF,double D);
int PFset_TCount(ParameterFrame PF,int TCount);
int PFset_NCh(ParameterFrame PF,int NCh);
int PFset_Ntot(ParameterFrame PF,int Ntot);
int PFset_TStep(ParameterFrame PF,double TStep);
int PFset_N_av(ParameterFrame PF,int N_av);
int PFset_Nstep_av(ParameterFrame PF,int Nstep_av);
int PFset_Time_check(ParameterFrame PF,int Time_check);
int PFset_Rcount(ParameterFrame PF,int Rcount);
int PFset_ClusterAlg(ParameterFrame PF,int ClusterAlg);
int PFset_CutOff(ParameterFrame PF,double CutOff);
int PFset_lambda(ParameterFrame PF,double lambda);
int PFset_SeedProt(ParameterFrame PF,int SeedProt);
int PFset_Attempts(ParameterFrame PF,int Attempts);
int PFset_FracSeed(ParameterFrame PF,double fracSeed);
int PFset_E0(ParameterFrame PF,double E0);
int PFset_sigma(ParameterFrame PF,double sigma);
int PFset_FracTrap(ParameterFrame PF,double fracTrap);
int PFset_Etrap(ParameterFrame PF,double Etrap);
int PFset_Tsigma(ParameterFrame PF,double Tsigma);
int PFset_TempStart(ParameterFrame PF,double TempStart);
int PFset_TempStep(ParameterFrame PF,int TempStep);
int PFset_TempInc(ParameterFrame PF,double TempInc);
int PFset_reOrg(ParameterFrame PF,double reOrg);
int PFset_AttemptToHop(ParameterFrame PF,double AttemptHop);
int PFset_gamma(ParameterFrame PF,double gamma);
int PFset_RelativePerm(ParameterFrame PF,double RelativePerm);
int PFset_MovieFrames(ParameterFrame PF, int MovieFrames);
double PFset_Tcv(ParameterFrame PF, double Tcv);

int PFget_method(ParameterFrame PF);
int PFget_Len(ParameterFrame PF);
int PFget_Wid(ParameterFrame PF);
int PFget_Hei(ParameterFrame PF);
int PFget_Px(ParameterFrame PF);
int PFget_Py(ParameterFrame PF);
int PFget_Pz(ParameterFrame PF);
int PFget_EndX(ParameterFrame PF);
int PFget_EndY(ParameterFrame PF);
int PFget_EndZ(ParameterFrame PF);
int PFget_XElecOn(ParameterFrame PF);
int PFget_YElecOn(ParameterFrame PF);
int PFget_ZElecOn(ParameterFrame PF);
double PFget_XFermiB(ParameterFrame PF);
double PFget_XFermiF(ParameterFrame PF);
double PFget_YFermiL(ParameterFrame PF);
double PFget_YFermiR(ParameterFrame PF);
double PFget_ZFermiB(ParameterFrame PF);
double PFget_ZFermiA(ParameterFrame PF);
double PFget_alphaxb(ParameterFrame PF);
double PFget_alphaxf(ParameterFrame PF);
double PFget_alphayl(ParameterFrame PF);
double PFget_alphayr(ParameterFrame PF);
double PFget_alphazb(ParameterFrame PF);
double PFget_alphaza(ParameterFrame PF);
double PFget_vX(ParameterFrame PF);
double PFget_vY(ParameterFrame PF);
double PFget_vZ(ParameterFrame PF);
double PFget_VoltageX(ParameterFrame PF);
double PFget_VoltageY(ParameterFrame PF);
double PFget_VoltageZ(ParameterFrame PF);
int PFget_VStepX(ParameterFrame PF);
int PFget_VStepY(ParameterFrame PF);
int PFget_VStepZ(ParameterFrame PF);
double PFget_VincX(ParameterFrame PF);
double PFget_VincY(ParameterFrame PF);
double PFget_VincZ(ParameterFrame PF);
double PFget_SiteDist(ParameterFrame PF);
double PFget_D(ParameterFrame PF);
int PFget_TCount(ParameterFrame PF);
int PFget_NCh(ParameterFrame PF);
int PFget_Ntot(ParameterFrame PF);
double PFget_TStep(ParameterFrame PF);
int PFget_N_av(ParameterFrame PF);
int PFget_Nstep_av(ParameterFrame PF);
int PFget_Time_check(ParameterFrame PF);
int PFget_Rcount(ParameterFrame PF);
int PFget_ClusterAlg(ParameterFrame PF);
double PFget_CutOff(ParameterFrame PF);
double PFget_lambda(ParameterFrame PF);
int PFget_SeedProt(ParameterFrame PF);
int PFget_Attempts(ParameterFrame PF);
double PFget_FracSeed(ParameterFrame PF);
double PFget_E0(ParameterFrame PF );
double PFget_sigma(ParameterFrame PF);
double PFget_FracTrap(ParameterFrame PF);
double PFget_Etrap(ParameterFrame PF);
double PFget_Tsigma(ParameterFrame PF);
double PFget_TempStart(ParameterFrame PF);
int PFget_TempStep(ParameterFrame PF);
double PFget_TempInc(ParameterFrame PF);
double PFget_reOrg(ParameterFrame PF);
double PFget_AttemptToHop(ParameterFrame PF);
double PFget_gamma(ParameterFrame PF);
double PFget_RelativePerm(ParameterFrame PF);
int PFget_MovieFrames(ParameterFrame PF);
int PFget_Tcv(ParameterFrame PF);
int PFget_Vcv(ParameterFrame PF);
#endif
