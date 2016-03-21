#include <stdlib.h>
#include <stdio.h>

#include "read.h"

int main(void){

	int method;

	int SLength;
	int SWidth;
	int SHeight;

	int PeriodicX;
	int PeriodicY;
	int PeriodicZ;

	int EndX;
	int EndY;
	int EndZ;

	int XElecOn;
	int YElecOn;
	int ZElecOn;

	double XFermiB;
	double XFermiF;
	double YFermiL;
	double YFermiR;
	double ZFermiB;
	double ZFermiA;

	double alphaxb;
	double alphaxf;
	double alphayl;
	double alphayr;
	double alphazb;
	double alphaza;

	double vX;
	double vY;
	double vZ;

	double VoltageX;
	double VoltageY;
	double VoltageZ;

	int VStepX;
	int VStepY;
	int VStepZ;

	double VincX;
	double VincY;
	double VincZ;

	double SiteDistance;
	double D;
	int TCount;
	int NCh;
	int Ntot;
	double TStep;
	int N_av;
	int Nstep_av;
	int Nstep_check;
	int Rcount;

	double CutOff;
	double lambda;

	int SeedProt;
	int Attempts;

	double fracSeed;
	double E0;
	double sigma;

	double fracTrap;
	double Etrap;
	double Tsigma;

	double TempStart;
	int TemperatureStep;
	double TemperatureInc;
	double reOrgEnergy;

	double AttemptToHop;
	double gamma;
	double RelativePerm;

	int MovieFrames;

	double Vcv;
	double Tcv;

	ParameterFrame PF = newParamFrame();

	printf("\nTesting: newParamFrame\n\n");
	printf("Testing:PFget_method\n");
	printf("method %d\n",PFget_method(PF));
	printf("Testing: PFget_Len\n");
	printf("SLength %d\n",PFget_Len(PF));
	printf("Testing: PFget_Wid\n");
	printf("SWidth %d\n",PFget_Wid(PF));
	printf("Testing: PFget_Hei\n");
	printf("SHeight %d\n",PFget_Hei(PF));
	printf("Testing: PFget_Px\n");
	printf("PeriodicX %d\n",PFget_Px(PF));
	printf("Testing: PFget_Py\n");
	printf("PeriodicY %d\n",PFget_Py(PF));
	printf("Testing: PFget_Pz\n");
	printf("PeriodicZ %d\n",PFget_Pz(PF));
	printf("Testing: PFget_EndX\n");
	printf("EndX %d\n",PFget_EndX(PF));
	printf("Testing: PFget_EndY\n");
	printf("EndY %d\n",PFget_EndY(PF));
	printf("Testing: PFget_EndZ\n");
	printf("EndZ %d\n",PFget_EndZ(PF));
	printf("Testing: PFget_XElecOn\n");
	printf("XElecOn %d\n",PFget_XElecOn(PF));
	printf("Testing: PFget_YElecOn\n");
	printf("YElecOn %d\n",PFget_YElecOn(PF));
	printf("Testing: PFget_ZElecOn\n");
	printf("ZElecOn %d\n",PFget_ZElecOn(PF));
	printf("Testing: PFget_XFermiB\n");
	printf("XFermiB %f\n",PFget_XFermiB(PF));
	printf("Testing: PFget_XFermiF\n");
	printf("XFermiF %f\n",PFget_XFermiF(PF));
	printf("Testing: PFget_YFermiL\n");
	printf("YFermiL %f\n",PFget_YFermiL(PF));
	printf("Testing: PFget_YFermiR\n");
	printf("YFermiR %f\n",PFget_YFermiR(PF));
	printf("Testing: PFget_ZFermiB\n");
	printf("ZFermiB %f\n",PFget_ZFermiB(PF));
	printf("Testing: PFget_ZFermiA\n");
	printf("ZFermiA %f\n",PFget_ZFermiA(PF));
	printf("Testing: PFget_alphaxb\n");
	printf("alphaxb %f\n",PFget_alphaxb(PF));
	printf("Testing: PFget_alphaxf\n");
	printf("alphaxf %f\n",PFget_alphaxf(PF));
	printf("Testing: PFget_alphayl\n");
	printf("alphayl %f\n",PFget_alphayl(PF));
	printf("Testing: PFget_alphayr\n");
	printf("alphayr %f\n",PFget_alphayr(PF));
	printf("Testing: PFget_alphazb\n");
	printf("alphazb %f\n",PFget_alphazb(PF));
	printf("Testing: PFget_alphaza\n");
	printf("alphaza %f\n",PFget_alphaza(PF));
	printf("Testing: PFget_vX\n");
	printf("vX %g\n",PFget_vX(PF));
	printf("Testing: PFget_vY\n");
	printf("vY %g\n",PFget_vY(PF));
	printf("Testing: PFget_vZ\n");
	printf("vZ %g\n",PFget_vZ(PF));
	printf("Testing: PFget_VoltageX\n");
	printf("VoltageX %g\n",PFget_VoltageX(PF));
	printf("Testing: PFget_VoltageY\n");
	printf("VoltageY %g\n",PFget_VoltageY(PF));
	printf("Testing: PFget_VoltageZ\n");
	printf("VoltageZ %g\n",PFget_VoltageZ(PF));
	printf("Testing: PFget_VStepX\n");
	printf("VStepX %d\n",PFget_VStepX(PF));
	printf("Testing: PFget_VStepY\n");
	printf("VStepY %d\n",PFget_VStepY(PF));
	printf("Testing: PFget_VStepZ\n");
	printf("VStepZ %d\n",PFget_VStepZ(PF));
	printf("Testing: PFget_VincX\n");
	printf("VincX %g\n",PFget_VincX(PF));
	printf("Testing: PFget_VincY\n");
	printf("VincY %g\n",PFget_VincY(PF));
	printf("Testing: PFget_VincZ\n");
	printf("VincZ %g\n",PFget_VincZ(PF));
	printf("Testing: PFget_SiteDist\n");
	printf("SiteDistance %g\n",PFget_SiteDist(PF));
	printf("Testing: PFget_D\n");
	printf("D %g\n",PFget_D(PF));
	printf("Testing: PFget_TCount\n");
	printf("TCount %d\n",PFget_TCount(PF));
	printf("Testing: PFget_NCh\n");
	printf("NCh %d\n",PFget_NCh(PF));
	printf("Testing: PFget_Ntot\n");
	printf("Ntot %d\n",PFget_Ntot(PF));
	printf("Testing: PFget_TStep\n");
	printf("TStep %g\n",PFget_TStep(PF));
	printf("Testing: PFget_N_av\n");
	printf("N_av %d\n",PFget_N_av(PF));
	printf("Testing: PFget_Nstep_av\n");
	printf("Nstep_av %d\n",PFget_Nstep_av(PF));
	printf("Testing: PFget_Time_check\n");
	printf("Time_check %d\n",PFget_Time_check(PF));
	printf("Testing: PFget_Rcount\n");
	printf("Rcount %d\n",PFget_Rcount(PF));
	printf("Testing: PFget_CutOff\n");
	printf("CutOff %g\n",PFget_CutOff(PF));
	printf("Testing: PFget_lambda\n");
	printf("lambda %g\n",PFget_lambda(PF));
	printf("Testing: PFget_SeedProt\n");
	printf("SeedProt %d\n",PFget_SeedProt(PF));
	printf("Testing: PFget_Attempts\n");
	printf("Attempts %d\n",PFget_Attempts(PF));
	printf("Testing: PFget_FracSeed\n");
	printf("FracSeed %f\n",PFget_FracSeed(PF));
	printf("Testing: PFget_E0\n");
	printf("E0 %f\n",PFget_E0(PF));
	printf("Testing: PFget_sigma\n");
	printf("sigma %f\n",PFget_sigma(PF));
	printf("Testing: PFget_FracTrap\n");
	printf("FracTrap %f\n",PFget_FracTrap(PF));
	printf("Testing: PFget_Etrap\n");
	printf("Etrap %f\n",PFget_Etrap(PF));
	printf("Testing: PFget_Tsigma\n");
	printf("Tsigma %f\n",PFget_Tsigma(PF));
	printf("Testing: PFget_TempStart\n");
	printf("TempStart %f\n",PFget_TempStart(PF));
	printf("Testing: PFget_TempStep\n");
	printf("TempStep %d\n",PFget_TempStep(PF));
	printf("Testing: PFget_TempInc\n");
	printf("TempInc %f\n",PFget_TempInc(PF));
	printf("Testing: PFget_reOrg\n");
	printf("reOrg %f\n",PFget_reOrg(PF));
	printf("Testing: PFget_AttemptToHop\n");
	printf("AttemptToHop %g\n",PFget_AttemptToHop(PF));
	printf("Testing: PFget_gamma\n");
	printf("gamma %f\n",PFget_gamma(PF));
	printf("Testing: PFget_RelativePerm\n");
	printf("RelativePerm %f\n",PFget_RelativePerm(PF));
	printf("MovieFrames %d\n",PFget_MovieFrames(PF));
	printf("Vcv %d\n",PFget_Vcv(PF));
	printf("Tcv %d\n",PFget_Tcv(PF));

	ReadParameter(&method,\
			&SLength, &SWidth, &SHeight,\
			&PeriodicX, &PeriodicY, &PeriodicZ,\
			&EndX, &EndY, &EndZ,\
			&XElecOn, &YElecOn, &ZElecOn,\
			&XFermiB, &XFermiF, &YFermiL,\
			&YFermiR, &ZFermiB, &ZFermiA,\
			&alphaxb, &alphaxf, &alphayl,\
			&alphayr, &alphazb, &alphaza,\
			&vX, &vY, &vZ,\
			&VoltageX, &VoltageY, &VoltageZ,\
			&VStepX, &VStepY, &VStepZ,\
			&VincX, &VincY, &VincZ,\
			&SiteDistance, &D, &TCount,\
			&NCh, &Ntot, &TStep, &N_av,\
			&Nstep_av, &Nstep_check, &Rcount, &CutOff,\
			&lambda, &SeedProt, &Attempts,\
			&fracSeed, &E0, &sigma,\
			&fracTrap, &Etrap, &Tsigma,\
			&TempStart, &TemperatureStep,\
			&TemperatureInc, &reOrgEnergy,\
			&AttemptToHop, &gamma,\
			&RelativePerm, &MovieFrames,\
			&Tcv, &Vcv);
	printf("\nTesting: ReadParameter\n\n");

	printf("Method %d\n",method);
	printf("SLength %d\n",SLength);
	printf("SWidth %d\n",SWidth);
	printf("SHeight %d\n",SHeight);
	printf("PeriodicX %d\n",PeriodicX);
	printf("PeriodicY %d\n",PeriodicY);
	printf("PeriodicZ %d\n",PeriodicZ);
	printf("EndX %d\n",EndX);
	printf("EndY %d\n",EndY);
	printf("EndZ %d\n",EndZ);
	printf("XElecOn %d\n",XElecOn);
	printf("YElecOn %d\n",YElecOn);
	printf("ZElecOn %d\n",ZElecOn);
	printf("XFermiB %f\n",XFermiB);
	printf("XFermiF %f\n",XFermiF);
	printf("YFermiL %f\n",YFermiL);
	printf("YFermiR %f\n",YFermiR);
	printf("ZFermiB %f\n",ZFermiB);
	printf("ZFermiA %f\n",ZFermiA);
	printf("alphaxb %f\n",alphaxb);
	printf("alphaxf %f\n",alphaxf);
	printf("alphayl %f\n",alphayl);
	printf("alphayr %f\n",alphayr);
	printf("alphazb %f\n",alphazb);
	printf("alphaza %f\n",alphaza);
	printf("vX %g\n",vX);
	printf("vY %g\n",vY);
	printf("vZ %g\n",vZ);
	printf("VoltageX %g\n",VoltageX);
	printf("VoltageY %g\n",VoltageY);
	printf("VoltageZ %g\n",VoltageZ);
	printf("VStepX %d\n",VStepX);
	printf("VStepY %d\n",VStepY);
	printf("VStepZ %d\n",VStepZ);
	printf("VincX %g\n",VincX);
	printf("VincY %g\n",VincY);
	printf("VincZ %g\n",VincZ);
	printf("SiteDistance %g\n",SiteDistance);
	printf("D %g\n",D);
	printf("TCount %d\n",TCount);
	printf("NCh %d\n",NCh);
	printf("Ntot %d\n",Ntot);
	printf("TStep %g\n",TStep);
	printf("N_av %d\n",N_av);
	printf("Nstep_av %d\n",Nstep_av);
	printf("Nstep_check %d\n",Nstep_check);
	printf("Rcount %d\n",Rcount);
	printf("CutOff %g\n",CutOff);
	printf("lambda %g\n",lambda);
	printf("SeedProt %d\n",SeedProt);
	printf("Attempts %d\n",Attempts);
	printf("FracSeed %f\n",fracSeed);
	printf("E0 %f\n",E0);
	printf("sigma %f\n",sigma);
	printf("FracTrap %f\n",fracTrap);
	printf("Etrap %f\n",Etrap);
	printf("Tsigma %f\n",Tsigma);
	printf("TempStart %f\n",TempStart);
	printf("TempStep %d\n",TemperatureStep);
	printf("TempInc %f\n",TemperatureInc);
	printf("reOrg %f\n",reOrgEnergy);
	printf("AttemptToHop %g\n",AttemptToHop);
	printf("gamma %f\n",gamma);
	printf("RelativePerm %f\n",RelativePerm);
	printf("Vcv %f\n",Vcv);
	printf("Tcv %f\n",Tcv);

	deleteParamFrame(&PF);

	ParameterFrame PF2 = newParamFrame_File();
	printf("\nTesting: newParamFrame_File\n\n");
	printf("Method %d\n",PFget_method(PF2));
	printf("SLength %d\n",PFget_Len(PF2));
	printf("SWidth %d\n",PFget_Wid(PF2));
	printf("SHeight %d\n",PFget_Hei(PF2));
	printf("PeriodicX %d\n",PFget_Px(PF2));
	printf("PeriodicY %d\n",PFget_Py(PF2));
	printf("PeriodicZ %d\n",PFget_Pz(PF2));
	printf("EndX %d\n",PFget_EndX(PF2));
	printf("EndY %d\n",PFget_EndY(PF2));
	printf("EndZ %d\n",PFget_EndZ(PF2));
	printf("XElecOn %d\n",PFget_XElecOn(PF2));
	printf("YElecOn %d\n",PFget_YElecOn(PF2));
	printf("ZElecOn %d\n",PFget_ZElecOn(PF2));
	printf("XFermiB %f\n",PFget_XFermiB(PF2));
	printf("XFermiF %f\n",PFget_XFermiF(PF2));
	printf("YFermiL %f\n",PFget_YFermiL(PF2));
	printf("YFermiR %f\n",PFget_YFermiR(PF2));
	printf("ZFermiB %f\n",PFget_ZFermiB(PF2));
	printf("ZFermiA %f\n",PFget_ZFermiA(PF2));
	printf("alphaxb %f\n",PFget_alphaxb(PF2));
	printf("alphaxf %f\n",PFget_alphaxf(PF2));
	printf("alphayl %f\n",PFget_alphayl(PF2));
	printf("alphayr %f\n",PFget_alphayr(PF2));
	printf("alphazb %f\n",PFget_alphazb(PF2));
	printf("alphaza %f\n",PFget_alphaza(PF2));
	printf("vX %g\n",PFget_vX(PF2));
	printf("vY %g\n",PFget_vY(PF2));
	printf("vZ %g\n",PFget_vZ(PF2));
	printf("VoltageX %g\n",PFget_VoltageX(PF2));
	printf("VoltageY %g\n",PFget_VoltageY(PF2));
	printf("VoltageZ %g\n",PFget_VoltageZ(PF2));
	printf("VStepX %d\n",PFget_VStepX(PF2));
	printf("VStepY %d\n",PFget_VStepY(PF2));
	printf("VStepZ %d\n",PFget_VStepZ(PF2));
	printf("VincX %g\n",PFget_VincX(PF2));
	printf("VincY %g\n",PFget_VincY(PF2));
	printf("VincZ %g\n",PFget_VincZ(PF2));
	printf("SiteDistance %g\n",PFget_SiteDist(PF2));
	printf("D %g\n",PFget_D(PF2));
	printf("TCount %d\n",PFget_TCount(PF2));
	printf("NCh %d\n",PFget_NCh(PF2));
	printf("Ntot %d\n",PFget_Ntot(PF2));
	printf("TStep %g\n",PFget_TStep(PF2));
	printf("N_av %d\n",PFget_N_av(PF2));
	printf("Nstep_av %d\n",PFget_Nstep_av(PF2));
	printf("Nstep_check %d\n",PFget_Time_check(PF2));
	printf("Rcount %d\n",PFget_Rcount(PF2));
	printf("CutOff %g\n",PFget_CutOff(PF2));
	printf("lambda %g\n",PFget_lambda(PF2));
	printf("SeedProt %d\n",PFget_SeedProt(PF2));
	printf("Attempts %d\n",PFget_Attempts(PF2));
	printf("FracSeed %f\n",PFget_FracSeed(PF2));
	printf("E0 %f\n",PFget_E0(PF2));
	printf("sigma %f\n",PFget_sigma(PF2));
	printf("FracTrap %f\n",PFget_FracTrap(PF2));
	printf("Etrap %f\n",PFget_Etrap(PF2));
	printf("Tsigma %f\n",PFget_Tsigma(PF2));
	printf("TempStart %f\n",PFget_TempStart(PF2));
	printf("TempStep %d\n",PFget_TempStep(PF2));
	printf("TempInc %f\n",PFget_TempInc(PF2));
	printf("reOrg %f\n",PFget_reOrg(PF2));
	printf("AttemptToHop %g\n",PFget_AttemptToHop(PF2));
	printf("gamma %f\n",PFget_gamma(PF2));
	printf("RelativePerm %f\n",PFget_RelativePerm(PF2));
	printf("MovieFrames %d\n",PFget_MovieFrames(PF2));
	printf("Vcv %d\n",PFget_Vcv(PF2));
	printf("Tcv %d\n",PFget_Tcv(PF2));

	deleteParamFrame(&PF2);
	return 0;
}
