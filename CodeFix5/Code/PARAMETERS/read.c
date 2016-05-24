#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

//#include "../MEM/mem.h"
#include "read.h"

struct _ParameterFrame{
	
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
	int ImageCharge;
	double IntrFermi;
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
	int Time_check;
	int Rcount;
	int ClusterAlg;
	double CutOff;
	double lambda;
	int ScaleAfterCorr;
	int SeedProt;
	int Attempts;
	double fracSeed;
	double E0;
	double sigma;
	double fracTrap;
	double Etrap;
	double Tsigma;
	double TempStart;
	int TempStep;
	double TempInc;
	double reOrg;
	double AttemptToHop;
	double gamma;
	double RelativePerm;
	int MovieFrames;
	double CutOffTime;
	double Tcv;
	double Vcv;
	double Tlag;
	int EndPtFile;
	int NumChargesTrack;
	int PathFile;
	int LogFile;
};

int deleteParamFrame(ParameterFrame * PF){

	if(PF==NULL){
		return -1;
	}

	free(*PF);

	return 0;
}

ParameterFrame newParamFrame(void){

	ParameterFrame PF = (ParameterFrame) malloc(sizeof(struct _ParameterFrame));

	if(PF==NULL){
		printf("ERROR unable to create new parameterframe\n");
		return NULL;
	}

	//initialize all variables to 0
	PF->method=0;
	PF->SLength=0;
	PF->SWidth=0;
	PF->SHeight=0;
	PF->PeriodicX=0;
	PF->PeriodicY=0;
	PF->PeriodicZ=0;
	PF->EndX=0;
	PF->EndY=0;
	PF->EndZ=0;
	PF->XElecOn=0;
	PF->YElecOn=0;
	PF->ZElecOn=0;
	PF->XFermiB=0;
	PF->XFermiF=0;
	PF->YFermiL=0;
	PF->YFermiR=0;
	PF->ZFermiB=0;
	PF->ZFermiA=0;
	PF->ImageCharge=0;
	PF->IntrFermi=0;
	PF->alphaxb=0;
	PF->alphaxf=0;
	PF->alphayl=0;
	PF->alphayr=0;
	PF->alphazb=0;
	PF->alphaza=0;
	PF->vX=0;
	PF->vY=0;
	PF->vZ=0;
	PF->VoltageX=0;
	PF->VoltageY=0;
	PF->VoltageZ=0;
	PF->VStepX=0;
	PF->VStepY=0;
	PF->VStepZ=0;
	PF->VincX=0;
	PF->VincY=0;
	PF->VincZ=0;
	PF->SiteDistance=0;
	PF->D=0;
	PF->TCount=0;
	PF->NCh=0;
	PF->Ntot=0;
	PF->TStep=0;
	PF->N_av=0;
	PF->Nstep_av=0;
	PF->Time_check=0;
	PF->Rcount=0;
	PF->ClusterAlg=0;
	PF->CutOff=0;
	PF->lambda=0;
	PF->ScaleAfterCorr=0;
	PF->SeedProt=0;
	PF->Attempts=0;
	PF->fracSeed=0;
	PF->E0=0;
	PF->sigma=0;
	PF->fracTrap=0;
	PF->Etrap=0;
	PF->Tsigma=0;
	PF->TempStart=0;
	PF->TempStep=0;
	PF->TempInc=0;
	PF->reOrg=0;
	PF->AttemptToHop=0;
	PF->gamma=0;
	PF->RelativePerm=0;
	PF->MovieFrames=0;
	PF->CutOffTime=0;
	PF->Tcv=0;
	PF->Vcv=0;
	PF->Tlag=0;
	PF->EndPtFile=0;
	PF->NumChargesTrack=0;
	PF->PathFile=0;
	PF->LogFile=0;
	return PF;
}


ParameterFrame newParamFrame_File(void){

	ParameterFrame PF = (ParameterFrame) malloc(sizeof(struct _ParameterFrame));

	if(PF==NULL)
		return NULL;


	int intval;
	double doubleval;

	char *buffer = NULL;
	int check;
	unsigned int position;
	int string_size,read_size;
	FILE *handler = fopen("../PARAMETERS/parameters.txt","r");

	if (handler)
	{
		printf("Successfully opened ../PARAMETERS/parameters.txt\n");
		//seek the last byte of the file
		fseek(handler,0,SEEK_END);
		//offset from the first to the last byte, or in other words, filesize
		string_size = ftell (handler);
		//go back to the start of the file
		rewind(handler);

		//allocate a string that can hold it all
		buffer = (char*) malloc (sizeof(char) * (string_size + 1) );
		//read it all in one operation
		read_size = fread(buffer,sizeof(char),string_size,handler);
		//fread doesnt set it so put a \0 in the last position
		//and buffer is now officialy a string
		buffer[string_size] = '\0';

		if (string_size != read_size) {
			//something went wrong, throw away the memory and set
			//the buffer to NULL
			free(buffer);
			buffer = NULL;
			fclose(handler);
		}
	}

	if(buffer){
		puts(buffer);

		check = match(buffer, "\nmethod");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("method %d\n",intval);
			PF->method = intval;
		}else{
			printf("ERROR when reading file can not find method!\n");
			exit(1);
		}
		check = match(buffer, "\nSLength");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("SLength %d\n",intval);
			PF->SLength = intval;
			if(intval<=1){
				printf("ERROR SLength set to a value less than 2!\n");
				printf("SLength must be at least 2 SiteNodes wide.\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find SLength!\n");
			exit(1);
		}

		check = match(buffer, "\nSWidth");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("SWidth %d\n",intval);
			PF->SWidth = intval;
			if(intval<=1){
				printf("ERROR SWidth set to a value less than 2!\n");
				printf("SWidth must be at least 2 SiteNodes wide.\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find SWidth!\n");
			exit(1);
		}

		check = match(buffer, "\nSHeight");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("SHeight %d\n",intval);
			PF->SHeight = intval;
			if(intval<=1){
				printf("ERROR SHeight set to a value less than 2!\n");
				printf("SHeight must be at least 2 SiteNodes wide.\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find SHeight!\n");
			exit(1);
		}

		check = match(buffer, "\nPeriodicX");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("PeriodicX %d\n",intval);
			PF->PeriodicX = intval;
		}else{
			printf("ERROR when reading file can not find PeriodicX!\n");
			exit(1);
		}

		check = match(buffer, "\nPeriodicY");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("PeriodicY %d\n",intval);
			PF->PeriodicY = intval;
		}else{
			printf("ERROR when reading file can not find PeriodicY!\n");
			exit(1);
		}

		check = match(buffer, "\nPeriodicZ");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("PeriodicZ %d\n",intval);
			PF->PeriodicZ = intval;
		}else{
			printf("ERROR when reading file can not find PeriodicZ!\n");
			exit(1);
		}

		check = match(buffer, "\nEndX");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("EndX %d\n",intval);
			PF->EndX = intval;
		}else{
			printf("ERROR when reading file can not find EndX!\n");
			exit(1);
		}

		check = match(buffer, "\nEndY");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("EndY %d\n",intval);
			PF->EndY = intval;
		}else{
			printf("ERROR when reading file can not find EndY!\n");
			exit(1);
		}

		check = match(buffer, "\nEndZ");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("EndZ %d\n",intval);
			PF->EndZ = intval;
		}else{
			printf("ERROR when reading file can not find EndZ!\n");
			exit(1);
		}

		check = match(buffer, "\nXElecOn");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("XElecOn %d\n",intval);
			PF->XElecOn = intval;
		}else{
			printf("ERROR when reading file can not find XElecOn!\n");
			exit(1);
		}

		check = match(buffer, "\nYElecOn");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("YElecOn %d\n",intval);
			PF->YElecOn = intval;
		}else{
			printf("ERROR when reading file can not find YElecOn!\n");
			exit(1);
		}

		check = match(buffer, "\nZElecOn");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("ZElecOn %d\n",intval);
			PF->ZElecOn = intval;
		}else{
			printf("ERROR when reading file can not find ZElecOn!\n");
			exit(1);
		}

		check = match(buffer, "\nXFermiB");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("XFermiB %g\n",doubleval);
			PF->XFermiB = doubleval;
		}else{
			printf("ERROR when reading file can not find XFermiB!\n");
			exit(1);
		}

		check = match(buffer, "\nXFermiF");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("XFermiF %g\n",doubleval);
			PF->XFermiF = doubleval;
		}else{
			printf("ERROR when reading file can not find XFermiF!\n");
			exit(1);
		}

		check = match(buffer, "\nYFermiL");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("YFermiL %g\n",doubleval);
			PF->YFermiL = doubleval;
		}else{
			printf("ERROR when reading file can not find YFermiL!\n");
			exit(1);
		}

		check = match(buffer, "\nYFermiR");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("YFermiR %g\n",doubleval);
			PF->YFermiR = doubleval;
		}else{
			printf("ERROR when reading file can not find YFermiR!\n");
			exit(1);
		}

		check = match(buffer, "\nZFermiB");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("ZFermiB %g\n",doubleval);
			PF->ZFermiB = doubleval;
		}else{
			printf("ERROR when reading file can not find ZFermiB!\n");
			exit(1);
		}

		check = match(buffer, "\nZFermiA");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("ZFermiA %g\n",doubleval);
			PF->ZFermiA = doubleval;
		}else{
			printf("ERROR when reading file can not find ZFermiA!\n");
			exit(1);
		}
		
		check = match(buffer, "\nImageCharge");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("ImageCharge %d\n",intval);
			PF->ImageCharge = intval;
			if (intval<0 || intval>1){
				printf("ERROR ImageCharge is set either above 1 or\n");
				printf("less than 0. It can either be 1 for on or \n");
				printf("0 for off\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find ImageCharge!\n");
			exit(1);
		}

		check = match(buffer, "\nIntrinsicFermi");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("IntrinsicFermi %g\n",doubleval);
			PF->IntrFermi = doubleval;
		}else{
			if(PF->method==1){
				printf("ERROR when reading file can not find IntrinsicFermi!\n");
				printf("This variable is needed to determine the number of\n");
				printf("thermally activaited carriers in the material when\n");
				printf("using the CELIV method\n");
				exit(1);
			}
		}

		//Alpha parameters must be converted from [1/nm] to [1/m]
		check = match(buffer, "\nalphaXb");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("alphaxb %g\n",doubleval);
			PF->alphaxb = doubleval*1E9;
		}else{
			printf("ERROR when reading file can not find alphaxb!\n");
			exit(1);
		}

		check = match(buffer, "\nalphaXf");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("alphaxf %g\n",doubleval);
			PF->alphaxf = doubleval*1E9;
		}else{
			printf("ERROR when reading file can not find alphaxf!\n");
			exit(1);
		}

		check = match(buffer, "\nalphaYl");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("alphayl %g\n",doubleval);
			PF->alphayl = doubleval*1E9;
		}else{
			printf("ERROR when reading file can not find alphayl!\n");
			exit(1);
		}

		check = match(buffer, "\nalphaYr");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("alphayr %g\n",doubleval);
			PF->alphayr = doubleval*1E9;
		}else{
			printf("ERROR when reading file can not find alphayr!\n");
			exit(1);
		}

		check = match(buffer, "\nalphaZb");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("alphazb %g\n",doubleval );
			PF->alphazb = doubleval*1E9; 
		}else{
			printf("ERROR when reading file can not find alphazb!\n");
			exit(1);
		}

		check = match(buffer, "\nalphaZa");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("alphaza %g\n",doubleval );
			PF->alphaza = doubleval*1E9; 
		}else{
			printf("ERROR when reading file can not find alphaza!\n");
			exit(1);
		}

		check = match(buffer, "\nRelativePermittivity");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("RelativePermittivity %g\n",doubleval );
			PF->RelativePerm = doubleval; 
		}else{
			printf("ERROR when reading file can not find RelativePermittivity!\n");
			exit(1);
		}

		check = match(buffer, "\nvX");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("vX %g\n",doubleval );
			PF->vX = doubleval; 
		}else{
			printf("ERROR when reading file can not find vX!\n");
			exit(1);
		}

		check = match(buffer, "\nvY");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("vY %g\n",doubleval );
			PF->vY = doubleval; 
		}else{
			printf("ERROR when reading file can not find vY!\n");
			exit(1);
		}

		check = match(buffer, "\nvZ");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("vZ %g\n",doubleval );
			PF->vZ = doubleval; 
		}else{
			printf("ERROR when reading file can not find vZ!\n");
			exit(1);
		}

		check = match(buffer, "\nVoltageX");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("VoltageX %g\n",doubleval );
			PF->VoltageX = doubleval; 

			if(PF->method == 1 && PF->VoltageX!=0){
				printf("ERROR method is CELIV but initial VoltageX is not 0\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find VoltageX!\n");
			exit(1);
		}

		check = match(buffer, "\nVoltageY");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval  = GrabDouble(position, &buffer[0] );
			printf("VoltageY %g\n",doubleval );
			PF->VoltageY =doubleval;
			
			if(PF->method == 1 && PF->VoltageY!=0){
				printf("ERROR method is CELIV but initial VoltageY is not 0\n");
				exit(1);
			} 
		}else{
			printf("ERROR when reading file can not find VoltageY!\n");
			exit(1);
		}

		check = match(buffer, "\nVoltageZ");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("VoltageZ %g\n",doubleval);
			PF->VoltageZ =doubleval;
			
			if(PF->method == 1 && PF->VoltageZ!=0){
				printf("ERROR method is CELIV but initial VoltageZ is not 0\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find VoltageZ!\n");
			exit(1);
		}

		check = match(buffer, "\nVStepX");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("VStepX %d\n",intval);
			PF->VStepX = intval;

			if(PF->method==1 && PF->VStepX!=0){
				printf("ERROR method is CELIV but VStepX is non 0\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find VStepX!\n");
			exit(1);
		}

		check = match(buffer, "\nVStepY");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("VStepY %d\n",intval);
			PF->VStepY = intval;
			
			if(PF->method==1 && PF->VStepY!=0){
				printf("ERROR method is CELIV but VStepY is non 0\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find VStepY!\n");
			exit(1);
		}

		check = match(buffer, "\nVStepZ");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("VStepZ %d\n",intval);
			PF->VStepZ = intval;
			
			if(PF->method==1 && PF->VStepZ!=0){
				printf("ERROR method is CELIV but VStepZ is non 0\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find VStepZ!\n");
			exit(1);
		}

		check = match(buffer, "\nVincX");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("VincX %g\n",doubleval);
			PF->VincX = doubleval;

		}else{
			printf("ERROR when reading file can not find VincX!\n");
			exit(1);
		}

		check = match(buffer, "\nVincY");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("VincY %g\n",doubleval);
			PF->VincY =doubleval;
			
		}else{
			printf("ERROR when reading file can not find VincY!\n");
			exit(1);
		}

		check = match(buffer, "\nVincZ");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("VincZ %g\n",doubleval);
			PF->VincZ = doubleval;
			
		}else{
			printf("ERROR when reading file can not find VincZ!\n");
			exit(1);
		}

		check = match(buffer, "\nSiteDistance");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("SiteDistance %g\n",doubleval);
			if(doubleval<0){
				printf("ERROR SiteDistance is negative\n");
				exit(1);
			}

			PF->SiteDistance = doubleval;
		}else{
			printf("ERROR when reading file can not find SiteDistance!\n");
			exit(1);
		}

		check = match(buffer, "\nD");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("D %g\n",doubleval);
			PF->D = doubleval;
		}else{
			printf("ERROR when reading file can not find D!\n");
			exit(1);
		}

		check = match(buffer, "\nTCount");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("TCount %d\n",intval);
			if(intval<0){
				printf("ERROR TCount is negative\n");
				exit(1);
			}
			PF->TCount = intval;
		}else{
			printf("ERROR when reading file can not find TCount!\n");
			exit(1);
		}

		check = match(buffer, "\nNCh");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("NCh %d\n",intval);
			if(intval<0){
				printf("ERROR NCh is negative\n");
				exit(1);
			}
			PF->NCh = intval;
			if(PF->method==1){
				if(PF->NCh>0){
					printf("ERROR number of charges should be set to 0\n");
					printf("in CELIV method because charges are generated\n");
					printf("thermally.\n");
					exit(1);
				}
			}

		}else{
			printf("ERROR when reading file can not find NCh!\n");
			exit(1);
		}

		PF->Ntot = (PF->TCount) * (PF->NCh);
		printf("Ntot %d\n",PF->Ntot);

		check = match(buffer, "\nTStep");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("TStep %g\n",doubleval);
			if(doubleval<0){
				printf("ERROR TStep is negative\n");
				exit(1);
			}
			PF->TStep = doubleval;
		}else{
			printf("ERROR when reading file can not find TStep!\n");
			exit(1);
		}

		check = match(buffer, "\nNstep_av");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("Nstep_av %d\n",intval);
			if(intval<0){
				printf("ERROR Nstep_av is negative\n");
				exit(1);
			}
			PF->Nstep_av = intval;
		}else{
			printf("ERROR when reading file can not find Nstep_av!\n");
			exit(1);
		}

		//Number of Averages recorded
		intval = (int) PF->TCount/ PF->Nstep_av;
		printf("N_av %d\n",intval);
		PF->N_av = intval;

		//Number of steps before checkpoint file
		check = match(buffer, "\nTime_check");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("Time_check %d\n",intval);
			if(intval<0){
				printf("ERROR Time_check is negative\n");
				exit(1);
			}
			PF->Time_check = intval;
		}else{
			printf("ERROR when reading file can not find Time_check!\n");
			exit(1);
		}

		check = match(buffer, "\nRcount");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("Rcount %d\n",intval);
			PF->Rcount = intval;
		}else{
			printf("ERROR when reading file can not find Rcount!\n");
			exit(1);
		}

		check = match(buffer, "\nClusterAlg");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("ClusterAlg %d\n",intval);
			PF->ClusterAlg = intval;
		}else{
			printf("ERROR when reading file can not find Rcount!\n");
			exit(1);
		}


		check = match(buffer, "\nCutOff");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("CutOff %g\n",doubleval);
			if(doubleval<0){
				printf("ERROR CutOff is negative\n");
				exit(1);
			}
			PF->CutOff = doubleval;
		}else{
			printf("ERROR when reading file can not find CutOff!\n");
			exit(1);
		}

		check = match(buffer, "\nlambda");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("lambda %g\n",doubleval);
			PF->lambda = doubleval;
		}else{
			printf("ERROR when reading file can not find lambda!\n");
			exit(1);
		}

		check = match(buffer, "\nScaleAfterCorr");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("ScaleAfterCorr %d\n",intval);
			if(intval<0){
				printf("ERROR ScaleAfterCorr is negative can only be 0-off or 1-on\n");
				exit(1);
			}else if(intval>1){
				printf("ERROR ScaleAfterCorr is greater than 1 can only be 0-off or 1-on\n");
				exit(1);
			}
			PF->ScaleAfterCorr = intval;
		}else{
			printf("ERROR when reading file can not find ScaleAfterCorr!\n");
			exit(1);
		}


		check = match(buffer, "\nSeedProt");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("SeedProt %d\n",intval);
			if(intval<0){
				printf("ERROR SeedProt is negative\n");
				exit(1);
			}else if(intval>2){
				printf("ERROR SeedProt is greater than 2\n");
				exit(1);
			}
			PF->SeedProt = intval;
		}else{
			printf("ERROR when reading file can not find SeedProt!\n");
			exit(1);
		}

		check = match(buffer, "\nAttempts");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("Attempts %d\n",intval);
			if(intval<0){
				printf("ERROR Attempts is negative\n");
				exit(1);
			}
			PF->Attempts = intval;
		}else{
			printf("ERROR when reading file can not find Attempts!\n");
			exit(1);
		}

		check = match(buffer, "\nfracSeed");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("fracSeed %g\n",doubleval);
			if(doubleval>1){
				printf("ERROR Seed fraction is greater than 1\n");
				exit(1);
			}else if(doubleval<0){
				printf("ERROR Seed fraction is negative\n");
				exit(1);
			}
			PF->fracSeed = doubleval;
		}else{
			printf("ERROR when reading file can not find fracSeed!\n");
			exit(1);
		}

		check = match(buffer, "\nE0");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("E0 %g\n",doubleval);
			PF->E0 = doubleval;
			if(PF->method==1){
				if(PF->E0<PF->IntrFermi){
					printf("ERROR mean of gaussian distribution for DOS\n");
					printf("is less than IntrFermi\n");
					exit(1);
				}
			}
		}else{
			printf("ERROR when reading file can not find E0!\n");
			exit(1);
		}

		check = match(buffer, "\nsigma");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("sigma %g\n",doubleval);
			PF->sigma = doubleval;
		}else{
			printf("ERROR when reading file can not find sigma!\n");
			exit(1);
		}

		check = match(buffer, "\nfracTrap");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("fracTrap %g\n",doubleval);
			PF->fracTrap = doubleval;
			if(doubleval>1){
				printf("ERROR Fraction of traps is greater than 1\n");
				exit(1);
			}else if(doubleval<0){
				printf("ERROR Fraction of traps is negative\n");
				exit(1);
			}else if(PF->ScaleAfterCorr!=0 && doubleval!=0){
				printf("ERROR Fraction of traps is non zero and ScaleAfterCorr is set to on\n");
				printf("At the moment the two options are not compatible\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find fracTrap!\n");
			exit(1);
		}

		check = match(buffer, "\nEtrap");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("Etrap %g\n",doubleval);
			PF->Etrap = doubleval;
			
			if(PF->fracTrap!=0){
				if(PF->method==1){
					if(PF->Etrap<(PF->IntrFermi)){
						printf("ERROR Etrap is less than intrinsic fermi level\n");
						exit(1);
					}
				}
			}
		}else{
			printf("ERROR when reading file can not find Etrap!\n");
			exit(1);
		}

		check = match(buffer, "\nTsigma");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("Tsigma %g\n",doubleval);
			PF->Tsigma = doubleval;
		}else{
			printf("ERROR when reading file can not find Tsigma!\n");
			exit(1);
		}

		check = match(buffer, "\nTempStart");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("TempStart %g\n",doubleval);
			PF->TempStart = doubleval;
			if(doubleval<0){
				printf("ERROR Starting temperature is negative\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find TempStart!\n");
			exit(1);
		}

		check = match(buffer, "\nTemperatureStep");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabInt(position, &buffer[0] );
			printf("TemperatureStep %g\n",doubleval);
			PF->TempStep = doubleval;
			if(doubleval<0){
				printf("ERROR Number of temperature steps is negative\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find TemperatureStep!\n");
			exit(1);
		}

		check = match(buffer, "\nTemperatureInc");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("TemperatureInc %g\n",doubleval);
			PF->TempInc =doubleval;
		}else{
			printf("ERROR when reading file can not find TemperatureInc!\n");
			exit(1);
		}

		check = match(buffer, "\nreOrgEnergy");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("reOrgEnergy %g\n",doubleval);
			PF->reOrg = doubleval;
		}else{
			printf("ERROR when reading file can not find reOrgEnergy!\n");
			exit(1);
		}

		check = match(buffer, "\nAttemptToHop");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("AttemptToHop %g\n",doubleval);
			if(doubleval<0){
				printf("ERROR Attempt to hop rate is negative\n");
				exit(1);
			}
			PF->AttemptToHop = doubleval;
		}else{
			printf("ERROR when reading file can not find AttemptToHop!\n");
			exit(1);
		}

		check = match(buffer, "\ngamma");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("gamma %g\n",doubleval);
			PF->gamma = doubleval*1E9;
		}else{
			printf("ERROR when reading file can not find gamma!\n");
			exit(1);
		}

		check = match(buffer,"\nMovieFrames");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0] );
			printf("MovieFrames %d\n",intval);
			PF->MovieFrames = intval;
		}else{
			printf("ERROR when reading file can not find MovieFrames!\n");
			exit(1);
		}

		check = match(buffer,"\nCutOffTime");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("CutOffTime %g\n",doubleval);
			PF->CutOffTime = doubleval;
			
			if(PF->CutOffTime<0){
				printf("ERROR CutOffTime assigned a value less than 0\n");
				printf("CutOffTime must be positive or 0");
				exit(1);
			}
			if((PF->CutOffTime)<(PF->TStep*PF->TCount)){
				printf("ERROR you have setup the simulation to end before\n");
				printf("all the charges have been inserted\n");
				printf("TCount*TStep is greater than CutOffTime\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find CutOffTime!\n");
			exit(1);
		}
		
		check = match(buffer,"\nVcv");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("Vcv %g\n",doubleval);
			PF->Vcv = doubleval;

			if(PF->method==1 && PF->Vcv<=0){
				printf("ERROR CELIV voltage ramp voltage set less or equal to 0\n");
				exit(1);
			}else if(PF->method==0 && PF->Vcv!=0){
				printf("ERROR Time of flight method specified but CELIV ramp\n");
				printf("voltage is non 0\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find Vcv\n");
			exit(1);
		}
		
		check = match(buffer,"\nTcv");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("Tcv %g\n",doubleval);
			PF->Tcv = doubleval;

			if(PF->method==1 && PF->Tcv<=0){
				printf("ERROR CELIV voltage ramp time set less or equal to 0\n");
				exit(1);
			}else if(PF->method==0 && PF->Tcv!=0){
				printf("ERROR Time of flight method specified but CELIV ramp\n");
				printf("time is non 0\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find Tcv!\n");
			exit(1);
		}
			
		check = match(buffer,"\nTlag");
		if(check!=-1){
			position = (unsigned int)check;
			doubleval = GrabDouble(position, &buffer[0] );
			printf("Tlag %g\n",doubleval);
			PF->Tlag = doubleval;
			if(PF->method==1 && PF->Tlag<0){
				printf("ERROR CELIV voltage ramp lag time set less than 0\n");
				exit(1);
			}else if(PF->method==0 && PF->Tlag!=0){
				printf("ERROR Time of flight method specified but CELIV ramp\n");
				printf("lag time is non 0\n");
				exit(1);
			}

		}else{
			printf("ERROR when reading file can not find Tlag!\n");
			exit(1);
		}			

		check = match(buffer,"\nEndPtFile");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0]);
			printf("EndPtFile %d\n",intval);
			PF->EndPtFile = intval;
			if(PF->EndPtFile<0 || PF->EndPtFile>1){
				printf("ERROR EndPtFile can only be set to 0 for off\n");
				printf("or set to 1 for on\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find EndPtFile!\n");
			exit(1);
		}

		check = match(buffer,"\nNumChargesTrack");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0]);
			printf("NumChargesTrack %d\n",intval);
			PF->NumChargesTrack = intval;
			if((PF->NumChargesTrack)<0 ){
				printf("ERROR NumChargesTrack can only be set to 0 or a positive integer\n");
				exit(1);
			}
		
		}else{
			printf("ERROR when reading file can not find NumChargesTrack!\n");
			exit(1);
		}

		check = match(buffer,"\nPathFile");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0]);
			printf("PathFile %d\n",intval);
			PF->PathFile = intval;
			if(PF->PathFile<0 || PF->PathFile>1){
				printf("ERROR PathFile can only be set to 0 for off\n");
				printf("or set to 1 for on\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find PathFile!\n");
			exit(1);
		}
		
		check = match(buffer,"\nLogFile");
		if(check!=-1){
			position = (unsigned int)check;
			intval = GrabInt(position, &buffer[0]);
			printf("LogFile %d\n",intval);
			PF->LogFile = intval;
			if(PF->LogFile<0 || PF->LogFile>1){
				printf("ERROR LogFile can only be set to 0 for off\n");
				printf("or set to 1 for on\n");
				exit(1);
			}
		}else{
			printf("ERROR when reading file can not find LogFile!\n");
			exit(1);
		}
		free(buffer);
		fclose(handler);
		return PF;
	}

	return NULL;
}

int ReadParameter(int * method,\
		int * SLength, int * SWidth, int * SHeight,\
		int * PeriodicX, int * PeriodicY, int * PeriodicZ,\
		int * EndX, int * EndY, int * EndZ,\
		int * XElecOn, int * YElecOn, int *ZElecOn,\
		double * XFermiB, double * XFermiF, double * YFermiL,\
		double * YFermiR, double * ZFermiB, double * ZFermiA,\
		int * ImageCharge, double * IntrFermi,\
		double * alphaxb, double * alphaxf, double * alphayl,\
		double * alphayr, double * alphazb, double * alphaza,\
		double * vX, double * vY, double * vZ,\
		double * VoltageX, double * VoltageY, double * VoltageZ,\
		int * VStepX, int * VStepY, int * VStepZ,\
		double * VincX, double * VincY, double * VincZ,\
		double * SiteDistance, double * D, int * TCount,\
		int * NCh, int * Ntot, double * TStep, int * N_av,\
		int * Nstep_av, int * Time_check, int * Rcount,int * ClusterAlg, double * CutOff,\
		double * lambda, int * ScaleAfterCorr,int * SeedProt, int * Attempts,\
		double * fracSeed, double * E0, double * sigma,\
		double * fracTrap, double * Etrap, double * Tsigma,\
		double * TempStart, int * TemperatureStep,\
		double * TemperatureInc, double * reOrgEnergy,\
		double * AttemptToHop, double * gamma,\
		double * RelativePerm, int * MovieFrames, double * CutOffTime,\
		double * Vcv, double * Tcv, double * Tlag, int * EndPtFile,\
		int * NumChargesTrack, int * PathFile, int * LogFile){

			char *buffer = NULL;
			int position;
			int string_size,read_size;
			FILE *handler = fopen("parameters.txt","r");
			FILE *handler2 = fopen("../PARAMETERS/parameters.txt","r");

			if(handler==NULL){
			   handler=handler2;
			}else{
			   fclose(handler2);
			}

			if (handler)
			{
				//seek the last byte of the file
				fseek(handler,0,SEEK_END);
				//offset from the first to the last byte, or in other words, filesize
				string_size = ftell (handler);
				//go back to the start of the file
				rewind(handler);

				//allocate a string that can hold it all
				buffer = (char*) malloc (sizeof(char) * (string_size + 1) );
				//read it all in one operation
				read_size = fread(buffer,sizeof(char),string_size,handler);
				//fread doesnt set it so put a \0 in the last position
				//and buffer is now officialy a string
				buffer[string_size] = '\0';

				if (string_size != read_size) {
					//something went wrong, throw away the memory and set
					//the buffer to NULL
					free(buffer);
					buffer = NULL;
				}
			}

			if(buffer!=NULL){
				puts(buffer);
					
				position = match(buffer, "\nmethod");
				*method = GrabInt(position, &buffer[0] );
				printf("method %d\n",*method);

				position = match(buffer, "\nSLength");
				*SLength = GrabInt(position, &buffer[0] );
				printf("SLength %d\n",*SLength);
				if(*SLength<=1){
					printf("ERROR SLength set to a value less than 2!\n");
					printf("SLength must be at least 2 SiteNodes wide.\n");
					exit(1);
				}

				position = match(buffer, "\nSWidth");
				*SWidth = GrabInt(position, &buffer[0] );
				printf("SWidth %d\n",*SWidth);
				if(*SWidth<=1){
					printf("ERROR SWidth set to a value less than 2!\n");
					printf("SWidth must be at least 2 SiteNodes wide.\n");
					exit(1);
				}

				position = match(buffer, "\nSHeight");
				*SHeight = GrabInt(position, &buffer[0] );
				printf("SHeight %d\n",*SHeight);
				if(*SHeight<=1){
					printf("ERROR SHeight set to a value less than 2!\n");
					printf("SHeight must be at least 2 SiteNodes wide.\n");
					exit(1);
				}

				position = match(buffer, "\nPeriodicX");
				*PeriodicX = GrabInt(position, &buffer[0] );
				printf("PeriodicX %d\n",*PeriodicX);

				position = match(buffer, "\nPeriodicY");
				*PeriodicY = GrabInt(position, &buffer[0] );
				printf("PeriodicY %d\n",*PeriodicY);

				position = match(buffer, "\nPeriodicZ");
				*PeriodicZ = GrabInt(position, &buffer[0] );
				printf("PeriodicZ %d\n",*PeriodicZ);

				position = match(buffer, "\nEndX");
				*EndX = GrabInt(position, &buffer[0] );
				printf("EndX %d\n",*EndX);

				position = match(buffer, "\nEndY");
				*EndY = GrabInt(position, &buffer[0] );
				printf("EndY %d\n",*EndY);

				position = match(buffer, "\nEndZ");
				*EndZ = GrabInt(position, &buffer[0] );
				printf("EndZ %d\n",*EndZ);

				position = match(buffer, "\nXElecOn");
				*XElecOn = GrabInt(position, &buffer[0] );
				printf("XElecOn %d\n",*XElecOn);

				position = match(buffer, "\nYElecOn");
				*YElecOn = GrabInt(position, &buffer[0] );
				printf("YElecOn %d\n",*YElecOn);

				position = match(buffer, "\nZElecOn");
				*ZElecOn = GrabInt(position, &buffer[0] );
				printf("ZElecOn %d\n",*ZElecOn);

				position = match(buffer, "\nXFermiB");
				*XFermiB = GrabDouble(position, &buffer[0] );
				printf("XFermiB %g\n",*XFermiB);

				position = match(buffer, "\nXFermiF");
				*XFermiF = GrabDouble(position, &buffer[0] );
				printf("XFermiF %g\n",*XFermiF);

				position = match(buffer, "\nYFermiL");
				*YFermiL = GrabDouble(position, &buffer[0] );
				printf("YFermiL %g\n",*YFermiL);

				position = match(buffer, "\nYFermiR");
				*YFermiR = GrabDouble(position, &buffer[0] );
				printf("YFermiR %g\n",*YFermiR);

				position = match(buffer, "\nZFermiB");
				*ZFermiB = GrabDouble(position, &buffer[0] );
				printf("ZFermiB %g\n",*ZFermiB);

				position = match(buffer, "\nZFermiA");
				*ZFermiA = GrabDouble(position, &buffer[0] );
				printf("ZFermiA %g\n",*ZFermiA);

				position = match(buffer, "\nImageCharge");
				*ImageCharge = GrabDouble(position, &buffer[0] );
				printf("ImageCharge %d\n",*ImageCharge);

				position = match(buffer, "\nIntrinsicFermi");
				*IntrFermi = GrabDouble(position, &buffer[0] );
				printf("IntrinsicFermi %g\n",*IntrFermi);

				position = match(buffer, "\nalphaXb");
				position = match(buffer, "\nalphaXb");
				*alphaxb = GrabDouble(position, &buffer[0] );
				printf("alphaXb %g\n",*alphaxb);

				position = match(buffer, "\nalphaXf");
				*alphaxf = GrabDouble(position, &buffer[0] );
				printf("alphaXf %g\n",*alphaxf);

				position = match(buffer, "\nalphaYl");
				*alphayl = GrabDouble(position, &buffer[0] );
				printf("alphaYl %g\n",*alphayl);

				position = match(buffer, "\nalphaYr");
				*alphayr = GrabDouble(position, &buffer[0] );
				printf("alphaYr %g\n",*alphayr);

				position = match(buffer, "\nalphaZb");
				*alphazb = GrabDouble(position, &buffer[0] );
				printf("alphaZb %g\n",*alphazb);

				position = match(buffer, "\nalphaZa");
				*alphaza = GrabDouble(position, &buffer[0] );
				printf("alphaZa %g\n",*alphaza);

				position = match(buffer, "\nRelativePermittivity");
				*RelativePerm = GrabDouble(position, &buffer[0] );
				printf("RelativePermittivity %g\n",*RelativePerm);

				position = match(buffer, "\nvX");
				*vX = GrabDouble(position, &buffer[0] );
				printf("vX %g\n",*vX);

				position = match(buffer, "\nvY");
				*vY = GrabDouble(position, &buffer[0] );
				printf("vY %g\n",*vY);

				position = match(buffer, "\nvZ");
				*vZ = GrabDouble(position, &buffer[0] );
				printf("vZ %g\n",*vZ);

				position = match(buffer, "\nVoltageX");
				*VoltageX = GrabDouble(position, &buffer[0] );
				printf("VoltageX %g\n",*VoltageX);
				if(*method == 1 && *VoltageX!=0){
					printf("ERROR method is CELIV but initial VoltageX is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVoltageY");
				*VoltageY = GrabDouble(position, &buffer[0] );
				printf("VoltageY %g\n",*VoltageY);
				if(*method == 1 && *VoltageY!=0){
					printf("ERROR method is CELIV but initial VoltageY is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVoltageZ");
				*VoltageZ = GrabDouble(position, &buffer[0] );
				printf("VoltageZ %g\n",*VoltageZ);
				if(*method == 1 && *VoltageZ!=0){
					printf("ERROR method is CELIV but initial VoltageZ is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVStepX");
				*VStepX = GrabInt(position, &buffer[0] );
				printf("VStepX %d\n",*VStepX);
				if(*method == 1 && *VStepX!=0){
					printf("ERROR method is CELIV but initial VStepX is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVStepY");
				*VStepY = GrabInt(position, &buffer[0] );
				printf("VStepY %d\n",*VStepY);
				if(*method == 1 && *VStepY!=0){
					printf("ERROR method is CELIV but initial VStepY is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVStepZ");
				*VStepZ = GrabInt(position, &buffer[0] );
				printf("VStepZ %d\n",*VStepZ);
				if(*method == 1 && *VStepZ!=0){
					printf("ERROR method is CELIV but initial VStepZ is not 0\n");
					exit(1);
				}

				position = match(buffer, "\nVincX");
				*VincX = GrabDouble(position, &buffer[0] );
				printf("VincX %g\n",*VincX);

				position = match(buffer, "\nVincY");
				*VincY = GrabDouble(position, &buffer[0] );
				printf("VincY %g\n",*VincY);

				position = match(buffer, "\nVincZ");
				*VincZ = GrabDouble(position, &buffer[0] );
				printf("VincZ %g\n",*VincZ);

				position = match(buffer, "\nSiteDistance");
				*SiteDistance = GrabDouble(position, &buffer[0] );
				printf("SiteDistance %g\n",*SiteDistance);

				position = match(buffer, "\nD");
				*D = GrabDouble(position, &buffer[0] );
				printf("D %g\n",*D);

				position = match(buffer, "\nTCount");
				*TCount = GrabInt(position, &buffer[0] );
				printf("TCount %d\n",*TCount);

				position = match(buffer, "\nNCh");
				*NCh = GrabInt(position, &buffer[0] );
				printf("NCh %d\n",*NCh);

				*Ntot = (*TCount) * (*NCh);
				printf("Ntot %d\n",*Ntot);

				position = match(buffer, "\nTStep");
				*TStep = GrabDouble(position, &buffer[0] );
				printf("TStep %g\n",*TStep);

				position = match(buffer, "\nNstep_av");
				*Nstep_av = GrabInt(position, &buffer[0] );
				printf("Nstep_av %d\n",*Nstep_av);

				//Number of Averages recorded
				*N_av = (int) *TCount/ *Nstep_av;
				printf("N_av %d\n",*N_av);

				//Number of steps before checkpoint file
				position = match(buffer, "\nTime_check");
				*Time_check = GrabInt(position, &buffer[0] );
				printf("Time_check %d\n",*Time_check);

				position = match(buffer, "\nRcount");
				*Rcount = GrabInt(position, &buffer[0] );
				printf("Rcount %d\n",*Rcount);

				position = match(buffer, "\nClusterAlg");
				*ClusterAlg = GrabInt(position, &buffer[0] );
				printf("ClusterAlg %d\n",*ClusterAlg);

				position = match(buffer, "\nCutOff");
				*CutOff = GrabDouble(position, &buffer[0] );
				printf("CutOff %g\n",*CutOff);

				position = match(buffer, "\nlambda");
				*lambda = GrabDouble(position, &buffer[0] );
				printf("lambda %g\n",*lambda);

				position = match(buffer, "\nScaleAfterCorr");
				*ScaleAfterCorr = GrabInt(position, &buffer[0] );
				printf("ScaleAfterCorr %d\n",*ScaleAfterCorr);

				position = match(buffer, "\nSeedProt");
				*SeedProt = GrabInt(position, &buffer[0] );
				printf("SeedProt %d\n",*SeedProt);

				position = match(buffer, "\nAttempts");
				*Attempts = GrabInt(position, &buffer[0] );
				printf("Attempts %d\n",*Attempts);

				position = match(buffer, "\nfracSeed");
				*fracSeed = GrabDouble(position, &buffer[0] );
				printf("fracSeed %g\n",*fracSeed);

				position = match(buffer, "\nE0");
				*E0 = GrabDouble(position, &buffer[0] );
				printf("E0 %g\n",*E0);

				position = match(buffer, "\nsigma");
				*sigma = GrabDouble(position, &buffer[0] );
				printf("sigma %g\n",*sigma);

				position = match(buffer, "\nfracTrap");
				*fracTrap = GrabDouble(position, &buffer[0] );
				printf("fracTrap %g\n",*fracTrap);

				position = match(buffer, "\nEtrap");
				*Etrap = GrabDouble(position, &buffer[0] );
				printf("Etrap %g\n",*Etrap);

				position = match(buffer, "\nTsigma");
				*Tsigma = GrabDouble(position, &buffer[0] );
				printf("Tsigma %g\n",*Tsigma);

				position = match(buffer, "\nTempStart");
				*TempStart = GrabDouble(position, &buffer[0] );
				printf("TempStart %g\n",*TempStart);

				position = match(buffer, "\nTemperatureStep");
				*TemperatureStep = GrabInt(position, &buffer[0] );
				printf("TemperatureStep %d\n",*TemperatureStep);

				position = match(buffer, "\nTemperatureInc");
				*TemperatureInc = GrabDouble(position, &buffer[0] );
				printf("TemperatureInc %g\n",*TemperatureInc);

				position = match(buffer, "\nreOrgEnergy");
				*reOrgEnergy = GrabDouble(position, &buffer[0] );
				printf("reOrgEnergy %g\n",*reOrgEnergy);

				position = match(buffer, "\nAttemptToHop");
				*AttemptToHop = GrabDouble(position, &buffer[0] );
				printf("AttemptToHop %g\n",*AttemptToHop);

				position = match(buffer, "\ngamma");
				*gamma = GrabDouble(position, &buffer[0] );
				printf("gamma %g\n",*gamma);

				position = match(buffer,"\nMovieFrames");
				*MovieFrames = GrabInt(position, & buffer[0]);
				printf("MovieFrames %d\n",*MovieFrames);

				position = match(buffer, "\nCutOffTime");
				*CutOffTime = GrabDouble(position, &buffer[0]);
				printf("CutOffTime %g\n",*CutOffTime);

				position = match(buffer, "\nTcv");
				*Tcv = GrabDouble(position, &buffer[0] );
				printf("Tcv %f\n",*Tcv);

				position = match(buffer, "\nVcv");
				*Vcv = GrabDouble(position, &buffer[0] );
				printf("Vcv %f\n",*Vcv);

				position = match(buffer, "\nTlag");
				*Tlag = GrabDouble(position, &buffer[0] );
				printf("Tlag %f\n",*Tcv);

				position = match(buffer, "\nEndPtFile");
				*EndPtFile = GrabInt(position, &buffer[0] );
				printf("EndPtFile %d\n",*EndPtFile);

				position = match(buffer, "\nNumChargesTrack");
				*NumChargesTrack = GrabInt(position, &buffer[0] );
				printf("NumChargesTrack %d\n",*NumChargesTrack);

				position = match(buffer, "\nPathFile");
				*PathFile = GrabInt(position, &buffer[0] );
				printf("PathFile %d\n",*PathFile);

				position = match(buffer, "\nLogFile");
				*LogFile = GrabInt(position, &buffer[0] );
				printf("LogFile %d\n",*LogFile);

				//converting gamma from [1/nm] to [1/m]
				(*gamma) = (*gamma)*1E9;
				//converting all alpha value from [1/nm\ to [1/m]
				(*alphaxb) = (*alphaxb)*1E9;
				(*alphaxf) = (*alphaxf)*1E9;
				(*alphayl) = (*alphayl)*1E9;
				(*alphayr) = (*alphayr)*1E9;
				(*alphazb) = (*alphazb)*1E9;
				(*alphaza) = (*alphaza)*1E9;

				fclose(handler);
				free(buffer);
					
		return 0;
			
		}

			
	return -1;
	
}

int match(char text[], char pattern[]) {
	int c, d, e, text_length, pattern_length, position = -1;

	text_length    = strlen(text);
	pattern_length = strlen(pattern);

	if (pattern_length > text_length) {
		return -1;
	}

	for (c = 0; c <= text_length - pattern_length; c++) {
		position = e = c;

		for (d = 0; d < pattern_length; d++) {
			if (pattern[d] == text[e]) {
				e++;
			}
			else {
				break;
			}
		}
		if (d == pattern_length) {
			return position;
		}
	}

	return -1;
}

int GrabInt(unsigned int position,char * buf ){

	char token[50];
	int count;
	unsigned int i;
	int j;
	int rv;

	//printf("Found at location %d\n", position + 1);
	count = 0;
	j = 0;

	for(i=0;i<30;i++){

		if(isspace(buf[position+i])){
			count++;
			if(count==3){
				break;
			}
		}

		if(count>1 && count<3){
			token[j]=buf[position+i];
			j++;
		}
	}
	token[j] = '\0';
	rv = (int) atoi(token);
	return rv;

}

double GrabDouble(unsigned int position,char * buf ){

	char token[50];
	char * ptr;
	int count;
	unsigned int i;
	int j;
	double rv;

	count = 0;
	j = 0;

	for(i=0;i<40;i++){

		if(isspace(buf[position+i])){
			count++;
			if(count==4){
				break;
			}
		}
		if(count>1 && count<3){
			token[j]=buf[position+i];
			j++;
		}
	}
	token[j] = '\0';
	rv = strtod(&token[0],&ptr);
	return rv;

}

int PFset_method(ParameterFrame PF, int method){
	if(PF == NULL || method !=0 || method !=1){
		return -1;
	}
	PF->method = method;
	return 0;
}

int PFset_Len(ParameterFrame PF,int SLength ){

	if(PF==NULL || SLength<0){
		return -1;
	}

	PF->SLength = SLength;

	return 0;
}

int PFset_Wid(ParameterFrame PF,int SWidth){

	if(PF==NULL || SWidth<0){
		return -1;
	}

	PF->SWidth = SWidth;

	return 0;
}

int PFset_Hei(ParameterFrame PF,int SHeight){

	if(PF==NULL || SHeight<0){
		return -1;
	}

	PF->SHeight = SHeight;

	return 0;
}

int PFset_Px(ParameterFrame PF,int PeriodicX){

	if(PF==NULL || PeriodicX<0 || PeriodicX>1){
		return -1;
	}

	PF->PeriodicX = PeriodicX;

	return 0;
}

int PFset_Py(ParameterFrame PF,int PeriodicY){

	if(PF==NULL || PeriodicY<0 || PeriodicY>1){
		return -1;
	}

	PF->PeriodicY = PeriodicY;

	return 0;
}

int PFset_Pz(ParameterFrame PF,int PeriodicZ){

	if(PF==NULL || PeriodicZ<0 || PeriodicZ>1){
		return -1;
	}

	PF->PeriodicZ = PeriodicZ;

	return 0;
}

int PFset_EndX(ParameterFrame PF,int EndX){

	if(PF==NULL || EndX<0 ){
		return -1;
	}

	PF->EndX = EndX;

	return 0;
}

int PFset_EndY(ParameterFrame PF,int EndY ){

	if(PF==NULL || EndY<0){
		return -1;
	}

	PF->EndY = EndY;

	return 0;
}

int PFset_EndZ(ParameterFrame PF,int EndZ){

	if(PF==NULL || EndZ<0){
		return -1;
	}

	PF->EndZ = EndZ;

	return 0;
}

int PFset_XElecOn(ParameterFrame PF,int XElecOn){

	if(PF==NULL || XElecOn<0 || XElecOn>1){
		return -1;
	}

	PF->XElecOn = XElecOn;

	return 0;
}

int PFset_YElecOn(ParameterFrame PF,int YElecOn ){

	if(PF==NULL || YElecOn<0 || YElecOn>1){
		return -1;
	}

	PF->YElecOn = YElecOn;

	return 0;
}

int PFset_ZElecOn(ParameterFrame PF,int ZElecOn ){

	if(PF==NULL || ZElecOn<0 || ZElecOn>1){
		return -1;
	}

	PF->ZElecOn = ZElecOn;

	return 0;
}

int PFset_XFermiB(ParameterFrame PF,double XFermiB ){

	if(PF==NULL){
		return -1;
	}

	PF->XFermiB = XFermiB;

	return 0;
}

int PFset_XFermiF(ParameterFrame PF,double XFermiF){

	if(PF==NULL){
		return -1;
	}

	PF->XFermiF = XFermiF;

	return 0;
}

int PFset_YFermiL(ParameterFrame PF,double YFermiL){

	if(PF==NULL){
		return -1;
	}

	PF->YFermiL = YFermiL;

	return 0;
}

int PFset_YFermiR(ParameterFrame PF,double YFermiR){

	if(PF==NULL){
		return -1;
	}

	PF->YFermiR = YFermiR;

	return 0;
}

int PFset_ZFermiB(ParameterFrame PF,double ZFermiB ){

	if(PF==NULL){
		return -1;
	}

	PF->ZFermiB = ZFermiB;

	return 0;
}

int PFset_ZFermiA(ParameterFrame PF,double ZFermiA ){

	if(PF==NULL){
		return -1;
	}

	PF->ZFermiA = ZFermiA;

	return 0;
}

int PFset_ImageCharge(ParameterFrame PF,int ImageCharge ){

	if(PF==NULL){
		return -1;
	}

	PF->ImageCharge = ImageCharge;

	return 0;
}

int PFset_IntrFermi(ParameterFrame PF,double IntrFermi ){

	if(PF==NULL){
		return -1;
	}

	PF->IntrFermi = IntrFermi;

	return 0;
}

int PFset_alphaxb(ParameterFrame PF,double alphaxb){

	if(PF==NULL){
		return -1;
	}

	PF->alphaxb = alphaxb;

	return 0;
}

int PFset_alphaxf(ParameterFrame PF,double alphaxf){

	if(PF==NULL){
		return -1;
	}

	PF->alphaxf = alphaxf;

	return 0;
}

int PFset_alphayl(ParameterFrame PF,double alphayl){

	if(PF==NULL){
		return -1;
	}

	PF->alphayl = alphayl;

	return 0;
}

int PFset_alphayr(ParameterFrame PF,double alphayr){

	if(PF==NULL){
		return -1;
	}

	PF->alphayr = alphayr;

	return 0;
}

int PFset_alphazb(ParameterFrame PF,double alphazb){

	if(PF==NULL){
		return -1;
	}

	PF->alphazb = alphazb;

	return 0;
}

int PFset_alphaza(ParameterFrame PF,double alphaza ){

	if(PF==NULL){
		return -1;
	}

	PF->alphaza = alphaza;

	return 0;
}

int PFset_vX(ParameterFrame PF,double vX ){

	if(PF==NULL || vX<0){
		return -1;
	}

	PF->vX = vX;

	return 0;
}

int PFset_vY(ParameterFrame PF,double vY ){

	if(PF==NULL || vY<0){
		return -1;
	}

	PF->vY = vY;

	return 0;
}

int PFset_vZ(ParameterFrame PF,double vZ ){

	if(PF==NULL || vZ<0){
		return -1;
	}

	PF->vZ = vZ;

	return 0;
}

int PFset_VoltageX(ParameterFrame PF,double VoltageX ){

	if(PF==NULL){
		return -1;
	}

	PF->VoltageX = VoltageX;

	return 0;
}

int PFset_VoltageY(ParameterFrame PF,double VoltageY ){

	if(PF==NULL){
		return -1;
	}

	PF->VoltageY = VoltageY;

	return 0;
}

int PFset_VoltageZ(ParameterFrame PF,double VoltageZ ){

	if(PF==NULL){
		return -1;
	}

	PF->VoltageZ = VoltageZ;

	return 0;
}

int PFset_VStepX(ParameterFrame PF,int VStepX ){

	if(PF==NULL){
		return -1;
	}

	PF->VStepX = VStepX;

	return 0;
}

int PFset_VStepY(ParameterFrame PF,int VStepY ){

	if(PF==NULL){
		return -1;
	}

	PF->VStepY = VStepY;

	return 0;
}

int PFset_VStepZ(ParameterFrame PF,int VStepZ ){

	if(PF==NULL){
		return -1;
	}

	PF->VStepZ = VStepZ;

	return 0;
}

int PFset_VincX(ParameterFrame PF,double VincX ){

	if(PF==NULL){
		return -1;
	}

	PF->VincX = VincX;

	return 0;
}

int PFset_VincY(ParameterFrame PF,double VincY ){

	if(PF==NULL){
		return -1;
	}

	PF->VincY = VincY;

	return 0;
}

int PFset_VincZ(ParameterFrame PF,double VincZ ){

	if(PF==NULL){
		return -1;
	}

	PF->VincZ = VincZ;

	return 0;
}

int PFset_SiteDist(ParameterFrame PF,double SiteDistance ){

	if(PF==NULL){
		return -1;
	}

	PF->SiteDistance = SiteDistance;

	return 0;
}

int PFset_D(ParameterFrame PF,double D ){

	if(PF==NULL){
		return -1;
	}

	PF->D = D;

	return 0;
}

int PFset_TCount(ParameterFrame PF,int TCount ){

	if(PF==NULL){
		return -1;
	}

	PF->TCount = TCount;

	return 0;
}

int PFset_NCh(ParameterFrame PF,int NCh ){

	if(PF==NULL){
		return -1;
	}

	PF->NCh = NCh;

	return 0;
}

int PFset_Ntot(ParameterFrame PF,int Ntot ){

	if(PF==NULL){
		return -1;
	}

	PF->Ntot = Ntot;

	return 0;
}

int PFset_TStep(ParameterFrame PF,double TStep ){

	if(PF==NULL){
		return -1;
	}

	PF->TStep = TStep;

	return 0;
}

int PFset_N_av(ParameterFrame PF,int N_av ){

	if(PF==NULL){
		return -1;
	}

	PF->N_av = N_av;

	return 0;
}

int PFset_Nstep_av(ParameterFrame PF,int Nstep_av ){

	if(PF==NULL){
		return -1;
	}

	PF->Nstep_av = Nstep_av;

	return 0;
}

int PFset_Time_check(ParameterFrame PF,int Time_check ){

	if(PF==NULL){
		return -1;
	}

	PF->Time_check = Time_check;

	return 0;
}

int PFset_Rcount(ParameterFrame PF,int Rcount ){

	if(PF==NULL){
		return -1;
	}

	PF->Rcount = Rcount;

	return 0;
}
int PFset_ClusterAlg(ParameterFrame PF,int ClusterAlg ){

	if(PF==NULL){
		return -1;
	}

	PF->ClusterAlg = ClusterAlg;

	return 0;
}


int PFset_CutOff(ParameterFrame PF,double CutOff ){

	if(PF==NULL){
		return -1;
	}

	PF->CutOff = CutOff;

	return 0;
}

int PFset_lambda(ParameterFrame PF,double lambda ){

	if(PF==NULL){
		return -1;
	}

	PF->lambda = lambda;

	return 0;
}

int PFset_ScaleAfterCorr(ParameterFrame PF,int ScaleAfterCorr ){

	if(PF==NULL){
		return -1;
	}

	PF->ScaleAfterCorr = ScaleAfterCorr;

	return 0;
}

int PFset_SeedProt(ParameterFrame PF,int SeedProt ){

	if(PF==NULL){
		return -1;
	}

	PF->SeedProt = SeedProt;

	return 0;
}

int PFset_Attempts(ParameterFrame PF,int Attempts ){

	if(PF==NULL){
		return -1;
	}

	PF->Attempts = Attempts;

	return 0;
}

int PFset_FracSeed(ParameterFrame PF,double fracSeed ){

	if(PF==NULL){
		return -1;
	}

	PF->fracSeed = fracSeed;

	return 0;
}

int PFset_E0(ParameterFrame PF,double E0 ){

	if(PF==NULL){
		return -1;
	}

	PF->E0 = E0;

	return 0;
}

int PFset_sigma(ParameterFrame PF,double sigma ){

	if(PF==NULL){
		return -1;
	}

	PF->sigma = sigma;

	return 0;
}

int PFset_FracTrap(ParameterFrame PF,double fracTrap ){

	if(PF==NULL){
		return -1;
	}

	PF->fracTrap = fracTrap;

	return 0;
}

int PFset_Etrap(ParameterFrame PF,double Etrap ){

	if(PF==NULL){
		return -1;
	}

	PF->Etrap = Etrap;

	return 0;
}

int PFset_Tsigma(ParameterFrame PF,double Tsigma ){

	if(PF==NULL){
		return -1;
	}

	PF->Tsigma = Tsigma;

	return 0;
}

int PFset_TempStart(ParameterFrame PF,double TempStart ){

	if(PF==NULL){
		return -1;
	}

	PF->TempStart = TempStart;

	return 0;
}

int PFset_TempStep(ParameterFrame PF,int TempStep ){

	if(PF==NULL){
		return -1;
	}

	PF->TempStep = TempStep;

	return 0;
}

int PFset_TempInc(ParameterFrame PF,double TempInc ){

	if(PF==NULL){
		return -1;
	}

	PF->TempInc = TempInc;

	return 0;
}

int PFset_reOrg(ParameterFrame PF,double reOrg ){

	if(PF==NULL){
		return -1;
	}

	PF->reOrg = reOrg;

	return 0;
}

int PFset_AttemptToHop(ParameterFrame PF,double AttemptHop ){

	if(PF==NULL){
		return -1;
	}

	PF->AttemptToHop = AttemptHop;

	return 0;
}

int PFset_gamma(ParameterFrame PF,double gamma ){

	if(PF==NULL){
		return -1;
	}

	PF->gamma = gamma;

	return 0;
}

int PFset_RelativePerm(ParameterFrame PF,double RelativePerm){

	if(PF==NULL){
		return -1;
	}

	PF->RelativePerm = RelativePerm;

	return 0;
}

int PFset_MovieFrames(ParameterFrame PF, int MovieFrames){

	if(PF==NULL){
		return -1;
	}

	PF->MovieFrames = MovieFrames;
	return 0;

}

int PFset_CutOffTime(ParameterFrame PF, double CutOffTime){
	if(PF==NULL || CutOffTime<0){
		return -1;
	}
	PF->CutOffTime = CutOffTime;
	return 0;
}

int PFset_Vcv(ParameterFrame PF, double Vcv){
	if(PF==NULL){
		return -1;
	}
	PF->Vcv=Vcv;
	return 0;
}

int PFset_Tcv(ParameterFrame PF, double Tcv){
	if(PF==NULL || Tcv<0){
		return -1;
	}
	PF->Tcv=Tcv;
	return 0;
}

int PFset_Tlag(ParameterFrame PF, double Tlag){
	if(PF==NULL || Tlag<0){
		return -1;
	}
	PF->Tlag=Tlag;
	return 0;
}

int PFset_EndPtFile(ParameterFrame PF, int EndPtFile){
	if(PF==NULL || EndPtFile<0 || EndPtFile>1){
		return -1;
	}
	PF->EndPtFile=EndPtFile;
	return 0;
}

int PFset_NumChargesTrack(ParameterFrame PF, int NumChargesTrack){
	if(PF==NULL || NumChargesTrack<0 ){
		return -1;
	}
	PF->NumChargesTrack=NumChargesTrack;
	return 0;
}

int PFset_PathFile(ParameterFrame PF, int PathFile){
	if(PF==NULL || PathFile<0 || PathFile>1){
		return -1;
	}
	PF->PathFile=PathFile;
	return 0;
}

int PFset_LogFile(ParameterFrame PF, int LogFile){
	if(PF==NULL || LogFile<0 || LogFile>1){
		return -1;
	}
	PF->LogFile=LogFile;
	return 0;
}

int PFget_method(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->method; 
}

int PFget_Len(ParameterFrame PF){

	if(PF==NULL ){
		return -1;
	}

	return	PF->SLength;

}

int PFget_Wid(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return	PF->SWidth;

}

int PFget_Hei(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->SHeight;
}

int PFget_Px(ParameterFrame PF){

	if(PF==NULL ){
		return -1;
	}

	return PF->PeriodicX;

}

int PFget_Py(ParameterFrame PF){

	if(PF==NULL ){
		return -1;
	}

	return PF->PeriodicY;
}

int PFget_Pz(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->PeriodicZ;
}

int PFget_EndX(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->EndX;
}

int PFget_EndY(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->EndY;
}

int PFget_EndZ(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->EndZ;
}

int PFget_XElecOn(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->XElecOn;

	return 0;
}

int PFget_YElecOn(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->YElecOn;
}

int PFget_ZElecOn(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->ZElecOn;
}

double PFget_XFermiB(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->XFermiB;
}

double PFget_XFermiF(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->XFermiF;
}

double PFget_YFermiL(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->YFermiL;
}

double PFget_YFermiR(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->YFermiR;
}

double PFget_ZFermiB(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->ZFermiB;
}

double PFget_ZFermiA(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->ZFermiA;
}

double PFget_ImageCharge(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->ImageCharge;
}

double PFget_IntrFermi(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->IntrFermi;
}

double PFget_alphaxb(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphaxb;
}

double PFget_alphaxf(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphaxf;
}

double PFget_alphayl(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphayl;
}

double PFget_alphayr(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphayr;
}

double PFget_alphazb(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphazb;
}

double PFget_alphaza(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->alphaza;
}

double PFget_vX(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->vX;
}

double PFget_vY(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->vY;
}

double PFget_vZ(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->vZ;
}

double PFget_VoltageX(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VoltageX;
}

double PFget_VoltageY(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VoltageY;
}

double PFget_VoltageZ(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VoltageZ;
}

int PFget_VStepX(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VStepX;
}

int PFget_VStepY(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VStepY;
}

int PFget_VStepZ(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VStepZ;
}

double PFget_VincX(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VincX;
}

double PFget_VincY(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VincY;
}

double PFget_VincZ(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->VincZ;
}

double PFget_SiteDist(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->SiteDistance;
}

double PFget_D(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->D;
}

int PFget_TCount(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->TCount;
}

int PFget_NCh(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->NCh;
}

int PFget_Ntot(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Ntot;
}

double PFget_TStep(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->TStep;
}

int PFget_N_av(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->N_av;
}

int PFget_Nstep_av(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Nstep_av;
}

int PFget_Time_check(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Time_check;
}

int PFget_Rcount(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->Rcount;
}

int PFget_ClusterAlg(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->ClusterAlg;
}

double PFget_CutOff(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->CutOff;
}

double PFget_lambda(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->lambda;
}

int PFget_ScaleAfterCorr(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->ScaleAfterCorr;
}

int PFget_SeedProt(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->SeedProt;
}

int PFget_Attempts(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Attempts;
}

double PFget_FracSeed(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->fracSeed;
}

double PFget_E0(ParameterFrame PF ){

	if(PF==NULL){
		return -1;
	}

	return PF->E0;
}

double PFget_sigma(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->sigma;
}

double PFget_FracTrap(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->fracTrap;
}

double PFget_Etrap(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Etrap;
}

double PFget_Tsigma(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->Tsigma;
}

double PFget_TempStart(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->TempStart;
}

int PFget_TempStep(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->TempStep;
}

double PFget_TempInc(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->TempInc;
}

double PFget_reOrg(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->reOrg;
}

double PFget_AttemptToHop(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->AttemptToHop;
}

double PFget_gamma(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->gamma;
}

double PFget_RelativePerm(ParameterFrame PF){

	if(PF==NULL){
		return -1;
	}

	return PF->RelativePerm;
}

int PFget_MovieFrames(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->MovieFrames;
}

double PFget_CutOffTime(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->CutOffTime;
}

double PFget_Vcv(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->Vcv;
}

double PFget_Tcv(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->Tcv;
}

double PFget_Tlag(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->Tlag;
}

int PFget_EndPtFile(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->EndPtFile;
}

int PFget_NumChargesTrack(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->NumChargesTrack;
}

int PFget_PathFile(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->PathFile;
}

int PFget_LogFile(ParameterFrame PF){
	if(PF==NULL){
		return -1;
	}
	return PF->LogFile;
}
