#include <stdlib.h>
#include <stdio.h>

#include "io.h"

/*FILE **openEndPtFile(char *FileName){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".endpt");
	FILE * EndPtOut = malloc(sizeof (FILE *) );
	if((EndPtOut=fopen(buf,"w"))==NULL){
		printf("ERROR! unable to write .endpt file!\n");
		exit(1);
	}
	return &(EndPtOut);
}
*/

int printToEndPtFile(char *FileName,int x, int y, int z, int ChargeId, long double GlobalTime){
	

	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".endpt");
	FILE * EndPtOut;

	if((EndPtOut=fopen(buf,"a"))==NULL){
		printf("ERROR unable to open file!\n");
		exit(1);
	}
	//printf("%d %d %d %d %Le\n",x,y,z,ChargeId,GlobalTime);
	fprintf(EndPtOut,"%d %d %d %d %Le\n",x,y,z,ChargeId,GlobalTime);
	fclose(EndPtOut);
	return 0;
}
/*
int closeEndPtFile(FILE *EndPtOut){

	fclose(EndPtOut);
	free(EndPtOut);
	return 0;
}
*/
/*
FILE *openPathFile(char *FileName){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".path");
	FILE * PathFile = malloc(sizeof (FILE *));
	if((PathFile = fopen(buf,"w"))==NULL){
		printf("ERROR! unable to write .path file!\n");
		exit(1);
	}
	return PathFile;

}
*/
int printToPathFile(char *FileName,int x, int y, int z, int ChargeId, double Time,long double GlobalTime){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".path");
	FILE * PathFile;
	if((PathFile=fopen(buf,"a"))==NULL){
		printf("ERROR! unable to write .path file!\n");
		exit(1);
	}
	fprintf(PathFile,"%d %d %d %d %g %Le\n",x,y,z,ChargeId,Time,GlobalTime);
	fclose(PathFile);
	return 0;
}
/*
int closePathFile(FILE *PathFile){
	
	fclose(PathFile);
	free(PathFile);
	return 0;
}
*/
/*
FILE *openLogFile(char *FileName){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".log");
	FILE * LogFile = malloc(sizeof(FILE *));
	if((LogFile = fopen(buf,"w"))==NULL){
		printf("ERROR! unable to write .log file!\n");
		exit(1);
	}
	return LogFile;

}

FILE *appendLogFile(char *FileName){

	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".log");
	FILE * LogFile = malloc(sizeof(FILE *));
	if((LogFile = fopen(buf,"a"))==NULL){
		printf("ERROR! unable to write .log file!\n");
		exit(1);
	}
	return LogFile;
}
*/
int LogFile_printHops( char *FileName, long int TotalHopAttempt, long int FailedHop){

	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".log");
	FILE * LogFile;
	if((LogFile = fopen(buf,"a"))==NULL){
		printf("ERROR! unable to open .log file!\n");
		exit(1);
	}
	fprintf(LogFile,"Total Hop Attempts  %ld\n",TotalHopAttempt);
	fprintf(LogFile,"Total Failed Hops   %ld\n",FailedHop);
	fclose(LogFile);
	return 0;
}

int LogFile_printTime( char *FileName, double time_spent, int seconds){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".log");
	FILE * LogFile;
	if((LogFile = fopen(buf,"a"))==NULL){
		printf("ERROR! unable to open .log file!\n");
		exit(1);
	}

	fprintf(LogFile,"Run Time %g seconds\n",time_spent);
	fprintf(LogFile,"user : %d secs\n", seconds);
	fclose(LogFile);
	return 0;
}
/*
int closeLogFile(FILE *LogFile){
	fclose(LogFile);
	free(LogFile);
	return 0;
}

*/
