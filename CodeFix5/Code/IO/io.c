#include <stdlib.h>
#include <stdio.h>

#include "io.h"

FILE *openEndPtFile(char *FileName){
	
	char buf[256];
	snprintf(buf,sizeof buf, "%s%s",FileName,".endpt");
	FILE * EndPtOut = malloc(sizeof (FILE *) );
	if((EndPtOut=fopen(buf,"w"))==NULL){
		printf("ERROR! unable to write .endpt file!\n");
		exit(1);
	}
	return (EndPtOut);
}

int printToEndPtFile(FILE *EndPtOut,int x, int y, int z, int ChargeId, long double GlobalTime){

	if(EndPtOut==NULL){
		printf("ERROR File identifier is NULL!\n");
		exit(1);
	}
	printf("%d %d %d %d %Le\n",x,y,z,ChargeId,GlobalTime);
	fprintf(EndPtOut,"%d %d %d %d %Le\n",x,y,z,ChargeId,GlobalTime);
	return 0;
}

int closeEndPtFile(FILE *EndPtOut){

	fclose(EndPtOut);
	free(EndPtOut);
	return 0;
}

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

int printToPathFile(FILE *PathFile,int x, int y, int z, int ChargeId, double Time,long double GlobalTime){
	
	fprintf(PathFile,"%d %d %d %d %g %Le\n",x,y,z,ChargeId,Time,GlobalTime);
	return 0;
}

int closePathFile(FILE *PathFile){
	
	fclose(PathFile);
	free(PathFile);
	return 0;
}

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

int LogFile_printHops( FILE *LogFile, long int TotalHopAttempt, long int FailedHop){

		fprintf(LogFile,"Total Hop Attempts  %ld\n",TotalHopAttempt);
		fprintf(LogFile,"Total Failed Hops   %ld\n",FailedHop);
		return 0;
}

int LogFile_printTime( FILE * LogFile, double time_spent, int seconds){
		fprintf(LogFile,"Run Time %g seconds\n",time_spent);
		fprintf(LogFile,"user : %d secs\n", seconds);
		return 0;
}

int closeLogFile(FILE *LogFile){
	fclose(LogFile);
	free(LogFile);
	return 0;
}


