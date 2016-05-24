#ifndef _MONTECARLO_IO_H_
#define _MONTECARLO_IO_H_

#include <stdio.h>

//Open the EndPtFile and pass out the file identifier
FILE *openEndPtFile(char *FileName);

//Print to the EndPtFile should contain the position the chargeid and the time
int printToEndPtFile(FILE *EndPtOut, int x, int y, int z, int ChargeId, long double GlobalTime);

//Close the EndPtFile
int closeEndPtFile(FILE *EndPtOut);

//Open the PathFile and pass out the file identifier
FILE *openPathFile(char *FileName);

//Print to the PathFile should contain the position the chargeid and the time on the site and the global time
int printToPathFile(FILE *EndPtOut, int x, int y, int z, int ChargeId, double time, long double GlobalTime);

//Close the PathFile
int closePathFile(FILE *EndPtOut);

//Open the LogFile and pass out the file identifier
FILE *openLogFile(char *FileName);

//Append to the LogFile and pass out the file identifier
FILE *appendLogFile(char *FileName);

//Print the total number of hops and the total failed hops 
int LogFile_printHops(FILE *LogFile, long int TotalHopAttempt, long int FailedHop);

//Print the total time the simulation took to run to the log file
int LogFile_printTime(FILE *LogFile, double time_spent, int seconds);

//Close the LogFile
int closeLogFile(FILE *LogFile);

#endif
