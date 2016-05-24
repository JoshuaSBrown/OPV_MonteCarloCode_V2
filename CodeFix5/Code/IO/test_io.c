#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../IO/io.h"

int main(){

	char FileName[] = "TestFile";
	printf("%s\n",FileName);
	
	FILE * EndPtFile;
	EndPtFile = openEndPtFile( FileName );
	printToEndPtFile( EndPtFile, 1, 1, 1, 4, 23.53);
	printToEndPtFile( EndPtFile, 2, 1, 1, 4, 53.83);
	closeEndPtFile( EndPtFile);
	
	FILE * PathFile; 
	PathFile = openPathFile( FileName );
	printToPathFile( PathFile, 1, 3, 5, 8, 243.01, 59094.1090);
	closePathFile( PathFile);

	return 0;
}
