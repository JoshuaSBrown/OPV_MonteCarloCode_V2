#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "analysis.h"

/* To run this simply type:
 * run_analysis filename.path
 */
int main(int argc, char *argv[]){
  
  if(argc<2){
    usage();
    exit(1);
  }
  char * FileName[64];
  printf("%s\n",argv[1]);
  *FileName = argv[1];
  readPath(*FileName);

  return 0;
}
