#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "analysis.h"

int main(int argc, char *argv[]){

  char * FileName[64];
  printf("%s\n",argv[1]);
  *FileName = argv[1];
  readPath(*FileName);

  return 0;
}
