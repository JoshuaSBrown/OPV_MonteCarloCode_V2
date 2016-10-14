
#include <stdlib.h>
#include "error.h"

int return_error_val(void){
  #ifdef _FORCE_HARD_CRASH_
  exit(1);
  #endif
  return -1;
}
