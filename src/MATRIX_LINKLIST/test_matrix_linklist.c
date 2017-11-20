#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "matrix_linklist.h"

int main(){ 
	
	printf("Starting\n");
  
  int rv;
  double rvd;
  double data[3];
  double data2[3];
  data[0] = -5.0;
  data[1] = -5.0;
  data[2] = -5.0;
	printf("Testing: newLL_MNode\n");
	LL_MNode LLrv=newLL_MNode(data, 3);
  assert(LLrv!=NULL);

  printf("Testing: newMatrix_LinkList\n");
  matrix_linklist LL = newMatrix_LinkList(data, 3);
  assert(LL!=NULL);

  printf("Testing: getMatrixLLNodeElem\n");
  rvd = getMatrixLLNodeElem(LL, 1, 1);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL, 1, 2);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL, 1, 3);
  assert(rvd==-5.0);

  printf("Testing: newBlankMatrixLinkList\n");
  matrix_linklist LL2 = newBlankMatrixLinkList();
  assert(LL2!=NULL);
  
  printf("Testing: deleteMatrixLL\n");
  rv = deleteMatrixLL(NULL);
  assert(rv==-1);
  rv = deleteMatrixLL(&LL);
  assert(rv==0);
  assert(LL==NULL);
   
  printf("Testing: deleteLL_MNode\n");
  rv = deleteLL_MNode(NULL);
  assert(rv==-1);
  rv = deleteLL_MNode(&LLrv);
  assert(rv==0);
  assert(LLrv==NULL);
  
  printf("Testing: addLL_MNode\n");
  fprintf(stderr,"Should trigger error message ");
  rv = addLL_MNode(NULL, data, 3);
  assert(rv==-1);
  fprintf(stderr,"Should trigger error message ");
  rv = addLL_MNode(LL2, data, 0);
  assert(rv==-1);
  rv = addLL_MNode(LL2, data, 3);
  data2[0]=4;
  data2[1]=1234;
  data2[2]=52.1;
  rv = addLL_MNode(LL2, data2, 3);
  assert(rv==0);
  rvd = getMatrixLLNodeElem(LL2, 1, 1);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 1, 2);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 1, 3);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 2, 1);
  assert(rvd==4);
  rvd = getMatrixLLNodeElem(LL2, 2, 2);
  assert(rvd==1234);
  rvd = getMatrixLLNodeElem(LL2, 2, 3);
  assert(rvd==52.1);

  printf("Testing: removeLL_MNode\n");
  rv = removeLL_MNode(LL2, 0);
  assert(rv==-1);
  rvd = getMatrixLLNodeElem(LL2, 1, 1);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 1, 2);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 1, 3);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 2, 1);
  assert(rvd==4);
  rvd = getMatrixLLNodeElem(LL2, 2, 2);
  assert(rvd==1234);
  rvd = getMatrixLLNodeElem(LL2, 2, 3);
  assert(rvd==52.1);
  rv = removeLL_MNode(NULL, 1);
  assert(rv==-1);
  rv = removeLL_MNode(LL2, 1);
  assert(rv==0);
  rvd = getMatrixLLNodeElem(LL2, 1, 1);
  assert(rvd==4);
  rvd = getMatrixLLNodeElem(LL2, 1, 2);
  assert(rvd==1234);
  rvd = getMatrixLLNodeElem(LL2, 1, 3);
  assert(rvd==52.1);

  printf("Testing: getMatrixLLlength\n");
  rv = addLL_MNode(LL2, data, 3);
  rv = addLL_MNode(LL2, data2, 3);
  rv = getMatrixLLlength(NULL);
  assert(rv==-1);
  rv = getMatrixLLlength(LL2);
  assert(rv==3);
  rv = addLL_MNode(LL2, data, 3);
  rv = getMatrixLLlength(LL2);
  assert(rv==4);

  printf("Testing: removeLL_MNodes\n");
  fprintf(stderr,"Should trigger Error message ");
  rv = removeLL_MNodes(NULL,2,3);
  assert(rv==-1);
  fprintf(stderr,"Should trigger Error message ");
  rv = removeLL_MNodes(LL2,3,2);
  assert(rv==-1);
  fprintf(stderr,"Should trigger Error message ");
  rv = removeLL_MNodes(LL2,0,3);
  assert(rv==-1);
  rv = getMatrixLLlength(LL2);
  assert(rv==4);
  rv = removeLL_MNodes(LL2,2,3);
  assert(rv==0);
  rv = getMatrixLLlength(LL2);
  assert(rv==2);
  rvd = getMatrixLLNodeElem(LL2, 1, 1);
  assert(rvd==4);
  rvd = getMatrixLLNodeElem(LL2, 1, 2);
  assert(rvd==1234);
  rvd = getMatrixLLNodeElem(LL2, 1, 3);
  assert(rvd==52.1);
  rvd = getMatrixLLNodeElem(LL2, 2, 1);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 2, 2);
  assert(rvd==-5.0);
  rvd = getMatrixLLNodeElem(LL2, 2, 3);
  assert(rvd==-5.0);
 
  printf("Testing: setMatrixLLValueAtRow\n");
  rv = setMatrixLLValueAtRow(LL2, 2, 1.0);
  assert(rv==0);
  rvd = getMatrixLLNodeElem(LL2, 1, 2);
  assert(rvd==1.0);
  rvd = getMatrixLLNodeElem(LL2, 2, 2);
  assert(rvd==1.0);

  
  rv = deleteMatrixLL(&LL2); 
  assert(rv==0);
  assert(LL2==NULL);

  return 0;
}
