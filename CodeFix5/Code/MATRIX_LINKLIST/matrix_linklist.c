
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "matrix_linklist.h"
#include "../MATRIX/matrix.h"
#include "../ERROR/error.h"

struct _matrix_linklist {
	int      length;
	LL_MNode start;
};

struct _LL_MNode {
	matrix   data;
	LL_MNode next;
};

LL_MNode newLL_MNode( double data[],int rows ){
	
	int i;
	
	if(rows <= 0){
		fprintf(stderr,"ERROR trying to initialize matrix_");	
		fprintf(stderr,"linklist node with no data\n");
		return NULL;	
	}

	LL_MNode node = (LL_MNode) malloc(sizeof(struct _LL_MNode));

	if(node==NULL){
		return NULL;
	}

	//Create Matrix
	node->data = newMatrix(rows,1);
	//Initialize the Matrix with the passed data
	for(i=0;i<rows;i++){
		setE(node->data,i+1,1,data[i]);
	}
	node->next = NULL;
	return node;
}

matrix_linklist newMatrix_LinkList( double  data[],int rows ){

	matrix_linklist LL = (matrix_linklist) malloc(sizeof(struct _matrix_linklist));

	if( LL == NULL){
		return NULL;
	}

	LL_MNode node = newLL_MNode(data,rows);
	
	if(node==NULL){
		//deleteMatrixLL( &LL);
		return NULL;
	}

	LL->length = 1;
	LL->start  = node;

	return LL;
}

matrix_linklist newBlankMatrixLinkList(){
	
	matrix_linklist LL = (matrix_linklist) malloc(sizeof(struct _matrix_linklist));

	if( LL == NULL){
		return NULL;
	}
	
	LL->length=0;
	LL->start = NULL;

	return LL;
}

//Delete entire matrix linklist
int deleteMatrixLL(matrix_linklist * LL){ 
	if(LL==NULL){
		return -1;
	}

	LL_MNode tempNode;
	LL_MNode tempNode2;

	tempNode = (*LL)->start;
	while (tempNode!=NULL){
		tempNode2 = tempNode;
		tempNode  = tempNode->next;
		deleteLL_MNode(&tempNode2);
	}

	free((*LL));

	*LL = NULL;
	return 0;
}

//deleting separate node, not part of linked list
int deleteLL_MNode( LL_MNode * node){ 
	if(node==NULL){
		return -1;
	}
  if((*node)==NULL){
    return -1;
  }
  matrix data;
  data = (*node)->data;
	deleteMatrix(&data);
  assert(data==NULL);
	free((*node));
  *node = NULL;
	return 0;
}

LL_MNode nextLL_MNode(LL_MNode node){
	if(node==NULL){
		return NULL;
	}

	return node->next;
}

int addLL_MNode(matrix_linklist * LL, double data[], int rows){

  #ifdef _ERROR_CHECKING_ON_	
	if(*LL==NULL ){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR addLL_MNode failed ");
		fprintf(stderr,"because *LL was NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif
	int row = sizeof(*data)/sizeof(double);
  #ifdef _ERROR_CHECKING_ON_	
	if(row<1){
    #ifdef _ERROR_
		fprintf(stderr,"ERROR addLL_MNode failed ");
		fprintf(stderr,"because data was less than 1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	//Cycle through matrix_linklist
	LL_MNode tempNode;
	LL_MNode oldNode;
	LL_MNode newNode;
	tempNode = (*LL)->start;
  //This means the first node in the LL has not
  //yet been assigned
  printf("addLL_MNode Step1\n");
  if (tempNode==NULL){
		newNode = newLL_MNode(data, rows);
		if(newNode==NULL){
  printf("addLL_MNode Step2\n");
			return -1;
		}
		(*LL)->length++;
		(*LL)->start = newNode;
  printf("addLL_MNode Step3\n");
		return 0;
	}

	while(tempNode!=NULL){
    //Searching for the last node in the LL
		oldNode  = tempNode;
		tempNode = tempNode->next;
	}
	
  //Appending a new node to the last position in the LL
	newNode = newLL_MNode(data, rows);
	if(newNode==NULL){
  printf("addLL_MNode Step4\n");
		return -1;
	}
	(*LL)->length++;
	oldNode->next = newNode;
  printf("addLL_MNode Step5 LL length %d\n",(*LL)->length);
	return 0;
}

int addLL_MNodeBegin(matrix_linklist * LL, double data[], int rows){
	
  #ifdef _ERROR_CHECKING_ON_
	if(*LL==NULL ){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR addLL_MNode failed ");
		fprintf(stderr,"because *LL was NULL\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif

	int row = sizeof(*data)/sizeof(double);
  #ifdef _ERROR_CHECKING_ON_
	if(row<1){
		#ifdef _ERROR_
		fprintf(stderr,"ERROR addLL_MNode failed ");
		fprintf(stderr,"because data was less than 1\n");
		#endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	
  LL_MNode tempNode;
	LL_MNode newNode;
	tempNode = (*LL)->start;
  //This means the first node in the LL has not
  //yet been assigned
  if (tempNode==NULL){
		newNode = newLL_MNode(data, rows);
		if(newNode==NULL){
			return -1;
		}
		(*LL)->length++;
		(*LL)->start = newNode;
		return 0;
	}

  //Appending a new node to the front of the LL
	newNode = newLL_MNode(data, rows);
	if(newNode==NULL){
		return -1;
	}
	(*LL)->length++;
  (*LL)->start     = newNode;
  newNode->next = tempNode;
	return 0;
}

//deleting one node from the linked list
//based on it's sequence in the LL
int removeLL_MNode(matrix_linklist * LL, int seq){ 	
  
  #ifdef _ERROR_CHECKING_ON_
	if(LL==NULL || seq<=0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR matrix linklist is NULL or ");
    fprintf(stderr,"seq is less than 1 in removeLL_MNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
    
	if(seq>(*LL)->length){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR seq is greater than the length");
    fprintf(stderr," of the matrix linklist\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
 	  return -1; 
	}
  #endif
	LL_MNode node    = (*LL)->start;
	LL_MNode oldnode = node;
	LL_MNode newnode;
 	int count;
	count = 1;
  //This case is triggered if the first node in the
  //LL is the one we are removing
	if (seq == count){
		newnode   = node->next;
		(*LL)->start = newnode;
		deleteLL_MNode(&node);
		(*LL)->length--;
		return 0;
	}

  //This is if we are searching for a node 
  //that is not the first in the LL, to remove. 	
	while(node!=NULL ){
	
		count++;	
		if(count==seq){
			newnode = node->next;
			oldnode->next = newnode;
			deleteLL_MNode(&node);
			(*LL)->length--;
			//Found node deleted it and 
			//closed the matrix_linklist
			return 0;
		}
		oldnode = node;
		node = node->next;
	}

	return -1;
}

int removeLL_MNodeEnd(matrix_linklist * LL){
  
  #ifdef _ERROR_CHECKING_ON_
  if((*LL)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR matrix linklist is NULL in removeLL_MNodeEnd function\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int rv;
  LL_MNode temp = (*LL)->start;
  LL_MNode begin;
  if((*LL)->length==1){
    rv = deleteLL_MNode(&temp);
    if (rv==-1){
      fprintf(stderr,"ERROR unable to delete LL_MNode in removeLL_MNodeEnd\n");
      return -1;
    }
    (*LL)->length--;
    (*LL)->start = NULL;
    return 0;
  }
  while(temp->next!=NULL){
    begin = temp;
    temp  = temp->next;
  }
  rv = deleteLL_MNode(&temp);
  if (rv==-1){
    fprintf(stderr,"ERROR unable to delete LL_MNode in removeLL_MNodeEnd\n");
    return -1;
  }
  (*LL)->length--;
  begin->next = NULL;
  return 0;
}

//Deleting nodes between and including those defined 
//by seq1 and seq2, seq1 must be less than seq2
int removeLL_MNodes(matrix_linklist * LL, int seq1, int seq2){

  #ifdef _ERROR_CHECKING_ON_
  if(*LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR matrix_Linklist is null in removeLL_MNodes\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  } 
  if(seq1>=seq2){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR seq1 is greater than or equal to seq2 in removeLL_MNodes\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(seq1<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR seq1 is less than or equal to 0 in removeLL_MNodes\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(seq2>(*LL)->length){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR seq2 is greater than the length of the matrix_linklist in removeLL_MNodes\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  LL_MNode node = (*LL)->start;
  LL_MNode oldnode = node;
  LL_MNode newnode;
  int count;
  count = 1;
  //This case is triggered if the first node in the
  //LL is the one we are removing
  if (seq1 == count){
    while(count<=seq2){
      newnode = node->next;
      deleteLL_MNode(&node);
      node = newnode;
      count++;
      (*LL)->length--;
    }
    (*LL)->start=newnode;
    return 0;
  }

  //This is if we are searching for nodes  
  //that are not the first in the LL, to remove. 	
  while(node!=NULL ){
    
    if(count==seq1){
      while(count<=seq2){
        newnode = node->next;
        deleteLL_MNode(&node);
        node = newnode;
        (*LL)->length--;
        //Found node deleted it and 
        //closed the matrix_linklist
        count++;
      }
      oldnode->next = newnode;
      return 0;
    }
    oldnode = node;
    node = node->next;
    count++;	
  }
  return -1;
}

int getMatrixLLlength(matrix_linklist LL){
  if(LL==NULL){
    return -1;
  }

  return LL->length;
}

double getMatrixLLstartAtRow(matrix_linklist LL, int R){
  if(LL==NULL){
    return -1;
  }
  
  if(LL->start==NULL){
    fprintf(stderr,"ERROR getMatrixLLstartAtRow has error because ");
    fprintf(stderr,"LL->start is NULL\n");
    return -1;
  }

  int rows;
  rows = getRows(LL->start->data);
  if(R>rows){
    fprintf(stderr,"ERROR getMatrixLLstartAtRow R value is larger ");
    fprintf(stderr,"than the number of rows\n");
    return -1;
  }
  return getE(LL->start->data,R,1);

}

int setMatrixLLElem(matrix_linklist LL, int seq, int R, double val){
  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"LL is NULL in setMatrixLLValueAtRow\n"); 
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  if(LL->start==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR setMatrixLLElem has error because ");
    fprintf(stderr,"LL->start is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  if(seq<=0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR setMatrixLLElem has error because ");
    fprintf(stderr,"seq<=0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  if(seq>getMatrixLLlength(LL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR setMatrixLLElem has error because ");
    fprintf(stderr,"seq>=getMatrixLLlength(LL)\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }

  #endif
  int rows;
  rows = getRows(LL->start->data);
  #ifdef _ERROR_CHECKING_ON_
  if(R>rows){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR setMatrixLLElem R value is larger ");
    fprintf(stderr,"than the number of rows\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  #endif

  int inc=1;
  LL_MNode temp = LL->start;
  while(inc<seq){
    temp = temp->next;
    inc++; 
  }
  setE(temp->data,R,1,val);
  
  return 0;
}

int setMatrixLLValueAtRow(matrix_linklist LL, int R, double val){
  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"LL is NULL in setMatrixLLValueAtRow\n"); 
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  if(LL->start==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getMatrixLLstartAtRow has error because ");
    fprintf(stderr,"LL->start is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  #endif
  int rows;
  rows = getRows(LL->start->data);
  #ifdef _ERROR_CHECKING_ON_
  if(R>rows){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR getMatrixLLstartAtRow R value is larger ");
    fprintf(stderr,"than the number of rows\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
     exit(1);
    #endif
    return -1;
  }
  #endif
  setE(LL->start->data,R,1,val);
  LL_MNode temp = LL->start->next;
  while(temp!=NULL){
    setE(temp->data,R,1,val); 
    temp = temp->next;
  }
  return 0;
}

int printMatrixLL(matrix_linklist LL){
  
  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR LL is NULL in printMatrixLL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif

  int rows;
  int i;
	printf("\nLength %d\n",LL->length);
	LL_MNode node = LL->start;
  int count;
  count = 1;
	while(node!=NULL){
    rows = getRows(node->data);
		printf("Sequence %d \t Rows of Data %d\n",count,rows);
    for(i=1;i<=rows;i++){
      printf("Row %d %g\n",i,getE(node->data,i,1));
    }
		node=node->next;
    count++;
	}
  printf("Finished with printMatrixLL\n");

	return 0;
}

int printMNode( LL_MNode node ){

  int rows;
  int i;
  int count;
  count = 1;
  while(node!=NULL){
    rows = getRows(node->data);
		printf("Sequence %d \t Rows of Data %d\n",count,rows);
    for(i=1;i<=rows;i++){
      printf("Row %d %g\n",i,getE(node->data,i,1));
    }
		node=node->next;
    count++;
	}
  return 0;
}

int getMatrixLLLastMatchAtRow(matrix_linklist LL, double match, int R){
  
  if(LL==NULL){
    fprintf(stderr,"ERROR LL is NULL in getMatrixLLLastMatch\n");
    return -1;
  }
  LL_MNode node = LL->start;
  int rows;
  int seq = 0;
  int lastseq = -1;
  while(node!=NULL){
    seq++;
    rows = getRows(node->data);
    
    if(R>rows){
      fprintf(stderr,"ERROR getMatrixLLLastMatch data data only\n");
      fprintf(stderr," has %d rows but you have requested\n",rows);
      fprintf(stderr," a match with row %d.\n",R);
      return -1;
    }

    if(getE(node->data,R,1)==match){
      lastseq = seq;
    }
    node=node->next;
  }
    
  return lastseq;
}

int getMatrixLLNumberMatchAtRow(matrix_linklist LL, double match, int R){

  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR LL is NULL in getMatrixLLNumberMatchAtRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif  
  
  int rows;
  LL_MNode node = LL->start;
  int count = 0;
  while(node!=NULL){

    rows = getRows(node->data);
    
    if(R>rows){
      fprintf(stderr,"ERROR getMatrixLLNumberMatchAtRow data only\n");
      fprintf(stderr," has %d rows but you have requested\n",rows);
      fprintf(stderr," a match with row %d.\n",R);
      return -1;
    }
    
    if(getE(node->data,R,1)==match){
      count++;
    }
    node=node->next;
  }
    
  return count;
}

matrix getMatrixLLMatchAtRow(matrix_linklist LL, double match, int R,  int * numMatches){

  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR LL is NULL in getMatrixLLMatchAtRow\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
  }
  #endif

  *numMatches =getMatrixLLNumberMatchAtRow( LL,(double)match,R);
 
  matrix matches; 
  #ifdef _ERROR_CHECKING_ON_
  if(*numMatches<0){
    matches=NULL;
    #ifdef _ERROR_
    fprintf(stderr,"numMatches in getMatrixLLMatchAtRow is less than 0\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return NULL;
  }
  #endif
  if(*numMatches==0){
    matches=NULL;
    return matches;
  }

  printf("LL print numMatches %d\n",*numMatches);
  matches = newMatrix(*numMatches,1);
  
  int index;
  int seq;
  seq = 1;
  index = 1;
  LL_MNode node = LL->start;
  while( node!=NULL){
    
    if(match==getE(node->data,R,1)){
      setE(matches,index,1,seq);
      printf("matchesRec %g\n",getE(matches,index,1));
      index++;
    }
    node=node->next;
    seq++;
  }
  
  return matches;
}

int getMatrixLLNumberElemGreaterThanMatchAtRow(matrix_linklist LL, double match, int R){

  if(LL==NULL){
    fprintf(stderr,"ERROR LL is NULL in getMatrixLLNumberElemGreaterThanMatchAtRow\n");
  }
  
  int rows;
  LL_MNode node = LL->start;
  int count = 0;
  while(node!=NULL){

    rows = getRows(node->data);
    
    if(R>rows){
      fprintf(stderr,"ERROR getMatrixLLNumberMatchAtRow data data only\n");
      fprintf(stderr," has %d rows but you have requested\n",rows);
      fprintf(stderr," a match with row %d.\n",R);
      return -1;
    }
    
    if(getE(node->data,R,1)>match){
      count++;
    }
    node=node->next;
  }
    
  return count;
}

double getMatrixLLNodeElem(matrix_linklist LL, int seq, int R){

  if(LL==NULL){
    fprintf(stderr,"ERROR LL is NULL in getMatrixLLNodeElem\n");
    exit(1);
  }
  if(seq<=0){
    fprintf(stderr,"ERROR seq is less than or equal to 0 in getMatrixLLNodeElem\n");
    exit(1);
  }
  if(seq>LL->length){
    fprintf(stderr,"ERROR seq is greater than the LL length in getMatrixLLNodeElem\n");
    exit(1);
  }
  if(R<1){
    fprintf(stderr,"ERROR R value is less than 0\n");
    exit(1);
  }
 
  int rows;
  LL_MNode node = LL->start;
  int sequence = 0;
  while(node!=NULL){
    sequence++;
    rows = getRows(node->data);
    if(sequence==seq){
      if(R>rows){
        fprintf(stderr,"ERROR getMatrixLLNodeElem R is larger than number of rows\n");
        fprintf(stderr,"total rows in matrix %d you have requested row %d\n",rows,R);
        exit(1);
      }
      return getE(node->data,R,1);
    }
    node=node->next;
  }
  return -1;    
}

int moveMatrixLLNodeToStart(matrix_linklist LL, int seq){
 
 #ifdef _ERROR_CHECKING_ON_
 if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToStart LL is NULL");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(seq<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToStart seq is less than 1\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(seq>getMatrixLLlength(LL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToStart seq is greater than ");
    fprintf(stderr,"length of the LL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(LL->start==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToStart LL->start is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  //If moving the first node there is nothing to do
  if(seq==1){
    return 0;
  }

  LL_MNode temp = NULL;
  LL_MNode begin = NULL;
  LL_MNode node = NULL;
  temp = LL->start;
  int i;
  i = 2;
  while(i<=getMatrixLLlength(LL)){
    begin = temp;
    temp = temp->next;
    if(i==seq){
      break;
    }
    i++;
  }
  if(i==getMatrixLLlength(LL)){
    begin->next = NULL; 
    node = temp;
    node->next = LL->start;
    LL->start=node;
  }else{
    /* Skipping center node (temp)     */
    begin->next = temp->next;
    /* Assigning the node the node will*/
    /* point to as the previous start  */
    temp->next = LL->start;
    /* Setting the beginning of the ll */
    /* to node                         */
    LL->start=temp;

  }
  return 0;
}

int moveMatrixLLNodeElemToEnd(matrix_linklist LL, int seq){
  if(LL==NULL){
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToEnd LL is NULL");
    return -1;
  }
  if(seq<1){
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToEnd seq is less than 1\n");
    return -1;
  }
  if(seq>getMatrixLLlength(LL)){
    fprintf(stderr,"ERROR moveMatrixLLNodeElemToEnd seq is greater than ");
    fprintf(stderr,"length of the LL\n");
    return -1;
  }
  //This means the node is already at the end so do nothing
  if(seq==getMatrixLLlength(LL)){
    return 0;
  }

  int i;
  LL_MNode node;
  LL_MNode temp = LL->start;
  //If the node at the front is being moved to the back
  if(seq==1){
    
    node = LL->start;
    LL->start = temp->next;
    while(temp->next!=NULL){
      temp = temp->next;
    }
    temp->next = node;
    node->next = NULL;
  }else{
    LL_MNode before;
    for(i=1;i<seq;i++){
      temp = temp->next;
    }
    before = temp;
    node = temp->next;
    before->next = node->next;
    while(temp->next!=NULL){
      temp = temp->next;
    }
    temp->next = node;
    node->next = NULL;
  }
  return 0;
}

int exchangeMatrixLLNodes(matrix_linklist LL, int seq1, int seq2){

  if(LL==NULL){
    fprintf(stderr,"ERROR exchangeMatrixLLNodes NULL LL\n");
    return -1;
  }
  if(seq1<1){
    fprintf(stderr,"ERROR exchangeMatrixLLNodes seq1 less than 1");
    return -1;
  }
  if(seq2<1){
    fprintf(stderr,"ERROR exchangeMatrixLLNodes seq2 less than 1");
    return -1;
  }
  if(seq1>LL->length){
    fprintf(stderr,"ERROR seq1 greater than size of LL\n");
    return -1;
  }
  if(seq2>LL->length){
    fprintf(stderr,"ERROR seq2 greater than size of LL\n");
    return -1;
  }

  //Nothing to do
  if(seq1==seq2){
    return 0;
  }
 
  LL_MNode node1;
  LL_MNode node2;
  LL_MNode temp;
  LL_MNode begin;
  LL_MNode after;
  LL_MNode begin1;
  LL_MNode after1;
  LL_MNode begin2;
  LL_MNode after2;
  int i;

  if(seq1==1){
    node1 = LL->start;
    temp  = LL->start;
    
    for(i=1;i<=getMatrixLLlength(LL);i++){
      if(i==seq2){
        node2 = temp;
        break;
      }
      begin = temp;
      temp = temp->next;
    }
    if(abs(seq1-seq2)>1){
      after = node2->next;
      LL->start = node2;
      node2->next = node1->next;
      begin->next = node1;
      node1->next = after;
    }else{
      LL->start = node2;
      node1->next = node2->next;
      node2->next = node1;
    }

  }else if(seq2==1){
    node2=LL->start;
    for(i=1;i<=getMatrixLLlength(LL);i++){
      if(i==seq2){
        node1 = temp;
        break;
      }
      begin = temp;
      temp  = temp->next;
    }
    if(abs(seq1-seq2)>1){
      after = node1->next;
      LL->start   = node1;
      node1->next = node2->next;
      begin->next = node2;
      node2->next = after;
    }else{
      LL->start = node1;
      node2->next = node1->next;
      node1->next = node2;
    }
  }else{
    
    temp  = LL->start;
    node1 = LL->start;
    node2 = LL->start;
    for(i=1;i<=getMatrixLLlength(LL);i++){
      if(i==seq1){
        begin1 = begin;
        node1  = temp;
      }else if(i==seq2){
        begin2 = begin;
        node2  = temp;
      }
      if(i>seq1 && i>seq2){
        break;
      }
      begin = temp;
      temp  = temp->next;
    }

    after1 = NULL;
    after2 = NULL;

    if(abs(seq1-seq2)>1){
      begin1->next = node2;
      begin2->next = node1;
      after1       = node2->next;
      after2       = node1->next;
      node1->next  = after2;
      node2->next  = after1;
    }else if (seq1<seq2){
      begin1->next = node2;
      node1->next  = after2;
      node2->next  = node1; 
    }else{
      begin2->next = node1;
      node2->next  = after1;
      node1->next  = node2;
    }
  }
 return 0; 
}

int incMatrixLLElem(matrix_linklist LL, int seq, int row){
  
  #ifdef _ERROR_CHECKING_ON_
  if(LL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR matrix_linklist in incMatrixLLElem is NULL.\n");
    #endif
    #ifdef _FORCE_CRASH_ON_
    exit(1);
    #endif
    return -1;
  }
  if(seq<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR seq is less than on in incMatrixLLElem.\n");
    #endif
    #ifdef _FORCE_CRASH_ON_
    exit(1);
    #endif
    return -1;
  }
  if(seq>getMatrixLLlength(LL)){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR seq is larger than matrix linklist in ");
    fprintf(stderr,"getMatrixLLlength\n");
    #endif
    #ifdef _FORCE_CRASH_ON_
    exit(1);
    #endif
    return -1;
  }
  #endif
  int i;
  LL_MNode temp;
  temp = LL->start;
  for(i=1;i<=getMatrixLLlength(LL);i++){
    if(i==seq){
      break;
    }
    temp=temp->next;
  }
  matrix mat;
  mat = temp->data;
  incrementE(&mat,row,1);
  return 0;
}
