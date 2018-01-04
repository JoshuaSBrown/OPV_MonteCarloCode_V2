
#include <stdio.h>
#include <stdlib.h>

#include "linklist2.h"

struct _linklist2 {
	int length;
	LLNode2 start;
};

struct _LLNode2 {
	int id;
	LLNode2 next;
};

LLNode2 newLLNode2( int id){
	
	LLNode2 node = (LLNode2) malloc(sizeof(struct _LLNode2));

	if(node==NULL){
		return NULL;
	}

	//Set id
	node->id = id;
	node->next = NULL;
	return node;
}

linklist2 newLinkList2( int id ){

	linklist2 LL = (linklist2) malloc(sizeof(struct _linklist2));

	if( LL == NULL){
		return NULL;
	}

	LLNode2 node = newLLNode2( id);
	
	if(node==NULL){
		deleteLL2( &LL);
		return NULL;
	}

	LL->length = 1;
	LL->start = node;

	return LL;
}

linklist2 newBlankLinkList2(){
	
	linklist2 LL = (linklist2) malloc(sizeof(struct _linklist2));

	if( LL == NULL){
		return NULL;
	}
	
	LL->length=0;
	LL->start = NULL;

	return LL;
}

//Delete entire linklist
int deleteLL2(linklist2 * LL){ 
	if(LL==NULL){
		return -1;
	}

	LLNode2 tempNode;
	LLNode2 tempNode2;

	tempNode = (*LL)->start;
	while (tempNode!=NULL){
		tempNode2=tempNode;
		tempNode = tempNode->next;
		deleteLLNode2(tempNode2);
	}

	free((*LL));

	*LL = NULL;
	return 0;
}

//deleting separate node, not part of linked list
int deleteLLNode2( LLNode2 node){ 
	if(node==NULL){
		return -1;
	}

	free(node);
	return 0;
}

LLNode2 nextLLNode2(LLNode2 node){
	if(node==NULL){
		return NULL;
	}

	return node->next;
}

int addLLNode2(linklist2 LL, int id){
	if(LL==NULL ){
		return -1;
	}

	//Cycle through linklist2
	LLNode2 tempNode;
	LLNode2 oldNode;
	LLNode2 newNode;
	tempNode = LL->start;
  //This means the first node in the LL has not
  //yet been assigned
  if (tempNode==NULL){
		newNode = newLLNode2(id);
		if(newNode==NULL){
			return -1;
		}
		LL->length++;
		LL->start = newNode;
		return 0;
	}

	while(tempNode!=NULL){
    //Searching for the last node in the LL
		oldNode = tempNode;
		tempNode = tempNode->next;
	}
	
  //Appending a new node to the last position in the LL
	newNode = newLLNode2(id);
	if(newNode==NULL){
		return -1;
	}
	LL->length++;
	oldNode->next = newNode;
	return 0;
}

//deleting one node from the linked list
//based on it's sequence in the LL
int removeLLNode2(linklist2 LL, int seq){ 	
  
  if(LL==NULL || seq<=0){
		return -1;
	}
    
  if(seq>LL->length){
    return -1;
  }
  
	LLNode2 node = LL->start;
	LLNode2 oldnode = node;
	LLNode2 newnode;
  int count;
  count = 1;
  //This case is triggered if the first node in the
  //LL is the one we are removing
	if (seq == count){
		newnode = node->next;
		LL->start=newnode;
		deleteLLNode2(node);
		LL->length--;
		return 0;
	}

  //This is if we are searching for a node 
  //that is not the first in the LL, to remove. 	
	while(node!=NULL ){
	
    count++;	
		if(count==seq){
			newnode = node->next;
			oldnode->next = newnode;
			deleteLLNode2(node);
			LL->length--;
			//Found node deleted it and 
			//closed the linklist2
			return 0;
		}
		oldnode = node;
		node = node->next;
	}

	return -1;
}

//Deleting nodes between and including those defined 
//by seq1 and seq2, seq1 must be less than seq2
int removeLLNodes2(linklist2 LL, int seq1, int seq2){

  if(LL==NULL){
    printf("ERROR Linklist2 is null in removeLLNodes2\n");
    return -1;
  } 
  if(seq1>=seq2){
    printf("ERROR seq1 is greater than or equal to seq2 in removeLLNodes2\n");
    return -1;
  }
  if(seq1<=0){
    printf("ERROR seq1 is less than or equal to 0 in removeLLNodes2\n");
    return -1;
  }
  if(seq2>LL->length){
    printf("ERROR seq2 is greater than the length of the linklist in removeLLNodes2\n");
    return -1;
  }

  LLNode2 node = LL->start;
  LLNode2 oldnode = node;
  LLNode2 newnode;
  int count;
  count = 1;
  //This case is triggered if the first node in the
  //LL is the one we are removing
  if (seq1 == count){
    while(count<=seq2){
      newnode = node->next;
      deleteLLNode2(node);
      node = newnode;
      count++;
      LL->length--;
    }
    LL->start=newnode;
    return 0;
  }

  //This is if we are searching for nodes  
  //that are not the first in the LL, to remove. 	
  while(node!=NULL ){
    
    if(count==seq1){
      while(count<=seq2){
        newnode = node->next;
        deleteLLNode2(node);
        node = newnode;
        LL->length--;
        //Found node deleted it and 
        //closed the linklist2
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

int getLLlength2(linklist2 LL){
  if(LL==NULL){
    return -1;
  }

  return LL->length;
}

int getLLstartID2(linklist2 LL){
  if(LL==NULL){
    return -1;
  }

  return LL->start->id;

}

int printLL2(linklist2 LL){
  if(LL==NULL){
    return -1;
  }

	printf("\nLength %d\n",LL->length);
	LLNode2 node = LL->start;
  int count;
  count = 1;
	while(node!=NULL){
		printf("Sequence %d \t Node %d\n",count,node->id);
		node=node->next;
    count++;
	}

	return 0;
}

int getLLLastMatch2(linklist2 LL, int match){
  
  if(LL==NULL){
    printf("ERROR LL is NULL in getLLLastMatch\n");
  }
  LLNode2 node = LL->start;
  int seq = 0;
  int lastseq = -1;
  while(node!=NULL){
    seq++;
    if(node->id==match){
      lastseq = seq;
    }
    node=node->next;
  }
    
  return lastseq;
}

int getLLNumberMatch2(linklist2 LL, int match){

  if(LL==NULL){
    printf("ERROR LL is NULL in getLLNumberMatch2\n");
  }
  
  LLNode2 node = LL->start;
  int count = 0;
  while(node!=NULL){
    if(node->id==match){
      count++;
    }
    node=node->next;
  }
    
  return count;
}

int getLLNodeID2(linklist2 LL, int seq){

  if(LL==NULL){
    printf("ERROR LL is NULL in getLLNodeID2\n");
    exit(1);
  }
  if(seq<=0){
    printf("ERROR seq is less than or equal to 0 in getLLNodeID2\n");
    exit(1);
  }
  if(seq>LL->length){
    printf("ERROR seq is greater than the LL length in getLLNodeID2\n");
    exit(1);
  }
 
  LLNode2 node = LL->start;
  int sequence = 0;
  while(node!=NULL){
    sequence++;
    if(sequence==seq){
      return node->id;
    }
    node=node->next;
  }
  return -1;    
}
