#include <stdio.h>
#include <stdlib.h>

#include "node.h"
#include "../ERROR/error.h"

struct _Node{
	int id;
	double p; 
	Node next;
	int flag[6];
};

////////////////////////////////////////////////////////////////////////
//Tools for accessing Nodes
Node newNode(int N_ID){

  #ifdef _ERROR_CHECKING_ON_
	if(N_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR N_ID is less than 0 in newNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif

	Node Nod = (Node) malloc(sizeof(struct _Node));

  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is Null when calling malloc in newNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	Nod->id=N_ID;
	Nod->p=0;
	Nod->flag[1]=0;
	Nod->flag[0]=0;
	Nod->flag[2]=0;
	Nod->flag[3]=0;
	Nod->flag[5]=0;
	Nod->flag[4]=0;
	Nod->next=NULL;
	return Nod;
}

int deleteNode(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL when calling deleteNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	free(Nod);
	return 0;
}

int deleteNodeAll(Node * Nod){
  #ifdef _ERROR_CHECKING_ON_
 	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL when calling deleteNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
  Node temp = *Nod;
  Node old;
  while(temp->next!=NULL){
    old = temp;
    temp = temp->next;
    deleteNode(old);
  }
  deleteNode(temp);
  *Nod=NULL;
  return 0;
}

int printNode(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Node in printNode is NULL\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	printf("Node id %d\n",Nod->id);
	printf("Pval %g\n",Nod->p);
	printf("Fro \t Beh \t Lef \t Rig \t Abo \t Bel\n");
	printf("%d \t %d \t %d \t %d \t %d \t %d\n",\
			Nod->flag[1], Nod->flag[0], Nod->flag[2], Nod->flag[3],\
			Nod->flag[5], Nod->flag[4]);
	return 0;
}

int getNode_id(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getNode_id\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->id;
}

int setNextNode(Node * Nod, Node Nod2){
  #ifdef _ERROR_CHECKING_ON_
  if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setNextNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if((*Nod)==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *Nod is NULL in setNextNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if(Nod2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod2 is NULL in setNextNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  if((*Nod)->next!=NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod->next is not NULL in setNextNode");
    fprintf(stderr," you would over write the current datastructure\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
  }
  #endif

  (*Nod)->next = Nod2;
  return 0;
}

Node getNextNode(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getNextNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
	return Nod->next;
}

int setNode_id(Node Nod, int ID){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setNode_id\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ID is less than 0 when calling setNode_id\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	Nod->id=ID;
	return 0;
}

double getNode_p(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getNode_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1.0;
	}
  #endif
	return Nod->p;
}

Node getNode_lastNode(Node Nod){

  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getNode_lastNode\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return NULL;
	}
  #endif
  
  while(Nod->next!=NULL){
    Nod=Nod->next;
  }
  return Nod;
}

int setNode_p(Node Nod, double pval){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"Nod is NULL in setNode_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(pval<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pval is less than 0.0 in setNode_p\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	Nod->p=pval;
	return 0;
}

int getFlagFro(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL is getFlagFro\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[1];
}

int getFlagBeh(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlagBeh\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[0];
}

int getFlagLef(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlagLef\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[2];
}

int getFlagRig(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlagRig\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[3];
}

int getFlagAbo(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlagAbo\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[5];
}

int getFlagBel(const_Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlagBel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	return Nod->flag[4];
}

int getFlag(const_Node Nod, int index){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in getFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(index<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR index is less than 0 in getFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(index>5){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR index is greater than 5 in getFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
  if(index>5){
    return -1;
  }
	return Nod->flag[index];
}

int setFlagFro(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Node is NULL in setFlagFro\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
	}
  #endif
	Nod->flag[1]=1;
	return 0;
}

int setFlagBeh(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlagBeh\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Nod->flag[0]=1;
	return 0;
}

int setFlagLef(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlagLef\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Nod->flag[2]=1;
	return 0;
}

int setFlagRig(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlagRig\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Nod->flag[3]=1;
	return 0;
}

int setFlagAbo(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlagAbo\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Nod->flag[5]=1;
	return 0;
}
int setFlagBel(Node Nod){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlagBel\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  #endif
	Nod->flag[4]=1;
	return 0;
}

int setFlag(Node Nod, int index){
  #ifdef _ERROR_CHECKING_ON_
	if(Nod==NULL){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR Nod is NULL in setFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
		return -1;
	}
  if(index<0){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR index is less than 0 in setFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  if(index>5){
		#ifdef _ERROR_
    fprintf(stderr,"ERROR index is greater than 5 in setFlag\n");
    #endif
    #ifdef _FORCE_HARD_CRASH_
    exit(1);
    #endif
    return -1;
  }
  #endif
	Nod->flag[index]=1;
	return 0;
}


