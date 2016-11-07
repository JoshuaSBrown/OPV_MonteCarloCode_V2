#include <stdio.h>
#include <stdlib.h>

#include "neighnode.h"
#include "../ERROR/error.h"

struct _NeighNode{
	int id;
	int hoplength;
	NeighNode next;
	Hop start;
};

struct _Hop{
  double p;
  double t;
  Hop next;
};

////////////////////////////////////////////////////////////////////////
//Tools for accessing Neigh Nodes
NeighNode newNeighNode(int N_ID){
	if(N_ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR N_ID is less than 0 in newNeighNode\n");
    #endif
		return NULL;
	}
	NeighNode NeighNod = (NeighNode) malloc(sizeof(struct _NeighNode)); 

	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc in newNeighNode returned NULL\n");
    #endif
		return NULL;
	}
	NeighNod->next=NULL;
	NeighNod->id=N_ID;
	NeighNod->start=NULL;
	NeighNod->hoplength=0;
	return NeighNod;
}

int deleteNeighNode(NeighNode NeighNod){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod in deleteNeighNode is NULL\n");
    #endif
		return -1;
	}

	if(NeighNod->start!=NULL){
		if(NeighNod->hoplength>0){
			printf("Hop length %d\n",NeighNod->hoplength);
			deleteAllHop(NeighNod->start);
			NeighNod->start=NULL;
			NeighNod->hoplength=0;
		}
	}

	free(NeighNod);
	return 0;
}

int printNeighNode(const_NeighNode NeighNod){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod in printNeighNode is NULL\n");
    #endif
		return -1;
	}
	printf("NeighNode id %d\n HopLength %d\n",NeighNod->id,NeighNod->hoplength);

	Hop h = NeighNod->start;
	while(h!=NULL){
		printHop(h);
		h = h->next;
	}
	return 0;
}

int getNeighNode_id(const_NeighNode NeighNod){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in getNeighNode_id\n");
    #endif
		return -1;
	}
	return NeighNod->id;
}

int setNextNeighNode(NeighNode * NeighNod,NeighNode * NeighNod2){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNextNeighNode\n");
    #endif
		return -1;
	}
  if(NeighNod2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod2 is NULL in setNextNeighNode\n");
    #endif
    return -1;
  }
	if(*NeighNod==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *NeighNod in setNextNeighNode is NULL \n");
    #endif
		return -1;
	}
  if(*NeighNod2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *NeighNod2 is NULL in setNextNeighNode\n");
    #endif
    return -1;
  }

	(*NeighNod)->next=(*NeighNod2);
 
	return 0;
}

int setNextNeighNodeToNULL(NeighNode * NeighNod){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNextNeighNodeToNULL\n");
    #endif
		return -1;
	}
	if(*NeighNod==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *NeighNod in setNextNeighNodeToNULL is NULL \n");
    #endif
		return -1;
	}

  (*NeighNod)->next=NULL;
  return 0;
}

NeighNode getNextNeigh(const_NeighNode NeighNod){
	if(NeighNod==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR *NeighNod is NULL in getNextNeigh\n");
    #endif
		return NULL;
	}
  //It is ok for the function to return NULL if the ->next Neigh
  //Node is NULL
	return NeighNod->next;
}

int setNeighNodeNew_p(NeighNode NeighNod, double pval){
	if (NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNeighNodeNew_p\n");
    #endif
		return -1;
	}
  if(pval<0.0){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pval is less than 0 in setNeighNodeNew_p\n");
    #endif
    return -1;
  }
	if (NeighNod->start==NULL){
		NeighNod->start = newHop();
		NeighNod->start->p=pval;
		NeighNod->hoplength++;
	}else{
		Hop h = NeighNod->start;
		while(h->next!=NULL){
			h = h->next;
		}
		h->next = newHop();
		h = h->next;
		h->p = pval;
		NeighNod->hoplength++;
	}
	return 0;
}

int setNeighNode_p(NeighNode NeighNod, double pval, int Elem){

	if (NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNeighNode_p\n");
    #endif
		return -1;
	}
  if(pval<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR pval is less than 0.0 in setNeighNode_p\n");
    #endif
    return -1;
  }
	if (NeighNod->start==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod->start is NULL in setNeighNode_p\n");
    #endif
		return -1;
	}
  if(Elem>NeighNod->hoplength){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is greater than NeighNod->hoplength in setNeighNode_p\n");
    #endif
    return -1;
  }
  if(Elem<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is less than 1 in setNeighNode_p\n");
    #endif
    return -1;
  }
	int inc=1;
	Hop h = NeighNod->start;
  if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setNeighNode_p\n");
    #endif
    return -1;
  }

	while (inc<Elem){
		inc++;
		h = h->next;
    if(h==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR h is NULL in setNeighNode_p\n");
      #endif
      return -1;
    }
	}
	h->p = pval;

	return 0;
}

int setNeighNode_t(NeighNode NeighNod, double time, int Elem){

	if (NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNeighNode_t\n");
    #endif
		return -1;
	}
  if(time<0.0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR time is less than 0.0 in setNeighNode_t\n");
    #endif
    return -1;
  }
	if (NeighNod->start==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod->start==NULL in setNeighNode_t\n");
    #endif
    return -1;
	}
  if(Elem>NeighNod->hoplength){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is greater than NeighNod->hoplength in setNeighNode_t\n");
    #endif
    return -1;
  }
  if(Elem<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is less than 1 in setNeighNode_t\n");
    #endif
    return -1;
  }
	int inc=1;
	Hop h = NeighNod->start;
  if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setNeighNode_t\n");
    #endif
    return -1;
  }
	while (inc<Elem){
		inc++;
		h = h->next;
    if(h==NULL){
      #ifdef _ERROR_
      fprintf(stderr,"ERROR h is NULL in setNeighNode_t\n");
      #endif
      return -1;
    }
	}
	h->t = time;

	return 0;
}

double getNeighNode_p(const_NeighNode NeighNod, int Elem){

	if (NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in getNeighNode_p\n");
    #endif
		return -1;
	}
  if (NeighNod->start==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod->start==NULL in getNeighNode_p\n");
    #endif
    return -1;
	}
  if(Elem>NeighNod->hoplength){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is greater than NeighNod->hoplength in getNeighNode_p\n");
    #endif
    return -1;
  }
  if(Elem<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is less than 1 in getNeighNode_p\n");
    #endif
    return -1;
  }
	int inc=1;
	Hop h = NeighNod->start;

	while (inc<Elem){
		inc++;
		h = h->next;
	}
	return h->p;
}

int getNeighNode_hoplength(NeighNode NeighNod){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in getNeighNode_hoplength\n");
    #endif
		return -1;
	}

	return NeighNod->hoplength;
}

double getNeighNode_t(NeighNode NeighNod, int Elem){

	if (NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in getNeighNode_t\n");
    #endif
		return -1;
	}
  if (NeighNod->start==NULL){ 
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod->start==NULL in getNeighNode_t\n");
    #endif
    return -1;
	}
  if(Elem>NeighNod->hoplength){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is greater than NeighNod->hoplength in getNeighNode_t\n");
    #endif
    return -1;
  }
  if(Elem<1){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Elem is less than 1 in getNeighNode_t\n");
    #endif
    return -1;
  }

	int inc=1;
	Hop h = NeighNod->start;
	while (inc<Elem){
		inc++;
		h = h->next;
	}
	return h->t;
}

int setNeighNode_id(NeighNode NeighNod, int ID){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod is NULL in setNeighNode_id\n");
    #endif
		return -1;
	}
  if(ID<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR ID is less than 0 in setNeighNode_id\n");
    #endif
    return -1;
  }
	NeighNod->id=ID;
	return 0;
}

int setNeighNode_hopstart(NeighNode NeighNod, Hop h){
	if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNOd is NULL in setNeighNode_hopstart\n");
    #endif
		return -1;
	}
  if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setNeeghNode_hopstart\n");
    #endif
    return -1;
  }
  if(NeighNod->start!=NULL){
     #ifdef _ERROR_
     fprintf(stderr,"ERROR start is already defined for NeighNod you");
     fprintf(stderr," must first delete the hops that are present ");
     fprintf(stderr,"before you can reinitiliaze it in setNeighNode_hopstart\n");
     #endif
     return -1;
  }
	NeighNod->start = h;
  NeighNod->hoplength++;
	return 0;
}

////////////////////////////////////////////////////////////////////////
//Tools for accessing Hop
Hop newHop(void){

	Hop h = (Hop) malloc(sizeof(struct _Hop));

	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in newHop\n");
    #endif
		return NULL;
	}

	h->t=0.0;
	h->p=0.0;
	h->next=NULL;
	return h;
}

int deleteAllHop(Hop h){//come back to deleteNeighNode

	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in deleteAllHop \n");
    #endif
		return -1;
	}
	Hop temp;
	while( h->next!=NULL){
		temp = h->next;
		h->next=NULL; //why necessary?
		free(h);
		h = temp;
	}
	h->next=NULL; //necessary condition, why restated here?
	free(h);

	return 0;
}

int printHop(Hop h){

	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in printHop\n");
    #endif
		return -1;
	}

	printf("Hop time %g pval %g\n",h->t,h->p);
	return 0;

}

int setHop_t(Hop h, double t){
	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setHop_t\n");
    #endif
		return -1;
	}
  if(t<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR t is less than 0 in setHop_t \n");
    #endif
    return -1;
  }
	h->t = t;
	return 0;
}

int setHop_p(Hop h, double p){
	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setHop_p\n");
    #endif
		return -1;
	}
  if(p<0){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR p is less than 0 in setHop_p\n");
    #endif
    return -1;
  }
	h->p = p;
	return 0;
}

int setHop_next(Hop h, Hop h2){
	if(h==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h is NULL in setHop_next\n");
    #endif
		return -1;
	}
  if(h2==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR h2 is NULL in setHop_next\n");
    #endif
    return -1;
  }
	h->next = h2;
	return 0;
}	

