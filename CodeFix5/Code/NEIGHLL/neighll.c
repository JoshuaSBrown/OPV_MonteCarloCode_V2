#include <stdio.h>
#include <stdlib.h>

#include "neighll.h"
#include "../NEIGHNODE/neighnode.h"
#include "../ERROR/error.h"

struct _NeighLL{
	int numNeigh;
	NeighNode start;
};

////////////////////////////////////////////////////////////////////////
NeighLL newNeighLL(void){

	NeighLL neiLL = (NeighLL) malloc(sizeof(struct _NeighLL));

	if(neiLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR malloc returned NULL in newNeighLL\n");
    #endif
		return NULL;
	}

	neiLL->numNeigh=0;
	neiLL->start=NULL;
	return neiLL;
}

int deleteNeighLL(NeighLL neighLL){
	if(neighLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR neighLL in deleteNeighLL is NULL\n");
    #endif
		return -1;
	}
	free(neighLL);
	return 0;
}

int deleteNeighLLAll(NeighLL neighLL){
	if(neighLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR neighLL in deleteNeighLLAll is NULL\n");
    #endif
		return -1;
	}

	NeighNode NeighNod;
	NeighNode tempNeighNod;
	NeighNod = neighLL->start;

	if(NeighNod!=NULL){
		while(getNextNeigh(NeighNod)!=NULL){
			tempNeighNod=getNextNeigh(NeighNod);
			setNextNeighNodeToNULL(&NeighNod);
			deleteNeighNode(NeighNod);
			NeighNod=tempNeighNod;
		}
		setNextNeighNodeToNULL(&NeighNod);
		deleteNeighNode(NeighNod);
		neighLL->start=NULL;
	}

	deleteNeighLL(neighLL);
	return 0;
}

int printNeighLL(const_NeighLL neighLL){

	if(neighLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR neighLL in printNeighLL is NULL\n");
    #endif
		return -1;
	}

	printf("NeighLL size %d\n",neighLL->numNeigh);

	NeighNode NeighNod;
	NeighNod = neighLL->start;

	if(NeighNod!=NULL){
		while(getNextNeigh(NeighNod)!=NULL){
			printNeighNode(NeighNod);
			NeighNod = getNextNeigh(NeighNod);
		}
		printNeighNode(NeighNod);
	}
	return 0;
}

int getNeighLL_numNeigh( NeighLL Nei){
	if(Nei==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nei in getNeighLL_numNeigh is NULL\n");
    #endif
		return -1;
	}
	return Nei->numNeigh;
}

int setNeighLL_start(NeighLL Nei, NeighNode NeighNod){
	if(Nei==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nei in setNeighLL_start is NULL\n");
    #endif
		return -1;
	}
  if(NeighNod==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeighNod in setNeighLL_start is NULL\n");
    #endif
    return -1;
  }
  if(Nei->start!=NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nei->start in setNeighLL_start is ");
    fprintf(stderr,"not NULL\n");
    #endif
    return -1;
  }
	Nei->start = NeighNod;
	Nei->numNeigh++;
  return 0;
}

NeighNode getNeighLL_start(NeighLL neighll){
	if(neighll==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR neighll in getNeighLL_start is NULL\n");
    #endif
		return NULL;
	}

  return neighll->start;
}

int setNeighLL_addNeighNode(NeighLL NeiLL, NeighNode Nei){
  if(Nei==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR Nei in setNeighLL_addNeighNode is NULL\n");
    #endif
    return -1;
  }
  if(NeiLL==NULL){
    #ifdef _ERROR_
    fprintf(stderr,"ERROR NeiLL in setNeighLL_addNeighNode is NULL\n");
    #endif
    return -1;
  }
 
  NeighNode temp;
  temp = NeiLL->start;
  if(temp==NULL){
    NeiLL->start=Nei; 
  }else{
    if(getNeighNode_id(temp)==getNeighNode_id(Nei)){
      return -2;
    }
    while(getNextNeigh(temp)!=NULL){
      temp = getNextNeigh(temp);
      if(getNeighNode_id(temp)==getNeighNode_id(Nei)){
        //Cannot add a node of the same id
        return -2;
      }
    }
    setNextNeighNode(&temp,&Nei);
  }
  NeiLL->numNeigh++;
  return 0;
}
