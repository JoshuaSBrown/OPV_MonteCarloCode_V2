
#include <stdio.h>
#include <stdlib.h>

#include "linklist.h"
//#include "../../../MEM/mem.h"

struct _linklist {
	int length;
	LLNode start;
};

struct _LLNode {
	int id;
  LLNode next;
};

LLNode newLLNode( int id){
	if( id<0 ){
		return NULL;
	}

	LLNode node = (LLNode) malloc(sizeof(struct _LLNode));

	if(node==NULL){
		return NULL;
	}

	//Set id
	node->id = id;
	node->next = NULL;

	return node;
}

linklist newLinkList( int id ){

	if(id<0){
		return NULL;
	}

	linklist LL = (linklist) malloc(sizeof(struct _linklist));

	if( LL == NULL){
		return NULL;
	}

	LLNode node = newLLNode( id);
	
	if(node==NULL){
		deleteLL( &LL);
		return NULL;
	}

	LL->length = 1;
	LL->start = node;

	return LL;
}

linklist newBlankLinkList(){
	
	linklist LL = (linklist) malloc(sizeof(struct _linklist));

	if( LL == NULL){
		return NULL;
	}
	
	LL->length=0;
	LL->start = NULL;

	return LL;
}

int deleteLL(linklist * LL){ //deleting entire linked list
	if(LL==NULL){
		return -1;
	}

	LLNode tempNode;
	LLNode tempNode2;

	tempNode = (*LL)->start;
	while (tempNode!=NULL){
		tempNode2=tempNode;
		tempNode = tempNode->next;
		deleteLLNode(tempNode2);
	}

	free((*LL));

	*LL = NULL;
	return 0;
}

int deleteLLNode( LLNode node){ //deleting separate node, not part of linked list
	if(node==NULL){
		return -1;
	}

	free(node);
	return 0;
}

LLNode nextLLNode(LLNode node){
	if(node==NULL){
		return NULL;
	}

	return node->next;
}

int addLLNode(linklist LL, int id){
	if(LL==NULL || id<0){
		return -1;
	}

	//Cycle through linklist
	LLNode tempNode;
	LLNode oldNode;
	LLNode newNode;
	tempNode = LL->start;
	if (tempNode==NULL){
		newNode = newLLNode(id);
		if(newNode==NULL){
			return -1;
		}
		LL->length++;
		LL->start = newNode;
		return 0;
	}

	while(tempNode!=NULL){

		if(tempNode->id==id){
			return -1;
			//This means a node with this id
			//is already in the link list not
			//allowed
		}
		oldNode = tempNode;
		tempNode = tempNode->next;
	}
	
	newNode = newLLNode(id);
	if(newNode==NULL){
		return -1;
	}
	LL->length++;
	oldNode->next = newNode;
	return 0;
}

int removeLLNode(linklist LL, int id){ //deleting one node from the linked list
	if(LL==NULL || id<0){
		return -1;
	}

	LLNode node = LL->start;
	LLNode oldnode = node;
	LLNode newnode;
	if (node->id==id){
		newnode = node->next;
		LL->start=newnode;
		deleteLLNode(node);
		LL->length--;
		return 0;
	}
	
	while(node!=NULL ){
		
		if(node->id==id){
			newnode = node->next;
			oldnode->next = newnode;
			deleteLLNode(node);
			LL->length--;
			//Found node deleted it and 
			//closed the linklist
			return 0;
		}
		oldnode = node;
		node = node->next;
	}

	return -1;
}

int getLLlength(linklist LL){
 	if(LL==NULL){
		return -1;
	}

	return LL->length;
}

int getLLstartID(linklist LL){
	if(LL==NULL){
		return -1;
	}

	return LL->start->id;

}

int printLL(linklist LL){
	if(LL==NULL){
		return -1;
	}

	printf("\nLength %d\n",LL->length);
	LLNode node = LL->start;

	while(node!=NULL){
		printf("Node %d\n",node->id);
		node=node->next;
	}

	return 0;

}
