
int readPath(char * FileName){

	char buf[256];
	int PFget_LogFile(ParameterFrame PF);

	snprintf(buf, sizeof buf,"%s%s",FileName,".path");
	
	FILE * PathIn;

	if((PathIn = fopen(buf,"r"))==NULL){
		if((PathIn = fopen(FileName,"r"))==NULL){
			printf("ERROR File %s and with extension %s do not exist!\n",FileName,buf);
			return -1;
		}else{
			
		}
	}else{

	}
}

int getPathCharge(int ID){

}

int getTrap(int ID){

}

int getPerc(){

}

int writePerc(){

}

int writeTrap(){

}
