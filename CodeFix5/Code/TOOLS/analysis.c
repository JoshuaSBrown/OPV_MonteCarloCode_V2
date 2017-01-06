#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../MATRIX/matrix.h"
#include "../LINKLIST2/linklist2.h"
#include "analysis.h"

int readPath(char * FileName){

	char buf[256];
	char bufRead[256];
  char FileName2[60];

  strncpy(FileName2,FileName,(strlen(FileName)-5));
  printf("FileName %s\nFileName2 %s\n",FileName,FileName2); 
  
  int CurrentCharge;  //The current charge being read from the
                      //.path file
  int inc;            //Increment used to determine the position 
                      //in the ChargePath matrix
  int rows;           //Track the number of rows in the ChargePath
                      //matrix so we know when to resize the matrix
  int FlagFinished;   //Determines when last charge has been read
  matrix ChargePath;

  //Properties loaded in from .path file
  int x, y, z;
  int max_y, max_z;
  int ChargeID;
  double dwelltime;
  double globaltime;
  int line;
  //Initializing variables
  CurrentCharge = 0;
  FlagFinished = 1;
  rows = 8;
  max_y = -1;
  max_z = -1;
 
	snprintf(buf, sizeof buf,"%s%s",FileName,".path");
	
  FILE * PathIn;

  if((PathIn = fopen(buf,"r"))==NULL){
    if((PathIn = fopen(FileName,"r"))==NULL){
      printf("ERROR File %s and with extension %s do not exist!\n",FileName,buf);
      return -1;
    }else{

      //File exists and we have opened it
      //We will proceed to read the information from the .path file
      //We will read one charge at a time starting with charge 0
      while(FlagFinished==1){
        //4 rows no need to store the chargeID that is stored in CurrentCharge
        //Also no need to store the global time
        inc = 0;
        ChargePath = newMatrix(rows,4);
        FlagFinished = 0;
        fgets(bufRead,256,PathIn);
        printf("Scanning in Charge %d\n",CurrentCharge);
        line =1;
        while(fscanf(PathIn,"%d %d %d %d %lf %lf",&x,&y,&z,&ChargeID,&dwelltime,&globaltime)!=EOF){
          if(y>max_y){
            max_y = y;
          }
          if(z>max_z){
            max_z = z;
          }
          printf("line %d\n",line);
          line++;
          if(ChargeID==CurrentCharge){

            FlagFinished = 1; 
            //Resize Matrix if needed
            if(inc>rows){
              rows = rows*2;
              printf("resizing Matrix\n");
              resizeRow(&ChargePath,rows);
            }

            //Store information
            setE(ChargePath,inc,1,(double)x);
            setE(ChargePath,inc,2,(double)y);
            setE(ChargePath,inc,3,(double)z);
            setE(ChargePath,inc,4,dwelltime);
            inc++;
          }
          fgets(bufRead,256,PathIn);
        }

        //resize the ChargePath matrix so there are no extra 0's
        //appended to the end
        printf("Resizing matrix part 2\n");
        resizeRow(&ChargePath,inc-1);
        
        //We have now collected all the data from the .path file
        //for the charge labeled CurrentCharge we can go ahead and 
        //analyize the data we have and append it in .trap and .perc
        //files so we are not using up to much memory
        printf("Sorting path into percolation and trap sites\n");
        if(FlagFinished ==1){
          getPercTrap(FileName2, ChargePath, CurrentCharge, max_y,max_z); 
        }
        //Return to the beginning of the file
        fseek(PathIn,0,SEEK_SET);

        deleteMatrix(&ChargePath);
        CurrentCharge++;
      }
    }
	}else{

	}
  return 0;
}

int getPercTrap(char * FileName, const_matrix ChargePath,int ChargeID, int y_max, int z_max){

  //First things first we will count up the total number of unique 
  //positions in the ChargePath and the total time spent on each
  //of them
  
  int inc;
  int R;
  int rows;
  int LLsize; 
  int i,j, k;
  double x,y,z;
  double x1,y1,z1;
  int flag;
  int flag2;
  double dwelltime;
  int SiteID;
  int SiteID2;
  int ID;
  int ID2;
  int seq;
	char buf[256];

  rows = getRows(ChargePath);
  inc = 1;
  R = 8;
  
/* Step 1: Each row of this matrix will track an individual site
 * that appears in the charge pathway file. The total time the 
 * charge spends on that site will be recorded for each site and
 * stored in the second column
 */
  matrix UniqueSites = newMatrix(R,5);

  for(i=1;i<=rows;i++){
    
    x = getE(ChargePath,i,1);
    y = getE(ChargePath,i,2);
    z = getE(ChargePath,i,3);
    SiteID = uniqueID(x,y,z,y_max,z_max); 

    flag = 0;
/* Before we proceed we will ensure that we have not already 
 * accounted for this site.
 */
    for(k=1;k<i;k++){
      x1 = getE(ChargePath,k,1);
      y1 = getE(ChargePath,k,2);
      z1 = getE(ChargePath,k,3);
      SiteID2 = uniqueID(x1,y1,z1,y_max,z_max);
      if(SiteID2==SiteID){
        flag = -1;
        break;
      }
    }

/* If the flag is 0 it means we have not yet accounted for the
 * site and we can proceed.
 */
    if(flag==0){ 
      dwelltime = getE(ChargePath,i,4);

      //Make sure that the matrix is large enough
      if(inc==R){
        R = R*2;
        resizeRow(&UniqueSites,R);
      } 
      setE(UniqueSites,inc,1,x);
      setE(UniqueSites,inc,2,y);
      setE(UniqueSites,inc,3,z);
      setE(UniqueSites,inc,4,dwelltime);
      setE(UniqueSites,inc,5,SiteID);
      for(j=(i+1);j<=rows;j++){
        x1 = getE(ChargePath,j,1);
        y1 = getE(ChargePath,j,2);
        z1 = getE(ChargePath,j,3);
        SiteID2 = uniqueID(x1,y1,z1,y_max,z_max);

        if(SiteID==SiteID2){
          dwelltime = dwelltime+getE(ChargePath,j,4);
          setE(UniqueSites,inc,4,dwelltime);
        }
      }
      inc++;
    }
  }
  
  //Finally resize the matrix so there are not any extra 0's 
  //appended to the bottom
  resizeRow(&UniqueSites,inc-1); 

/* Step 2: Now that we have a matrix with the total times on 
 * each site we will create a linklist of the charge pathways.
 * We will use linklist because it is less intensive to remove
 * elements from a linklist than from a matrix. 
 */ 
  
  linklist2 ChargePercLL;
  linklist2 ChargeTrapLL;
  ConvertMatrixLinkList(ChargePath, &ChargePercLL, y_max, z_max);
  ConvertMatrixLinkList(ChargePath, &ChargeTrapLL, y_max, z_max);
  R = 8;

/* Step 3: Here we will remove from the linklist any sites that
 * are not part of the percolation pathway. This is done by 
 * observing the reappearance of sites. 
 */
  LLsize = getLLlength2(ChargePercLL);
  for(i=1;i<=LLsize;i++){
    ID = getLLNodeID2(ChargePercLL,i); 
    flag2 = getLLNumberMatch2(ChargePercLL, ID);
    if(flag2>1){
      seq = getLLLastMatch2(ChargePercLL, ID);
      removeLLNodes2(ChargePercLL, (i+1), seq);
      LLsize = getLLlength2(ChargePercLL);
    }else if(flag2==0){
      printf("ERROR flag2 should be at lest 1 in getPerc\n");
    }
  }  

/* Step 4: we can now append this information to the .perc file
 */
  rows = getRows(ChargePath);
	snprintf(buf, sizeof buf,"%s%s",FileName,".perc");
  
  printLL2(ChargePercLL);

  FILE * PercOut;
  if(ChargeID==0){
    if((PercOut = fopen(buf,"w"))!=NULL){
        
        LLsize = getLLlength2(ChargePercLL);
        fprintf(PercOut,"%d %d\n",ChargeID,LLsize);
        for(i=1;i<=LLsize;i++){
          ID = getLLNodeID2(ChargePercLL,i);
          printf("ID of node %d\n",ID);
          //Cycle through the matrix and print the information for the
          //site if there is a match.
          for(j=1;j<=rows;j++){
            x = getE(UniqueSites,j,1);
            y = getE(UniqueSites,j,2);
            z = getE(UniqueSites,j,3);
            ID2 = getE(UniqueSites,j,5);
            printf("ID of perc site %d\n",ID2);
            if(ID==ID2){
              dwelltime = getE(UniqueSites,j,4);
              fprintf(PercOut,"%d %d %d %g\n",(int)x,(int)y,(int)z,dwelltime);
              break;
            }
          } 
        }
        fclose(PercOut);
    }else{
      printf("ERROR unable to write to file %s with extension %s!\n",FileName,buf);
      return -1;
    }
  }else{
    if((PercOut = fopen(buf,"a"))!=NULL){
        LLsize = getLLlength2(ChargePercLL);
        fprintf(PercOut,"%d %d\n",ChargeID,LLsize);
        for(i=1;i<=LLsize;i++){
          ID = getLLNodeID2(ChargePercLL, i);
          //Cycle through the matrix and print the information for the
          //site if there is a match.
          for(j=1;j<=rows;j++){
            x = getE(UniqueSites,j,1);
            y = getE(UniqueSites,j,2);
            z = getE(UniqueSites,j,3);
            ID2 = getE(UniqueSites,j,5);
            if(ID==ID2){
              dwelltime = getE(UniqueSites,j,4);
              fprintf(PercOut,"%d %d %d %g\n",(int)x,(int)y,(int)z,dwelltime);
              break;
            }
          } 
        }
        fclose(PercOut);
    }else{
      printf("ERROR unable to append to file %s with extension %s!\n",FileName,buf);
      return -1;
    }
  }

/* Step 5: Now we need to find the sites that make up the traps
 * this can be done by removing all sites in the linlist that 
 * only appear once and are also in the percolation pathway
 */
  LLsize = getLLlength2(ChargePercLL);
  for(i=1;i<=LLsize;i++){
    ID = getLLNodeID2(ChargePercLL,i); 
    flag2 = getLLNumberMatch2(ChargeTrapLL, ID);
    if(flag2==1){
    //This means it is not part of a trap so we will remove it from 
    //the trap LL
      seq = getLLLastMatch2(ChargeTrapLL, ID);
      removeLLNode2(ChargeTrapLL, seq);
    }
  }

/* Step 6: Now we are going to print the .trap file in the same
 * format as the .perc file. 
 */
  rows = getRows(ChargePath);
	snprintf(buf, sizeof buf,"%s%s",FileName,".trap");
	
  FILE * TrapOut;
  if(ChargeID==0){
    if((TrapOut = fopen(buf,"w"))!=NULL){

      LLsize = getLLlength2(ChargeTrapLL);
      fprintf(TrapOut,"%d %d\n",ChargeID,LLsize);
      for(i=1;i<=LLsize;i++){
        ID = getLLNodeID2(ChargeTrapLL,i);
        //Cycle through the matrix and print the information for the
        //site if there is a match.
        for(j=1;j<=rows;j++){
          x = getE(UniqueSites,j,1);
          y = getE(UniqueSites,j,2);
          z = getE(UniqueSites,j,3);
          ID2 = getE(UniqueSites,j,5);
          if(ID==ID2){
            dwelltime = getE(UniqueSites,j,4);
            fprintf(TrapOut,"%d %d %d %g\n",(int)x,(int)y,(int)z,dwelltime);
            break;
          }
        } 
      }
      fclose(TrapOut);
    }else{
      printf("ERROR unable to write to file %s with extension %s!\n",FileName,buf);
      return -1;
    }
  }else{
    if((TrapOut = fopen(buf,"a"))!=NULL){
        LLsize = getLLlength2(ChargeTrapLL);
        fprintf(TrapOut,"%d %d\n",ChargeID,LLsize);
        for(i=1;i<=LLsize;i++){
          ID = getLLNodeID2(ChargeTrapLL,i);
          //Cycle through the matrix and print the information for the
          //site if there is a match.
          for(j=1;j<=rows;j++){
            x = getE(UniqueSites,j,1);
            y = getE(UniqueSites,j,2);
            z = getE(UniqueSites,j,3);
            ID2 = getE(UniqueSites,j,5);
            if(ID==ID2){
              dwelltime = getE(UniqueSites,j,4);
              fprintf(TrapOut,"%d %d %d %g\n",(int)x,(int)y,(int)z,dwelltime);
              break;
            }
          } 
        }
        fclose(TrapOut);
    }else{
      printf("ERROR unable to append to file %s with extension %s!\n",FileName,buf);
      return -1;
    }
  }
  
/* Step 8: Free matrices and linklist from memory
 */
  deleteMatrix(&UniqueSites);
  deleteLL2(&ChargePercLL);
  deleteLL2(&ChargeTrapLL);
  
  return 0;
}

int uniqueID(int x, int y, int z, int wid, int hei){ 
  return x*hei*wid+hei*y+z;
}

int ConvertMatrixLinkList(const_matrix ChargePath, linklist2 * ChargePathLL,int wid, int hei){
  
  int SiteID;
  int i;
  int rows;
  double x, y, z;
  rows = getRows(ChargePath);
  
  for(i=1;i<=rows;i++){
    x = getE(ChargePath,i,1);
    y = getE(ChargePath,i,2);
    z = getE(ChargePath,i,3);
    SiteID = uniqueID((int)x,(int)y,(int)z,wid,hei);
 
    if(i==1){
      *ChargePathLL = newLinkList2(SiteID);
    }else{
      addLLNode2(*ChargePathLL,SiteID);
    }
  }

  return 0;
}

