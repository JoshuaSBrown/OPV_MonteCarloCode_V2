

//Read the .path file store the path of each charge in an array
int readPath(char * FileName);

//Given .path file data find the pathway of a single charge and
//store in an array or matrix
int getPathCharge(int ID);

//Given the pathway of a single charge find the trap sites where 
//the charge does not progress throught the sample but goes in circles
int getTrap();

//Gitven the pathway of a single charge find the percolation path
//this excludes the site where the charge did not progress through
//the sample
int getPerc();

//Write a .perc file storing all the sites that charges visited that
//were essential for crossing the system
int writePerc(char * FileName);

//Write a .trap file storing all the sites that charges visited that
//were not essential for crossing the system
int writeTrap(char * FileName);

