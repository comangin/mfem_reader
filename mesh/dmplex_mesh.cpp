#include <iostream>
#include <cstdio>                                                                                                    
#include <vector>                                                                                                    
#include <algorithm>                                                                                                 
#include <map>


#include <petsc.h>
#include <petscdmplex.h>                                                                      
#include <petscviewerhdf5.h>
  


using namespace std;


void petsc_test() 
{

  int argc = 1;
  char *argv[] = {(char*)"program_name"};
  char **argv_ptr = argv;
  
  cout << "Compile with petsc" << endl;
     
  PetscInitialize(&argc,&argv_ptr,PETSC_NULL,PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD,"Hello World\n");
  PetscFinalize();
	
}


void 
