#include "mfem.hpp"                                                                                                                            
#include "../common/mfem-common.hpp"                                                                                                           
#include <fstream>                                                                                                                             
#include <limits>                                                                                                                              
#include <cstdlib>                                                                                                                             
#include <chrono>                                                                                                                              
#include <sys/resource.h>      
#include <petsc.h>
#include <iostream>


using namespace mfem;                                                                                                                            
using namespace std;   

int main(int argc, char *argv[]) {                                                                                                               
    // Initialisation de PETSc                                                                                                                 
    PetscInitialize(&argc, &argv, nullptr, nullptr);                                                                                           
                                                                                                                                               
    // Récupération des options                                                                                                                
    PetscBool flagM, flagO;                                                                                                                    
    char filename[256], objectname[256];                                                                                                       
    PetscOptionsGetString(nullptr, nullptr, "-m", filename, sizeof(filename), &flagM);                                                         
    PetscOptionsGetString(nullptr, nullptr, "-o", objectname, sizeof(objectname), &flagO);                                                     
                                                                                                                                               
    if (flagM) {                                                                                                                           
        // Lire le maillage depuis le fichier                                                                                                  
        Mesh mesh(filename, 1, 1, false);                                                                                                         
                                                                                                                                               
        // Vérifier quelques informations de base sur le maillage                                                                                  
        mesh.PrintInfo();                                                                                                                          
                                                                                                                                               
        std::cout << "Nombre de sommets : " << mesh.GetNV() << std::endl;                                                                          
        std::cout << "Nombre d'arêtes : " << mesh.GetNEdges() << std::endl;                                                                        
        std::cout << "Nombre de faces : " << mesh.GetNFaces() << std::endl;                                                                        
        std::cout << "Nombre d'éléments : " << mesh.GetNE() << std::endl;                                                                          
    } else {
        std::cout << "No filename provided." << std::endl;                                                                                  
    }

    if (flagO) {
        std::cout << "Object name: " << objectname << std::endl;
    } else {
        std::cout << "No object name provided." << std::endl;
    }
                                                                                                                                               
    // Faire quelque chose avec le maillage                                                                                                    
    // ... (par exemple, traiter ou visualiser le maillage)                                                                                    
                                                                                                                                               
    // Finaliser PETSc                                                                                                                         
    PetscFinalize();                                                                                                                           
    return 0;                                                                                                                                  
}
