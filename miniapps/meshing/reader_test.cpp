#include <petsc.h>
#include <sys/resource.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include "../common/mfem-common.hpp"
#include "mfem.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  PetscBool flagM = PETSC_FALSE, flagO = PETSC_FALSE, flagF = PETSC_FALSE,
            flagD = PETSC_FALSE;
  char filename[256] = {0}, objectname[256] = {0}, outname[256] = {0};

  PetscOptionsGetString(nullptr, nullptr, "-m", filename, sizeof(filename),
                        &flagM);
  PetscOptionsGetString(nullptr, nullptr, "-o", objectname, sizeof(objectname),
                        &flagO);
  PetscOptionsGetString(nullptr, nullptr, "-file", outname, sizeof(outname),
                        &flagF);
  PetscOptionsGetString(nullptr, nullptr, "-d", outname, sizeof(outname),
                        &flagD);

  Mesh mesh(filename);

  if (flagM)
  {
    if (flagD)
    {
      // Debug Mode
      for (int i = 0; i < mesh.GetNV(); i++)
      {
        const double *node = mesh.GetVertex(i);
        cout << "Noeud " << i << " : (" << node[0] << ", " << node[1] << ", "
             << node[2] << ")" << endl;
      }

      for (int i = 0; i < mesh.GetNE(); i++)
      {
        Element *el = mesh.GetElement(i);
        const int *vertices = el->GetVertices();
        int num_vertices = el->GetNVertices();

        cout << "Element " << i << " : Sommets = [ ";
        for (int j = 0; j < num_vertices; j++)
        {
          cout << vertices[j] << " ";
        }
        cout << "]" << endl;
      }
    }

    mesh.PrintInfo();
    cout << "Nombre de sommets : " << mesh.GetNV() << endl;
    cout << "Nombre d'arêtes : " << mesh.GetNEdges() << endl;
    cout << "Nombre de faces : " << mesh.GetNFaces() << endl;
    cout << "Nombre d'éléments : " << mesh.GetNE() << endl;
  }

  if (flagF)
  {
    std::ofstream ofs(outname);
    mesh.Print(ofs);
  }

  PetscFinalize();
  return 0;
}
