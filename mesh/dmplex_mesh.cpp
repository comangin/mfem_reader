#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>

#include "petscdm.h"
#include "petscdmlabel.h"
#include "petscds.h"
#include "dmplex_mesh.hpp"
#include <petsc.h>
#include <petscdmplex.h>
#include <petscviewerhdf5.h>
#include "mesh_headers.hpp"
#include "../fem/fem.hpp"

using namespace std;

void swap(int tab[], int i, int j) {
    int temp = tab[i];
    tab[i] = tab[j];
    tab[j] = temp;
}

namespace mfem
{
  // Load an .h5 file with PETSc
  PetscErrorCode Mesh::LoadMeshHDF5fromfile(const std::string &filename, bool &is_dmplex)
  {
    PetscErrorCode ierr;
    PetscViewer viewer;

    struct _n_DMPlexStorageVersion version = {3, 0, 0};

    if (filename.length() < 3 || filename.substr(filename.size() - 3) != ".h5")
    {
      // std::cerr << "The file '" << filename << "' is not a DMplex" << std::endl;
      is_dmplex = false;
      return 0;
    }

    PetscInitialize(nullptr,nullptr, nullptr, nullptr);
    string objectname = filename;
    size_t pos = objectname.rfind(".h5");
    objectname = objectname.substr(0, pos);
    PetscBool flg;
    cout << "DEBUG HERE" << endl;
    cout << filename << endl; 
    PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename.c_str(), objectname.c_str(), PETSC_TRUE, &dm));
    PetscCall(PetscObjectSetName((PetscObject)dm, objectname.c_str()));
    PetscCall(DMSetOptionsPrefix(dm, "loaded_"));
    PetscCall(DMViewFromOptions(dm, NULL, "-dm_view"));
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5SetDMPlexStorageVersionReading(viewer, &version));
    is_dmplex = true;
    return ierr;
  }

  PetscErrorCode Mesh::LoaderHDF5(int generate_edges, const std::string &parse_tag)
  {
    int curved = 0, read_gf = 1;
    bool finalize_topo = true;

    ReadDmplex(curved, read_gf);
    FinalizeTopology();
    CheckElementOrientation(true);
    Finalize();
    return 0;
  }

  PetscErrorCode FinalizeHDF5(bool refine, bool fix_orientation)
  {
    // Implementation of FinalizeHDF5
    return 0;
  }

  PetscErrorCode Mesh::LoadDmplex(int generate_edges, int refine, bool fix_orientation = true)
  {
    std ::string tag_parse = "";
    LoaderHDF5(generate_edges, tag_parse);
    Finalize(refine, fix_orientation);
    return 0;
  }

  PetscErrorCode Mesh::ReadDmplex(int curved, int read_gf)
  {
    Vec coordinates;
    PetscErrorCode ierr;
    IS globalVertexNumbers = NULL;
    PetscInt dim, coordDim, nValues, numCellsStart, numCellsEnd;
    PetscReal *coords;
    PetscInt *cones, numElements, dm_dim;
    PetscMPIInt rank;

    map<int, int> vertices_map;

    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    coordDim = 3;

    // VERTICES //
    PetscCall(DMGetCoordinates(dm, &coordinates));
    PetscCall(DMGetDimension(dm, &dm_dim));
    spaceDim = dm_dim; // Need modification here
    Dim = dm_dim;
    PetscCall(VecGetLocalSize(coordinates, &nValues));
    PetscCall(VecGetArray(coordinates, &coords));
    NumOfVertices = nValues / coordDim;
    vertices.SetSize(NumOfVertices);
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    real_t coordLocal[3];

    for (PetscInt i = 0; i < NumOfVertices; ++i)
    {
      for (int d = 0; d < coordDim; ++d)
      {
        coordLocal[d] = coords[i * coordDim + d];
      }
      vertices[i] = Vertex(coordLocal, coordDim);
      PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    }

    // ELEMENTS //
    DMPolytopeType celltype;
    PetscBool hasLabel;

    DMPlexGetHeightStratum(dm, 0, &numCellsStart, &numCellsEnd);
    NumOfElements = numCellsEnd - numCellsStart;
    elements.SetSize(NumOfElements);

    for (PetscInt i = numCellsStart; i < numCellsEnd; ++i)
    {
      const PetscInt *closure = NULL;
      PetscInt closureSize;

      DMPlexGetTransitiveClosure(dm, i, PETSC_TRUE, &closureSize, (PetscInt **)&closure);
      CHKERRQ(ierr);
      PetscInt vertex_tetra[closureSize];
      PetscInt Nv = 0;
      PetscInt vStart, vEnd;

      ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
      CHKERRQ(ierr);

      for (PetscInt cl = 0; cl < closureSize * 2; cl += 2)
      {
        PetscInt vertex = closure[cl];

        if (vertex >= vStart && vertex < vEnd)
        {
          vertex_tetra[Nv++] = vertex;
        }
      }

      for (PetscInt j = 0; j < Nv; ++j)
      {
        vertex_tetra[j] = vertex_tetra[j] - numCellsEnd;
      }

      // PHYSICAL GROUP //
      PetscBool hasLabel;
      PetscInt groupID;
      DMHasLabel(dm, "Cell Sets", &hasLabel);

      if (hasLabel)
      {
        DMGetLabelValue(dm, "Cell Sets", i, &groupID);
      }

      if (groupID == 0)
        groupID = 1;

      // TYPE ELEMENT //
      DMPlexGetCellType(dm, i, &celltype);
      int tag; 
      switch (celltype)
      {
      case 0:
        elements[i] = new Point(&vertex_tetra[0], groupID);
        break;
      case 1:
        elements[i] = new Segment(&vertex_tetra[0], groupID);
        break;
      case 3:
        elements[i] = new Triangle(&vertex_tetra[0], groupID);
        break;

      case 4:
        elements[i] = new Quadrilateral(&vertex_tetra[0], groupID);
        break;

      case 5:
                 // TODO: Handle SEG_PRIM
        break;

      case 6:
        swap(vertex_tetra,0,1); // Inversion par PETsc 
        elements[i] = new Tetrahedron(&vertex_tetra[0], groupID);
        break;

      case 7:
        swap(vertex_tetra,1,3); // Inversion par PETsc 
        elements[i] = new Hexahedron(&vertex_tetra[0], groupID);
        break;

      case 8:
        swap(vertex_tetra,1,2); // Inversion par PETsc 
        elements[i] = new Wedge(&vertex_tetra[0], groupID);
        break;

      case 9:
        // TODO: Handle TRI_PRISM_TENSOR
        break;

      case 10:
        // TODO: Handle QUAD_PRISM_TENSOR
        break;

      case 11:
        swap(vertex_tetra,1,3); // Inversion par PETsc 
        elements[i] = new Pyramid(&vertex_tetra[0], groupID);
        break;

      case 12:
        // TODO: Handle PV_GHOST
        break;

      case 13:
        // TODO: Handle INTERIOR_GHOST
        break;

      case 14:
      case 15:
      case 16:
      case 17:
        std::cerr << "Unknown or unsupported cell type: " << celltype << std::endl;
        break;

      default:
        std::cerr << "Unhandled cell type: " << celltype << std::endl;
        break;
      }

      DMPlexRestoreTransitiveClosure(dm, i, PETSC_TRUE, &closureSize, (PetscInt **)&closure);
      CHKERRQ(ierr);
    }
    return 0;
  }

}
