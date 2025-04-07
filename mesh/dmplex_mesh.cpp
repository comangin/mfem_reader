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



namespace mfem
{
  // Load an .h5 file with PETSc
  PetscErrorCode Mesh::LoadMeshHDF5fromfile(const std::string &filename,bool &is_dmplex) {
    PetscErrorCode ierr;
    PetscViewer viewer;

    struct _n_DMPlexStorageVersion version = {3, 0, 0};

    if (filename.length() < 3 || filename.substr(filename.size() - 3) != ".h5") {
      std::cerr << "The file '" << filename << "' is not a DMplex" << std::endl;
      is_dmplex = false;
      return 0;
    }

    char objectname[PETSC_MAX_PATH_LEN];
    PetscBool flg;
    PetscOptionsGetString(nullptr, nullptr, "-o", objectname, sizeof(objectname), &flg);

    std::cout << "File name: " << filename << " and mesh name: " << objectname << std::endl;

    PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename.c_str(),objectname, PETSC_TRUE, &dm));
    PetscCall(PetscObjectSetName((PetscObject)dm, objectname));
    PetscCall(DMSetOptionsPrefix(dm, "loaded_"));
    PetscCall(DMViewFromOptions(dm, NULL, "-dm_view"));
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5SetDMPlexStorageVersionReading(viewer, &version));
    is_dmplex=true;
    return ierr;
  }

  
  PetscErrorCode Mesh::LoaderHDF5(int generate_edges, const std::string &parse_tag) {
    int curved = 0, read_gf = 1;
    bool finalize_topo = true;

    ReadDmplex(curved,read_gf);
    FinalizeTopology();
    PrintCharacteristics();
    
    // if (curved && read_gf)
    // {
    //    Nodes = new GridFunction(this, input);

    //    own_nodes = 1;
    //    spaceDim = Nodes->VectorDim();
    //    if (ncmesh) { ncmesh->spaceDim = spaceDim; }

    //    // Set vertex coordinates from the 'Nodes'
    //    SetVerticesFromNodes(Nodes);
    // }

    // // If a parse tag was supplied, keep reading the stream until the tag is
    // // encountered.
    // if (mfem_version >= 12)
    // {
    //    string line;
    //    do
    //    {
    //       skip_comment_lines(input, '#');
    //       MFEM_VERIFY(input.good(), "Required mesh-end tag not found");
    //       getline(input, line);
    //       filter_dos(line);
    //       // mfem v1.2 may not have parse_tag in it, e.g. if trying to read a
    //       // serial mfem v1.2 mesh as parallel with "mfem_serial_mesh_end" as
    //       // parse_tag. That's why, regardless of parse_tag, we stop reading if
    //       // we find "mfem_mesh_end" which is required by mfem v1.2 format.
    //       if (line == "mfem_mesh_end") { break; }
    //    }
    //    while (line != parse_tag);
    // }
    // else if (mfem_nc_version >= 10)
    // {
    //    string ident;
    //    skip_comment_lines(input, '#');
    //    input >> ident;
    //    MFEM_VERIFY(ident == "mfem_mesh_end",
    //                "invalid mesh: end of file tag not found");
    // }

    // // Finalize(...) should be called after this, if needed.
    return 0;
  }

  PetscErrorCode FinalizeHDF5(bool refine, bool fix_orientation) {
    // Implementation of FinalizeHDF5
    return 0;
  }

  PetscErrorCode Mesh::LoadDmplex(int generate_edges,int refine, bool fix_orientation = true)
  {
    std :: string tag_parse = "";
    LoaderHDF5(generate_edges,tag_parse);
    //Finalize(refine, fix_orientation);
    
    return 0;
  }


  PetscErrorCode Mesh::ReadDmplex(int curved, int read_gf) {
    Vec coordinates;
    PetscErrorCode ierr;
    IS              globalVertexNumbers = NULL;
    PetscInt dim, coordDim, nValues, numCellsStart, numCellsEnd;
    PetscReal *coords;
    PetscInt *cones, numElements,dm_dim;
    PetscMPIInt rank;

    map<int, int> vertices_map;
     
    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    coordDim = 3;

    // VERTICES // 
    PetscCall(DMGetCoordinates(dm, &coordinates));
    PetscCall(DMGetDimension(dm,&dm_dim));
    spaceDim=coordDim; 
    Dim=dm_dim;
    PetscCall(VecGetLocalSize(coordinates, &nValues));
    PetscCall(VecGetArray(coordinates, &coords));
    cout << "Size vec : " << nValues << endl; 
    NumOfVertices = nValues /coordDim;
    vertices.SetSize(NumOfVertices);
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Process %d : Sommets locaux = %d\n", rank, NumOfVertices));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    real_t coordLocal[3];
    
    for (PetscInt i = 0; i < NumOfVertices; ++i) {
      for (int d = 0; d < coordDim; ++d) {
       	coordLocal[d] = coords[i * coordDim + d];
      } 
      vertices[i]=Vertex(coordLocal,coordDim);
      cout << "Points " << i << " : " << coordLocal[0] << " " << coordLocal[1] << " "  << coordLocal[2] << endl;
      PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    }


    // ELEMENTS // 
    DMPolytopeType celltype;
    PetscBool hasLabel;
  
    DMPlexGetHeightStratum(dm, 0, &numCellsStart, &numCellsEnd);
    NumOfElements=numCellsEnd-numCellsStart;
    elements.SetSize(NumOfElements);
    cout << "Num Of Elements : " << NumOfElements << endl; 
    
    for (PetscInt i = numCellsStart; i < numCellsEnd; ++i) {
      const PetscInt *closure = NULL;
      PetscInt closureSize;

      DMPlexGetTransitiveClosure(dm, i, PETSC_TRUE, &closureSize, (PetscInt**)&closure);CHKERRQ(ierr);

      PetscPrintf(PETSC_COMM_WORLD, "Element %d :  ", i - numCellsStart);

      PetscInt vertex_tetra[closureSize];
      PetscInt Nv = 0;
      PetscInt vStart, vEnd;

      ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);

      for (PetscInt cl = 0; cl < closureSize * 2; cl += 2) {
    	PetscInt vertex = closure[cl];

    	if (vertex >= vStart && vertex < vEnd) {
    	  vertex_tetra[Nv++] = vertex;
    	}
      }

      for (PetscInt j = 0; j < Nv; ++j) {
    	cout << vertex_tetra[j]-numCellsEnd+1 << " ";
	    vertex_tetra[j]=vertex_tetra[j]-numCellsEnd;
      }

      // PHYSICAL GROUP //
      DMHasLabel(dm, "celltype", &hasLabel);
      PetscInt groupID;
      DMGetLabelValue(dm, "Cell Sets", i, &groupID);
      if (groupID != -1) {
	    PetscPrintf(PETSC_COMM_WORLD, "belongs to the physical group %d\n", groupID);
      } else {
	      PetscPrintf(PETSC_COMM_WORLD, "does not belong to any physical group\n");
      }


      // TYPE ELEMENT //
      DMPlexGetCellType(dm,i,&celltype);
      cout << "Type : " << celltype << endl;

      switch (celltype)
      	{
      	case 0:
	        elements[i]=new Point(&vertex_tetra[0],groupID);
      	  break;

      	case 1:
      	  // Handle SEGMENTS case
      	  break;

      	case 2:
      	  // Handle POINT_PRISM case
      	  break;

      	case 3:
      	  elements[i]=new Triangle(&vertex_tetra[0],groupID);
      	  break;

      	case 4:
      	  elements[i]=new Quadrilateral(&vertex_tetra[0],groupID);
      	  break;

      	case 5:
      	  // Handle SEG_PRIM case
      	  break;

      	case 6:
	        elements[i]=new Tetrahedron(&vertex_tetra[0],groupID);
      	  break;

      	case 7:
      	  elements[i]=new Hexahedron(&vertex_tetra[0],groupID);
      	  break;

      	case 8:
      	  // Handle TRI_PRISM case
      	  break;

      	case 9:
      	  // Handle TRI_PRISM_TENSOR case
      	  break;

      	case 10:
      	  // Handle QUAD_PRISM_TENSOR case
      	  break;

      	case 11:
      	  // Handle PYRAMID case
      	  break;

      	case 12:
      	  // Handle PV_GHOST case
      	  break;

      	case 13:
      	  // Handle INTERIOR_GHOST case
      	  break;

      	case 14:
      	  // Handle UNKNOWN case
      	  break;

      	case 15:
      	  // Handle UNKNOWN_CELL case
      	  break;

      	case 16:
      	  // Handle UNKNOWN_FACE case
      	  break;

      	case 17:
      	  // Handle POLYTOPES case
      	  break;

      	default:
      	  // Handle unexpected cell types
      	  break;
      	}

      DMPlexRestoreTransitiveClosure(dm, i, PETSC_TRUE, &closureSize, (PetscInt**)&closure);CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD, "\n");
    }

    DMView(dm, PETSC_VIEWER_STDOUT_WORLD);

    return 0;

 }


}
