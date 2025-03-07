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

PetscErrorCode Mesh::LoadMeshHDF5fromfile(const std::string &filename,bool &is_dmplex) {
    PetscErrorCode ierr;
    PetscViewer viewer;

    struct _n_DMPlexStorageVersion version = {3, 0, 0};

    // Vérifier si l'extension du fichier est .h5
    if (filename.length() < 3 || filename.substr(filename.size() - 3) != ".h5") {
      std::cerr << "Erreur : Le fichier '" << filename << "' n'a pas une extension .h5 valide." << std::endl;
      is_dmplex = false;
      return 0;
    }

    char objectname[PETSC_MAX_PATH_LEN];
    PetscBool flg;

    PetscOptionsGetString(nullptr, nullptr, "-o", objectname, sizeof(objectname), &flg);


    std::cout << "Nom du fichier : " << filename << " et nom du maillage :" << objectname <<  std::endl;

    PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename.c_str(),objectname, PETSC_TRUE, &dm));
    PetscCall(PetscObjectSetName((PetscObject)dm, objectname));
    PetscCall(DMSetOptionsPrefix(dm, "loaded_"));
    PetscCall(DMViewFromOptions(dm, NULL, "-dm_view"));

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5SetDMPlexStorageVersionReading(viewer, &version));

    std::cout << "Avec succes" << std::endl;
    is_dmplex=true;


    return ierr;
  }


PetscErrorCode Mesh::LoaderHDF5(int generate_edges, const std::string &parse_tag) {

  int curved = 0, read_gf = 1;
  bool finalize_topo = true;

  const char     *name;

  PetscObjectGetName((PetscObject)dm, &name);
  std::cout << "Le nom de l'objet maillage est : " << name << std::endl;

  ReadDmplex(curved,read_gf);

    DMDestroy(&dm);
    PetscFinalize();
    exit(0);

  // if (finalize_topo)
  // {
  //    // don't generate any boundary elements, especially in parallel
  //    bool generate_bdr = false;

  //    FinalizeTopology(generate_bdr);
  // }

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
    cout << "hello" << endl;
    LoaderHDF5(generate_edges,tag_parse);
    //  FinalizeHDF5(refine, fix_orientation);
    return 0;
  }


PetscErrorCode Mesh::ReadDmplex(int curved, int read_gf) {
    Vec coordinates;
    PetscInt dim, coordDim, nValues;
    PetscReal *coords;
    PetscInt nVertices;
    PetscMPIInt rank;

    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    coordDim = 3;
    dim = 3;

    const char *name;
    PetscObjectGetName((PetscObject)dm, &name);
    std::cout << "Le nom de l'objet maillage est : " << name << std::endl;

    PetscCall(DMGetCoordinates(dm, &coordinates));
    PetscCall(VecGetLocalSize(coordinates, &nValues));
    PetscCall(VecGetArray(coordinates, &coords));
    NumOfVertices = nValues / coordDim;
    vertices.SetSize(NumOfVertices);
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Process %d : Sommets locaux = %d\n", rank, nVertices));
    real_t coordLocal[3];
    for (PetscInt i = 0; i < nVertices; ++i) {
      for (int d = 0; d < coordDim; ++d) {
       	coordLocal[d] = coords[i * coordDim + d];
	  }
      vertices[i]=Vertex(coordLocal,3);
      cout << "Ok pour " << i << " coordonnée : " << coordLocal[0] << " " << coordLocal[1] << " "  << coordLocal[2] << endl;
    }
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    PetscCall(VecRestoreArray(coordinates, &coords));
    //vertices[ver] = Vertex(coord, gmsh_dim);
    // PetscViewer viewer;
    // PetscViewerVTKOpen(PETSC_COMM_WORLD, "output_mesh_ascii.vtu", FILE_MODE_WRITE, &viewer);
    // PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
    // DMView(dm, viewer);
    // PetscViewerDestroy(&viewer);


    PrintInfo();

    return 0;
 }
}
