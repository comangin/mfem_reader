#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>

#include "petscdm.h"
#include "petscdmlabel.h"
#include "petscds.h"
#include "dmplex_mesh.hpp"
#include <petsc.h>
#include <petscdmplex.h>
#include <petscviewerhdf5.h>
#include "mesh_headers.hpp"
#include "../fem/fem.hpp"

void swap(int tab[], int i, int j)
{
   int temp = tab[i];
   tab[i] = tab[j];
   tab[j] = temp;
}

namespace mfem
{
// Load an .h5 file with PETSc
PetscErrorCode Mesh::LoadMeshHDF5fromfile(const std::string &filename,
                                          bool &is_dmplex)
{
   PetscErrorCode ierr;
   PetscViewer viewer;

   struct _n_DMPlexStorageVersion version = {3, 0, 0};

   if (filename.length() < 3 || filename.substr(filename.size() - 3) != ".h5")
   {
      is_dmplex = false;
      return 0;
   }

   PetscInitialize(nullptr, nullptr, nullptr, nullptr);
   std::string objectname = filename;
   size_t pos1 = objectname.find_last_of("/");
   size_t pos2 = objectname.rfind(".h5");
   objectname = objectname.substr(pos1+1, pos2-pos1-1);
   PetscBool flg;
   PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename.c_str(),
                                  objectname.c_str(), PETSC_TRUE, &dm));
   PetscCall(PetscObjectSetName((PetscObject)dm, objectname.c_str()));
   // PetscCall(DMView(dm, PETSC_VIEWER_STDOUT_WORLD));
   PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(),
                                 FILE_MODE_READ, &viewer));
   PetscCall(PetscViewerHDF5SetDMPlexStorageVersionReading(viewer, &version));
   is_dmplex = true;
   return ierr;
}

PetscErrorCode Mesh::LoaderHDF5(int generate_edges)
{
   int curved = 0, read_gf = 1;

   ReadDmplex(curved, read_gf);
   FinalizeTopology();
   return 0;
}

PetscErrorCode FinalizeHDF5(bool refine, bool fix_orientation)
{
   // Implementation of FinalizeHDF5
   return 0;
}

PetscErrorCode Mesh::LoadDmplex(int generate_edges, int refine,
                                bool fix_orientation = true)
{
   LoaderHDF5(generate_edges);
   Finalize(refine, fix_orientation);

   CheckElementOrientation(true);
   return 0;
}

PetscErrorCode Mesh::ReadDmplex(int curved, int read_gf)
{
   Vec coordinates;
   PetscErrorCode ierr;
   IS globalVertexNumbers = NULL;
   PetscInt dim, coordDim, nValues, numCellsStart, numCellsEnd;
   PetscReal *coords;
   PetscInt *cones, numElements, dm_dim, sp_dim;
   PetscMPIInt rank;

   std::map<int, int> vertices_map;

   PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

   // VERTICES //
   PetscCall(DMGetCoordinates(dm, &coordinates));
   PetscCall(DMGetDimension(dm, &dm_dim));
   PetscCall(DMGetCoordinateDim(dm, &sp_dim));
   spaceDim = sp_dim;
   Dim = dm_dim;
   PetscCall(VecGetLocalSize(coordinates, &nValues));
   PetscCall(VecGetArray(coordinates, &coords));
   NumOfVertices = nValues / spaceDim;
   vertices.SetSize(NumOfVertices);
   PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
   real_t coordLocal[3];

   for (PetscInt i = 0; i < NumOfVertices; ++i)
   {
      for (int d = 0; d < spaceDim; ++d)
      {
         coordLocal[d] = coords[i * spaceDim + d];
      }
      vertices[i] = Vertex(coordLocal, spaceDim);
      PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
   }

   // ELEMENTS IN CELL SETS //
   DMPlexGetHeightStratum(dm, 0, &numCellsStart, &numCellsEnd);
   NumOfElements = numCellsEnd - numCellsStart;
   elements.SetSize(NumOfElements);

   for (PetscInt i = numCellsStart; i < numCellsEnd; ++i)
   {
      PetscInt *closure = NULL;
      PetscInt closureSize;

      PetscCall(DMPlexGetTransitiveClosure(dm, i, PETSC_TRUE,
                                           &closureSize, reinterpret_cast<PetscInt **>(&closure)));
      PetscInt vertex_tetra[closureSize];
      PetscInt Nv = 0;
      PetscInt vStart, vEnd;

      PetscCall(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));

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

      // ELEMENTS PHYSICAL GROUPS //
      DMPolytopeType celltype;
      PetscBool hasLabel, stratumHasPoint;
      PetscInt groupID = -1, numLabelValues;
      DMLabel cellSetslabel;
      DMHasLabel(dm, "Cell Sets", &hasLabel);

      if (hasLabel)
      {
         DMGetLabel(dm, "Cell Sets", &cellSetslabel);
         DMLabelGetNumValues(cellSetslabel, &numLabelValues);
         for (int j = 1; j <= numLabelValues; j++)
         {
            DMLabelStratumHasPoint(cellSetslabel, j, i, &stratumHasPoint);
            if (stratumHasPoint) { groupID = j; }
         }
      }

      if (groupID == -1)
      {
         groupID = 1;
      }

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
            swap(vertex_tetra, 0, 1);  // Inversion par PETsc
            elements[i] = new Tetrahedron(&vertex_tetra[0], groupID);
            break;

         case 7:
            swap(vertex_tetra, 1, 3);  // Inversion par PETsc
            elements[i] = new Hexahedron(&vertex_tetra[0], groupID);
            break;

         case 8:
            swap(vertex_tetra, 1, 2);  // Inversion par PETsc
            elements[i] = new Wedge(&vertex_tetra[0], groupID);
            break;

         case 9:
            // TODO: Handle TRI_PRISM_TENSOR
            break;

         case 10:
            // TODO: Handle QUAD_PRISM_TENSOR
            break;

         case 11:
            swap(vertex_tetra, 1, 3);  // Inversion par PETsc
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
            mfem::err << "Unknown/unsupported cell type: " << celltype << std::endl;
            break;

         default:
            mfem::err << "Unhandled cell type: " << celltype << std::endl;
            break;
      }

      PetscCall(DMPlexRestoreTransitiveClosure(dm, i, PETSC_TRUE,
                                               &closureSize, reinterpret_cast<PetscInt **>(&closure)));
   }

   // FACETS IN FACE SETS -> BOUNDARY ELEMENTS //
   DMLabel faceSetslabel;
   PetscBool hasLabel;
   PetscInt numLabelValues, NumOfFacets;
   IS is;
   PetscInt p;
   const PetscInt *points;
   DMHasLabel(dm, "Face Sets", &hasLabel);
   if (hasLabel)
   {
      DMGetLabel(dm, "Face Sets", &faceSetslabel);
      DMLabelGetNumValues(faceSetslabel, &numLabelValues);
      DMLabelGetStratumSize(faceSetslabel, 0, &NumOfFacets);
      NumOfBdrElements = NumOfFacets;
      boundary.SetSize(NumOfBdrElements);
      DMLabelGetStratumIS(faceSetslabel, 0, &is);
      if (is)
      {
         ISGetIndices(is, &points);
         for (PetscInt i = 0; i < NumOfFacets; ++i)
         {
            PetscInt *closure = NULL;
            PetscInt closureSize;

            PetscCall(DMPlexGetTransitiveClosure(dm, points[i], PETSC_TRUE,
                                                 &closureSize, reinterpret_cast<PetscInt **>(&closure)));
            PetscInt vertex_tetra[closureSize];
            PetscInt Nv = 0;
            PetscInt vStart, vEnd;

            PetscCall(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));

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

            DMPolytopeType celltype;
            PetscBool stratumHasPoint;
            PetscInt groupID = -1;

            for (int j = 1; j <= numLabelValues; j++)
            {
               DMLabelStratumHasPoint(faceSetslabel, j, points[i], &stratumHasPoint);
               if (stratumHasPoint)
               {
                  groupID = j;
               }
            }

            if (groupID == -1)
            {
               groupID = 1;
            }

            // TYPE ELEMENT //
            DMPlexGetCellType(dm, points[i], &celltype);
            int tag;
            switch (celltype)
            {
               case 0:
                  boundary[i] = new Point(&vertex_tetra[0], groupID);
                  break;
               case 1:
                  boundary[i] = new Segment(&vertex_tetra[0], groupID);
                  break;
               case 3:
                  boundary[i] = new Triangle(&vertex_tetra[0], groupID);
                  break;
               case 4:
                  boundary[i] = new Quadrilateral(&vertex_tetra[0], groupID);
                  break;
               default:
                  mfem::err << "Unhandled cell type: " << celltype << std::endl;
                  break;
            }

            PetscCall(DMPlexRestoreTransitiveClosure(dm, points[i], PETSC_TRUE,
                                                     &closureSize, reinterpret_cast<PetscInt **>(&closure)));
            ISRestoreIndices(is, &points);
         }
      }
   }

   this->RemoveUnusedVertices();
   this->RemoveInternalBoundaries();
   display_mesh();
   return 0;
}

}
