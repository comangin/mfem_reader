#ifndef MFEM_DMPLEX_MESH
#define MFEM_DMPLEX_MESH

#include <petsc.h>
#include <petscdmplex.h>
#include <petscviewerhdf5.h>
#include "mesh.hpp"

PetscErrorCode FinalizeHDF5(bool refine = false, bool fix_orientation = false);


#endif
