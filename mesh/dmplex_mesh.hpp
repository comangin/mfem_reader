#ifndef MFEM_DMPLEX_MESH
#define MFEM_DMPLEX_MESH

#include <petsc.h>
#include <petscdmplex.h>
#include <petscviewerhdf5.h>

PetscErrorCode LoadMeshHDF5fromfile(const std :: string &filename, bool &is_dmplex);


PetscErrorCode LoaderHDF5(int generate_edges = 0,
	    std::string parse_tag = "");

PetscErrorCode FinalizeHDF5(bool refine = false, bool fix_orientation = false);

PetscErrorCode LoadDmplex(int generate_edges = 0,                                                                    
			  int refine = 1, bool fix_orientation = true);

PetscErrorCode ConversionDmplextoMfem(int curved, int read_gf);

#endif
