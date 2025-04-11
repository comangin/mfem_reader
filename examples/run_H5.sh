#!/bin/bash
cd ../
cd build || exit 1

echo "Please compile Mfem with PETSc support"

if [ ! -f Makefile ]; then
    echo "Makefile not found. Stopping."
    exit 1
fi

make -j 10 reader_test

exec=$PWD/miniapps/meshing/reader_test

## Phase test 
cd ../data/hdf5

echo "=============================="
echo "== Meshes without periodic test =="
echo "=============================="

echo "Mesh : elementary_mesh"
$exec -m elementary_mesh.msh > res_msh
$exec -m elementary_mesh.h5 -o elementary > res_h5
diff -s res_msh res_h5

echo "Mesh : geometrie"
$exec -m geometrie.msh > res_msh
$exec -m geometrie.h5 -o geometrie > res_h5
diff -s res_msh res_h5

echo "Mesh : rectangle_mesh"
$exec -m rectangle_mesh.msh > res_msh
$exec -m rectangle_mesh.h5 -o rectangle > res_h5
diff -s res_msh res_h5

echo "Mesh : square_mesh"
$exec -m square_mesh.msh > res_msh
$exec -m square_mesh.h5 -o square > res_h5
diff -s res_msh res_h5

echo "Mesh : boolean"
$exec -m boolean.msh > res_msh
$exec -m boolean.h5 -o boolean > res_h5
diff -s res_msh res_h5

echo "Mesh : bspline_bezier_patches"
$exec -m bspline_bezier_patches.msh > res_msh
$exec -m bspline_bezier_patches.h5 -o bspline_bezier_patches > res_h5

echo "Mesh : bspline_filling"
$exec -m bspline_filling.msh
$exec -m bspline_filling.h5 -o bspline_filling > res_h5

echo ""
echo "==========================="
echo "== Meshes with periodic test =="
echo "==========================="

echo "Mesh : gmsh-3d"
$exec -m gmsh-3d.msh > res_msh
$exec -m gmsh-3d.h5 -o gmsh > res_h5
diff -s res_msh res_h5

echo "Mesh : mesh-periodic"
$exec -m mesh-periodic.msh > res_msh
$exec -m mesh-periodic.h5 -o mesh > res_h5
diff -s res_msh res_h5

echo "Mesh : periodic-rectangle"
$exec -m periodic-rectangle.msh > res_msh
$exec -m periodic-rectangle.h5 -o periodic > res_h5
diff -s res_msh res_h5

echo "Mesh : periodic-square"
$exec -m periodic-square.msh > res_msh
$exec -m periodic-square.h5 -o periodic > res_h5
diff -s res_msh res_h5