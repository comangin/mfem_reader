#!/bin/bash

# Fonction pour comparer uniquement les statistiques des fichiers de sortie
compare_stats() {
    grep -E 'Dimension|Space dimension|Number of vertices|Number of edges|Number of elements|Number of bdr elem|Euler Number|h_min|h_max|kappa_min|kappa_max' res_msh > stats_msh
    grep -E 'Dimension|Space dimension|Number of vertices|Number of edges|Number of elements|Number of bdr elem|Euler Number|h_min|h_max|kappa_min|kappa_max' res_h5 > stats_h5

    echo "=> Comparing statistics:"
    diff_found=0
    while IFS= read -r line_msh && IFS= read -r line_h5 <&3; do
        key_msh=$(echo "$line_msh" | cut -d':' -f1 | xargs)
        val_msh=$(echo "$line_msh" | cut -d':' -f2- | xargs)
        val_h5=$(echo "$line_h5" | cut -d':' -f2- | xargs)
        if [ "$val_msh" != "$val_h5" ]; then
            echo "   ❌ $key_msh: msh = $val_msh | h5 = $val_h5"
            diff_found=1
        fi
    done < stats_msh 3< stats_h5

    if [ "$diff_found" -eq 0 ]; then
        echo "   ✅ All statistics match."
    fi

    rm -f stats_msh stats_h5
}


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
$exec -m elementary_mesh.h5  > res_h5
compare_stats

echo "Mesh : geometrie"
$exec -m geometrie.msh > res_msh
$exec -m geometrie.h5 > res_h5
compare_stats

echo "Mesh : rectangle_mesh"
$exec -m rectangle_mesh.msh > res_msh
$exec -m rectangle_mesh.h5  > res_h5
compare_stats

echo "Mesh : square_mesh"
$exec -m square_mesh.msh > res_msh
$exec -m square_mesh.h5 > res_h5
compare_stats

echo "Mesh : boolean"
$exec -m boolean.msh > res_msh
$exec -m boolean.h5  > res_h5
compare_stats

echo "Mesh : bspline_bezier_patches"
$exec -m bspline_bezier_patches.msh > res_msh
$exec -m bspline_bezier_patches.h5  > res_h5
compare_stats

echo "Mesh : bspline_filling"
$exec -m bspline_filling.msh > res_msh
$exec -m bspline_filling.h5  > res_h5
compare_stats

echo ""
echo "==========================="
echo "== Meshes with periodic test =="
echo "==========================="

echo "Mesh : gmsh-3d"
$exec -m gmsh-3d.msh > res_msh
$exec -m gmsh-3d.h5 > res_h5
compare_stats

echo "Mesh : mesh-periodic"
$exec -m mesh-periodic.msh > res_msh
$exec -m mesh-periodic.h5 > res_h5
compare_stats

echo "Mesh : periodic-rectangle"
$exec -m periodic-rectangle.msh > res_msh
$exec -m periodic-rectangle.h5 > res_h5
compare_stats

echo "Mesh : periodic-square"
$exec -m periodic-square.msh > res_msh
$exec -m periodic-square.h5 > res_h5
compare_stats

rm -f res_msh res_h5
