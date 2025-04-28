#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
MFEM_DIR=${SCRIPT_DIR}/..
(cd ${MFEM_DIR}/data/gmsh;
 gmsh -3 -refine -setnumber partitioner 1 periodic-cube.geo
 gmsh -3 -refine -setnumber partitioner 0 periodic-cube.geo
 gmsh -3 -refine -setnumber partitioner 1 cube-periodic.geo
 gmsh -3 -refine -setnumber partitioner 0 cube-periodic.geo
 cp cube-per*msh periodic-cube*msh ../../build/examples
 )
(cd ${MFEM_DIR}/build/examples;
 make -j 12 ex8p;
 mpirun  -np 4 ./ex8p -no-vis -no-gp -m periodic-cube.msh >out_ex8a.no-gp
 mpirun  -np 4 ./ex8p -no-vis -gp -m    periodic-cube_    >out_ex8a.gp
 mpirun  -np 4 ./ex8p -no-vis -no-gp -m cube-periodic.msh >out_ex8b.no-gp
 mpirun  -np 4 ./ex8p -no-vis -gp -m    cube-periodic_    >out_ex8b.gp
 for i in out_ex8*; do grep Total $i> r_$i; done
 echo " "
 echo "Result of the diff 8a :"
 diff -s r_out_ex8a.gp r_out_ex8a.no-gp
 echo "Result of the diff 8b :"
 diff -s r_out_ex8b.gp r_out_ex8b.no-gp
 )
(cd ${MFEM_DIR}/build/examples;
 make -j 12 ex4p;
 mpirun  -np 4 ./ex4p -no-vis -no-gp -m periodic-cube.msh >out_ex4a.no-gp
 mpirun  -np 4 ./ex4p -no-vis -gp -m    periodic-cube_    >out_ex4a.gp
 mpirun  -np 4 ./ex4p -no-vis -no-gp -m cube-periodic.msh >out_ex4b.no-gp
 mpirun  -np 4 ./ex4p -no-vis -gp -m    cube-periodic_    >out_ex4b.gp
 echo " "
 echo "Result of the diff 4a :"
 diff -s out_ex4a.gp out_ex4a.no-gp
 echo "Result of the diff 4b :"
 diff -s out_ex8b.gp out_ex8b.no-gp
 )
