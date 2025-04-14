#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
MFEM_DIR=${SCRIPT_DIR}/..
(cd ${MFEM_DIR}/data/gmsh;
 gmsh -3 -refine -setnumber partitioner 1 periodic-cube.geo
 gmsh -3 -refine -setnumber partitioner 0 periodic-cube.geo)
(cd ${MFEM_DIR}/build/examples;
 make -j 12 ex8p;
 mpirun  -np 4 ./ex8p -no-gm -m periodic-cube.msh >out_ex8.no-gm
 mpirun  -np 4 ./ex8p -gm -m periodic-cube_ >out_ex8.gm
 for i in out_ex8*; do grep Total $i> r_$i; done
 echo " "
 echo "Result of the diff :"
 diff -s r_out_ex8.gm r_out_ex8.no-gm
 )
