// ~/Softs/gmsh/bin/gmsh -3 -setnumber partitioner 0 -refine  cube-periodic.geo // One single msh file
// ~/Softs/gmsh/bin/gmsh -3 -setnumber partitioner 1 -refine  cube-periodic.geo // Several partitions and files with METIS partitioner
periodic = 0;
ep = 1.;
mp = 0.5*ep;
cl = mp/2.1;


Point(1) = { 0,   0,  0, cl} ; 
Point(2) = { 0,   0, ep, cl} ;
Point(3) = { 0,  ep, ep, cl} ;
Point(4) = { 0,  ep,  0, cl} ;

Point(5) = { mp,  0,  0, cl} ; 
Point(6) = { mp,  0, ep, cl} ;
Point(7) = { mp, ep, ep, cl} ;
Point(8) = { mp, ep,  0, cl} ;

Point(9)  = { ep,  0,  0, cl} ; 
Point(10) = { ep,  0, ep, cl} ;
Point(11) = { ep, ep, ep, cl} ;
Point(12) = { ep, ep,  0, cl} ;

Line(1) = {1,2}; 
Line(2) = {2,3};
Line(3) = {3,4}; 
Line(4) = {4,1};
Line Loop(1) = {  1, 2, 3, 4 } ;
Plane Surface(1) = {1};

Line(5) = {5,6}; 
Line(6) = {6,7};
Line(7) = {7,8}; 
Line(8) = {8,5};
Line Loop(2) = {  5, 6, 7, 8 } ;
Plane Surface(2) = {2};

Line(9)  = {9 ,10}; 
Line(10) = {10,11};
Line(11) = {11,12}; 
Line(12) = {12, 9};
Line Loop(3) = {  9, 10, 11, 12 } ;
Plane Surface(3) = {3};

Line(13) = {2, 6};
Line(14) = {1, 5};
Line(15) = {3, 7};
Line(16) = {4, 8};
Line(17) = {6, 10};
Line(18) = {5, 9};
Line(19) = {7, 11};
Line(20) = {8, 12};

Line Loop(4) = { 1, 13, -5, -14 };
Plane Surface(4) = {4};

Line Loop(5) = { 3, 16, -7, -15 };
Plane Surface(5) = {5};

Line Loop(6) = { 5, 17, -9, -18 };
Plane Surface(6) = {6};
Line Loop(7) = { 7, 20, -11, -19 };
Plane Surface(7) = {7};

Line Loop(8) = { 8, 18, -12, -20 };
Plane Surface(8) = {8};

Line Loop(9) = { 6, 19, -10, -17 };
Plane Surface(9) = {9};

Line Loop(10) = { 8, -14, -4, 16 };
Plane Surface(10) = {10};

Line Loop(11) = { 6, -15, -2, 13 };
Plane Surface(11) = {11};

If (periodic > 0)
  Periodic Surface {3} = {1} Translate {ep,0,0};
  Periodic Surface {5} = {4} Translate {0,ep,0};
  Periodic Surface {7} = {6} Translate {0,ep,0};
  Periodic Surface {9} = {8} Translate {0,0,ep};
  Periodic Surface {11} = {10} Translate {0,0,ep};
EndIf

Surface Loop (1) = {1, 2, 4, 5, 10, 11};
Volume (1) = {1};

Surface Loop (2) = {2, 3, 6, 7, 8, 9};
Volume (2) = {2};

Physical Volume(1)  = {1, 2};

//#For vo In {1:2}
//#    Physical Volume(vo)  = {vo}; 
//#EndFor
For su In {1:11}
    Physical Surface(su)  = {su}; 
EndFor
//For li In {1:100}
//    Physical Line(li)  = {li}; 
//EndFor
//For pt In {1:100}
//    Physical Point(pt)  = {pt}; 
//EndFor

Mesh 3;

// We now define several constants to fine-tune how the mesh will be partitioned
// We now define several constants to fine-tune how the mesh will be partitioned
DefineConstant[
  partitioner = {1, Choices{0="None", 1="Metis", 2="SimplePartition"},
    Name "Parameters/0Mesh partitioner"}
  N = {4, Min 1, Max 256, Step 1,
    Name "Parameters/1Number of partitions"}
  topology = {1, Choices{0, 1},
    Name "Parameters/2Create partition topology (BRep)?"}
  ghosts = {0, Choices{0, 1},
    Name "Parameters/3Create ghost cells?"}
  physicals = {0, Choices{0, 1},
    Name "Parameters/3Create new physical groups?"}
  write = {1, Choices {0, 1},
    Name "Parameters/3Write file to disk?"}
  split = {1, Choices {0, 1},
    Name "Parameters/4Write one file per partition?"}
];

If (partitioner > 0)
  // Should we create the boundary representation of the partition entities?
  Mesh.PartitionCreateTopology = topology;
  
  // Should we create ghost cells?
  Mesh.PartitionCreateGhostCells = ghosts;
  
  // Should we automatically create new physical groups on the partition entities?
  Mesh.PartitionCreatePhysicals = physicals;
  
  // Should we keep backward compatibility with pre-Gmsh 4, e.g. to save the mesh
  // in MSH2 format?
  Mesh.PartitionOldStyleMsh2 = 0;
  
  // Should we save one mesh file per partition?
  Mesh.PartitionSplitMeshFiles = split;
EndIf

If (partitioner == 1)
  // Use Metis to create N partitions
  PartitionMesh N;
  // Several options can be set to control Metis: `Mesh.MetisAlgorithm' (1:
  // Recursive, 2: K-way), `Mesh.MetisObjective' (1: min. edge-cut, 2:
  // min. communication volume), `Mesh.PartitionTriWeight' (weight of
  // triangles), `Mesh.PartitionQuadWeight' (weight of quads), ...
EndIf

If (partitioner == 2)
  // Use the `SimplePartition' plugin to create chessboard-like partitions
  Plugin(SimplePartition).NumSlicesX = N;
  Plugin(SimplePartition).NumSlicesY = 1;
  Plugin(SimplePartition).NumSlicesZ = 1;
  Plugin(SimplePartition).Run;
EndIf

