SetFactory("OpenCASCADE");
// 0 for tetrahedra, 1 for hexahedra
tet_or_hex = 0;
// We now define several constants to fine-tune how the mesh will be partitioned
DefineConstant[
  partitioner = {0, Choices{0="None", 1="Metis", 2="SimplePartition"},
    Name "Parameters/0Mesh partitioner"}
  N = {2, Min 1, Max 256, Step 1,
    Name "Parameters/1Number of partitions"}
  topology = {1, Choices{0, 1},
    Name "Parameters/2Create partition topology (BRep)?"}
  ghosts = {0, Choices{0, 1},
    Name "Parameters/3Create ghost cells?"}
  physicals = {1, Choices{0, 1},
    Name "Parameters/3Create new physical groups?"}
  write = {1, Choices {0, 1},
    Name "Parameters/3Write file to disk?"}
  split = {1, Choices {0, 1},
    Name "Parameters/4Write one file per partition?"}
  periodic = {0, Choices{0="NotPeriodic", 1="Periodic"},
    Name "Parameters/5Mesh periodicity"}
];

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Characteristic Length {:} = .5;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

If (periodic > 0)
  Periodic Curve {1} = {-3};
  Periodic Curve {2} = {-4};
EndIf

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};

If (tet_or_hex > 0)
   Recombine Surface {1};
   out[] = Extrude {0, 0, 1} { Surface{1}; Layers{4}; Recombine; };
Else
   out[] = Extrude {0, 0, 1} { Surface{1}; Layers{2};  };
EndIf

Physical Volume(1) = {1}; 
//Physical Surface(1) = {1,out[0],out[2],out[3],out[4],out[5]};

Mesh 3;

If (periodic >0)
    Periodic Surface {out[0]} = {1} Translate {0, 0, 1};
    Periodic Surface {out[4]} = {out[2]} Translate {0, 1, 0};
    Periodic Surface {out[3]} = {out[5]} Translate {1, 0, 0};
EndIf

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
  Plugin(SimplePartition).NumSlicesX = 1;
  Plugin(SimplePartition).NumSlicesY = 1;
  Plugin(SimplePartition).NumSlicesZ = N;
  Plugin(SimplePartition).Run;
EndIf


