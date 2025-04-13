// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "gmsh.hpp"
#include "vtk.hpp"

namespace mfem
{

int BarycentricToGmshTet(int *b, int ref)
{
   int i = b[0];
   int j = b[1];
   int k = b[2];
   int l = b[3];
   bool ibdr = (i == 0);
   bool jbdr = (j == 0);
   bool kbdr = (k == 0);
   bool lbdr = (l == 0);
   if (ibdr && jbdr && kbdr)
   {
      return 0;
   }
   else if (jbdr && kbdr && lbdr)
   {
      return 1;
   }
   else if (ibdr && kbdr && lbdr)
   {
      return 2;
   }
   else if (ibdr && jbdr && lbdr)
   {
      return 3;
   }
   int offset = 4;
   if (jbdr && kbdr) // Edge DOF on j == 0 and k == 0
   {
      return offset + i - 1;
   }
   else if (kbdr && lbdr) // Edge DOF on k == 0 and l == 0
   {
      return offset + ref - 1 + j - 1;
   }
   else if (ibdr && kbdr) // Edge DOF on i == 0 and k == 0
   {
      return offset + 2 * (ref - 1) + ref - j - 1;
   }
   else if (ibdr && jbdr) // Edge DOF on i == 0 and j == 0
   {
      return offset + 3 * (ref - 1) + ref - k - 1;
   }
   else if (ibdr && lbdr) // Edge DOF on i == 0 and l == 0
   {
      return offset + 4 * (ref - 1) + ref - k - 1;
   }
   else if (jbdr && lbdr) // Edge DOF on j == 0 and l == 0
   {
      return offset + 5 * (ref - 1) + ref - k - 1;
   }

   // Recursive numbering for the faces
   offset += 6 * (ref - 1);
   if (kbdr)
   {
      int b_out[3];
      b_out[0] = j-1;
      b_out[1] = i-1;
      b_out[2] = ref - i - j - 1;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   else if (jbdr)
   {
      int b_out[3];
      b_out[0] = i-1;
      b_out[1] = k-1;
      b_out[2] = ref - i - k - 1;
      offset += (ref - 1) * (ref - 2) / 2;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   else if (ibdr)
   {
      int b_out[3];
      b_out[0] = k-1;
      b_out[1] = j-1;
      b_out[2] = ref - j - k - 1;
      offset += (ref - 1) * (ref - 2);
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   else if (lbdr)
   {
      int b_out[3];
      b_out[0] = ref-j-k-1;
      b_out[1] = j-1;
      b_out[2] = k-1;
      offset += 3 * (ref - 1) * (ref - 2) / 2;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }

   // Recursive numbering for interior
   {
      int b_out[4];
      b_out[0] = i-1;
      b_out[1] = j-1;
      b_out[2] = k-1;
      b_out[3] = ref - i - j - k - 1;
      offset += 2 * (ref - 1) * (ref - 2);
      return offset + BarycentricToGmshTet(b_out, ref-4);
   }
}

int CartesianToGmshQuad(int idx_in[], int ref)
{
   int i = idx_in[0];
   int j = idx_in[1];
   // Do we lie on any of the edges
   bool ibdr = (i == 0 || i == ref);
   bool jbdr = (j == 0 || j == ref);
   if (ibdr && jbdr) // Vertex DOF
   {
      return (i ? (j ? 2 : 1) : (j ? 3 : 0));
   }
   int offset = 4;
   if (jbdr) // Edge DOF on j==0 or j==ref
   {
      return offset + (j ? 3*ref - 3 - i : i - 1);
   }
   else if (ibdr) // Edge DOF on i==0 or i==ref
   {
      return offset + (i ? ref - 1 + j - 1 : 4*ref - 4 - j);
   }
   else // Recursive numbering for interior
   {
      int idx_out[2];
      idx_out[0] = i-1;
      idx_out[1] = j-1;
      offset += 4 * (ref - 1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
}

int CartesianToGmshHex(int idx_in[], int ref)
{
   int i = idx_in[0];
   int j = idx_in[1];
   int k = idx_in[2];
   // Do we lie on any of the edges
   bool ibdr = (i == 0 || i == ref);
   bool jbdr = (j == 0 || j == ref);
   bool kbdr = (k == 0 || k == ref);
   if (ibdr && jbdr && kbdr) // Vertex DOF
   {
      return (i ? (j ? (k ? 6 : 2) : (k ? 5 : 1)) :
              (j ? (k ? 7 : 3) : (k ? 4 : 0)));
   }
   int offset = 8;
   if (jbdr && kbdr) // Edge DOF on x-directed edge
   {
      return offset + (j ? (k ? 12*ref-12-i: 6*ref-6-i) :
                       (k ? 8*ref-9+i: i-1));
   }
   else if (ibdr && kbdr) // Edge DOF on y-directed edge
   {
      return offset + (k ? (i ? 10*ref-11+j: 9*ref-10+j) :
                       (i ? 3*ref-4+j: ref-2+j));
   }
   else if (ibdr && jbdr) // Edge DOF on z-directed edge
   {
      return offset + (i ? (j ? 6*ref-7+k: 4*ref-5+k) :
                       (j ? 7*ref-8+k: 2*ref-3+k));
   }
   else if (ibdr) // Face DOF on x-directed face
   {
      int idx_out[2];
      idx_out[0] = i ? j-1 : k-1;
      idx_out[1] = i ? k-1 : j-1;
      offset += (12 + (i ? 3 : 2) * (ref - 1)) * (ref - 1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   else if (jbdr) // Face DOF on y-directed face
   {
      int idx_out[2];
      idx_out[0] = j ? ref-i-1 : i-1;
      idx_out[1] = j ? k-1 : k-1;
      offset += (12 + (j ? 4 : 1) * (ref - 1)) * (ref - 1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   else if (kbdr) // Face DOF on z-directed face
   {
      int idx_out[2];
      idx_out[0] = k ? i-1 : j-1;
      idx_out[1] = k ? j-1 : i-1;
      offset += (12 + (k ? 5 : 0) * (ref - 1)) * (ref - 1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   else // Recursive numbering for interior
   {
      int idx_out[3];
      idx_out[0] = i-1;
      idx_out[1] = j-1;
      idx_out[2] = k-1;

      offset += (12 + 6 * (ref - 1)) * (ref - 1);
      return offset + CartesianToGmshHex(idx_out, ref-2);
   }
}

int WedgeToGmshPri(int idx_in[], int ref)
{
   int i = idx_in[0];
   int j = idx_in[1];
   int k = idx_in[2];
   int l = ref - i -j;
   bool ibdr = (i == 0);
   bool jbdr = (j == 0);
   bool kbdr = (k == 0 || k == ref);
   bool lbdr = (l == 0);
   if (ibdr && jbdr && kbdr)
   {
      return k ? 3 : 0;
   }
   else if (jbdr && lbdr && kbdr)
   {
      return k ? 4 : 1;
   }
   else if (ibdr && lbdr && kbdr)
   {
      return k ? 5 : 2;
   }
   int offset = 6;
   if (jbdr && kbdr)
   {
      return offset + (k ? 6 * (ref - 1) + i - 1: i - 1);
   }
   else if (ibdr && kbdr)
   {
      return offset + (k ? 7 * (ref -1) + j-1 : ref - 1 + j - 1);
   }
   else if (ibdr && jbdr)
   {
      return offset + 2 * (ref - 1) + k - 1;
   }
   else if (lbdr && kbdr)
   {
      return offset + (k ? 8 * (ref -1) + j - 1 : 3 * (ref - 1) + j - 1);
   }
   else if (jbdr && lbdr)
   {
      return offset + 4 * (ref - 1) + k - 1;
   }
   else if (ibdr && lbdr)
   {
      return offset + 5 * (ref - 1) + k - 1;
   }
   offset += 9 * (ref-1);
   if (kbdr) // Triangular faces at k=0 and k=ref
   {
      int b_out[3];
      b_out[0] = k ? i-1 : j-1;
      b_out[1] = k ? j-1 : i-1;
      b_out[2] = ref - i - j - 1;
      offset += k ? (ref-1)*(ref-2) / 2: 0;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   offset += (ref-1)*(ref-2);
   if (jbdr) // Quadrilateral face at j=0
   {
      int idx_out[2];
      idx_out[0] = i-1;
      idx_out[1] = k-1;
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   else if (ibdr) // Quadrilateral face at i=0
   {
      int idx_out[2];
      idx_out[0] = k-1;
      idx_out[1] = j-1;
      offset += (ref-1)*(ref-1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   else if (lbdr) // Quadrilateral face at l=ref-i-j=0
   {
      int idx_out[2];
      idx_out[0] = j-1;
      idx_out[1] = k-1;
      offset += 2*(ref-1)*(ref-1);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   offset += 3*(ref-1)*(ref-1);
   // The Gmsh Prism interiors are a tensor product of segments of order ref-2
   // and triangles of order ref-3
   {
      int b_out[3];
      b_out[0] = i-1;
      b_out[1] = j-1;
      b_out[2] = ref - i - j - 1;
      int ot = BarycentricToVTKTriangle(b_out, ref-3);
      int os = (k==1) ? 0 : (k == ref-1 ? 1 : k);
      return offset + (ref-1) * ot + os;
   }
}

int CartesianToGmshPyramid(int idx_in[], int ref)
{
   int i = idx_in[0];
   int j = idx_in[1];
   int k = idx_in[2];
   // Do we lie on any of the edges
   bool ibdr = (i == 0 || i == ref-k);
   bool jbdr = (j == 0 || j == ref-k);
   bool kbdr = (k == 0);
   if (ibdr && jbdr && kbdr)
   {
      return i ? (j ? 2 : 1): (j ? 3 : 0);
   }
   else if (k == ref)
   {
      return 4;
   }
   int offset = 5;
   if (jbdr && kbdr)
   {
      return offset + (j ? (6 * ref - 6 - i) : (i - 1));
   }
   else if (ibdr && kbdr)
   {
      return offset + (i ? (3 * ref - 4 + j) : (ref - 2 + j));
   }
   else if (ibdr && jbdr)
   {
      return offset + (i ? (j ? 6 : 4) : (j ? 7 : 2 )) * (ref-1) + k - 1;
   }
   offset += 8*(ref-1);
   if (jbdr)
   {
      int b_out[3];
      b_out[0] = j ? ref - i - k - 1 : i - 1;
      b_out[1] = k - 1;
      b_out[2] = (j ? i - 1 : ref - i - k - 1);
      offset += (j ? 3 : 0) * (ref - 1) * (ref - 2) / 2;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   else if (ibdr)
   {
      int b_out[3];
      b_out[0] = i ? j - 1: ref - j - k - 1;
      b_out[1] = k - 1;
      b_out[2] = (i ? ref - j - k - 1: j - 1);
      offset += (i ? 2 : 1) * (ref - 1) * (ref - 2) / 2;
      return offset + BarycentricToVTKTriangle(b_out, ref-3);
   }
   else if (kbdr)
   {
      int idx_out[2];
      idx_out[0] = k ? i-1 : j-1;
      idx_out[1] = k ? j-1 : i-1;
      offset += 2 * (ref - 1) * (ref - 2);
      return offset + CartesianToGmshQuad(idx_out, ref-2);
   }
   offset += (2 * (ref - 2) + (ref - 1)) * (ref - 1) ;
   {
      int idx_out[3];
      idx_out[0] = i-1;
      idx_out[1] = j-1;
      idx_out[2] = k-1;
      return offset + CartesianToGmshPyramid(idx_out, ref-3);
   }
}

void GmshHOSegmentMapping(int order, int *map)
{
   map[0] = 0;
   map[order] = 1;
   for (int i=1; i<order; i++)
   {
      map[i] = i + 1;
   }
}

void GmshHOTriangleMapping(int order, int *map)
{
   int b[3];
   int o = 0;
   for (b[1]=0; b[1]<=order; ++b[1])
   {
      for (b[0]=0; b[0]<=order-b[1]; ++b[0])
      {
         b[2] = order - b[0] - b[1];
         map[o] = BarycentricToVTKTriangle(b, order);
         o++;
      }
   }
}

void GmshHOQuadrilateralMapping(int order, int *map)
{
   int b[2];
   int o = 0;
   for (b[1]=0; b[1]<=order; b[1]++)
   {
      for (b[0]=0; b[0]<=order; b[0]++)
      {
         map[o] = CartesianToGmshQuad(b, order);
         o++;
      }
   }
}

void GmshHOTetrahedronMapping(int order, int *map)
{
   int b[4];
   int o = 0;
   for (b[2]=0; b[2]<=order; ++b[2])
   {

      for (b[1]=0; b[1]<=order-b[2]; ++b[1])
      {
         for (b[0]=0; b[0]<=order-b[1]-b[2]; ++b[0])
         {
            b[3] = order - b[0] - b[1] - b[2];
            map[o] = BarycentricToGmshTet(b, order);
            o++;
         }
      }
   }
}

void GmshHOHexahedronMapping(int order, int *map)
{
   int b[3];
   int o = 0;
   for (b[2]=0; b[2]<=order; b[2]++)
   {
      for (b[1]=0; b[1]<=order; b[1]++)
      {
         for (b[0]=0; b[0]<=order; b[0]++)
         {
            map[o] = CartesianToGmshHex(b, order);
            o++;
         }
      }
   }
}

void GmshHOWedgeMapping(int order, int *map)
{
   int b[3];
   int o = 0;
   for (b[2]=0; b[2]<=order; b[2]++)
   {
      for (b[1]=0; b[1]<=order; b[1]++)
      {
         for (b[0]=0; b[0]<=order - b[1]; b[0]++)
         {
            map[o] = WedgeToGmshPri(b, order);
            o++;
         }
      }
   }
}

void GmshHOPyramidMapping(int order, int *map)
{
   int b[3];
   int o = 0;
   for (b[2]=0; b[2]<=order; b[2]++)
   {
      for (b[1]=0; b[1]<=order - b[2]; b[1]++)
      {
         for (b[0]=0; b[0]<=order - b[2]; b[0]++)
         {
            map[o] = CartesianToGmshPyramid(b, order);
            o++;
         }
      }
   }
}

#ifdef MFEM_USE_MPI

typedef Triple<int,int,const std::vector<int>*> Tr_iivi;
bool mytriple_compare (Tr_iivi a, Tr_iivi b)
{
  return (a.one < b.one);
}
  
// ParGmshMesh implementation
// This function loads a parallel GMSH mesh (that has been read previously through
// Mesh::ReadGmshMesh and returns the parallel MFEM mesh corresponding to it.
ParGmshMesh::ParGmshMesh(MPI_Comm comm, std::string gmsh_file,
                         int refine, int generate_edges, bool fix_orientation)
{
   // Set the communicator for gtopo
   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
   gtopo.SetComm(comm);

   MyComm = comm;
   MPI_Comm_size(MyComm, &NRanks);
   MPI_Comm_rank(MyComm, &MyRank);

   std::ifstream ifs(gmsh_file);
   MFEM_VERIFY(ifs.good(), "Mesh file " << gmsh_file << " not found.");
   std::string mesh_type;
   //TODO : curved and high-order
   int curved = 0, read_gf=1;
   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
   ifs >> std::ws;
   getline(ifs, mesh_type);
   Mesh::ReadGmshMesh(ifs, curved, read_gf, true);

   
   ListOfIntegerSets  groups;
   IntegerSet         group;

   // The first group is the local one
   group.Recreate(1, &MyRank);
   groups.Insert(group);

   MFEM_ASSERT(Dim >= 3 || Dim < 1 || GetNFaces() == 0,
               "[proc " << MyRank << "]: invalid state");
   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
   PairIntVectMap *gmshE = gmesh->gmshE;
   VerMap &vinfo = gmesh->vertices_info;
   EltMap &einfo = gmesh->elts_info;
   Array<int> eleRanks;

   // Determine shared faces
   std::vector<int> sface_group, sfaces;
   if (Dim > 2)
   {
     for (auto const& elt : einfo) {
       const int no_elt = elt.first;
       const uint64_t DimEntity = elt.second.one[0];
       const uint64_t TagEntity = elt.second.one[1];
       const uint64_t elt_type  = elt.second.one[2];
       if (DimEntity == 2) {
	 const std::vector<int> &vlist = elt.second.two;
	 const auto &myPair = gmshE[DimEntity][TagEntity];
	 const std::vector<int> &procs = myPair.one;       
	 const std::vector<int> &phys  = myPair.two;       
	 MFEM_VERIFY (procs.size() > 0 || phys.size() > 0, "GMSH reader internal error");
	 bool contains = std::binary_search(procs.begin(),
					    procs.end(), MyRank);
	 int shared_psize = procs.size();
	 if (contains && shared_psize > 1) {
	   std::cerr << "shared_face " << no_elt << " elt_type " << elt_type << " tage " << \
	     TagEntity << " psize " << shared_psize << " contains " << contains << std::endl;
	   eleRanks.SetSize(shared_psize);
	   for (int i=0; i<shared_psize; i++) eleRanks[i]=procs[i];
	   group.Recreate(shared_psize, eleRanks);
	   sfaces.push_back(no_elt); //TODO: change tag here for a local id
	   sface_group.push_back(groups.Insert(group) - 1);
	 }
       }
     }
       
   }
   
   // Determine shared vertices
   std::vector<Tr_iivi> sverts;

   for (auto const& vit : vinfo)
     {
       const int gmsh_vindex = vit.first;
       const std::array<uint64_t,4> &myv = vit.second;
       const int ver = myv[0];
       const int DimEntity = myv[1];
       const int TagEntity = myv[2];
       const int data = myv[3];
       PairIntVectMap &myDimMap = gmshE[DimEntity];
       PairIntVectMap::const_iterator it = myDimMap.find(TagEntity);
       MFEM_VERIFY(it != myDimMap.end(), "Error reading GMSH file");
       // If more than one proc sharing this vertex, do stuff
       const std::vector<int> &shared_procs = (it->second.one);
       const int shared_psize = shared_procs.size();
       if (shared_psize > 1) {
	 sverts.push_back(Tr_iivi(gmsh_vindex,ver,&(it->second.one)));
//TOREMOVE	 if (MyRank == 3) std::cout << "list_sh_v1 gmsh_idx " << \
//TOREMOVE			    gmsh_vindex << " ver " << ver << " " << shared_psize << std::endl;
       }
     }
   // Sort sverts based on value of GMSH vertex numbering 
   std::sort (sverts.begin(), sverts.end(), mytriple_compare);

   // Fill svert_group and svert_list
   Array<int> svert_group(sverts.size());
   Array<int> svert_list(sverts.size());

   int j = 0;
   for (auto const &tripl : sverts)
     {
       const int ver = tripl.two;
       const std::vector<int> *shared_procs = tripl.three;
       const int shared_psize = shared_procs->size();
       
       eleRanks.SetSize(shared_psize);
       for (int i=0; i<shared_psize; i++) eleRanks[i]=(*shared_procs)[i];
       group.Recreate(shared_psize, eleRanks);
       svert_list[j] = ver;
       svert_group[j] = (groups.Insert(group) - 1);
       j++;
       //TOREMOVE       if (MyRank == 3) std::cout << "list_sh_v2 " << ver << " " << shared_psize << std::endl;
     }
   // Erase all data within svert
   sverts.resize(0);
   
   group_stria.MakeI(groups.Size()-1);
   group_squad.MakeI(groups.Size()-1);
   group_stria.MakeJ();
   group_squad.MakeJ();
   group_stria.ShiftUpI();
   group_squad.ShiftUpI();

   group_sedge.MakeI(groups.Size()-1);
   group_sedge.MakeJ();
   group_sedge.ShiftUpI();

   // Build group_svert
   group_svert.MakeI(groups.Size()-1);
   for (int i = 0; i < svert_group.Size(); i++)
   {
      group_svert.AddAColumnInRow(svert_group[i]);
   }
   group_svert.MakeJ();
   for (int i = 0; i < svert_group.Size(); i++)
   {
      group_svert.AddConnection(svert_group[i], i);
   }
   group_svert.ShiftUpI();
 
   shared_edges.SetSize(0);  
   sedge_ledge. SetSize(0);

   svert_lvert.SetSize(svert_list.Size());
   for (int i = 0; i < svert_list.Size(); i++)
   {
      svert_lvert[i] = svert_list[i];
   }

   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
   MPI_Barrier(MyComm);
   // Build the group communication topology
   gtopo.Create(groups, 822);

   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
   // Determine sedge_ledge and sface_lface
   FinalizeParTopo();

      // Set nodes for higher order mesh
   if (curved) // curved mesh
   {
     //???
   }
   Finalize(refine, fix_orientation);
   std::cerr  << __FILE__ << " log " << __LINE__ << std::endl;
}

#endif  // MFEM_USE_MPI


} // namespace mfem
