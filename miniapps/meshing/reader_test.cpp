#include <sys/resource.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include "../common/mfem-common.hpp"
#include "mfem.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
    const char *mesh_file = nullptr;
    const char *mesh_output_file = nullptr;
    const char *vtk_output_file = nullptr;
    bool debug_mode = false;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&mesh_output_file, "-mesh", "--mesh-output", "Output mesh file.");
    args.AddOption(&vtk_output_file, "-vtk", "--vtk-output", "Output VTK file.");
    args.AddOption(&debug_mode, "-dbg", "--debug", "-no-dbg", "--no-debug", "Enable or disable debug output.");
    args.Parse();

    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }

    args.PrintOptions(std::cout);

    if (!mesh_file)
    {
        fprintf(stderr, "Erreur : Aucun fichier de maillage spécifié avec -m.\n");
        return 1;
    }

    Mesh mesh(mesh_file,1,1);
    mesh.PrintInfo();
    if (debug_mode)
    {
        printf("DEBUG MODE ACTIVE \n");
        for (int i = 0; i < mesh.GetNV(); i++)
        {
            const double *node = mesh.GetVertex(i);
            printf("Noeud %d : (%.6f, %.6f, %.6f)\n", i, node[0], node[1], node[2]);
        }

        for (int i = 0; i < mesh.GetNE(); i++)
        {
            Element *el = mesh.GetElement(i);
            const int *vertices = el->GetVertices();
            int num_vertices = el->GetNVertices();

            printf("Élément %d : Sommets = [", i);
            for (int j = 0; j < num_vertices; j++)
            {
                printf(" %d", vertices[j]);
            }
            printf(" ]\n");
        }
    }

    
    printf("Nombre de sommets : %d\n", mesh.GetNV());
    printf("Nombre d'arêtes : %d\n", mesh.GetNEdges());
    printf("Nombre de faces : %d\n", mesh.GetNFaces());
    printf("Nombre d'éléments : %d\n", mesh.GetNE());

    if (mesh_output_file)
    {
        ofstream mesh_out(mesh_output_file);
        mesh.Print(mesh_out);
        mesh_out.close();
        printf("Maillage exporté dans : %s\n", mesh_output_file);
    }

    if (vtk_output_file)
    {
        ofstream vtk_out(vtk_output_file);
        mesh.PrintVTK(vtk_out);
        vtk_out.close();
        printf("Maillage exporté au format VTK dans : %s\n", vtk_output_file);
    }
    return 0;
}
