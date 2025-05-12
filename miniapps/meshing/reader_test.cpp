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
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&mesh_output_file, "-mesh", "--mesh-output",
                   "Output mesh file.");
    args.AddOption(&vtk_output_file, "-vtk", "--vtk-output",
                   "Output VTK file.");
    args.AddOption(&debug_mode, "-dbg", "--debug", "-no-dbg", "--no-debug",
                   "Enable or disable debug output.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    if (!mesh_file)
    {
        cerr << "Erreur : Aucun fichier de maillage spécifié avec -m." << endl;
        return 1;
    }

    Mesh mesh(mesh_file);

    if (debug_mode)
    {
        for (int i = 0; i < mesh.GetNV(); i++)
        {
            const double *node = mesh.GetVertex(i);
            cout << "Noeud " << i << " : (" << node[0] << ", " << node[1] << ", "
                 << node[2] << ")" << endl;
        }

        for (int i = 0; i < mesh.GetNE(); i++)
        {
            Element *el = mesh.GetElement(i);
            const int *vertices = el->GetVertices();
            int num_vertices = el->GetNVertices();

            cout << "Element " << i << " : Sommets = [ ";
            for (int j = 0; j < num_vertices; j++)
            {
                cout << vertices[j] << " ";
            }
            cout << "]" << endl;
        }
    }

    mesh.PrintInfo();
    cout << "Nombre de sommets : " << mesh.GetNV() << endl;
    cout << "Nombre d'arêtes : " << mesh.GetNEdges() << endl;
    cout << "Nombre de faces : " << mesh.GetNFaces() << endl;
    cout << "Nombre d'éléments : " << mesh.GetNE() << endl;

    if (mesh_output_file)
    {
        ofstream mesh_out(mesh_output_file);
        mesh.Print(mesh_out);
        mesh_out.close();
        cout << "Maillage exporté dans : " << mesh_output_file << endl;
    }

    if (vtk_output_file)
    {
        ofstream vtk_out(vtk_output_file);
        mesh.PrintVTK(vtk_out);
        vtk_out.close();
        cout << "Maillage exporté au format VTK dans : " << vtk_output_file << endl;
    }

    return 0;
}
