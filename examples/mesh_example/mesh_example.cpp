
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>

#include "bo/mesh.hpp"
#include "bo/io/mesh_io.hpp"

using namespace boost::filesystem;
using namespace bo::io;


// Directory where test data is stored.
std::string DataDirectory;

// Paths to data used in this example.
struct Paths
{
    void SetUp()
    {
        InDir = path(DataDirectory);
        PlyMeshPath01 = InDir / "mesh01.ply";
        PlyMeshPath02 = InDir / "mesh02.ply";

        OutDir = InDir / "results";
        create_directory(OutDir);
        PlyMeshOutPath = OutDir / "mesh_result.ply";
    }

    path InDir;
    path PlyMeshPath01;
    path PlyMeshPath02;

    path OutDir;
    path PlyMeshOutPath;
};

Paths paths;


void AssertPathExists(const path& filepath)
{
    BOOST_ASSERT(exists(filepath));
    BO_UNUSED(filepath);
}

void Join()
{
    typedef bo::Mesh<float> Mesh;

    AssertPathExists(paths.PlyMeshPath01);
    AssertPathExists(paths.PlyMeshPath02);
    Mesh mesh01 = mesh_from_ply(paths.PlyMeshPath01.string());
    Mesh mesh02 = mesh_from_ply(paths.PlyMeshPath02.string());

    mesh01.join(mesh02);
    mesh_to_ply(mesh01, paths.PlyMeshOutPath.string());
}


int main(int argc, char* argv[])
{
    // Extract directory with data for tests from command-line, or apply default value.
    if (argc > 1)
        DataDirectory.assign(argv[1]);
    else if (argc == 1)
        // Apply default value which is "./data/mesh" directory.
        DataDirectory.assign((initial_path() /= "data/mesh").string());
    else
        DataDirectory.assign("");

    paths.SetUp();

    Join();

    //bo::Mesh::Normal normal = mesh.get_vertex_normal(1);

    ////mesh.to_ply(std::string("..\\debug\\pts2.ply"));
}
