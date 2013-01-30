
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>

#include "bo/methods/parallel_planes_tiling.hpp"
#include "bo/methods/triangulation.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

#include "bo/extended_std.hpp"
#include <boost/format.hpp>

using namespace boost::filesystem;
using namespace bo::methods::surfaces;
using namespace bo::io;

typedef MinSpanPropagation<float> TilingAlgo;
typedef Triangulation<float> TriangAlgo;

// Directory where test data is stored.
std::string DataDirectory;

// Paths to data used in this example.
struct Paths
{
    void SetUp()
    {
        path InDir = path(DataDirectory);

        RawClosedPath = InDir / "2d_contour4_512x512_8bit.raw";
        PlyClosedPath01 = InDir / "2d_contour4_512x512_8bit.ply";
        PlyClosedPath02 = InDir / "2d_contour4_512x512_8bit_transformed.ply";
        PlySheepPath = InDir / "sheep_plane.ply";
        FemurInDir = InDir / "femur";
        PlyFemurPath01 = FemurInDir / "femur_plane01.ply";
        PlyFemurPath02 = FemurInDir / "femur_plane02.ply";

        OutDir = InDir / "results";
        create_directory(OutDir);
        PlyClosedOutPath = OutDir / "2d_contour4_512x512_8bit_result.ply";
        PlyClosedOutMeshPath = OutDir / "2d_contour4_512x512_8bit_mesh.ply";
        PlySheepOutPath = OutDir / "sheep_result.ply";
        FemurOutDir = OutDir / "femur";
        create_directory(FemurOutDir);
        PlyFemurOutPath01 = OutDir / "femur_plane01_result.ply";
        PlyFemurOutPath02 = OutDir / "femur_plane02_result.ply";
        PlyFemurOutPath0102 = OutDir / "femur_result_0102.ply";
        PlyFemurOutMeshPath = OutDir / "femur_result.ply";
    }

    path RawClosedPath;
    path PlyClosedPath01;
    path PlyClosedPath02;
    path PlySheepPath;
    path FemurInDir;
    path PlyFemurPath01;
    path PlyFemurPath02;

    path OutDir;
    path PlyClosedOutPath;
    path PlyClosedOutMeshPath;
    path PlySheepOutPath;
    path FemurOutDir;
    path PlyFemurOutPath01;
    path PlyFemurOutPath02;
    path PlyFemurOutPath0102;
    path PlyFemurOutMeshPath;
};

Paths paths;


void AssertPathExists(const path& filepath)
{
    BOOST_ASSERT(exists(filepath));
}

void PropagateClosed()
{
    TilingAlgo tiling;

    AssertPathExists(paths.RawClosedPath);
    TilingAlgo::Image2D test_image = load_raw_image_8bpps<float>(
                paths.RawClosedPath.string(), 512, 512);
    TilingAlgo::ParallelPlanePtr plane_data = tiling.load_plane(test_image);

    TilingAlgo::PropagationResult contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = TilingAlgo::Mesh::from_vertices(contour.points.get());
    mesh_to_ply(mesh, paths.PlyClosedOutPath.string());
}

void PropagateFemur01()
{
    TilingAlgo tiling;

    AssertPathExists(paths.PlyFemurPath01);
    TilingAlgo::Mesh test_mesh = mesh_from_ply(paths.PlyFemurPath01.string());
    TilingAlgo::ParallelPlanePtr plane_data =
        boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh.get_all_vertices());

    TilingAlgo::PropagationResult contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = TilingAlgo::Mesh::from_vertices(contour.points.get());
    mesh_to_ply(mesh, paths.PlyFemurOutPath01.string());
}

void PropagateSheep()
{
    TilingAlgo tiling;

    AssertPathExists(paths.PlySheepPath);
    TilingAlgo::Mesh test_mesh = mesh_from_ply(paths.PlySheepPath.string());
    TilingAlgo::ParallelPlanePtr plane_data =
            boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh.get_all_vertices());

    TilingAlgo::PropagationResult contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = TilingAlgo::Mesh::from_vertices(contour.points.get());
    mesh_to_ply(mesh, paths.PlySheepOutPath.string());
}

void ChrisitiansenFemur()
{
    TilingAlgo tiling;
    TriangAlgo triang;

    AssertPathExists(paths.PlyFemurPath01);
    TilingAlgo::Mesh test_mesh1 = mesh_from_ply(paths.PlyFemurPath01.string());
    TilingAlgo::ParallelPlanePtr plane_data1 =
            boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh1.get_all_vertices());

    AssertPathExists(paths.PlyFemurPath02);
    TilingAlgo::Mesh test_mesh2 = mesh_from_ply(paths.PlyFemurPath02.string());
    TilingAlgo::ParallelPlanePtr plane_data2 =
            boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh2.get_all_vertices());

    TilingAlgo::PropagationResult contour1 = tiling.propagate(plane_data1, 0.5f);
    TilingAlgo::PropagationResult contour2 = tiling.propagate(plane_data2, 0.2f);

    TilingAlgo::Mesh mesh = triang.christiansen_triangulation(contour1.points,
        !contour1.has_hole, contour2.points, !contour2.has_hole);
    mesh_to_ply(mesh, paths.PlyFemurOutPath0102.string());
}

void ChrisitiansenClosed()
{
    TilingAlgo tiling;
    TriangAlgo triang;

    AssertPathExists(paths.PlyClosedPath01);
    TilingAlgo::Mesh test_mesh1 = mesh_from_ply(paths.PlyClosedPath01.string());
    TilingAlgo::ParallelPlanePtr plane_data1 =
            boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh1.get_all_vertices());

    AssertPathExists(paths.PlyClosedPath02);
    TilingAlgo::Mesh test_mesh2 = mesh_from_ply(paths.PlyClosedPath02.string());
    TilingAlgo::ParallelPlanePtr plane_data2 =
            boost::make_shared<TilingAlgo::ParallelPlane>(test_mesh2.get_all_vertices());

    TilingAlgo::PropagationResult contour1 = tiling.propagate(plane_data1, 0.5f);
    TilingAlgo::PropagationResult contour2 = tiling.propagate(plane_data2, 0.2f);

    TilingAlgo::Mesh mesh = triang.christiansen_triangulation(contour1.points,
        !(contour1.has_hole), contour2.points, !(contour2.has_hole));
    mesh_to_ply(mesh, paths.PlyClosedOutMeshPath.string());
}

void ChrisitiansenFemurFull()
{
    typedef std::vector<path> ContourData;
    typedef std::vector<TilingAlgo::PropagationResult> Contours;

    TilingAlgo tiling;
    TriangAlgo triang;
    ContourData contour_data;
    Contours contours;
    TilingAlgo::Mesh result_mesh;

    // Load planes paths.
    AssertPathExists(paths.FemurInDir);
    for (directory_iterator it(paths.FemurInDir); it != directory_iterator(); ++it)
    {
        if ((is_regular_file(*it)) &&
            (it->path().filename().string().substr(0, 11) == "femur_plane"))
            contour_data.push_back(it->path());
    }

    // Sort paths in lexicographic order, since directory iteration is not ordered
    // on some file systems.
    std::sort(contour_data.begin(), contour_data.end());

    // Extract contours from contour data.
    for (ContourData::const_iterator it = contour_data.begin(); it != contour_data.end(); ++it)
    {
        AssertPathExists(*it);
        TilingAlgo::Mesh mesh = mesh_from_ply(it->string());
        TilingAlgo::ParallelPlanePtr plane_data =
                boost::make_shared<TilingAlgo::ParallelPlane>(mesh.get_all_vertices());
        TilingAlgo::PropagationResult contour = tiling.propagate(plane_data, 0.5f);

        contours.push_back(contour);
    }

    // Tile pair of contours and join it with the result mesh.
    for (Contours::const_iterator it = contours.begin() + 1; it != contours.end(); ++it)
    {
        TilingAlgo::Mesh mesh = triang.christiansen_triangulation((it - 1)->points,
            !(it - 1)->has_hole, it->points, !it->has_hole);
        result_mesh.join(mesh);
    }

    mesh_to_ply(result_mesh, paths.PlyFemurOutMeshPath.string());
}


int main(int argc, char* argv[])
{
    // Extract directory with data for tests from command-line, or apply default value.
    if (argc > 1)
        DataDirectory.assign(argv[1]);
    else if (argc == 1)
        // Apply default value which is "./data/tiling" directory.
        DataDirectory.assign((initial_path() /= "data/tiling").string());
    else
        DataDirectory.assign("");

    paths.SetUp();

    PropagateSheep();
    ChrisitiansenClosed();
    ChrisitiansenFemur();
    ChrisitiansenFemurFull();
}
