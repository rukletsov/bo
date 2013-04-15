
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>

#include "bo/methods/complex_propagation.hpp"
#include "bo/methods/triangulation.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

#include "bo/extended_std.hpp"
#include <boost/format.hpp>

using namespace boost::filesystem;
using namespace bo::methods::surfaces;
using namespace bo::io;

typedef bo::Mesh<float> Mesh3D;
typedef ComplexPropagation<float> PropagAlgo;
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
        SheepInDir = InDir / "sheep";
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
    path SheepInDir;
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
    BO_UNUSED(filepath);
}

void PropagateClosed()
{
    AssertPathExists(paths.RawClosedPath);
    PropagAlgo::Image2D test_image = load_raw_image_8bpps<float>(
                paths.RawClosedPath.string(), 512, 512);

    PropagAlgo::Ptr tiling_ptr = PropagAlgo::from_raw_image(test_image, 3.f, 10.f, 0.5f, 20.f);

    tiling_ptr->propagate();
    Mesh3D mesh = Mesh3D::from_vertices(tiling_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlyClosedOutPath.string());
}

void PropagateFemur01()
{
    AssertPathExists(paths.PlyFemurPath01);
    Mesh3D test_mesh = mesh_from_ply(paths.PlyFemurPath01.string());
    PropagAlgo::Ptr tiling_ptr = PropagAlgo::from_mesh(test_mesh, 3.f, 7.f, 0.5f, 20.f);

    tiling_ptr->propagate();
    Mesh3D mesh = Mesh3D::from_vertices(tiling_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlyFemurOutPath01.string());
}

void PropagateSheep()
{
    AssertPathExists(paths.PlySheepPath);
    Mesh3D test_mesh = mesh_from_ply(paths.PlySheepPath.string());
    PropagAlgo::Ptr tiling_ptr = PropagAlgo::from_mesh(test_mesh, 5.f, 10.f, 0.5f, 20.f);

    tiling_ptr->propagate();
    Mesh3D mesh = Mesh3D::from_vertices(tiling_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlySheepOutPath.string());
}

void ChrisitiansenFemur()
{
    AssertPathExists(paths.PlyFemurPath01);
    Mesh3D test_mesh1 = mesh_from_ply(paths.PlyFemurPath01.string());
    PropagAlgo::Ptr tiling_ptr1 = PropagAlgo::from_mesh(test_mesh1, 3.f, 7.f, 0.5f, 20.f);

    AssertPathExists(paths.PlyFemurPath02);
    Mesh3D test_mesh2 = mesh_from_ply(paths.PlyFemurPath02.string());
    PropagAlgo::Ptr tiling_ptr2 = PropagAlgo::from_mesh(test_mesh2, 3.f, 7.f, 0.5f, 20.f);

    tiling_ptr1->propagate();
    tiling_ptr2->propagate();

    TriangAlgo triang(tiling_ptr1->contour(), !tiling_ptr1->has_hole(),
                      tiling_ptr2->contour(), !tiling_ptr2->has_hole());
    TriangAlgo::MeshPtr mesh_ptr = triang.christiansen();
    mesh_to_ply(*mesh_ptr, paths.PlyFemurOutPath0102.string());
}

void ChrisitiansenClosed()
{
    AssertPathExists(paths.PlyClosedPath01);
    Mesh3D test_mesh1 = mesh_from_ply(paths.PlyClosedPath01.string());
    PropagAlgo::Ptr tiling_ptr1 = PropagAlgo::from_mesh(test_mesh1, 3.f, 10.f, 0.5f, 20.f);

    AssertPathExists(paths.PlyClosedPath02);
    Mesh3D test_mesh2 = mesh_from_ply(paths.PlyClosedPath02.string());
    PropagAlgo::Ptr tiling_ptr2 = PropagAlgo::from_mesh(test_mesh2, 3.f, 10.f, 0.5f, 20.f);

    tiling_ptr1->propagate();
    tiling_ptr2->propagate();

    TriangAlgo triang(tiling_ptr1->contour(), !tiling_ptr1->has_hole(),
                      tiling_ptr2->contour(), !tiling_ptr2->has_hole());
    TriangAlgo::MeshPtr mesh_ptr = triang.christiansen();
    mesh_to_ply(*mesh_ptr, paths.PlyClosedOutMeshPath.string());
}

void ChrisitiansenFemurFull()
{
    typedef std::vector<path> ContourData;
    typedef PropagAlgo::Points3DConstPtrs PlaneData;

    ContourData contour_data;
    PlaneData plane_data;
    Mesh3D result_mesh;

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
    BOOST_FOREACH (path contour_path, contour_data)
    {
        AssertPathExists(contour_path);
        Mesh3D mesh = mesh_from_ply(contour_path.string());
        plane_data.push_back(boost::make_shared<PropagAlgo::Points3D>(mesh.get_all_vertices()));
    }

    // Prepare weights for all slices.
    PropagAlgo::Weights weights;
    weights.push_back(0.8f);
    weights.push_back(0.5f);
    weights.push_back(0.3f);
    weights.push_back(0.1f);

    // Run propagation for every plane.
    PropagAlgo::Ptrs props = PropagAlgo::create(plane_data, weights, 3.f, 7.f, 0.5f, 20.f);
    BOOST_FOREACH (PropagAlgo::Ptr propagator, props)
    { propagator->propagate(); }

    // Tile pair of contours and join it with the result mesh.
    for (PropagAlgo::Ptrs::const_iterator it = props.begin() + 1; it != props.end(); ++it)
    {
        TriangAlgo triang((*(it - 1))->contour(), !((*(it - 1))->has_hole()),
                          (*it)->contour(), !((*it)->has_hole()));
        TriangAlgo::MeshPtr mesh_ptr = triang.christiansen();
        result_mesh.join(*mesh_ptr);
    }

    mesh_to_ply(result_mesh, paths.PlyFemurOutMeshPath.string());
}

void ChrisitiansenSheepFull()
{
    typedef std::vector<path> ContourData;
    typedef PropagAlgo::Points3DConstPtrs PlaneData;

    ContourData contour_data;
    PlaneData plane_data;
    Mesh3D result_mesh;

    // Load planes paths.
    AssertPathExists(paths.SheepInDir);
    for (directory_iterator it(paths.SheepInDir); it != directory_iterator(); ++it)
    {
        if (is_regular_file(*it))
            contour_data.push_back(it->path());
    }

    // Sort paths in lsexicographic order, since directory iteration is not ordered
    // on some file systems.
    std::sort(contour_data.begin(), contour_data.end());

    // Extract contours from contour data.
    BOOST_FOREACH (path contour_path, contour_data)
    {
        AssertPathExists(contour_path);
        Mesh3D mesh = mesh_from_ply(contour_path.string());
        plane_data.push_back(boost::make_shared<PropagAlgo::Points3D>(mesh.get_all_vertices()));
    }

    // Prepare weights for all slices.
    PropagAlgo::Weights weights;
    weights.push_back(0.8f);
    weights.push_back(0.5f);

    // Run propagation for every plane.
    PropagAlgo::Ptrs props = PropagAlgo::create(plane_data, weights, 2.f, 5.f, 0.4f, 15.f);
    BOOST_FOREACH (PropagAlgo::Ptr propagator, props)
    { propagator->propagate(); }

    // Tile pair of contours and join it with the result mesh.
    for (PropagAlgo::Ptrs::const_iterator it = props.begin() + 1; it != props.end(); ++it)
    {
        PropagAlgo::Ptr cur_contour = *it;
        PropagAlgo::Ptr prev_contour = *(it - 1);

        TriangAlgo triang(prev_contour->contour(), !(prev_contour->has_hole()),
                          cur_contour->contour(), !(cur_contour->has_hole()));
        Mesh3D mesh = *triang.christiansen();
        result_mesh.join(mesh);
    }

    mesh_to_ply(result_mesh, paths.PlySheepOutPath.string());
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
    PropagateClosed();
    ChrisitiansenClosed();
    ChrisitiansenFemur();
    ChrisitiansenFemurFull();
    ChrisitiansenSheepFull();
}
