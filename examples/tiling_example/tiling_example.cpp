
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>
#include <boost/chrono.hpp>

#include "bo/config.hpp"
#include "bo/surfaces/complex_propagation.hpp"
#include "bo/surfaces/triangulation.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

using namespace boost::filesystem;
using namespace bo::surfaces;
using namespace bo::io;

typedef bo::Mesh<float> FloatMesh;
typedef ComplexPropagation<float> Propagator;
typedef Triangulation<float> Triangulator;

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
    Propagator::Image2D test_image = load_raw_image_8bpps<float>(
                paths.RawClosedPath.string(), 512, 512);

    Propagator::Ptr contour_descriptor_ptr = Propagator::from_raw_image(test_image,
            3.f, 10.f, 0.3f, 0.4f, 20.f);
    contour_descriptor_ptr->propagate();

    FloatMesh mesh = FloatMesh::from_vertices(contour_descriptor_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlyClosedOutPath.string());
}

void PropagateFemur01()
{
    AssertPathExists(paths.PlyFemurPath01);
    FloatMesh test_mesh = mesh_from_ply(paths.PlyFemurPath01.string());

    Propagator::Ptr contour_descriptor_ptr = Propagator::from_mesh(test_mesh,
            3.f, 7.f, 0.3f, 0.4f, 20.f);
    contour_descriptor_ptr->propagate();

    FloatMesh mesh = FloatMesh::from_vertices(contour_descriptor_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlyFemurOutPath01.string());
}

void PropagateSheep()
{
    AssertPathExists(paths.PlySheepPath);
    FloatMesh test_mesh = mesh_from_ply(paths.PlySheepPath.string());

    Propagator::Ptr contour_descriptor_ptr = Propagator::from_mesh(test_mesh,
            5.f, 10.f, 0.3f, 0.4f, 20.f);
    contour_descriptor_ptr->propagate();

    FloatMesh mesh = FloatMesh::from_vertices(contour_descriptor_ptr->contour().get());
    mesh_to_ply(mesh, paths.PlySheepOutPath.string());
}

void ChrisitiansenFemur()
{
    AssertPathExists(paths.PlyFemurPath01);
    FloatMesh test_mesh1 = mesh_from_ply(paths.PlyFemurPath01.string());
    Propagator::Ptr contour_descriptor_ptr1 = Propagator::from_mesh(test_mesh1,
            3.f, 7.f, 0.3f, 0.4f, 20.f);

    AssertPathExists(paths.PlyFemurPath02);
    FloatMesh test_mesh2 = mesh_from_ply(paths.PlyFemurPath02.string());
    Propagator::Ptr contour_descriptor_ptr2 = Propagator::from_mesh(test_mesh2,
            3.f, 7.f, 0.3f, 0.4f, 20.f);

    contour_descriptor_ptr1->propagate();
    contour_descriptor_ptr2->propagate();

    FloatMesh result_mesh = Triangulator::christiansen(contour_descriptor_ptr1,
                                                       contour_descriptor_ptr2);
    mesh_to_ply(result_mesh, paths.PlyFemurOutPath0102.string());
}

void ChrisitiansenClosed()
{
    AssertPathExists(paths.PlyClosedPath01);
    FloatMesh test_mesh1 = mesh_from_ply(paths.PlyClosedPath01.string());
    Propagator::Ptr contour_descriptor_ptr1 = Propagator::from_mesh(test_mesh1,
            3.f, 10.f, 0.3f, 0.4f, 20.f);

    AssertPathExists(paths.PlyClosedPath02);
    FloatMesh test_mesh2 = mesh_from_ply(paths.PlyClosedPath02.string());
    Propagator::Ptr contour_descriptor_ptr2 = Propagator::from_mesh(test_mesh2,
            3.f, 10.f, 0.3f, 0.4f, 20.f);

    contour_descriptor_ptr1->propagate();
    contour_descriptor_ptr2->propagate();

    FloatMesh result_mesh = Triangulator::christiansen(contour_descriptor_ptr1,
                                                       contour_descriptor_ptr2);
    mesh_to_ply(result_mesh, paths.PlyClosedOutMeshPath.string());
}

void ChrisitiansenFemurFull()
{
    typedef std::vector<path> PlanesPaths;
    typedef Propagator::Points3DConstPtrs PlanesData;

    PlanesPaths planes_paths;
    PlanesData planes_data;

    // Load planes paths.
    AssertPathExists(paths.FemurInDir);
    for (directory_iterator it(paths.FemurInDir); it != directory_iterator(); ++it)
    {
        if ((is_regular_file(*it)) &&
            (it->path().filename().string().substr(0, 11) == "femur_plane"))
            planes_paths.push_back(it->path());
    }

    // Sort paths in lexicographic order, since directory iteration is not ordered
    // on some file systems.
    std::sort(planes_paths.begin(), planes_paths.end());

    // Extract contours from contour data.
    BOOST_FOREACH (path contour_path, planes_paths)
    {
        AssertPathExists(contour_path);
        FloatMesh mesh = mesh_from_ply(contour_path.string());
        planes_data.push_back(boost::make_shared<Propagator::Points3D>(mesh.get_all_vertices()));
    }

    // Prepare weights for all slices.
    Propagator::Weights weights;
    weights.push_back(0.8f);
    weights.push_back(0.5f);
    weights.push_back(0.3f);
    weights.push_back(0.1f);

    // Run propagation for every plane.
    Propagator::Ptrs contour_descriptors = Propagator::create(planes_data, weights,
            3.f, 7.f, 0.3f, 0.4f, 20.f);
    BOOST_FOREACH (Propagator::Ptr contour_descriptor, contour_descriptors)
    { contour_descriptor->propagate(); }

    // Triangulate the contours and save the outpur mesh.
    FloatMesh result_mesh = Triangulator::christiansen(contour_descriptors);
    mesh_to_ply(result_mesh, paths.PlyFemurOutMeshPath.string());
}

void ChrisitiansenSheepFull()
{
    typedef std::vector<path> PlanesPaths;
    typedef Propagator::Points3DConstPtrs PlanesData;

    PlanesPaths planes_paths;
    PlanesData planes_data;

    // Load planes paths.
    AssertPathExists(paths.SheepInDir);
    for (directory_iterator it(paths.SheepInDir); it != directory_iterator(); ++it)
    {
        if (is_regular_file(*it) && (it->path().extension() == path(".ply")))
            planes_paths.push_back(it->path());
    }

    // Sort paths in lsexicographic order, since directory iteration is not ordered
    // on some file systems.
    std::sort(planes_paths.begin(), planes_paths.end());

    // Extract contours from contour data.
    BOOST_FOREACH (path contour_path, planes_paths)
    {
        AssertPathExists(contour_path);
        FloatMesh mesh = mesh_from_ply(contour_path.string());
        planes_data.push_back(boost::make_shared<Propagator::Points3D>(mesh.get_all_vertices()));
    }

    typedef boost::chrono::steady_clock BoostTimer;
    BoostTimer::time_point start = BoostTimer::now();

    // 1. Prepare propagators for all slices. Use no neighbouring information.
//    Propagator::Ptrs contour_descriptors = Propagator::create(planes_data,
//            2.f, 5.f, 0.3f, 0.4f, 15.f);

    // 2. Prepare weights and propagators for all slices.
    Propagator::Weights weights;
    weights.push_back(0.8f);
    weights.push_back(0.5f);
    Propagator::Ptrs contour_descriptors = Propagator::create(planes_data, weights,
            2.f, 5.f, 0.3f, 0.4f, 15.f);

    BoostTimer::time_point after_construction = BoostTimer::now();

    // Run propagation for every plane.
    BOOST_FOREACH (Propagator::Ptr contour_descriptor, contour_descriptors)
    { contour_descriptor->propagate(); }
    BoostTimer::time_point after_propagation = BoostTimer::now();

    // Triangulate the contours and save the outpur mesh.
    FloatMesh result_mesh = Triangulator::christiansen(contour_descriptors);
    BoostTimer::time_point after_triangulation = BoostTimer::now();

    std::cout << "Sheep propagators construction time: "
              << boost::chrono::duration<double>(after_construction - start).count()
              << "s " << std::endl;
    std::cout << "Sheep propagation time: "
              << boost::chrono::duration<double>(after_propagation - after_construction).count()
              << "s " << std::endl;
    std::cout << "Sheep triangulation time: "
              << boost::chrono::duration<double>(after_triangulation - after_propagation).count()
              << "s " << std::endl;
    std::cout << "Sheep total time: "
              << boost::chrono::duration<double>(after_triangulation - start).count()
              << "s " << std::endl;

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
