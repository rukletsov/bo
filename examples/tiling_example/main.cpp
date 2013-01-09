
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>

#include "bo/methods/parallel_planes_tiling.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

// Directory where test data is stored.
std::string DataDirectory;

// Paths to data used in the example.
static const std::string raw_test_filename = "2d_contour4_512x512_8bit.raw";
static const std::string ply_femur_filename = "femur/femur_plane01.ply";
static const std::string ply_sheep_filename = "sheep_plane.ply";

using namespace bo::methods::surfaces;
using namespace bo::io;

boost::filesystem3::path raw_test_filepath_;
boost::filesystem3::path ply_femur_filepath_;
boost::filesystem3::path ply_sheep_filepath_;

void AssertFileExists(const boost::filesystem3::path& filepath)
{
    BOOST_ASSERT(boost::filesystem3::exists(boost::filesystem3::path(filepath)));
}

void SetUp()
{
    // Get data directory name. Input directory can have or not have trail slashes.
    raw_test_filepath_ = boost::filesystem3::path(DataDirectory) /= raw_test_filename;
    ply_femur_filepath_ = boost::filesystem3::path(DataDirectory) /= ply_femur_filename;
    ply_sheep_filepath_ = boost::filesystem3::path(DataDirectory) /= ply_sheep_filename;
}

void RunPropagation();
void RunChrisitiansen();


int main(int argc, char* argv[])
{
    // Extract directory with data for tests from command-line, or apply default value.
    if (argc > 1)
        DataDirectory.assign(argv[1]);
    else if (argc == 1)
        // Apply default value which is "./data/tiling" directory.
        DataDirectory.assign((boost::filesystem3::initial_path() /= "data/tiling").string());
    else
        DataDirectory.assign("");

    SetUp();

    RunChrisitiansen();
}

void RunPropagation()
{
    typedef MinSpanPropagation<float> TilingAlgo;

    MinSpanPropagation<float> tiling;

//    ASSERT_TRUE(IsFileAvailable(raw_test_filepath_));
//    TilingAlgo::Image2D test_image = load_raw_image_8bpps<float>(raw_test_filepath_.string(), 512, 512);
//    TilingAlgo::ParallelPlanePtr plane_data = tiling.load_plane(test_image);

//    ASSERT_TRUE(IsFileAvailable(ply_femur_filepath_));
//    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_femur_filepath_.string());
//    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    AssertFileExists(ply_sheep_filepath_);
    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_sheep_filepath_.string());
    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    TilingAlgo::ParallelPlanePtr contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = tiling.to_mesh(contour);
    mesh_to_ply(mesh, (boost::filesystem3::path(DataDirectory) /= "result_contour.ply").string());
}

void RunChrisitiansen()
{
    typedef MinSpanPropagation<float> TilingAlgo;

    MinSpanPropagation<float> tiling;

//    ASSERT_TRUE(IsFileAvailable(raw_test_filepath_));
//    TilingAlgo::Image2D test_image = load_raw_image_8bpps<float>(raw_test_filepath_.string(), 512, 512);
//    TilingAlgo::ParallelPlanePtr plane_data = tiling.load_plane(test_image);

    AssertFileExists(ply_femur_filepath_);
    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_femur_filepath_.string());
    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

//    ASSERT_TRUE(IsFileAvailable(ply_sheep_filepath_));
//    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_sheep_filepath_.string());
//    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    TilingAlgo::ParallelPlanePtr contour1 = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::ParallelPlanePtr contour2 = tiling.propagate(plane_data, 0.2f);

    for (TilingAlgo::ParallelPlane::iterator it = contour2->begin();
         it != contour2->end(); ++it)
        it->z() = it->z() + 3;

    TilingAlgo::Mesh mesh = tiling.christiansen_triangulation(contour1, contour2);
    mesh_to_ply(mesh, (boost::filesystem3::path(DataDirectory) /= "result_mesh.ply").string());
}
