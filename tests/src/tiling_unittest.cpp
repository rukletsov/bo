
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "bo/methods/parallel_planes_tiling.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

extern std::string DataDirectory;
static const std::string raw_test_filename = "2d_contour4_512x512_8bit.raw";
static const std::string ply_femur_filename = "femur_plane.ply";
static const std::string ply_sheep_filename = "sheep_plane.ply";

using namespace bo::methods::surfaces;
using namespace bo::io;


// Create a so-called "text fixture" using base class form GTEST.
class TilingTest : public testing::Test
{
protected:
    virtual void SetUp()
    {
        // Get data directory name. Input directory can have or not have trail slashes.
        raw_test_filepath_ = boost::filesystem3::path(DataDirectory) /= raw_test_filename;
        ply_femur_filepath_ = boost::filesystem3::path(DataDirectory) /= ply_femur_filename;
        ply_sheep_filepath_ = boost::filesystem3::path(DataDirectory) /= ply_sheep_filename;
    }

    ::testing::AssertionResult IsFileAvailable(boost::filesystem3::path filepath)
    {
        return
            (boost::filesystem3::exists(boost::filesystem3::path(filepath))
             ? ::testing::AssertionSuccess() << "file \"" << filepath.string() << "\" found"
             : ::testing::AssertionFailure() << "file \"" << filepath.string() << "\" not found");
    }

    boost::filesystem3::path raw_test_filepath_;
    boost::filesystem3::path ply_femur_filepath_;
    boost::filesystem3::path ply_sheep_filepath_;
};

TEST_F(TilingTest, TempPropagate)
{
    typedef MinSpanPropagation<float> TilingAlgo;

    MinSpanPropagation<float> tiling;

//    ASSERT_TRUE(IsFileAvailable(raw_test_filepath_));
//    TilingAlgo::Image2D test_image = load_raw_image_8bpps<float>(raw_test_filepath_.string(), 512, 512);
//    TilingAlgo::ParallelPlanePtr plane_data = tiling.load_plane(test_image);

//    ASSERT_TRUE(IsFileAvailable(ply_femur_filepath_));
//    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_femur_filepath_.string());
//    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    ASSERT_TRUE(IsFileAvailable(ply_sheep_filepath_));
    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_sheep_filepath_.string());
    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    TilingAlgo::ParallelPlanePtr contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = tiling.to_mesh(contour);
    mesh_to_ply(mesh, (boost::filesystem3::path(DataDirectory) /= "result.ply").string());
}
