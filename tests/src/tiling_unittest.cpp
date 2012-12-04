
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "bo/methods/parallel_planes_tiling.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

extern std::string DataDirectory;
static const std::string test_filename = "2d_contour4_512x512_8bit.raw";

using namespace bo::methods::surfaces;
using namespace bo::io;


// Create a so-called "text fixture" using base class form GTEST.
class TilingTest : public testing::Test
{
protected:
    virtual void SetUp()
    {
        // Get data directory name. Input directory can have or not have trail slashes.
        test_filepath = boost::filesystem3::path(DataDirectory) /= test_filename;
    }

    ::testing::AssertionResult IsTestFileAvailable(boost::filesystem3::path filepath)
    {
        return
            (boost::filesystem3::exists(boost::filesystem3::path(filepath))
             ? ::testing::AssertionSuccess() << "file \"" << filepath.string() << "\" found"
             : ::testing::AssertionFailure() << "file \"" << filepath.string() << "\" not found");
    }

    boost::filesystem3::path test_filepath;
};

TEST_F(TilingTest, Temp)
{
    typedef MinSpanPropagation<float> TilingAlgo;

    MinSpanPropagation<float> tiling;

    TilingAlgo::Image2D test_image = load_raw_image_8bpps<float>(test_filepath.string(), 512, 512);
//    std::cout << test_image(112, 132) << " ";
//    std::cout << test_image(113, 132) << " ";
//    std::cout << test_image(114, 132) << " ";
//    std::cout << test_image(115, 132) << " ";
//    std::cout << test_image(116, 132) << " ";

    TilingAlgo::ParallelPlanePtr plane_data = tiling.load_plane(test_image);
    TilingAlgo::Mesh mesh = tiling.propagate(plane_data);
    mesh_to_ply(mesh, (boost::filesystem3::path(DataDirectory) /= "result.ply").string());
}
