
#define BOOST_FILESYSTEM_VERSION 3

#include <cstdio>
#include <string>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/assert.hpp>

#include "bo/methods/parallel_planes_tiling.hpp"
#include "bo/io/raw_image_2d_io.hpp"
#include "bo/io/mesh_io.hpp"

#include "bo/extended_std.hpp"
#include <boost/format.hpp>

// Directory where test data is stored.
std::string DataDirectory;

// Paths to data used in the example.
static const std::string raw_test_filename = "2d_contour4_512x512_8bit.raw";
static const std::string ply_femur_dirname = "femur";
static const std::string ply_femur_filename1 = "femur/femur_plane01.ply";
static const std::string ply_femur_filename2 = "femur/femur_plane02.ply";
static const std::string ply_sheep_filename = "sheep_plane.ply";

using namespace boost::filesystem;
using namespace bo::methods::surfaces;
using namespace bo::io;

path raw_test_filepath_;
path ply_femur_dirpath_;
path ply_femur_filepath1_;
path ply_femur_filepath2_;
path ply_sheep_filepath_;

void AssertPathExists(const path& filepath)
{
    BOOST_ASSERT(exists(filepath));
}

void SetUp()
{
    // Get data directory name. Input directory can have or not have trail slashes.
    raw_test_filepath_ = path(DataDirectory) /= raw_test_filename;
    ply_femur_dirpath_ = path(DataDirectory) /= ply_femur_dirname;
    ply_femur_filepath1_ = path(DataDirectory) /= ply_femur_filename1;
    ply_femur_filepath2_ = path(DataDirectory) /= ply_femur_filename2;
    ply_sheep_filepath_ = path(DataDirectory) /= ply_sheep_filename;
}

void RunPropagation();
void RunChrisitiansen(const path &contour1_path, const path &contour2_path);
void RunFemurChrisitiansen();


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

    SetUp();

    RunChrisitiansen(ply_femur_filepath1_, ply_femur_filepath2_);
    RunFemurChrisitiansen();
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

    AssertPathExists(ply_sheep_filepath_);
    TilingAlgo::Mesh test_mesh = mesh_from_ply(ply_sheep_filepath_.string());
    TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(test_mesh.get_all_vertices()));

    TilingAlgo::ParallelPlanePtr contour = tiling.propagate(plane_data, 0.5f);
    TilingAlgo::Mesh mesh = tiling.to_mesh(contour);
    mesh_to_ply(mesh, (path(DataDirectory) /= "result_contour.ply").string());
}

void RunChrisitiansen(const path& contour1_path, const path& contour2_path)
{
    typedef MinSpanPropagation<float> TilingAlgo;

    MinSpanPropagation<float> tiling;

    AssertPathExists(contour1_path);
    TilingAlgo::Mesh test_mesh1 = mesh_from_ply(contour1_path.string());
    TilingAlgo::ParallelPlanePtr plane_data1(new TilingAlgo::ParallelPlane(test_mesh1.get_all_vertices()));

    AssertPathExists(contour2_path);
    TilingAlgo::Mesh test_mesh2 = mesh_from_ply(contour2_path.string());
    TilingAlgo::ParallelPlanePtr plane_data2(new TilingAlgo::ParallelPlane(test_mesh2.get_all_vertices()));

    TilingAlgo::ParallelPlanePtr contour1 = tiling.propagate(plane_data1, 0.5f);
    TilingAlgo::ParallelPlanePtr contour2 = tiling.propagate(plane_data2, 0.2f);

    TilingAlgo::Mesh mesh = tiling.christiansen_triangulation(contour1, contour2);
    mesh_to_ply(mesh, (path(DataDirectory) /= "result_mesh.ply").string());
}

void RunFemurChrisitiansen()
{
    typedef MinSpanPropagation<float> TilingAlgo;
    typedef std::vector<path> ContourData;
    typedef std::vector<TilingAlgo::ParallelPlanePtr> Contours;
    typedef std::vector<TilingAlgo::Mesh> Meshes;

    MinSpanPropagation<float> tiling;
    ContourData contour_data;
    Contours contours;
    Meshes meshes;

    // Load planes paths.
    AssertPathExists(ply_femur_dirpath_);
    for (directory_iterator it(ply_femur_dirpath_); it != directory_iterator(); ++it)
    {
        if ((is_regular_file(*it)) &&
            (it->path().filename().string().substr(0, 11) == "femur_plane"))
            contour_data.push_back(it->path());
    }

    // Sort paths in lexicographic order, since directory iteration is not ordered
    // on some file systems.
    std::sort(contour_data.begin(), contour_data.end());
//    std::cout << "List of loaded contours: " << contour_data;

    // Extract contours from contour data.
    for (ContourData::const_iterator it = contour_data.begin(); it != contour_data.end(); ++it)
    {
        AssertPathExists(*it);
        TilingAlgo::Mesh mesh = mesh_from_ply(it->string());
        TilingAlgo::ParallelPlanePtr plane_data(new TilingAlgo::ParallelPlane(mesh.get_all_vertices()));
        TilingAlgo::ParallelPlanePtr contour = tiling.propagate(plane_data, 0.5f);

        contours.push_back(contour);
    }

    // Tile pair of contours.
    int i = 1;
    create_directory(path(DataDirectory) /= "result");
    for (Contours::const_iterator it = contours.begin() + 1; it != contours.end(); ++it)
    {
        TilingAlgo::Mesh mesh = tiling.christiansen_triangulation(*(it - 1), *it);
        std::string outpath = ((path(DataDirectory) /= "result") /=
                (boost::format("%1%.ply") % i).str()).string();
        std::cout << outpath << std::endl;
        mesh_to_ply(mesh, outpath);

        //meshes.push_back(tiling.christiansen_triangulation(*(it - 1), *it));
        ++i;
    }

    // Stitch stripes.
}
