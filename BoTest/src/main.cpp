
#include "gtest/gtest.h"

#include "debug_alloc.hpp"

// Directory where test data is stored.
std::string DataDirectory;


int main(int argc, char* argv[])
{
    // Enable MSVC's debug heap for detecting memory leaks. This includes changing
    // new operator to one with more info about the leaked block, tracking of all
    // allocations and redirecting output to the stderr.
#if defined(_MSC_VER) && defined(_DEBUG)
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);

    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);

    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);

    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif // defined(_MSC_VER) && defined(_DEBUG)


    // Extract GTest's command-line arguments and prepare test environment.
    testing::InitGoogleTest(&argc, argv);

    // Extract directory with data for tests from command-line, or apply default value.
    if (argc > 1)
        DataDirectory.assign(argv[1]);
    else if (argc == 1)
        // Apply default value which is "./data" directory.
        DataDirectory.assign("./data");
    else
        DataDirectory.assign("");

    // Run all defined tests.
    return RUN_ALL_TESTS();



    //        |
    // TODO:  |  move this to mesh_unittest.
    //        V

    //// Load mesh.
    //bo::Mesh mesh = bo::Mesh::from_ply(std::string("..\\data\\mesh.ply"));

    //// Load point cloud.
    //

    //std::cout << mesh << std::endl;

    //bo::Mesh::Normal normal = mesh.get_vertex_normal(1);

    ////mesh.to_ply(std::string("..\\debug\\pts2.ply"));
}
