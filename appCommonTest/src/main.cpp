
#include "gtest/gtest.h"

#if defined(_MSC_VER) && defined(_DEBUG)
#   define _CRTDBG_MAP_ALLOC
#   include <cstdlib>
#   include <crtdbg.h>
#endif // defined(_MSC_VER) && defined(_DEBUG)

// Directory where test data is stored.
std::string DataDirectory;


int main(int argc, char* argv[])
{
    // Enable MSVC's debug heap for detecting memory leaks. This includes changing
    // new operator to one with more info about the leaked block, tracking of all
    // allocations and redirecting output to the stderr.
#if defined(_MSC_VER) && defined(_DEBUG)
#   define new new(_CLIENT_BLOCK, __FILE__, __LINE__)

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





    //// Load mesh.
    //common::Mesh mesh = common::Mesh::from_ply(std::string("..\\data\\mesh.ply"));

    //// Load point cloud.
    //

    //std::cout << mesh << std::endl;

    //common::Mesh::Normal normal = mesh.get_vertex_normal(1);

    ////mesh.to_ply(std::string("..\\debug\\pts2.ply"));

	return 0;
}
