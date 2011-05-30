
#include "gtest/gtest.h"

#if defined(_MSC_VER) && defined(_DEBUG)
#   define _CRTDBG_MAP_ALLOC
#   include <stdlib.h>
#   include <crtdbg.h>
#endif


int main(int argc, char* argv[])
{
    // Dump detected memory leaks into the stderr for debug mode.
#ifdef _MSC_VER
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);

    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);

    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);

    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif


    // Run all defined tests.
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();



    //// operator[] test.
    //double temp1 = vec1[1];
    //vec1[1] = 10.;
    //// Should call assertion.
    ////vec4[3] = 10.;



    //// Other functions test.
    //vec3.size();
    //vec1.swap(vec2);
    //vec4.assign(f, 2);
    //std::cout << vec1 << std::endl;






    //// Load mesh.
    //common::Mesh mesh = common::Mesh::from_ply(std::string("..\\data\\mesh.ply"));

    //// Load point cloud.
    //

    //std::cout << mesh << std::endl;

    //common::Mesh::Normal normal = mesh.get_vertex_normal(1);

    ////mesh.to_ply(std::string("..\\debug\\pts2.ply"));

	return 0;
}
