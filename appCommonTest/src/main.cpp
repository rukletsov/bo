
#include <iostream>
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

    // Run all declared tests.
    std::cout << "Running all tests." << std::endl;

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

    //common::Vector<int, 3> vec01(1, 2, 3);
    //double d = vec01.eucl_norm();
    //std::cout << vec01[3];


    //// Constructor tests.
    //common::Vector<double, 3> vec1;
    //common::Vector<double, 3> vec2(5.f);

    //float f[2] = {1.f, 2.2f};
    //common::Vector<int, 3> vec3(f, 2);

    //int f2[5] = {3, 4, 10, 45, 3};
    //common::Vector<double, 3> vec4(f2, 5);

    //common::Vector<double, 3> vec5(vec2);


    //// operator[] test.
    //double temp1 = vec1[1];
    //vec1[1] = 10.;
    //// Should call assertion.
    ////vec4[3] = 10.;


    //// Operators + boost generated operators.
    //common::Vector<double, 3> vec6 = vec2 + vec4;
    //vec6 *= 3;
    //vec6 - vec4;


    //// Partially available accessors.
    //common::Vector<double, 4> vec20(0., 4., 56., -1.);
    //vec1.x();
    //vec20.w() = 34.;


    //// Simple functions.
    //common::Vector<double, 3> vec21 = vec1.cross_product(vec2);

    //std::size_t temp2 = vec1.min_index();
    //double temp3 = vec4.min();

    //std::size_t temp4 = vec3.max_index();
    //int temp5 = vec3.max();

    //bool temp6 = vec3[temp4] == temp5;

    //double temp7 = vec2.sum();
    //double temp8 = vec4.product();
    //double temp9 = vec1.avg();


    //// Norm and normalization test.
    //double temp10 = vec1.eucl_norm();
    //float temp11;
    //vec1.eucl_norm(temp11);
    //
    //double temp12 = vec4.eucl_norm();
    //int temp13;
    //vec4.eucl_norm(temp13);

    //common::Vector<double, 3> vec7 = vec4.normalized();


    //// Other functions test.
    //vec3.size();
    //vec1.swap(vec2);
    //vec4.assign(f, 2);
    //std::cout << vec1 << std::endl;


    //// Speed test.
    //// Vector
    //boost::int64_t freq = common::get_proc_freq();
    //boost::int64_t start = common::get_proc_ticks();

    //common::Vector<double, 3000> vec10(10.);
    //common::Vector<double, 1000> res;
    //for (int i = 0; i < 1000; ++i)
    //    res[i] = vec10.sum();

    //boost::int64_t total = common::get_proc_ticks() - start;

    //std::cout << "result: " << vec10.normalized() << std::endl 
    //    << "Vector<> ticks: " << total << std::endl;




    //// Load mesh.
    //common::Mesh mesh = common::Mesh::from_ply(std::string("..\\data\\mesh.ply"));

    //// Load point cloud.
    //

    //std::cout << mesh << std::endl;

    //common::Mesh::Normal normal = mesh.get_vertex_normal(1);

    ////mesh.to_ply(std::string("..\\debug\\pts2.ply"));

	return 0;
}
