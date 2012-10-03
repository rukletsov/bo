
#include <iostream>
#include "bo/performance.hpp"
#include "bo/vector.hpp"

#define VEC_SIZE 10

int main(int argc, char* argv[])
{
    typedef bo::Vector<int, VEC_SIZE> TestVec;
    const int repetitions_count = 100000000;

    TestVec* vec1 = new TestVec(-569);
    TestVec* vec2 = new TestVec(9899);
    TestVec* vec3 = new TestVec(98999);

    // Vector += Scalar
    bo::Timer timer;
    for (int i = 0; i < repetitions_count; ++i)
        (*vec1) += 879556;
    std::cout << "Vector += Scalar took " << timer.elapsed() << std::endl;

    // Vector + Scalar
    timer.restart();
    for (int i = 0; i < repetitions_count; ++i)
        (*vec2) = *vec1 + 879556;
    std::cout << "Vector + Scalar took " << timer.elapsed() << std::endl;

    // Vector /= Scalar
    timer.restart();
    for (int i = 0; i < repetitions_count; ++i)
        (*vec1) /= 29;
    std::cout << "Vector /= Scalar took " << timer.elapsed() << std::endl;

    // Vector / Scalar
    timer.restart();
    for (int i = 0; i < repetitions_count; ++i)
        (*vec2) = *vec1 / 29;
    std::cout << "Vector / Scalar took " << timer.elapsed() << std::endl;

    // Vector -= Vector
    timer.restart();
    for (int i = 0; i < repetitions_count; ++i)
        *vec3 -= *vec1;
    std::cout << "Vector -= Vector took " << timer.elapsed() << std::endl;

    // Vector - Vector
    timer.restart();
    for (int i = 0; i < repetitions_count; ++i)
        *vec3 = *vec1 - *vec2;
    std::cout << "Vector - Vector took " << timer.elapsed() << std::endl;

    delete vec1;
    delete vec2;
    delete vec3;

    return 0;
}
