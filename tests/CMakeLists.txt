
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

PROJECT (BoTest)

SET (BoTest_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vector_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/timer_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/raw_image_2d_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ini_reader_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dpg_hough_transform_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/d25_active_contours_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/topology_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/convex_hull_3d_unittest.cpp
)

SET (BoTest_HDRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/debug_alloc.hpp
)

# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (./include
                     $ENV{GTEST_ROOT}/include
                     ${BO_ROOT}
)
     
# For additional dependencies check the root CMakeLists.
LINK_DIRECTORIES ($ENV{GTEST_ROOT}/stage/lib
                  ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
                  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)

# Compile the application.
ADD_EXECUTABLE (BoTest
                ${BoTest_SRCS}
                ${BoTest_HDRS}
)

# Add dependencies. MSVC distinguishes between debug and release runtime, therefore both
# debug and release versions of gtest should be available when using MSVC compiler. Also
# add boost dependencies if compiler is not MSVC (boost has auto-link support for it).
IF (MSVC)
    TARGET_LINK_LIBRARIES (BoTest
                           bo
                           debug gtestd optimized gtest
    )
ELSE ()
    TARGET_LINK_LIBRARIES (BoTest
                           bo
                           gtest
                           boost_system
                           boost_filesystem
    )
ENDIF (MSVC)

ADD_CUSTOM_COMMAND(TARGET BoTest
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${TESTS_PATH}/data"
        "$<TARGET_FILE_DIR:BoTest>/data/"
        COMMENT "Copying test data to the runtime directory."
)

# GTest depends on pthread library on platforms, where it is availabe (e.g. GNUCXX).
FIND_PACKAGE (Threads)
IF (${CMAKE_USE_PTHREADS_INIT})
    TARGET_LINK_LIBRARIES (BoTest
                           ${CMAKE_THREAD_LIBS_INIT}
    )
ENDIF (${CMAKE_USE_PTHREADS_INIT})

# Specify the name suffix for debug build.
SET_TARGET_PROPERTIES (BoTest PROPERTIES DEBUG_POSTFIX "d")

# Compiler-specific features and options.
IF (MSVC)
    # Disable C4245 warning. Conversions from negatibe integers to size_t are used
    # widely within the tests. It's done on purpose in order to test incorrect
    # indexation.
    ADD_DEFINITIONS ("/wd4245")

    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (BoTest PROPERTIES COMPILE_FLAGS
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)

# Create tests from the current binary and pass a directory where test data is stored
# as a command-line argument.
FUNCTION (TestCase test_name)
    ADD_TEST (NAME ${test_name}Test
              COMMAND BoTest --gtest_filter=${test_name}*
                                    "${TESTS_PATH}/data"
    )
ENDFUNCTION (TestCase test_name)

TestCase(Vector)
TestCase(Mesh)
TestCase(Timer)
TestCase(RawImage2D)
TestCase(IniReader)
TestCase(Math)
TestCase(HoughTransform)
TestCase(Geometry)
TestCase(Topology)
TestCase(ConvexHull3D)
