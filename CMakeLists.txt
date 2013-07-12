
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

# User-defined options for the project.
OPTION (BuildTests "Build all of Bo's tests" OFF)
OPTION (BuildExamples "Build all of Bo's examples" ON)
OPTION (BUILD_SHARED_LIBS "Build shared libraries (DLLs)" OFF)

IF (MSVC)
    OPTION (UsePch "Enable use of precompiled headers" OFF)
ENDIF (MSVC)


PROJECT (Bo)

# Cache some paths.
SET (BO_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET (TESTS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests)
SET (EXAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Check that one of the officially supported compilers is used.
IF (MSVC)
    MESSAGE (STATUS "Supported compiler used: MSVC")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE (STATUS "Supported compiler used: GCC")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    MESSAGE (STATUS "Supported compiler used: Clang")
ELSE ()
    MESSAGE ("Unsupported compiler used: ${CMAKE_CXX_COMPILER_ID}")
ENDIF (MSVC)

# Locate and enable Boost.
FIND_PACKAGE (Boost COMPONENTS regex system filesystem chrono REQUIRED)
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})

# Increase default warning level for sources in the current directory and below.
# Supports MSVC, GCC and Clang compilers.
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # Remove default warning level from .cpp and .c sources to avoid conflicts with
    # desired warning, which is manually set later.
    IF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        # Replace exisiting warning level with the empty one.
        STRING (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ENDIF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")

    IF (CMAKE_C_FLAGS MATCHES "/W[0-4]")
        # Replace exisiting warning level with the empty one.
        STRING (REGEX REPLACE "/W[0-4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    ENDIF (CMAKE_C_FLAGS MATCHES "/W[0-4]")

    # Set the desired warning level. If it has been already set and differs from the
    # provided one, a D9025 warning is raised.
    ADD_DEFINITIONS ("/W4")

ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    ADD_DEFINITIONS ("-Wall -Wextra -Wno-long-long -Wno-variadic-macros -pedantic")

ENDIF ()

# Specify project's C++ source files.
SET (BO_SRCS
    # pch.cpp is added later on demand.
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/mesh_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/ini_reader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/surfaces/d25_active_contours.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/surfaces/mbutterfly_mesh_subdiv.cpp
)

# A special section for C sources, which should be compiled using C compiler, not C++
SET (BO_CSRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.1.2/rply/rply.c
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.c
)

SET (BO_HDRS
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/pch.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/config.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/vector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/mesh.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/triangle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/raw_image_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/image_operations.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/linear_filtering_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/core/kdtree.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/distances/distances_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/distances/mesh_closest_face_sr.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/ini_reader.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/mesh_io.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/raw_image_2d_io.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/config_parser/basic_configuration.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/config_parser/uon_basic_configuration.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/blas_extensions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/blas_conversions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/mean.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/pca.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/topology.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/math/transformation_3d.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/mrf.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/likelihood_functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/mrf_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/node_types.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/optimization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/prior_functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf/type_values.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/recognition/dpg_hough_transform_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/recognition/hough_transform.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/registration/icp_3d.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/christiansen_tiling.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/complex_propagator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/convex_hull_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/d25_active_contours.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/mbutterfly_mesh_subdiv.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/triangulation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/arched_strip.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/container_traversers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/indexed_tstrip.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/propagation_result.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/inertial_propagation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/centrifugal_propagation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/tangential_propagation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/surfaces/detail/total_propagation.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/utils/logging.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/utils/performance.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/utils/std_extensions.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/kdtree.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/allocator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/function.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/iterator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/node.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/region.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.1.2/rply/rply.h
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.h
)

# Enable precompiled headers if requested by the user.
IF (UsePch)

    # Set preprocessor directive BO_USE_PCH so we can use it in code sources.
    ADD_DEFINITIONS (-DBO_USE_PCH)

    # MSVC-specific handling of precompiled headers.
    IF (MSVC)
        # Print info message.
        MESSAGE (STATUS "Precompiled header enabled through \"pch.h\" and \"pch.cpp\"")

        # Set "use precompiled header" option for every ".cpp" file.
        FOREACH (src_file ${BO_SRCS})
            SET_SOURCE_FILES_PROPERTIES (${src_file} PROPERTIES COMPILE_FLAGS "/Yupch.h")
        ENDFOREACH (src_file ${BO_SRCS})

        # Set "create precompiled header" option only for "pch.cpp" file.
        SET_SOURCE_FILES_PROPERTIES (${CMAKE_CURRENT_SOURCE_DIR}/src/pch.cpp
                                     PROPERTIES COMPILE_FLAGS "/Ycpch.h")

        # Add "pch.cpp" file to the list of C++ sources.
        SET (BO_SRCS
            ${CMAKE_CURRENT_SOURCE_DIR}/src/pch.cpp
            ${BO_SRCS})

    ENDIF (MSVC)

    # Enable usage for precompiled header for non-MSVC compilers (NOT TESTED, can
    # miss some actions. Should be expanded separately for every supported compiler.

ENDIF (UsePch)

# Compile the library.
ADD_LIBRARY (bo
             ${BO_HDRS}
             ${BO_SRCS}
             ${BO_CSRCS}
)

INCLUDE_DIRECTORIES (${BO_ROOT})

# Let the sources know, that this is our own source.
SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_DEFINITIONS BO_SOURCE)

# Add shared library support if requested by a user.
IF (BUILD_SHARED_LIBS)
    SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_DEFINITIONS BO_SHARED_LIBRARY)
ENDIF (BUILD_SHARED_LIBS)

# Add boost dependencies.
TARGET_LINK_LIBRARIES (bo
                       ${Boost_LIBRARIES}
)

# Specify the name suffix for debug build.
SET_TARGET_PROPERTIES (bo PROPERTIES DEBUG_POSTFIX "d")

# Compiler-specific features and options.
IF (MSVC)
    # If MSVC2010 compiler (or newer) is used, target Windows 7. This means that
    # the newest Windows API will be used but it can be partially unavailable on
    # previous versions of Windows. For older compilers a default value will be
    # applied in the "pch.h" file. For more information see
    #     http://msdn.microsoft.com/en-us/library/aa383745%28v=vs.85%29.aspx
    IF ((MSVC10) OR ((MSVC_VERSION) EQUAL 1600) OR ((MSVC_VERSION) GREATER 1600))
        MESSAGE (STATUS "Targeting Windows 7 API: _WIN32_WINNT = 0x0601")
        SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_DEFINITIONS _WIN32_WINNT=0x0601)
    ENDIF ((MSVC10) OR ((MSVC_VERSION) EQUAL 1600) OR ((MSVC_VERSION) GREATER 1600))

    # Disable C4512 warning. There are a some boost classes, for which assignment
    # operator cannot be generated. The warning seems not to be important.
    SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_FLAGS "/wd4512")

    # Disable C4127 warning. Constant conditional expressions should not always
    # trigger a warning (e.g. for while(1)). But MSVC does, which is annoying.
    SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_FLAGS "/wd4127")

    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_FLAGS
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")

    # Supress "no public symbols found" MSVC's linker warning (this can occur for
    # the "pch.obj" file when using precompiled headers.
    SET_TARGET_PROPERTIES (bo PROPERTIES STATIC_LIBRARY_FLAGS "/IGNORE:4221")

ENDIF (MSVC)


# Add test project only if requested.
IF (BuildTests)
    ENABLE_TESTING()
    FIND_PACKAGE (GTest)

    # Do not enable tests if GTest is not found.
    IF (GTEST_FOUND)
        ADD_SUBDIRECTORY (${TESTS_PATH})
    ELSE ()
        MESSAGE ("GTest not found, tests disabled.")
    ENDIF (GTEST_FOUND)

ENDIF (BuildTests)

# Add examples if requested.
IF (BuildExamples)
    ADD_SUBDIRECTORY (${EXAMPLES_PATH})
ENDIF (BuildExamples)
