
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

PROJECT (Bo)

# Specify project's C++ source files.
SET (BO_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/mesh_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/ini_reader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/methods/d25_active_contours.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/methods/mbutterfly_mesh_subdiv.cpp
)

# A special section for C sources, which should be compiled using C compiler, not C++
SET (BO_CSRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.1.2/rply/rply.c
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.c
)

SET (BO_HDRS
    ${CMAKE_CURRENT_SOURCE_DIR}/pch.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/config.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/bo/vector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mesh.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/triangle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/logging.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/performance.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/extended_math.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/extended_std.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/raw_image_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/transformation_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/mrf.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/image_operations.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/kdtree.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/topology.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/blas/blas.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/blas/conversions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/blas/pca.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/d25_active_contours.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/distances_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/icp_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/linear_filtering_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/mbutterfly_mesh_subdiv.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/mesh_closest_face_sr.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/dpg_hough_transform_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/complex_propagation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/triangulation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/methods/convex_hull_3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/mrf_2d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/likelihood_functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/node_types.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/optimization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/prior_functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/mrf/type_values.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/surfaces/container_traversers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/internal/surfaces/arched_strip.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/mesh_io.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/ini_reader.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/raw_image_2d_io.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/config_parser/basic_configuration.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bo/io/config_parser/uon_basic_configuration.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.1.2/rply/rply.h

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.h

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/kdtree.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/allocator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/function.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/iterator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/node.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/region.hpp
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


# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (./
                     ./3rdparty
                     ./3rdparty/rply-1.1.2
)

# Let the sources know, that this is our own source.
ADD_DEFINITIONS (-DBO_SOURCE)

# Add shared library support if requested by a user.
IF (BUILD_SHARED_LIBS)
    ADD_DEFINITIONS (-DBO_SHARED_LIBRARY)
ENDIF (BUILD_SHARED_LIBS)

# Compile the library.
ADD_LIBRARY (bo
             ${BO_HDRS}
             ${BO_SRCS}
             ${BO_CSRCS}
)

# Add boost dependencies if compiled not with msvc (boost has auto-link support for it).
IF (NOT MSVC)
    TARGET_LINK_LIBRARIES (bo
                           boost_regex)
ENDIF (NOT MSVC)


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
        ADD_DEFINITIONS (-D_WIN32_WINNT=0x0601)
    ENDIF ((MSVC10) OR ((MSVC_VERSION) EQUAL 1600) OR ((MSVC_VERSION) GREATER 1600))

    # Disable C4512 warning. There are a some boost classes, for which assignment
    # operator cannot be generated. The warning seems not to be important.
    ADD_DEFINITIONS ("/wd4512")

    # Disable C4127 warning. Constant conditional expressions should not always
    # trigger a warning (e.g. for while(1)). But MSVC does, which is annoying.
    ADD_DEFINITIONS ("/wd4127")

    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (bo PROPERTIES COMPILE_FLAGS
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")

    # Supress "no public symbols found" MSVC's linker warning (this can occur for
    # the "pch.obj" file when using precompiled headers.
    SET_TARGET_PROPERTIES (bo PROPERTIES STATIC_LIBRARY_FLAGS "/IGNORE:4221")

ENDIF (MSVC)
