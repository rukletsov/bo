
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# User-defined options for the library.
IF (MSVC)
    OPTION (UsePch "Enable use of precompiled headers" OFF)
ENDIF (MSVC)


PROJECT (LibCommon)

# Specify project's C++ source files.
SET (LibCommon_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/d25_active_contours.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ini_reader.cpp
)

# A special section for C sources, which should be compiled using C compiler, not C++
SET (LibCommon_CSRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.01/rply/rply.c
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.c
)

SET (LibCommon_HDRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pch.h

    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/vector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/mesh.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/triangle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/errors.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/performance.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common/raw_image.hpp	
    ${CMAKE_CURRENT_SOURCE_DIR}/include/blas/blas.hpp	
    ${CMAKE_CURRENT_SOURCE_DIR}/include/methods/d25_active_contours.hpp	
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io/ini_reader.hpp	
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io/config_parser/basic_configuration.hpp	
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io/config_parser/uon_basic_configuration.hpp	
	
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.01/rply/rply.h

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svd/svd.h

    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/kdtree.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/allocator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/function.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/iterator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/node.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libkdtree++_0.7.0/kdtree++/region.hpp
)


# Enable usage for precompiled header for MSVC compiler for C++ sources.
IF ((MSVC) AND (UsePch))
    # Set preprocessor directive USE_PCH so we can use it code sources.
    ADD_DEFINITIONS (-DUSE_PCH)

    # Set "use precompiled header" option for every ".cpp" file.
    FOREACH (src_file ${LibCommon_SRCS})
        SET_SOURCE_FILES_PROPERTIES (${src_file} PROPERTIES COMPILE_FLAGS "/Yupch.h")
    ENDFOREACH (src_file ${LibCommon_SRCS})

    # Set "create precompiled header" option only for "pch.cpp" file.
    SET_SOURCE_FILES_PROPERTIES (${CMAKE_CURRENT_SOURCE_DIR}/src/pch.cpp
                                 PROPERTIES COMPILE_FLAGS "/Ycpch.h")

    # Add "pch.cpp" file to the list of C++ sources.
    SET (LibCommon_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/pch.cpp
        ${LibCommon_SRCS})

ENDIF ((MSVC) AND (UsePch))


# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (./include
                     ./include/blas
                     ./include/common
                     ./include/methods
                     ./include/io
                     ./include/io/config_parser
                     ./3rdparty
                     ./3rdparty/rply-1.01
                     ./3rdparty/libkdtree++_0.7.0
)

# Compile our Library.
ADD_LIBRARY (LibCommon
             ${LibCommon_HDRS}
             ${LibCommon_SRCS}
             ${LibCommon_CSRCS}
)

# Specify the name suffix for debug build.
SET_TARGET_PROPERTIES (LibCommon PROPERTIES DEBUG_POSTFIX "d")


# Supress MSVC-specific compiler and linker warnings.
IF (MSVC)
    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (LibCommon PROPERTIES COMPILE_FLAGS 
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")                         
  
    # Supress "no public symbols found" MSVC's linker warning (this can occur for
    # the "pch.obj" file when using precompiled headers.
    SET_TARGET_PROPERTIES (LibCommon PROPERTIES STATIC_LIBRARY_FLAGS "/IGNORE:4221")
ENDIF (MSVC)
