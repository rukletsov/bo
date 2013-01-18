
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

# User-defined options for the project.
OPTION (BuildTests "Build all of Bo's tests" OFF)
OPTION (BuildExamples "Build all of Bo's examples" ON)
OPTION (BUILD_SHARED_LIBS "Build shared libraries (DLLs)" OFF)

IF (MSVC)
    OPTION (UsePch "Enable use of precompiled headers" OFF)
ENDIF (MSVC)


PROJECT (BoProject)

SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET (TESTS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# Gloabl dependencies (for all projects).
INCLUDE_DIRECTORIES ($ENV{BOOST_ROOT})

LINK_DIRECTORIES ($ENV{BOOST_ROOT}/stage/lib)

# Projects list with root paths.
SET (BO_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/Bo)
SET (BOTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/tests)
SET (BOEXAMPLES_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/examples)

ADD_SUBDIRECTORY (${BO_ROOT})

# Add test project only if requested.
IF (BuildTests)
    ENABLE_TESTING()
    ADD_SUBDIRECTORY (${BOTEST_ROOT})
ENDIF (BuildTests)

# Add examples if requested.
IF (BuildExamples)
    ADD_SUBDIRECTORY (${BOEXAMPLES_ROOT})
ENDIF (BuildExamples)
