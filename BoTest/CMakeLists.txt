
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

PROJECT (BoTest)

SET (BoTest_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vector_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/timer_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/raw_image_2d_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ini_reader_unittest.cpp
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
# debug and release versions of gtest should be available when using MSVC compiler.
IF (MSVC)
    TARGET_LINK_LIBRARIES (BoTest
                           bo
                           debug gtestd optimized gtest
    )
ELSE ()
    TARGET_LINK_LIBRARIES (BoTest
                           bo
                           gtest
    )
ENDIF (MSVC)

# Add boost dependencies if compiled not with msvc (boost has auto-link support for it).
IF (NOT MSVC)
    TARGET_LINK_LIBRARIES (BoTest
                           boost_filesystem
                           boost_system)
ENDIF (NOT MSVC)

ADD_CUSTOM_COMMAND(TARGET BoTest
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${TESTS_OUTPUT_PATH}/data"
        "$<TARGET_FILE_DIR:BoTest>/data/"
        COMMENT "Copy test data to the runtime directory."
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

# Suppress some warnings and set warning level for MSVC and GNUCXX.
IF (MSVC)
    # Remove default warning level from .cpp and .c sources to avoid conflicts with
    # desired warning which is manually set later.
    IF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        # Replace exisiting warning level with a desired one.
        STRING (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ENDIF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")

    IF (CMAKE_C_FLAGS MATCHES "/W[0-4]")
        # Replace exisiting warning level with a desired one.
        STRING (REGEX REPLACE "/W[0-4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    ENDIF (CMAKE_C_FLAGS MATCHES "/W[0-4]")

    # Set a desired warning level. If it has been already set and differs from the
    # provided one, a D9025 warning is raised.
    ADD_DEFINITIONS ("/W4")

    # Disable C4245 warning. Conversions from negatibe integers to size_t are used
    # widely within the tests. It's done on purpose in order to test incorrect
    # indexation.
    ADD_DEFINITIONS ("/wd4245")

ELSEIF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    SET_TARGET_PROPERTIES (BoTest PROPERTIES COMPILE_FLAGS
                           "-Wall -Wextra -Wno-long-long -Wno-variadic-macros -pedantic")
ENDIF (MSVC)

# Create tests from the current binary and pass a directory where test data is stored
# as a command-line argument.
FUNCTION (TestCase test_name)
    ADD_TEST (NAME ${test_name}Test
              COMMAND BoTest --gtest_filter=${test_name}*
                                    "${TESTS_OUTPUT_PATH}/data"
    )
ENDFUNCTION (TestCase test_name)

TestCase(Vector)
TestCase(Mesh)
TestCase(Timer)
TestCase(RawImage2D)
TestCase(IniReader)


# Supress MSVC-specific compiler and linker warnings.
IF (MSVC)
    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (BoTest PROPERTIES COMPILE_FLAGS 
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)
