
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (AppCommonTest)

SET (AppCommonTest_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vector_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/timer_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/raw_image_2d_unittest.cpp
)

SET (AppCommonTest_HDRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/debug_alloc.hpp
)

# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (./include
                     $ENV{GTEST_ROOT}/include
                     ${LIBCOMMON_ROOT}
)
     
# For additional dependencies check the root CMakeLists.
LINK_DIRECTORIES ($ENV{GTEST_ROOT}/stage/lib
                  ${LIBRARY_OUTPUT_PATH}
)

# Compile our Library.
ADD_EXECUTABLE (AppCommonTest
                ${AppCommonTest_SRCS}
                ${AppCommonTest_HDRS}
)

TARGET_LINK_LIBRARIES (AppCommonTest
                       LibCommon
                       debug gtestd optimized gtest
)

# On Linux GTest depends on pthread library (so it should be included).
IF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    TARGET_LINK_LIBRARIES (AppCommonTest
                           pthread
    )
ENDIF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

# Specify the name suffix for debug build.
SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES DEBUG_POSTFIX "d")

# Rise warning level to maximum for MSVC and GNUCXX.
IF (MSVC)
    IF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        # Replace exisiting warning level with a desired one.
        STRING (REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ELSE ()
        # Set a desired warning level if there was no set previously.
        ADD_DEFINITIONS ("/W4")
    ENDIF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")

    # Disable C4512 warning. There are a some boost classes, for which assignment
    # operator cannot be generated. The warning seems not to be important.
    ADD_DEFINITIONS ("/wd4512")

ELSEIF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES COMPILE_FLAGS
                           "-Wall -Wextra -pedantic")
ENDIF (MSVC)

# Copy binary to the directory with test data.
INSTALL (TARGETS AppCommonTest DESTINATION ${TESTS_OUTPUT_PATH})

# Create tests from the current binary and pass a directory where test data is stored
# as a command-line argument.
FUNCTION (TestCase test_name)
    ADD_TEST (NAME ${test_name}Test
              COMMAND AppCommonTest --gtest_filter=${test_name}*
                                    "${TESTS_OUTPUT_PATH}/data"
    )
ENDFUNCTION (TestCase test_name)

TestCase(Vector)
TestCase(Mesh)
TestCase(Timer)
TestCase(RawImage2D)


# Supress MSVC-specific compiler and linker warnings.
IF (MSVC)
    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES COMPILE_FLAGS 
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)
