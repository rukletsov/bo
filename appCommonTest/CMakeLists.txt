
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (AppCommonTest)

SET (AppCommonTest_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vector_unittest.cpp
)

# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES ($ENV{GTEST_ROOT}/include
                     ${LIBCOMMON_ROOT}/include
)
     
# For additional dependencies check the root CMakeLists.
LINK_DIRECTORIES ($ENV{GTEST_ROOT}/stage/lib
                  ${LIBRARY_OUTPUT_PATH}
)

# Compile our Library.
ADD_EXECUTABLE (AppCommonTest
                ${AppCommonTest_SRCS}
)

TARGET_LINK_LIBRARIES (AppCommonTest
                       LibCommon
                       gtest
)

IF (MSVC)
    # Specify the name suffix for debug build.
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES DEBUG_POSTFIX "-d")
    
    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES COMPILE_FLAGS 
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)
