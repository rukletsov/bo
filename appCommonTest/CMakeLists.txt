
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (AppCommonTest)

SET (AppCommonTest_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
)

# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (${LIBCOMMON_ROOT}/include
)
     
# For additional dependencies check the root CMakeLists.
LINK_DIRECTORIES (${LIBRARY_OUTPUT_PATH}
)

# Compile our Library.
ADD_EXECUTABLE (AppCommonTest
                ${AppCommonTest_SRCS}
)

TARGET_LINK_LIBRARIES (AppCommonTest
                       debug LibCommon-d optimized LibCommon
)

IF (MSVC)
    # Specify output path and name suffix for debug build.
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES PREFIX "../")
    SET_TARGET_PROPERTIES (AppCommonTest PROPERTIES DEBUG_POSTFIX "-d")
ENDIF (MSVC)
