
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (Common)

SET (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Gloabl dependencies (for all projects).
INCLUDE_DIRECTORIES ($ENV{BOOST_ROOT}
)

LINK_DIRECTORIES ($ENV{BOOST_ROOT}/stage/lib
)

# Projects list with root paths.
SET (LIBCOMMON_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/libCommon)
SET (APPCOMMONTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/appCommonTest)

ADD_SUBDIRECTORY (${LIBCOMMON_ROOT})
ADD_SUBDIRECTORY (${APPCOMMONTEST_ROOT})