
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (LibCommon)

SET (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET (BIN_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

SET (LibCommon_SRCS
	${CMAKE_CURRENT_SOURCE_DIR}/src/stdafx.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/point.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.cpp
)

SET (LibCommon_HDRS
	${CMAKE_CURRENT_SOURCE_DIR}/include/stdafx.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/point.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/mesh.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/triangle.hpp
)

INCLUDE_DIRECTORIES (BEFORE
                    ./include
                    $ENV{BOOST_ROOT}
)
                    
LINK_DIRECTORIES ($ENV{BOOST_ROOT}/stage/lib
)

# Compile our Library.
ADD_LIBRARY (LibCommon
			${LibCommon_HDRS}
			${LibCommon_SRCS}
)

IF (MSVC)
    # Enable usage for precompiled header for MSVC compiler.
	SET_TARGET_PROPERTIES (LibCommon PROPERTIES COMPILE_FLAGS "/Yustdafx.h")    
    SET_SOURCE_FILES_PROPERTIES (${CMAKE_CURRENT_SOURCE_DIR}/src/stdafx.cpp 
								PROPERTIES COMPILE_FLAGS "/Ycstdafx.h")
    
	SET_TARGET_PROPERTIES (LibCommon PROPERTIES PREFIX "../")
	SET_TARGET_PROPERTIES (LibCommon PROPERTIES DEBUG_POSTFIX "_debug")
ENDIF (MSVC)
