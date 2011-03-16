
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

PROJECT (LibCommon)

SET (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
SET (BIN_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

SET (LibCommon_SRCS
	${CMAKE_CURRENT_SOURCE_DIR}/src/stdafx.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/point.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.cpp
)

# A special section for C sources, which should be compilde using C compiler, not C++
SET (LibCommon_CSRCS
	${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.01/rply.c
)

SET (LibCommon_HDRS
	${CMAKE_CURRENT_SOURCE_DIR}/include/stdafx.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/point.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/mesh.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/triangle.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/errors.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/performance.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/common/raw_image.hpp
	
	${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/rply-1.01/rply.h
)

INCLUDE_DIRECTORIES (BEFORE
                    ./include
                    ./include/common
                    ./3rdparty/rply-1.01
                    $ENV{BOOST_ROOT}
)
                    
LINK_DIRECTORIES ($ENV{BOOST_ROOT}/stage/lib
)

# Compile our Library.
ADD_LIBRARY (LibCommon
			${LibCommon_HDRS}
			${LibCommon_SRCS}
			${LibCommon_CSRCS}
)

IF (MSVC)
  # Enable usage for precompiled header for MSVC compiler for C++ sources.
	# SET_TARGET_PROPERTIES (LibCommon PROPERTIES COMPILE_FLAGS "/Yustdafx.h")    
	FOREACH (src_file ${LibCommon_SRCS})
    SET_SOURCE_FILES_PROPERTIES (${src_file} PROPERTIES COMPILE_FLAGS "/Yustdafx.h")
  ENDFOREACH (src_file ${LibCommon_SRCS})

  SET_SOURCE_FILES_PROPERTIES (${CMAKE_CURRENT_SOURCE_DIR}/src/stdafx.cpp 
	                					  	PROPERTIES COMPILE_FLAGS "/Ycstdafx.h")
    
	SET_TARGET_PROPERTIES (LibCommon PROPERTIES PREFIX "../")
	SET_TARGET_PROPERTIES (LibCommon PROPERTIES DEBUG_POSTFIX "_debug")
ENDIF (MSVC)
