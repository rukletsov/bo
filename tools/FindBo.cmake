#
# Locate the Bö library.
#
# Defines the following variables:
#
#   Bo_FOUND - Found the Bö library.
#   Bo_INCLUDE_DIRS - Include directories.
#   Bo_LIBRARIES - Release [and debug] versions of libbo.
#
# Accepts the following variables as input:
#
#   BO_ROOT - (as an environment variable) The root directory of Bö.
#
# Example Usage:
#
#    # Don't forget to provide CMake with the path to this file.
#    # For example, if you put this file into <your_project_root>/cmake-modules,
#    # the following command will make it visible to CMake.
#    SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)
#
#    FIND_PACKAGE (Bo REQUIRED)
#    INCLUDE_DIRECTORIES (${Bo_INCLUDE_DIRS})
#
#    ADD_EXECUTABLE (BoUser bo_user.cpp)
#    TARGET_LINK_LIBRARIES (BoUser ${Bo_LIBRARIES})
#


# Finds particular version (debug | release, .a | .so) of the Bö library.
FUNCTION (FindBoLibrary libvar libname)
    FIND_LIBRARY (${libvar}
        NAMES ${libname} 
        HINTS $ENV{BO_ROOT} ${Bo_INCLUDE_DIR}
        PATH_SUFFIXES lib bin
    )

    # This shows the variable only in the advanced mode of cmake-gui.
    MARK_AS_ADVANCED (${libvar})
ENDFUNCTION (FindBoLibrary libvar libname)


# Find the directory, that is parent to the "bo" directory, where library headers are.
# Show this variable only in the advanced mode of cmake-gui.
FIND_PATH (Bo_INCLUDE_DIR bo/config.hpp
    HINTS $ENV{BO_ROOT}
)
MARK_AS_ADVANCED (BO_INCLUDE_DIR)

# Locate release and debug versions of the library.
FindBoLibrary(Bo_LIBRARY bo)
FindBoLibrary(Bo_LIBRARY_DEBUG bod)

# Use the standard CMake tool to handle FIND_PACKAGE() options and set the BO_FOUND
# variable. Note that the debug version is not required to be present.
INCLUDE (${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (Bo DEFAULT_MSG Bo_INCLUDE_DIR Bo_LIBRARY)

# Define variable Bo_FOUND additionally to BO_FOUND for consistency.
SET (Bo_FOUND BO_FOUND)

# In case of success define more variables using CMake de facto naming conventions.
IF (Bo_FOUND)
    SET (Bo_INCLUDE_DIRS ${Bo_INCLUDE_DIR})
    
    IF (Bo_LIBRARY AND Bo_LIBRARY_DEBUG)
        SET (Bo_LIBRARIES debug ${Bo_LIBRARY_DEBUG} optimized ${Bo_LIBRARY})
    ELSE ()
        SET (Bo_LIBRARIES ${Bo_LIBRARY})
    ENDIF (Bo_LIBRARY AND Bo_LIBRARY_DEBUG)
ENDIF (Bo_FOUND)
