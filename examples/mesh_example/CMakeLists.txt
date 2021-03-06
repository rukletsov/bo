
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.4)

SET (MeshExample_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_example.cpp
)

# For additional dependencies check the root CMakeLists.
INCLUDE_DIRECTORIES (${BO_ROOT}
)
     
# For additional dependencies check the root CMakeLists.
LINK_DIRECTORIES (${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
                  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)

# Compile the application.
ADD_EXECUTABLE (MeshExample
                ${MeshExample_SRCS}
)

# Add dependencies.
TARGET_LINK_LIBRARIES (MeshExample
                       ${Boost_LIBRARIES}
                       bo
)

ADD_CUSTOM_COMMAND(TARGET MeshExample
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${EXAMPLES_PATH}/data/mesh"
        "$<TARGET_FILE_DIR:MeshExample>/data/mesh"
        COMMENT "Copying Tiling's data to the runtime directory."
)

# Specify the name suffix for debug build.
SET_TARGET_PROPERTIES (MeshExample PROPERTIES DEBUG_POSTFIX "d")

# Compiler-specific features and options.
IF (MSVC)
    # Supress secure warnings (C4996) for MSVC compiler.
    SET_TARGET_PROPERTIES (MeshExample PROPERTIES COMPILE_FLAGS
                           "/D _CRT_SECURE_NO_DEPRECATE /D _SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)
