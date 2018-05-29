project(func)

cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(POLICY CMP0043)
    cmake_policy( SET CMP0043  OLD )
endif()

find_package(CGAL COMPONENTS)
include( ${CGAL_USE_FILE} )

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Add specific Find.cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules )


set(HDRS triangulation.h)

set(SRCS main.cpp)

    include(AddFileDependencies)

    add_file_dependencies( main.cpp triangulation.h )

    # The executable itself.
    add_executable( func ${SRCS} ${HDRS})

#    target_link_libraries( func ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries( func ${CGAL_LIBRARIES} )
