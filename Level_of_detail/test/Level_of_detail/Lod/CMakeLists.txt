# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( lod )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED COMPONENTS filesystem)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# 3rd party libraries
find_package(GDAL REQUIRED)
if (NOT GDAL_FOUND)
  message(STATUS "This project requires the GDAL library, and will not be compiled.")
  return()
endif()
  
find_package(OpenCV 3.2 REQUIRED)
if (NOT OpenCV_FOUND)
  message(STATUS "This project requires the OpenCV 3.2 library, and will not be compiled.")
  return()
endif()

find_package(Eigen3 3.1.0 REQUIRED)
if (NOT EIGEN3_FOUND)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

include( ${EIGEN3_USE_FILE} )
include_directories(${OpenCV_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  
# 3rd parties to clean
include_directories("/usr/local/include/ooqp")
include_directories("/usr/local/include/kinetic")
set(ma27_LIBRARY ma27)
set(OOQP_LIBS ooqpbase ooqpdense ooqpgensparse ooqpgondzio ooqpsparse)
set(KINETIC_LIBRARY kinetic)

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

set(CMAKE_CXX_STANDARD 11)

set(sources main.cpp)

add_executable(lod ${sources})

target_link_libraries(lod ${OpenCV_LIBS} ${GDAL_LIBRARY} ${Boost_LIBRARIES} ${ma27_LIBRARY} ${OOQP_LIBS} ${KINETIC_LIBRARY})



