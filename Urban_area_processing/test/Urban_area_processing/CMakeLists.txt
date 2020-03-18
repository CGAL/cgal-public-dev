# Created by the script cgal_create_CMakeLists This is the CMake script for
# compiling a set of CGAL applications.
project(Urban_area_processing_Tests)

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

# CGAL and its components.
find_package(CGAL QUIET COMPONENTS Core)
if(NOT CGAL_FOUND)
  message(
    STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# Include a helper file.
include(${CGAL_USE_FILE})

# Eigen and its components.
find_package(EIGEN3 REQUIRED)
if(NOT EIGEN3_FOUND)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

# Boost and its components.
find_package(Boost REQUIRED COMPONENTS filesystem)
if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

include_directories(BEFORE include)
include_directories(${EIGEN3_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS})

set(sources
test_boundary_extraction_LIDAR.cpp)

# Creating entries for all C++ files with "main" routine.
include(CGAL_CreateSingleSourceCGALProgram)
add_executable(test_boundary_extraction_LIDAR ${sources})

target_link_libraries(test_boundary_extraction_LIDAR
                      ${Boost_LIBRARIES})
