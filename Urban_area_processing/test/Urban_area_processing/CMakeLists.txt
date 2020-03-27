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

# Creating entries for all C++ files with "main" routine.
include(CGAL_CreateSingleSourceCGALProgram)

set(targets 
  test_boundary_extraction_LIDAR
  test_boundary_extraction_with_repair
  test_extract_boundary_with_holes
  test_merge_and_orient_segments
  test_merge_and_orient_contours
  test_filter_points
  test_find_holes
)

foreach(target ${targets})
  create_single_source_cgal_program("${target}.cpp")
  if(TARGET ${target})
    target_link_libraries(${target} PUBLIC ${Boost_LIBRARIES})
  endif()
endforeach()
