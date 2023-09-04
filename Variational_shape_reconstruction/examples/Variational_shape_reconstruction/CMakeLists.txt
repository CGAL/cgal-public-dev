# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Variational_shape_reconstruction_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires Eigen 3.1 (or greater), and will not be compiled.")
  return()
endif()
#set(CMAKE_BUILD_TYPE Debug)
find_package(SCIP QUIET)
include(CGAL_SCIP_support)
if(NOT TARGET CGAL::SCIP_support)
  find_package(GLPK QUIET)
  include(CGAL_GLPK_support)
  if(NOT TARGET CGAL::GLPK_support)
    message("NOTICE: This project requires either SCIP or GLPK, and will not be compiled.")
    return()
  endif()
endif()



# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
    create_single_source_cgal_program("${cppfile}")
    get_filename_component(target ${cppfile} NAME_WE)
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
    if(TARGET CGAL::SCIP_support)
      target_link_libraries(${target} PUBLIC CGAL::SCIP_support)
    else()
      target_link_libraries(${target} PUBLIC CGAL::GLPK_support)
    endif()
endforeach()
