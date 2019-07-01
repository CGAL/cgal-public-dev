# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Shape_regularization_Examples)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Use Eigen.
  find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
  # Find OSQP library and headers
  find_package(osqp REQUIRED)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    set(sources segment_regularization.cpp)
    add_executable(Shape_regularization_Examples ${sources})
    # Link the OSQP shared library
    target_link_libraries(Shape_regularization_Examples osqp::osqp)

    #create_single_source_cgal_program(
    # "segment_regularization.cpp")
  endif()


else()
  message(WARNING 
  "This program requires the CGAL library, and will not be compiled.")
endif()
