# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project( Triangle_mesh_sampling )

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

create_single_source_cgal_program( "sample_elimination.cpp" )

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(sample_elimination PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Tests will use sequential code.")
endif()



find_package(Eigen3 3.3.0 REQUIRED)
include(CGAL_Eigen3_support)


target_link_libraries(sample_elimination PUBLIC CGAL::Eigen3_support)


