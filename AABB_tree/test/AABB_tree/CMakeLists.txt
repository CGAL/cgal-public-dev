# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( AABB_tree_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

  add_executable(aabb_distance_triangle_test-NOKDTREE aabb_distance_triangle_test.cpp)
  target_compile_definitions(aabb_distance_triangle_test-NOKDTREE PRIVATE -DNO_KDTREE)
  target_link_libraries(aabb_distance_triangle_test-NOKDTREE CGAL::CGAL)
  add_test(NAME aabb_distance_triangle_test-NOKDTREE
    COMMAND aabb_distance_triangle_test-NOKDTREE
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

