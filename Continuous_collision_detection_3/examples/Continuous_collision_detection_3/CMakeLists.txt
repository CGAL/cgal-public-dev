# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.1...3.23)
project(Collisions_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# create the executable of the application
# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})

  create_single_source_cgal_program("${cppfile}")

  get_filename_component(current_program ${cppfile} NAME_WLE)
  
  target_include_directories(${current_program} PRIVATE ../../include/CGAL PRIVATE ../../include/CGAL/Continuous_collision_detection_3/internal)

  if(CGAL_Qt5_FOUND)
    target_link_libraries(${current_program} PUBLIC CGAL::CGAL_Basic_viewer)
  else()
    message(FATAL_ERROR Couldn't find qt5)
  endif()

endforeach()