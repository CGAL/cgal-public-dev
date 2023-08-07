# # Created by the script cgal_create_cmake_script
# # This is the CMake script for compiling a CGAL application.

# cmake_minimum_required(VERSION 3.25)
# project( sandbox CXX )

# find_package(CGAL REQUIRED)

# # create a target per cppfile
# file(
#   GLOB cppfiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
# )
# foreach(cppfile ${cppfiles})
#   create_single_source_cgal_program("${cppfile}")
# endforeach()
cmake_minimum_required(VERSION 3.1...3.23)
project(Collisions_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

#create the executable of the application
create_single_source_cgal_program("collision_candidates.cpp")
create_single_source_cgal_program("update_functor.cpp")
target_include_directories(collision_candidates PRIVATE ../include)
target_include_directories(update_functor PRIVATE ../include)
if(CGAL_Qt5_FOUND)
  #link it with the required CGAL libraries
  target_link_libraries(collision_candidates PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(update_functor PUBLIC CGAL::CGAL_Basic_viewer)
else()
  message(FATAL_ERROR Couldn't find qt5)
endif()