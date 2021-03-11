# Created by the script cgal_create_CMakeLists This is the CMake script for
# compiling a set of CGAL applications.
project(Levels_of_detail_Tests)

# kinetic 3
cmake_policy(SET CMP0042 NEW)

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

# CGAL and its components.
find_package(CGAL QUIET COMPONENTS Core)
if(NOT CGAL_FOUND)
  message(
    STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# Required by kinetic3.
# set(CMAKE_CXX_FLAGS
# "${CMAXE_CXX_FLAGS} -DBOOST_NO_CXX11_SCOPED_ENUMS -fdeclspec")

set(CMAKE_CXX_FLAGS
"${CMAXE_CXX_FLAGS} -DBOOST_NO_CXX11_SCOPED_ENUMS -fdeclspec -Wno-format -Wno-ignored-attributes -Wno-parentheses -Wno-return-type -Wno-macro-redefined -Wno-switch -Wno-logical-op-parentheses -Wno-deprecated-declarations -Wno-user-defined-warnings -Wno-unused-result")

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

# 3rd party libraries.

# OpenCV.
find_package(OpenCV 3.4 REQUIRED)
if(NOT OpenCV_FOUND)
  message(
    STATUS
      "This project requires the OpenCV library, and will not be compiled.")
  return()
endif()

# GDAL.
find_package(GDAL REQUIRED)
if(NOT GDAL_FOUND)
  message(
    STATUS "This project requires the GDAL library, and will not be compiled.")
  return()
endif()

include_directories(BEFORE include)
include_directories(${EIGEN3_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${OpenCV_INCLUDE_DIRS}
                    ${GDAL_INCLUDE_DIR})

# Kinetic 2.
include_directories("/usr/local/include")
add_library(kinetic2 STATIC IMPORTED)
set_property(TARGET kinetic2
             PROPERTY IMPORTED_LOCATION /usr/local/lib/libkinetic2.a)
set(KINETIC_LIBRARY2 kinetic2)

# Kinetic 3.
include_directories("/usr/local/include")
add_library(kinetic3 STATIC IMPORTED)
set_property(TARGET kinetic3
             PROPERTY IMPORTED_LOCATION /usr/local/lib/libkinetic3.a)
set(KINETIC_LIBRARY3 kinetic3)

set(sources
test_lod_basic.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/event_queue.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/event.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/intersection_line.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/octree_base_vertex.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/octree_base.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/oriented_bounding_box_2.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/oriented_bounding_box_3.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/oriented_bounding_box.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/parameters.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_edge.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_facet.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_multiple.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_polyhedron.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_simple.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_vertex_octree.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition_vertex.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/partition.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/planar_segment.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/ply_in.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/point_cloud_vertex.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_directions.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_edge.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_group.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_link.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_node.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_segment_r.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_segment_s.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_segment.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_set.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_tree.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_vertex_octree_data.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_vertex_octree.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_vertex_r.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_vertex_s.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon_vertex.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/polygon.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/propagation_multiple.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/propagation_simple.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/propagation.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/segment_translation.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/segment.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/stats.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_construct.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_draw.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_event_handlers.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_object.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_segment_initializers.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane_stopping_condition.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/support_plane.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/universe.cpp
../../include/CGAL/Levels_of_detail/internal/Partitioning/Kinetic/src/vars.cpp
)

# OSQP.
find_package(osqp REQUIRED)

# FLANN.
add_library(flann SHARED IMPORTED)
set_property(TARGET flann
             PROPERTY IMPORTED_LOCATION /usr/local/Cellar/flann/1.9.1_10/lib/libflann.dylib)
set(FLANN_LIBRARY flann)

# Creating entries for all C++ files with "main" routine.
include(CGAL_CreateSingleSourceCGALProgram)
add_executable(test_lod_basic ${sources})

target_link_libraries(test_lod_basic
                      ${Boost_LIBRARIES}
                      ${OpenCV_LIBS}
                      ${GDAL_LIBRARY}
                      ${KINETIC_LIBRARY2}
                      ${KINETIC_LIBRARY3}
                      ${FLANN_LIBRARY}
                      osqp::osqp)
