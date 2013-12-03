if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(NAME Mesh_3_examples)

include_directories(../../include)
include_directories(../../../Triangulation_3/include)
include_directories(../../../STL_Extension/include)
include_directories(../../../AABB_tree/include)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if (MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

CGAL_example(SOURCE mesh_implicit_sphere.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE mesh_implicit_sphere_variable_size.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE mesh_two_implicit_spheres_with_balls.cpp REQUIRED_DEPENDENCIES CGAL)
# create_single_source_cgal_program( "mesh_implicit_domains.cpp" "implicit_functions.cpp" )
CGAL_example(SOURCE mesh_polyhedral_domain.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE mesh_polyhedral_domain_with_features.cpp REQUIRED_DEPENDENCIES CGAL)

CGAL_example(SOURCE mesh_optimization_example.cpp REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO)
CGAL_example(SOURCE mesh_optimization_lloyd_example.cpp REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO)

if(CGAL_ImageIO_USE_ZLIB)
  CGAL_example(SOURCE mesh_3D_image.cpp REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO)
  CGAL_example(SOURCE mesh_3D_image_variable_size.cpp REQUIRED_DEPENDENCIES CGAL CGAL_ImageIO)
else()
  message( STATUS "The examples mesh_3D_image.cpp and mesh_3D_image_variable_size.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
endif()
