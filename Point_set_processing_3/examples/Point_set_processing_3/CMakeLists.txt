if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(NAME Point_set_processing_3_examples)

# VisualC++ optimization for applications dealing with large data
if (MSVC)
  # Use /EHa option to catch stack overflows.
  # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
  string(REGEX REPLACE "/EH[asc]*" "/EHa" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

  # Use /FR to turn on IntelliSense
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FR")

  # Allow Windows applications to use up to 3GB of RAM
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

  # Turn off stupid VC++ warnings
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")

  # Prints new compilation options
  message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
  message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
  message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
  message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
endif()

# Executables that do *not* require EIGEN or LAPACK
CGAL_example(SOURCE average_spacing_example.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE grid_simplification_example.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE normals_example.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE property_map.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE random_simplification_example.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE read_write_xyz_point_set_example.cpp REQUIRED_DEPENDENCIES CGAL)
CGAL_example(SOURCE remove_outliers_example.cpp REQUIRED_DEPENDENCIES CGAL)

CGAL_example(SOURCE jet_smoothing_example.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(SOURCE normal_estimation.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(SOURCE jet_smoothing_example.cpp
  TARGET jet_smoothing_example_LAPACK
  REQUIRED_DEPENDENCIES CGAL LAPACK)
CGAL_example(SOURCE normal_estimation.cpp
  TARGET normal_estimation_LAPACK
  REQUIRED_DEPENDENCIES CGAL LAPACK)
