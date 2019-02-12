# Created by the script cgal_create_CMakeLists This is the CMake script for
# compiling a set of CGAL applications.

project(Levels_of_detail_Test)

cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_STANDARD 11)

# CGAL and its components.
find_package(CGAL QUIET COMPONENTS)

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

# 3rd party libraries.

# OpenCV.
find_package(OpenCV 3.2 REQUIRED)
if(NOT OpenCV_FOUND)

  message(
    STATUS
      "This project requires the OpenCV 3.2 library, and will not be compiled.")
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

# Kinetic.
include_directories("/usr/local/include/kinetic")
add_library(kinetic STATIC IMPORTED)
set_property(TARGET kinetic
             PROPERTY IMPORTED_LOCATION /usr/local/lib/libkinetic.a)
set(KINETIC_LIBRARY kinetic)

# Creating entries for all C++ files with "main" routine.
include(CGAL_CreateSingleSourceCGALProgram)
add_executable(test_lod_basic test_lod_basic.cpp)

target_link_libraries(test_lod_basic
                      ${Boost_LIBRARIES}
                      ${OpenCV_LIBS}
                      ${GDAL_LIBRARY}
                      ${KINETIC_LIBRARY})
