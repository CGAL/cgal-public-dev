# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)

project( demo )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

find_package(Eigen3 3.3.0 REQUIRED)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: These examples require the Eigen library (3.3 or greater), and will not be compiled.")
  return()
endif()


# Check for ccache
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory
set(demo_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(demo_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADER_GLOB ${demo_INCLUDE_DIR}/*.hpp)
file(GLOB SOURCES_GLOB ${demo_SOURCE_DIR}/*.cpp)

set(demo_HEADERS ${HEADER_GLOB})
set(demo_SOURCES ${SOURCES_GLOB})

# include for local package


# Creating entries for target: demo
# ############################

add_executable(demo ${demo_SOURCES} ${demo_HEADERS})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS demo )

# Link the executable to CGAL and third-party libraries
target_link_libraries(demo PRIVATE CGAL::CGAL CGAL::Eigen3_support)

target_include_directories(
  demo PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                         $<BUILD_INTERFACE:${demo_INCLUDE_DIR}>)
