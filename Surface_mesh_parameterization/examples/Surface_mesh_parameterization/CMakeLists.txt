if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(NAME Surface_mesh_parameterization_examples)

# Require packages new or improved since CGAL 3.4
include_directories (BEFORE ../../../Installation/include/)

# Include this package's headers first
include_directories (BEFORE . include ../../include)

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
  
  # Print new compilation options
  message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
  message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
  message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
  message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
endif()

# g++ 4.1 creates an infinite loop in Surface_mesh_parameterization
# if strict-aliasing optimization is on (implied by -O2, -O3 and -Os).
if (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-strict-aliasing")
  SET (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -fno-strict-aliasing")

  # Print new compilation options
  string( TOUPPER "${CMAKE_BUILD_TYPE}" CGAL_BUILD_TYPE_UPPER )
  message( STATUS "USING CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
  message( STATUS "USING EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
endif()

CGAL_example(SOURCE Authalic_parameterization.cpp REQUIRED_DEPENDENCIES CGAL 
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(SOURCE Mesh_cutting_parameterization.cpp REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(SOURCE Simple_parameterization.cpp REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES Boost_program_options)
CGAL_example(SOURCE Square_border_parameterization.cpp REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES Boost_program_options)

CGAL_example(SOURCE Complete_parameterization_example.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(SOURCE Eigen_parameterization.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(SOURCE Taucs_parameterization.cpp REQUIRED_DEPENDENCIES CGAL TAUCS)
CGAL_example(SOURCE polyhedron_ex_parameterization.cpp REQUIRED_DEPENDENCIES CGAL 
  OPTIONAL_DEPENDENCIES TAUCS)
