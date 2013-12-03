if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Surface_reconstruction_points_3
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES)

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

# Temporary debugging stuff
ADD_DEFINITIONS( "-DDEBUG_TRACE" ) # turn on traces

CGAL_example(SOURCE poisson_reconstruction_example.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(SOURCE poisson_reconstruction.cpp REQUIRED_DEPENDENCIES CGAL Eigen3)
CGAL_example(TARGET poisson_reconstruction_example_taucs 
  SOURCE poisson_reconstruction_example.cpp REQUIRED_DEPENDENCIES CGAL TAUCS)
CGAL_example(TARGET poisson_reconstruction_example_taucs 
  SOURCE poisson_reconstruction.cpp REQUIRED_DEPENDENCIES CGAL TAUCS)
