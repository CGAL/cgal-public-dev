include(CGAL_Macros)

include(CMakeParseArguments)


function(CGAL_external_library)
  set(options NO_OPTION) # no options
  set(oneValueArgs NAME DEFAULT_ENABLED PREFIX VERSION)
  set(multiValueArgs) # no multi-value args
  cmake_parse_arguments(CGAL_external_library "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  list(APPEND CGAL_EXTERNAL_LIBRARIES ${CGAL_external_library_NAME})
  set(CGAL_EXTERNAL_LIBRARIES ${CGAL_EXTERNAL_LIBRARIES} PARENT_SCOPE)
  
  if(${CGAL_external_library_NO_OPTION})
    set(WITH_${CGAL_external_library_NAME} ON PARENT_SCOPE)
  else()
    option(WITH_${CGAL_external_library_NAME}
      "Enable support for the external library ${CGAL_external_library_NAME}"
      ${CGAL_external_library_DEFAULT_ENABLED})
  endif()

  if("${CGAL_external_library_PREFIX}" STREQUAL "")
    # The name is the prefix
    set(CGAL_${CGAL_external_library_NAME}_PREFIX ${CGAL_external_library_NAME} PARENT_SCOPE)
  else()
    set(CGAL_${CGAL_external_library_NAME}_PREFIX ${CGAL_external_library_PREFIX} PARENT_SCOPE)
  endif()

  if(NOT "${CGAL_external_library_VERSION}" STREQUAL "")
    set(CGAL_${CGAL_external_library_NAME}_VERSION ${CGAL_external_library_VERSION} PARENT_SCOPE)
  endif()
endfunction()

set(CGAL_EXTERNAL_LIBRARIES)
# This is the place to tell which external libs are supported.
CGAL_external_library(NAME GMP
  DEFAULT_ENABLED ON)
CGAL_external_library(NAME MPFR
  DEFAULT_ENABLED ON)
CGAL_external_library(NAME Qt3
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME QT
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME ZLIB
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME OpenGL
  DEFAULT_ENABLED ON
  PREFIX OPENGL)
CGAL_external_library(NAME LEDA
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME MPFI
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME RS
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME RS3
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME OpenNL
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME TAUCS
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME Eigen3
  DEFAULT_ENABLED OFF
  PREFIX EIGEN3)
CGAL_external_library(NAME BLAS
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME LAPACK
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME QGLVIEWER
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME ESBTL
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME Coin3D
  DEFAULT_ENABLED OFF
  PREFIX COIN3D)
CGAL_external_library(NAME NTL
  DEFAULT_ENABLED OFF)
CGAL_external_library(NAME IPE
  DEFAULT_ENABLED OFF)
# Boost is non-negotiable
CGAL_external_library(NAME Boost
  DEFAULT_ENABLED OFF
  NO_OPTION
  VERSION 1.33.1)
#  There exists FindF2C, FindMKL, but they are only used to
#  support supporting libs.

if(NOT WIN32)
  # GMPXX is not supported on WIN32 machines
  CGAL_external_library(NAME GMPXX
    DEFAULT_ENABLED OFF)
endif()

foreach(lib ${CGAL_EXTERNAL_LIBRARIES})
  set(vlib ${CGAL_${lib}_PREFIX})

  if(WITH_${lib})
    find_package(${lib} ${CGAL_${lib}_VERSION} QUIET)
    if(${vlib}_FOUND)
      message(STATUS "${lib} has been found:") 
      message(STATUS "  Use${lib}-file:      ${${vlib}_USE_FILE}") 
      message(STATUS "  ${lib} include:      ${${vlib}_INCLUDE_DIR}")
      message(STATUS "  ${lib} libraries:    ${${vlib}_LIBRARIES}")
      message(STATUS "  ${lib} definitions:  ${${vlib}_DEFINITIONS}")
    else()
      # Never allow erroneous configurations. Should we rather use
      # find_package(... REQUIRED)?
      message(FATAL_ERROR "The external library ${lib} has been requested, but could not be found.")
    endif()
  endif()
endforeach()
unset(vlib)

include(CGAL_TweakFindBoost)

# Special handling still required.
#
# if (${lib} STREQUAL "GMP") 
#   get_dependency_version(GMP)
# endif()

# if (${lib} STREQUAL "MPFR") 
#   set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
#   set( MPFR_DEPENDENCY_LIBRARIES   ${GMP_LIBRARIES} )
#   get_dependency_version(MPFR)
# endif()

# if (${lib} STREQUAL "LEDA") 
#   # special case for LEDA - add a flag
#   message( STATUS "$LEDA cxx flags:   ${LEDA_CXX_FLAGS}" )
#   uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
# endif()
