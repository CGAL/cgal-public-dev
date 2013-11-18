include(CGAL_Macros)

# This is the place to tell which external libs are supported.

# TODO This is a list of pairs, the first is the library prefix, the second
# is the status of the dependency. It can have three possible values:
# - ALWAYS: this library is always required
# - ON    : this library is optional but searched for by default
# - OFF   : this library is optional but searched for by default
#
# Remarks: 
#  Coin is used in KDS, but no FindCoin or FindCOIN exists
#
#  There exists FindF2C, FindIPE, FindMKL, but they are only used to
#  support supporting libs.
set(CGAL_EXTERNAL_LIBRARIES GMP MPFR Qt3 QT ZLIB OpenGL LEDA MPFI RS RS3
  OpenNL TAUCS EIGEN3 BLAS LAPACK QGLVIEWER ESBTL COIN3D NTL IPE
  CACHE INTERNAL "List of supported extenal libraries" FORCE)

if(NOT WIN32)
  # GMPXX is not supported on WIN32 machines
  list(INSERT CGAL_EXTERNAL_LIBRARIES 1 GMPXX)
endif()

foreach(lib ${CGAL_EXTERNAL_LIBRARIES}) 
  option(WITH_${lib} "Enable support for the external library ${lib}" FALSE)

  if(WITH_${lib})
    find_package(${lib})
    if(${lib}_FOUND)
      message(STATUS "${lib} has been found:") 
      message(STATUS "  Use${lib}-file:      ${${lib}_USE_FILE}") 
      message(STATUS "  ${lib} include:      ${${lib}_INCLUDE_DIR}")
      message(STATUS "  ${lib} libraries:    ${${lib}_LIBRARIES}")
      message(STATUS "  ${lib} definitions:  ${${lib}_DEFINITIONS}")
    else()
      # Never allow erroneous configurations. Should we rather use
      # find_package(... REQUIRED)?
      message(ERROR "The external library ${lib} has been requested, but could not be found.")
    endif()
  endif()
endforeach()

include(CGAL_TweakFindBoost)
# In the documentation, we say we require Boost-1.39, but technically we
# require 1.33.1. Some packages may require more recent versions, though.
find_package(Boost 1.33.1 REQUIRED thread system)
message(STATUS "Using BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'")

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

# finally setup Boost

