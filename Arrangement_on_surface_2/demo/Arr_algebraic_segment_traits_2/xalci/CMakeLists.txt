# This is the CMake script for compiling a CGAL application.

project( Arrangement_2 )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Core Qt3 )
include( ${CGAL_USE_FILE} )

find_package( Qt3 )

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND)
  if ( APPLE )
    add_definitions( -Wno-c++11-extensions )
  endif ( APPLE )
  # include( ${QT_USE_FILE} )
  include_directories( ./ )
  # Arrangement package includes
  include_directories( BEFORE ../../include )

  add_executable( arrangement_2 misc.cpp  poly2ntl.cpp  xalci.cpp )

  target_link_libraries( arrangement_2
    ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
    ${QT_LIBRARIES}
  )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
