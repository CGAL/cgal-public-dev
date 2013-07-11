# This is the CMake script for compiling a CGAL application.

project( arr2lcc ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Core )

find_package( Boost COMPONENTS program_options REQUIRED )
include( ${CGAL_USE_FILE} )


if ( CGAL_FOUND AND Boost_FOUND)
  
  include_directories( ./ )
  # Arrangement package includes
  include_directories( BEFORE ../../include )

  include_directories( ${Boost_INCLUDE_DIR} )

   add_executable( arr2lcc
	arr2lcc.cpp
        convertors.h
      )   
   

  target_link_libraries( arr2lcc
    ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
    ${QT_LIBRARIES}
    ${Boost_LIBRARIES}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Boost, and will not be compiled.")

endif()

# add_subdirectory( tests )
