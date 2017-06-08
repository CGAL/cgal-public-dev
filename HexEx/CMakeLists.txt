
project( HexEx )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

OPTION(CGAL_USE_LCC_VIEWER "Use LCC viewer." OFF)

SET ( CMAKE_CXX_FLAGS "-std=c++11" )

if (CGAL_USE_LCC_VIEWER)
  find_package(CGAL COMPONENTS Qt5)
  include("CMakeLCCViewerQt.inc")
else (CGAL_USE_LCC_VIEWER)
  find_package(CGAL QUIET)
endif (CGAL_USE_LCC_VIEWER)
  
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../include)


  add_executable(HexExtr HexExtr.cpp)
  target_link_libraries(HexExtr ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
    ${MAP_VIEWER_LIBRARIES})

#if (CGAL_USE_LCC_VIEWER)
 # qt5_use_modules(triangulation_to_LCC ${MAP_VIEWER_MODULES})
#endif (CGAL_USE_LCC_VIEWER)
  
else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

