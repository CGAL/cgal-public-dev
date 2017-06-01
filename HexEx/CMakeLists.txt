
project( HexEx )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL QUIET)


if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../include)


  add_executable(triangulation_to_LCC triangulation_to_LCC.cpp)
  target_link_libraries(triangulation_to_LCC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${MAP_VIEWER_LIBRARIES})
  
else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

