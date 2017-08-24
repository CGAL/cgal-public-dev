
project( HexEx )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

OPTION(CGAL_USE_LCC_VIEWER "Use LCC viewer." OFF)

SET ( CMAKE_CXX_FLAGS "-std=c++11 -ggdb" )

include_directories(/usr/include/eigen3)

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

  add_executable(display-gmap display-gmap.cpp)
  target_link_libraries(display-gmap ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${MAP_VIEWER_LIBRARIES})

  add_executable(main hexextr.cpp post_processing.cpp connection_extraction.cpp hexahedron_extraction.cpp sanitization.cpp preprocessing.cpp main.cpp)
  target_link_libraries(main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${MAP_VIEWER_LIBRARIES})

  #add_executable(test3 hexextr.cpp post_processing.cpp connection_extraction.cpp hexahedron_extraction.cpp sanitization.cpp preprocessing.cpp test3.cpp)
  #target_link_libraries(test3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${MAP_VIEWER_LIBRARIES})

  add_executable(maintest maintest.cpp)
  target_link_libraries(maintest ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${MAP_VIEWER_LIBRARIES})
  



if (CGAL_USE_LCC_VIEWER)
  qt5_use_modules(main ${MAP_VIEWER_MODULES})
  qt5_use_modules(test3 ${MAP_VIEWER_MODULES})
  qt5_use_modules(maintest ${MAP_VIEWER_MODULES})
  qt5_use_modules(display-gmap ${MAP_VIEWER_MODULES})
endif (CGAL_USE_LCC_VIEWER)
  
else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

