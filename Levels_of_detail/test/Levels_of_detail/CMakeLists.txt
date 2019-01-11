# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Levels_of_detail_Tests )

cmake_minimum_required( VERSION 2.8.11 )

set( CMAKE_CXX_STANDARD 11 )

# CGAL and its components.
find_package( CGAL QUIET COMPONENTS )

if ( NOT CGAL_FOUND )

  message( STATUS "This project requires the CGAL library, and will not be compiled." )
  return()  

endif()

# Include a helper file.
include( ${CGAL_USE_FILE} )

# Boost and its components.
find_package( Boost REQUIRED COMPONENTS filesystem )

if ( NOT Boost_FOUND )

  message( STATUS "This project requires the Boost library, and will not be compiled." )
  return()  

endif()

# 3rd party libraries.
find_package( OpenCV 3.2 REQUIRED )
if ( NOT OpenCV_FOUND )
    
  message( STATUS "This project requires the OpenCV 3.2 library, and will not be compiled." )
  return()

endif()

find_package( GDAL REQUIRED )
if ( NOT GDAL_FOUND )
    
  message( STATUS "This project requires the GDAL library, and will not be compiled." )
  return()

endif()

include_directories( ${Boost_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} )
  
# 3rd parties to clean.
include_directories( "/usr/local/include/ooqp" )
include_directories( "/usr/local/include/kinetic" )

add_library( ma27 STATIC IMPORTED )
set_property( TARGET ma27 PROPERTY IMPORTED_LOCATION /usr/local/lib/libma27.a )

add_library( qpbase STATIC IMPORTED )
set_property( TARGET qpbase PROPERTY IMPORTED_LOCATION /usr/local/lib/libooqpbase.a )

add_library( qpdense STATIC IMPORTED )
set_property( TARGET qpdense PROPERTY IMPORTED_LOCATION /usr/local/lib/libooqpdense.a )

add_library( qpgensparse STATIC IMPORTED )
set_property( TARGET qpgensparse PROPERTY IMPORTED_LOCATION /usr/local/lib/libooqpgensparse.a )

add_library( qpgondzio STATIC IMPORTED )
set_property( TARGET qpgondzio PROPERTY IMPORTED_LOCATION /usr/local/lib/libooqpgondzio.a )

add_library( qpsparse STATIC IMPORTED )
set_property( TARGET qpsparse PROPERTY IMPORTED_LOCATION /usr/local/lib/libooqpsparse.a )

add_library( kinetic STATIC IMPORTED )
set_property( TARGET kinetic PROPERTY IMPORTED_LOCATION /usr/local/lib/libkinetic.a )

set( ma27_LIBRARY ma27 )
set( OOQP_LIBS qpbase qpdense qpgensparse qpgondzio qpsparse )
set( KINETIC_LIBRARY kinetic )

# Creating entries for all C++ files with "main" routine.
include( CGAL_CreateSingleSourceCGALProgram )
add_executable( test_lod_basic test_lod_basic.cpp )

target_link_libraries(test_lod_basic 
  ${Boost_LIBRARIES} ${OpenCV_LIBS} ${GDAL_LIBRARY} 
  -L/usr/local/Cellar/gcc/8.2.0/lib/gcc/8/ -lgfortran -L/usr/local/Cellar/openblas/0.3.3/lib/ -lblas -llapack ${ma27_LIBRARY} ${OOQP_LIBS} ${KINETIC_LIBRARY} )
