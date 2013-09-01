# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( cyclide )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core Qt3 )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )

# include helper file
include( ${CGAL_USE_FILE} )

# Qt3
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".
find_package(Qt3-patched QUIET )

if ( NOT QT3_FOUND )

  message(STATUS "This project requires the Qt3 library, and will not be compiled.")
  return()  

endif()

if ( CGAL_Qt3_FOUND AND QT3_FOUND )
  
  include( Qt3Macros-patched )

endif()

find_package(GLUT)

if ( NOT GLUT_FOUND )

  message(STATUS "This project requires the Glut library, and will not be compiled.")
  return()  

endif()

# Creating entries for target: cyclide
# ############################

if ( CGAL_Qt3_FOUND AND QT3_FOUND )

qt3_automoc(  cyclide.cpp mainwnd.cpp math.cpp moc_mainwnd.cpp moc_ogl_view.cpp ogl_exts.cpp ogl_view.cpp render.cpp zz_arrangements.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories( BEFORE ${CMAKE_CURRENT_BINARY_DIR} )
 
  include_directories( ${QT3_INCLUDE_DIR} )

endif()

add_executable( cyclide  cyclide.cpp mainwnd.cpp math.cpp moc_mainwnd.cpp moc_ogl_view.cpp ogl_exts.cpp ogl_view.cpp render.cpp zz_arrangements.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS cyclide )

# Link the executable to CGAL and third-party libraries
target_link_libraries(cyclide  ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GLUT_LIBRARIES} )

