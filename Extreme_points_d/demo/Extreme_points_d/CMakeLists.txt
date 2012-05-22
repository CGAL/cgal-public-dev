# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Extreme_points_Demo)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(Qt4)

# TODO: remove last part..
include_directories (BEFORE ../../include /home/chrigu/cgal_svn/trunk/GraphicsView/include/)

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

  include(${QT_USE_FILE})
  #--------------------------------
  # Demo: Extreme_points_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt4_wrap_ui( DT_UI_FILES Extreme_points_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( DT_RESOURCE_FILES ./Extreme_points_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt4_generate_moc( Extreme_points_2.cpp Extreme_points_2.moc )

  # The executable itself.
  add_executable  ( Extreme_points_2 Extreme_points_2.cpp Extreme_points_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Extreme_points_2 )

  # Link with Qt libraries
  target_link_libraries( Extreme_points_2 ${QT_LIBRARIES} )
  # And with CGAL libraries
  target_link_libraries( Extreme_points_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, and Qt4, and will not be compiled.")

endif()
