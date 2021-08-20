# This is the CMake script for compiling the Minimal_angle_remeshing_Demo.

project( Minimal_angle_remeshing_Demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Gui Svg)

include( ${CGAL_USE_FILE} )


if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)

  qt5_wrap_ui( UI_FILES Main_window.ui isotropic_parameters.ui minimal_angle_parameters.ui)

  include(AddFileDependencies)

  qt5_generate_moc( "main_window.h" "${CMAKE_CURRENT_BINARY_DIR}/main_window_moc.cpp" )
  add_file_dependencies( main_window_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/main_window.h" )

  qt5_generate_moc( "viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/viewer_moc.cpp" )
  add_file_dependencies( viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/viewer.h" )

  qt5_generate_moc( "scene.h" "${CMAKE_CURRENT_BINARY_DIR}/scene_moc.cpp" )
  add_file_dependencies( scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/scene.h" )

  qt5_generate_moc( "isotropic_parameters.h" "${CMAKE_CURRENT_BINARY_DIR}/isotropic_parameters_moc.cpp" )
  add_file_dependencies( isotropic_parameters_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/isotropic_parameters.h" )

  qt5_generate_moc( "minimal_angle_parameters.h" "${CMAKE_CURRENT_BINARY_DIR}/minimal_angle_parameters_moc.cpp" )
  add_file_dependencies( minimal_angle_parameters_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/minimal_angle_parameters.h" )

  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES minimal_angle_remeshing_demo.qrc )

  add_file_dependencies( minimal_angle_remeshing_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/main_window_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/viewer_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/scene_moc.cpp" 
	                       "${CMAKE_CURRENT_BINARY_DIR}/isotropic_parameters_moc.cpp"
		       "${CMAKE_CURRENT_BINARY_DIR}/minimal_angle_parameters_moc.cpp")

  add_executable  ( Minimal_angle_remeshing_Demo minimal_angle_remeshing_demo.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} 
    #${CGAL_Qt5_MOC_FILES}
    )
  # Link with Qt libraries
  target_link_libraries( Minimal_angle_remeshing_Demo PRIVATE
    Qt5::OpenGL Qt5::Gui Qt5::Xml
    CGAL::CGAL
    CGAL::CGAL_Qt5
    )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Minimal_angle_remeshing_Demo )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Minimal_angle_remeshing_Demo)

else (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(RR_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(RR_MISSING_DEPS "the CGAL library, ${RR_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt5_FOUND)
    set(RR_MISSING_DEPS "CGAL_Qt5, ${RR_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(RR_MISSING_DEPS "Qt5, ${RR_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${RR_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )
