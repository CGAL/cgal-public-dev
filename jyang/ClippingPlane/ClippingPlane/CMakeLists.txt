# cmake ../ -DCMAKE_BUILD_TYPE=Debug

cmake_minimum_required(VERSION 3.1...3.15)
project (Clipping_plane_test)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

##################
find_package(CGAL REQUIRED COMPONENTS Qt5)
find_package(OpenGL REQUIRED)

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt5 REQUIRED COMPONENTS Xml Script OpenGL Svg)

add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)

# ui file, created wih Qt Designer
qt5_wrap_ui(uis clippingplane.ui)

# qrc files (resources files, that contain icons, at least)
# qt5_add_resources (CGAL_Qt5_RESOURCE_FILES file.qrc)

add_executable(test_clipping_plane main.cpp clippingplane.cpp widget.cpp
  ${uis} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

target_link_libraries(test_clipping_plane PUBLIC
  CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui Qt5::OpenGL OpenGL::GL)
