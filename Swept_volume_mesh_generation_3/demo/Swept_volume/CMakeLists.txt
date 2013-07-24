# This is the CMake script for compiling the CGAL Mesh_3 demo. 


project( Mesh_3_demo )

cmake_minimum_required(VERSION 2.4.5)

set(MESH_3_ROOT $ENV{MESH_3_ROOT})

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

link_directories(${MESH_3_ROOT}/demo/Mesh_3/)
include_directories(BEFORE ${MESH_3_ROOT}/demo/Mesh_3/)
include_directories(BEFORE ${MESH_3_ROOT}/demo/Mesh_3/include/)
include_directories(BEFORE ${MESH_3_ROOT}/include/)



include_directories (../../src)
add_subdirectory (../../src ../../src)
set (EXTRA_LIBS ${EXTRA_LIBS} swept_volume_tools)


if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4 ImageIO)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find GLEW (for OpenGL-1.5 and OpenGL extensions)
find_package(GLEW)


# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer)
endif(QT4_FOUND)


if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  if(GLEW_FOUND)
    include_directories ( ${GLEW_INCLUDE_DIR} )
    add_definitions(-DSCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE)
  else(GLEW_FOUND)
    message(STATUS "NOTICE: GLEW library is not found. 3D images rendering will not be available.")
  endif(GLEW_FOUND)


  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt4_wrap_ui( UI_FILES ${MESH_3_ROOT}/demo/Mesh_3/ui_files/MainWindow.ui )
  qt4_wrap_ui( meshingUI_FILES  ${MESH_3_ROOT}/demo/Mesh_3/ui_files/Meshing_dialog.ui)
  qt4_wrap_ui( optimUI_FILES ${MESH_3_ROOT}/demo/Mesh_3/ui_files/Smoother_dialog.ui ui_files/LocalOptim_dialog.ui)
  qt4_wrap_ui( ribUI_FILES  ${MESH_3_ROOT}/demo/Mesh_3/ui_files/Rib_dialog.ui)
  qt4_wrap_ui( funcUI_FILES  ${MESH_3_ROOT}/demo/Mesh_3/ui_files/Function_dialog.ui)

  include(AddFileDependencies)

  # Shared demo files directory
  set(DEMO_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/CGAL_demo")
  set(DEMO_H_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/CGAL_demo")

  qt4_automoc( Scene_swept_volume_3_item.cpp )

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  #set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})
  # prefix for Mesh_3 libs to avoid conflicts with Polyhedron demo
  set(MESH_3_LIB_PREFIX "mesh_3_demo_")

  set(SCENE_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_item")
  set(SCENE_C3T3_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_c3t3_item")


  add_library(scene_swept_volume_3_item SHARED
    Scene_swept_volume_3_item.cpp Scene_swept_volume_3_item.moc)
  target_link_libraries(scene_swept_volume_3_item ${SCENE_ITEM_LIB} ${EXTRA_LIBS})

  add_definitions(-DUSE_FORWARD_DECL)

  ###########
  # PLUGINS #
  ###########

  # polyhedron_demo_plugin is a macro copied from demo/Polyhedron/CMakeLists.txt
  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      set(other_sources ${ARGN})
      set(option "")
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} ${EXTRA_LIBS})
    # Link with scene_item
    target_link_libraries( ${plugin_name} ${SCENE_ITEM_LIB} ${EXTRA_LIBS} )
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${EXTRA_LIBS})
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(io_swept_volume_3_plugin Io_swept_volume_3_plugin)
  target_link_libraries(io_swept_volume_3_plugin scene_swept_volume_3_item ${EXTRA_LIBS})

  polyhedron_demo_plugin(swept_volume_3_plugin Swept_volume_3_plugin
                         ${MESH_3_ROOT}/demo/Mesh_3/Meshing_thread.cpp
                         ${MESH_3_ROOT}/demo/Mesh_3/Scene_c3t3_item.moc  
                         ${meshingUI_FILES})  

  target_link_libraries(swept_volume_3_plugin
                        ${SCENE_C3T3_ITEM_LIB}
                        scene_swept_volume_3_item
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY}
                        ${EXTRA_LIBS})

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt4 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(MESH_3_MISSING_DEPS "Qt4, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(MESH_3_MISSING_DEPS "OpenGL, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(MESH_3_MISSING_DEPS "QGLViewer, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)


