# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Swept_volume_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories (../../src)
add_subdirectory (../../src ../../src)
set (EXTRA_LIBS ${EXTRA_LIBS} swept_volume_tools)

include_directories(BEFORE ../../demo/Swept_volume/)
include_directories(BEFORE ../../include/)
include_directories(BEFORE ./include/)

set(MESH_3_ROOT $ENV{MESH_3_ROOT})
link_directories(${MESH_3_ROOT}/demo/Mesh_3/)
include_directories(BEFORE ${MESH_3_ROOT}/demo/Mesh_3/)
include_directories(BEFORE ${MESH_3_ROOT}/demo/Mesh_3/include/)
include_directories(BEFORE ${MESH_3_ROOT}/include/)

#include_directories(BEFORE $ENV{CGAL_TRUNK}/Triangulation_3/include/)


set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()


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
# find_package(GLEW)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

if(CUDA_FOUND)
  create_single_source_cgal_program( "Swept_volume_with_vhull_3.cpp")
  target_link_libraries (Swept_volume_with_vhull_3 ${EXTRA_LIBS} )

else()  
#  create_single_source_cgal_program( "Octree_3.cpp" )
#  target_link_libraries (Octree_3 ${EXTRA_LIBS})
#  create_single_source_cgal_program( "HashSet_Voxel.cpp" )
#  target_link_libraries (HashSet_Voxel ${EXTRA_LIBS})
#  create_single_source_cgal_program( "Voxelization_3.cpp")
#  target_link_libraries (Voxelization_3 ${EXTRA_LIBS})
  create_single_source_cgal_program( "Swept_volume_with_vhull_3.cpp")
  target_link_libraries (Swept_volume_with_vhull_3 ${EXTRA_LIBS})
#  create_single_source_cgal_program( "Swept_volume_with_surface_mesh_3.cpp")
#  target_link_libraries (Swept_volume_with_surface_mesh_3 ${EXTRA_LIBS})
#  create_single_source_cgal_program( "Gen_triangle_soup.cpp")
#  target_link_libraries (Gen_triangle_soup ${EXTRA_LIBS})
endif()


#target_link_libraries ( Cuda_Sweep_Voxelizer_3 Cuda_Sweep_Voxelizer_3_kernels ${EXTRA_LIBS})

else()

	message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

