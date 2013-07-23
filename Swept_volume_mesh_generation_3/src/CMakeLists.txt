cmake_minimum_required(VERSION 2.6)

include_directories(BEFORE ../include/)
include_directories(BEFORE ../test/Swept_volume/include/)



# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4 ImageIO)
include( ${CGAL_USE_FILE} )


add_library(swept_volume_tools SHARED Triangle_voxelizer.cpp Voxelization_3.cpp)

