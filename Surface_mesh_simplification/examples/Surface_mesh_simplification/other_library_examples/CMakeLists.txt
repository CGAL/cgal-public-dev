project ( other_library_examples )

cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

add_executable(
openmesh_trimesh_simplification
openmesh_trimesh_simplification.cpp
)

target_link_libraries(
openmesh_trimesh_simplification
OpenMeshCore
OpenMeshTools
)


add_executable(
openmesh_polymesh_simplification
openmesh_polymesh_simplification.cpp
)

target_link_libraries(
openmesh_polymesh_simplification
OpenMeshCore
OpenMeshTools
)
