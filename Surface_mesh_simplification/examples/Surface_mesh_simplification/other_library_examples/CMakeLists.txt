project ( other_library_examples )

cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

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

find_package(Threads)

add_executable(
  libigl_simplification
  libigl_simplification.cpp
)

target_link_libraries(
  libigl_simplification
  ${CMAKE_THREAD_LIBS_INIT}
)

find_package(pmp)

add_executable(
  pmp_simplification
  pmp_simplification.cpp
)

target_link_libraries(
  pmp_simplification
  pmp
)
