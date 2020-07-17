cmake_minimum_required(VERSION 3.1...3.14)
project( AABB_tree_benchmark )

find_package(CGAL  QUIET )
FIND_PACKAGE(embree 3.0 QUIET)

include_directories(../../include)


if ( CGAL_FOUND AND embree_FOUND )

    ADD_EXECUTABLE(EmbreeTutorial EmbreeTutorial.cpp)
    TARGET_LINK_LIBRARIES(EmbreeTutorial CGAL embree)

    ADD_EXECUTABLE(TestEmbreeAABB TestEmbreeAABB.cpp)
    TARGET_LINK_LIBRARIES(TestEmbreeAABB CGAL embree)

    ADD_EXECUTABLE(MultipleRayEmbreeAabb MultipleRayEmbreeAabb.cpp)
    TARGET_LINK_LIBRARIES(MultipleRayEmbreeAabb CGAL embree)

    ADD_EXECUTABLE(PolyhedronTest PolyhedronTest.cpp)
    TARGET_LINK_LIBRARIES(PolyhedronTest CGAL embree)

else()

    message(STATUS "This program requires the CGAL and the Embree library, and will not be compiled.")

endif()
