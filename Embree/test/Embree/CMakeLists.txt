cmake_minimum_required(VERSION 3.1...3.14)
project( AABB_tree_benchmark )

find_package(CGAL  QUIET )
FIND_PACKAGE(embree 3.0 QUIET)

include_directories(../../include)


if ( CGAL_FOUND AND embree_FOUND )

    ADD_EXECUTABLE(aabb_any_all_benchmark aabb_any_all_benchmark.cpp)
    TARGET_LINK_LIBRARIES(aabb_any_all_benchmark CGAL embree)

    # ADD_EXECUTABLE(closest_point closest_point.cpp)
    # TARGET_LINK_LIBRARIES(closest_point CGAL embree)


else()

    message(STATUS "This program requires the CGAL and the Embree library, and will not be compiled.")

endif()
