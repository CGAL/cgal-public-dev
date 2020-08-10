cmake_minimum_required(VERSION 3.1...3.14)
project( Embree_test )

find_package(CGAL  QUIET )
FIND_PACKAGE(embree 3.0 QUIET)

include_directories(../../include)


if ( CGAL_FOUND AND embree_FOUND )

    ADD_EXECUTABLE(aabb_any_all_benchmark aabb_any_all_benchmark.cpp)
    TARGET_LINK_LIBRARIES(aabb_any_all_benchmark CGAL::CGAL embree)

    ADD_EXECUTABLE(aabb_correctness_triangle_test aabb_correctness_triangle_test.cpp)
    TARGET_LINK_LIBRARIES(aabb_correctness_triangle_test CGAL::CGAL embree)

    ADD_EXECUTABLE(aabb_test_empty_tree aabb_test_empty_tree.cpp)
    TARGET_LINK_LIBRARIES(aabb_test_empty_tree CGAL::CGAL embree)

    ADD_EXECUTABLE(aabb_test_ray_intersection aabb_test_ray_intersection.cpp)
    TARGET_LINK_LIBRARIES(aabb_test_ray_intersection CGAL::CGAL embree)
    
    ADD_EXECUTABLE(aabb_distance_edge_test aabb_distance_edge_test.cpp)
    TARGET_LINK_LIBRARIES(aabb_distance_edge_test CGAL::CGAL embree)

    ADD_EXECUTABLE(aabb_intersection_triangle_test aabb_intersection_triangle_test.cpp)
    TARGET_LINK_LIBRARIES(aabb_intersection_triangle_test CGAL::CGAL embree)    

else()

    message(STATUS "This program requires the CGAL and the Embree library, and will not be compiled.")

endif()
