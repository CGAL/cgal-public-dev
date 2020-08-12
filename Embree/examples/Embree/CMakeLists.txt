cmake_minimum_required(VERSION 3.1...3.14)
project( Embree_examples )

find_package(CGAL  QUIET )
FIND_PACKAGE(embree 3.0 QUIET)

include_directories(../../include)


if ( CGAL_FOUND AND embree_FOUND )

    ADD_EXECUTABLE(surface_mesh_embree surface_mesh_embree.cpp)
    TARGET_LINK_LIBRARIES(surface_mesh_embree CGAL::CGAL embree)
    
    ADD_EXECUTABLE(polyhedron_embree polyhedron_embree.cpp)
    TARGET_LINK_LIBRARIES(polyhedron_embree CGAL::CGAL embree)

    ADD_EXECUTABLE(closest_point closest_point.cpp)
    TARGET_LINK_LIBRARIES(closest_point CGAL::CGAL embree)

    ADD_EXECUTABLE(ray_shooting_embree ray_shooting_embree.cpp)
    TARGET_LINK_LIBRARIES(ray_shooting_embree CGAL::CGAL embree)

    ADD_EXECUTABLE(segment_shooting_embree segment_shooting_embree.cpp)
    TARGET_LINK_LIBRARIES(segment_shooting_embree CGAL::CGAL embree)

else()

    message(STATUS "This program requires the CGAL and the Embree library, and will not be compiled.")

endif()
