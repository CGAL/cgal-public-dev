#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Polyhedron_3.h>


#include <cassert>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::Embree::Triangle_mesh_geometry<Polyhedron, K> TriangleMesh;
typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(a, b, c, d);

    Tree tree;

    // Configure the ray.
    Point rayOrigin(0.1f, 0.2f, -1.0f);
    Vector rayDirection(0.0f, 0.0f, 1.0f); /*Direction need not be normalized.*/
    Ray ray(rayOrigin, rayDirection);


    // Test calls to all functions but those who have `!empty()` as
    // precondition.
    CGAL::Emptyset_iterator devnull;
    tree.all_intersections(ray, devnull);
    tree.all_intersected_primitives(ray, devnull);
    assert(!tree.any_intersected_primitive(ray));
    assert(!tree.any_intersection(ray));
    //Cannot call tree.bbox();
    // tree.build();
    tree.clear();
    //Cannot call tree.closest_*(...)
    assert(tree.do_intersect(ray) == false);
    assert(tree.empty());
    //Do not call tree.insert(...)
    assert(tree.number_of_intersected_primitives(ray) == 0);
    // Cannot call tree.rebuild(..)
    assert(tree.size() == 0);
    // Cannot call tree.squared_distance(..)

    return EXIT_SUCCESS;
}
