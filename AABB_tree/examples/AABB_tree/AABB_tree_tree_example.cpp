// Author(s) : Camille Wormser, Pierre Alliez

#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;
typedef boost::optional< Tree::Intersection_and_primitive_id<Plane>::Type > Plane_intersection;
typedef Tree::Primitive_id Primitive_id;

int main()
{
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(p, q, r, s);

    Point p2(2.0, 0.0, 0.0);
    Point q2(0.0, 2.0, 0.0);
    Point r2(0.0, 0.0, 2.0);
    Point s2(0.1, 0.1, 0.1);
    Polyhedron polyhedron2;
    polyhedron2.make_tetrahedron(p2, q2, r2, s2);

    // constructs AABB trees
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end(),polyhedron);
    Tree tree2(polyhedron2.facets_begin(),polyhedron2.facets_end(),polyhedron2);

    std::cout << tree.number_of_intersected_primitives(tree2) << " intersections" << std::endl;

    return EXIT_SUCCESS;
}
