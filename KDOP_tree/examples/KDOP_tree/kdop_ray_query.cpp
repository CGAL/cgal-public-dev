/*
 * An example of ray query to compute first intersections.
 *
 */

#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>

// prescribe the number of directions in the k-dop
const unsigned int NUM_DIRECTIONS = 14;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Ray_3 Ray;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

typedef boost::optional< Tree_kdop::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Need mesh file!" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;

  // create a ray
  Point p1(1., 1., 1.);
  Point p2(1./3., 1./6., 2./3.);

  Ray ray(p1, p2);

  // read the mesh into the k-dop tree
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );

  // build the tree, including splitting primitives
  // and computing k-dops with pre-defined directions
  tree_kdop.build();

  // ray query to get first intersections
  Ray_intersection intersection = tree_kdop.first_intersection(ray);

  if (intersection) {
    const Point* p_kdop = boost::get<Point>( &(intersection->first) );
    std::cout << *p_kdop << std::endl;
  }
  else {
    std::cout << "No intersection" << std::endl;
  }

  return 0;
}
