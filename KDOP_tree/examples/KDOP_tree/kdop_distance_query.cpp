/*
 * An example of distance query to compute closest points.
 *
 */

#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>

// prescribed number of directions for the k-dop
const unsigned int NUM_DIRECTIONS = 14;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

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

  const char* pointsFile = argv[2];
  std::ifstream pointsf(pointsFile);

  // point query
  Point point(1./3., 2./3., 1./6.);

  // read the mesh into the k-dop tree
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );

  // build the k-dop tree, including splitting primitives and computing k-dops
  tree_kdop.build();

  // distance query to get closest points
  Point closest_point = tree_kdop.closest_point(point);
  std::cout << closest_point << std::endl;

  return 0;
}
