/*
 * kdop_test_do_intersect.cpp
 *
 *  Created on: 18 Jun 2019
 *      Author: xx791
 */

#include <iostream>
#include <fstream>
#include <list>

#include <boost/timer.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/KDOP_tree/KDOP_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef K::Segment_3 Segment;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef CGAL::KDOP_tree::KDOP_face_graph_triangle_primitive<Mesh> Primitive;

const unsigned int NUM_DIRECTION = 6;

typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTION, K, Primitive> Traits;
typedef CGAL::KDOP_tree::KDOP_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tetrahedron.off";

  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;

  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  // user-defined directions for k-dops
  std::vector< Point > kdop_directions;

  for (int i = 0; i < 3; ++i) {
    std::vector<double> direction(3);
    direction[0] = 0., direction[1] = 0., direction[2] = 0.;

    direction[i] = 1.;

    Point direction1(direction[0], direction[1], direction[2]);
    kdop_directions.push_back(direction1);

    direction[i] = -1.;

    Point direction2(direction[0], direction[1], direction[2]);
    kdop_directions.push_back(direction2);
  }

  // input k-dop directions to the tree
  tree.set_kdop_directions(kdop_directions);

  // build the tree, including splitting primitives and computing k-dops
  tree.build();

  // ray intersection
  // \todo implemented as a line segment intersection at the moment, need to generalise it.
  Point p1(1., 0., 0.);
  Point p2(1., 1., 1.);

  Ray ray_query(p1, p2);

  bool is_intersect = tree.do_intersect(ray_query);

  if (is_intersect == true) std::cout << "intersected" << std::endl;
  else std::cout << "not intersected" << std::endl;

  return 0;
}


