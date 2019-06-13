/*
 * kdop_test_kdop_tree.cpp
 *
 *  Created on: 12 Jun 2019
 *      Author: xx791
 */

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
typedef CGAL::KDOP_tree::KDOP_traits<K, Primitive> Traits;
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
  std::vector< std::vector<double> > kdop_directions;

  std::vector<double> kdop_direction(3);

  kdop_direction[0] = 1., kdop_direction[1] = 1., kdop_direction[2] = 1.;
  kdop_directions.push_back(kdop_direction);

  kdop_direction[0] = 1., kdop_direction[1] = 2., kdop_direction[2] = 3.;
  kdop_directions.push_back(kdop_direction);

  // input k-dop directions to the tree
  tree.set_kdop_directions(kdop_directions);

  // build the tree, including splitting primitives and computing k-dops
  tree.build();


  return 0;
}
