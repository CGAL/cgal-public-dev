/*
 * kdop_test_aabb_kdop_tree.cpp
 *
 *  Created on: 20 Jun 2019
 *      Author: xx791
 */

#include <iostream>
#include <fstream>
#include <list>

#include <boost/timer.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

// AABB tree includes
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// KDOP tree includes
#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>
#include <CGAL/KDOP_tree/KDOP_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef K::Segment_3 Segment;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

// AABB tree type definitions
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_aabb;
typedef CGAL::AABB_traits<K, Primitive_aabb> Traits_aabb;
typedef CGAL::AABB_tree<Traits_aabb> Tree_aabb;

// KDOP tree type definitions
const unsigned int NUM_DIRECTIONS = 14;

typedef CGAL::KDOP_tree::KDOP_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tetrahedron.off";

  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;

  // create a set of random rays, centred at p0.
  std::vector< Ray > rays;

  const Point p0(0., 50., 10.); // the centre of the rays
  const double radius = 20.; // the radius of the rays

  const int num_alpha = 10;
  const int num_beta = 10;
  for (int i = 0; i < num_alpha; ++i) {
    double alpha = i*(2.*M_PI/num_alpha);
    for (int j = 0; j < num_beta; ++j) {
      double beta = -M_PI/2. + j*(M_PI/num_beta);

      double x = p0.x() + radius*std::cos(beta)*std::cos(alpha);
      double y = p0.y() + radius*std::cos(beta)*std::sin(alpha);
      double z = p0.z() + radius*std::sin(beta);

      const Point p(x, y, z);
      Ray ray(p0, p);
      rays.push_back(ray);
    }
  }

#ifdef WRITE_FILE
  // write rays to file
  std::string rayFile("bunny_ray_test.obj");
  std::ofstream rayf(rayFile.c_str());

  for (int i = 0; i < rays.size(); ++i) {
    Ray ray = rays[i];

    Point source = ray.source();
    Point target = ray.second_point();

    rayf << "v " << source.x() << " " << source.y() << " " << source.z() << std::endl;
    rayf << "v " << target.x() << " " << target.y() << " " << target.z() << std::endl;
  }

  for (int i = 0; i < rays.size(); ++i) {
    rayf << "l " << 2*i + 1 << " " << 2*i + 2 << std::endl;
  }
#endif

  //===========================================================================
  // AABB tree build
  //===========================================================================
  Tree_aabb tree_aabb( faces(mesh).first, faces(mesh).second, mesh );

  tree_aabb.build();

  //===========================================================================
  // KDOP tree build
  //===========================================================================
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );

  // user-defined directions for k-dops
  // (number of directions = NUM_DIRECTIONS)
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

  kdop_directions.push_back(Point(1., 1., 1.));
  kdop_directions.push_back(Point(-1., -1., -1.));

  kdop_directions.push_back(Point(-1., 1., 1.));
  kdop_directions.push_back(Point(1., -1., -1.));

  kdop_directions.push_back(Point(-1., -1., 1.));
  kdop_directions.push_back(Point(1., 1., -1.));

  kdop_directions.push_back(Point(1., -1., 1.));
  kdop_directions.push_back(Point(-1., 1., -1.));

  // input k-dop directions to the tree
  tree_kdop.set_kdop_directions(kdop_directions);

  // build the tree, including splitting primitives and computing k-dops
  tree_kdop.build();

  //===========================================================================
  // Ray intersection check using AABB tree and KDOP tree
  //===========================================================================
  int num_error = 0;
  for (int i = 0; i < rays.size(); ++i) {
    Ray ray_query = rays[i];

    // AABB tree
    bool is_intersect_aabb = tree_aabb.do_intersect(ray_query);

    // KDOP tree
    bool is_intersect_kdop = tree_kdop.do_intersect(ray_query);

    if (is_intersect_aabb != is_intersect_kdop) {
      std::cout << "ERROR!" << std::endl;
      num_error += 1;
      break;
    }

  }

  if (num_error == 0) std::cout << "The do_intersect result of KDOP is the same as AABB." << std::endl;

  return 0;
}


