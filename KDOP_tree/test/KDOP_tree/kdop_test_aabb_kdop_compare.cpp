/*
 * kdop_test_aabb_kdop_compare.cpp
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
//#define DEBUG_
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

void read_points(std::ifstream& pointsf, std::vector<Point>& points);

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Need mesh file and points file!" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;

  const char* pointsFile = argv[2];
  std::ifstream pointsf(pointsFile);

  std::vector<Point> points;
  read_points(pointsf, points);

  // create a set of random rays, centred at points read from the file.
  std::vector< Ray > rays;

  const double radius = 20.; // the radius of the rays
  const int num_alpha = 10;
  const int num_beta = 10;

  for (int i = 0; i < points.size(); ++i) {
    Point p0 = points[i];

    for (int ii = 0; ii < num_alpha; ++ii) {
      double alpha = ii*(2.*CGAL_PI/num_alpha);
      for (int jj = 0; jj < num_beta; ++jj) {
        double beta = -CGAL_PI/2. + jj*(CGAL_PI/num_beta);

        double x = p0.x() + radius*std::cos(beta)*std::cos(alpha);
        double y = p0.y() + radius*std::cos(beta)*std::sin(alpha);
        double z = p0.z() + radius*std::sin(beta);

        const Point p(x, y, z);
        Ray ray(p0, p);
        rays.push_back(ray);
      }
    }
  }

#ifdef WRITE_FILE
  // write rays to file
  std::string rayFile("ray_file_test.obj");
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

  boost::timer timer_aabb;

  // build AABB tree
  tree_aabb.build();

  // check intersections with rays
  for (int i = 0; i < rays.size(); ++i) {
    Ray ray_query = rays[i];
    bool is_intersect_aabb = tree_aabb.do_intersect(ray_query);
  }

  std::cout << timer_aabb.elapsed() << std::endl;

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

  boost::timer timer_kdop;

  // build KDOP tree, including splitting primitives and computing k-dops
  tree_kdop.build();

  // check intersections with rays
  for (int i = 0; i < rays.size(); ++i) {
    Ray ray_query = rays[i];
    bool is_intersect_kdop = tree_kdop.do_intersect(ray_query);
  }

  std::cout << timer_kdop.elapsed() << std::endl;

  return 0;
}

void read_points(std::ifstream& pointsf, std::vector<Point>& points)
{
  std::string line;
  while ( std::getline(pointsf, line) ) {
    std::stringstream line_stream;
    line_stream.str(line);
    double x, y, z;
    line_stream >> x >> y >> z;

    Point p(x, y, z);

    points.push_back(p);
  }
}
