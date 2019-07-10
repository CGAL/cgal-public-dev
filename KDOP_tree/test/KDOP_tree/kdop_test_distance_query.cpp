/*
 * kdop_test_distance_query.cpp
 *
 *  Created on: 9 Jul 2019
 *      Author: xx791
 */
int COUNTER_AABB = 0;
int COUNTER_KDOP = 0;
int COUNTER_TRIANGLES_AABB = 0;
int COUNTER_TRIANGLES_KDOP = 0;

//#define CHECK_CORRECTNESS

#define AABB_TIMING
#define KDOP_TIMING

#include <iostream>
#include <fstream>
#include <list>

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

#include <CGAL/Timer.h>

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

// AABB tree type definitions
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_aabb;
typedef CGAL::AABB_traits<K, Primitive_aabb> Traits_aabb;
typedef CGAL::AABB_tree<Traits_aabb> Tree_aabb;

// KDOP tree type definitions
const unsigned int NUM_DIRECTIONS = 6;

typedef CGAL::KDOP_tree::KDOP_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

typedef CGAL::Timer Timer;

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

  std::cout << "read points from file" << std::endl;
  std::vector<Point> points;
  read_points(pointsf, points);

  Timer t;

#ifdef AABB_TIMING
  //===========================================================================
  // AABB tree build
  //===========================================================================
  t.start();
  Tree_aabb tree_aabb( faces(mesh).first, faces(mesh).second, mesh );
  tree_aabb.build();
  t.stop();
  std::cout << "Build time AABB tree: " << t.time() << " sec."<< std::endl;
#endif

#ifdef KDOP_TIMING
  //===========================================================================
  // KDOP tree build
  //===========================================================================
  t.reset();
  t.start();
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );
  tree_kdop.build();
  t.stop();
  std::cout << "Build time " << NUM_DIRECTIONS << "-DOP tree: " << t.time() << " sec."<< std::endl << std::endl;
#endif

  // distance query
#ifdef CHECK_CORRECTNESS
  int num_error = 0;
  for (int i = 0; i < points.size(); ++i) {
    const Point& point = points[i];

    // AABB tree
    Point closest_point_aabb = tree_aabb.closest_point(point);

    // KDOP tree
    Point closest_point_kdop = tree_kdop.closest_point(point);

    bool is_same = K().equal_3_object()(closest_point_aabb, closest_point_kdop);

    if (is_same == false) {
      std::cout << "ERROR!" << std::endl;
      num_error += 1;
    }

    if (num_error == 0) {
      std::cout << "The closest_point query result of KDOP is the same as AABB." << std::endl;
    }
    else {
      std::cout << num_error << " differences for " << points.size() << " queries." << std::endl;
      return -1;
    }

  }
#endif

#ifdef AABB_TIMING
  t.reset();
  t.start();
  for (int i = 0; i < points.size(); ++i) {
    const Point& point = points[i];
    Point closest_point = tree_aabb.closest_point(point);
  }
  t.stop();
  std::cout << t.time() << " sec. for "   << points.size() << " closest_point queries with an AABB tree" << std::endl;
#endif

#ifdef KDOP_TIMING
  t.reset();
  t.start();
  for (int i = 0; i < points.size(); ++i) {
    const Point& point = points[i];
    Point closest_point = tree_kdop.closest_point(point);
  }
  t.stop();
  std::cout << t.time() << " sec. for "  << points.size() << " closest_point queries with a " << NUM_DIRECTIONS << "-DOP tree" << std::endl;
#endif

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
