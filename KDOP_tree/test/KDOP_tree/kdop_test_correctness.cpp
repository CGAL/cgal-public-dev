/*
 * kdop_test_correctness.cpp
 * Check the correctness of k-DOP tree by comparing the ray and distance queries
 * with AABB tree.
 */

#define CHECK_RAY_QUERY
#ifdef CHECK_RAY_QUERY
//#define CHECK_DO_INTERSECT
//#define CHECK_NUMBER_OF_INTERSECTED_PRIMITIVES
//#define CHECK_ALL_INTERSECTED_PRIMITIVES
//#define CHECK_ANY_INTERSECTED_PRIMITIVE
//#define CHECK_ALL_INTERSECTIONS
//#define CHECK_ANY_INTERSECTION
#define CHECK_FIRST_INTERSECTION
//#define CHECK_FIRST_INTERSECTED_PRIMITIVE
#endif

//#define CHECK_DISTANCE_QUERY
#ifdef CHECK_DISTANCE_QUERY
//#define CHECK_CLOSEST_POINT
//#define CHECK_SQUARED_DISTANCE
#define CHECK_CLOSEST_POINT_AND_PRIMITIVE
#endif

//#define WRITE_FILE

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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

typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
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

typedef typename Tree_aabb::Primitive_id Primitive_id_aabb;

typedef boost::optional< Tree_aabb::Intersection_and_primitive_id<Ray>::Type > Ray_intersection_aabb;

// KDOP tree type definitions
const unsigned int NUM_DIRECTIONS = 14;

//typedef CGAL::KDOP_tree::KDOP_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

typedef typename Tree_kdop::Primitive_id Primitive_id_kdop;

typedef boost::optional< Tree_kdop::Intersection_and_primitive_id<Ray>::Type > Ray_intersection_kdop;

void read_points(std::ifstream& pointsf, std::vector<Point>& points);

struct Skip {
  face_descriptor fd;

  Skip(const face_descriptor fd)
    : fd(fd)
  {}

  bool operator()(const face_descriptor& t) const
  { if(t == fd){
      std::cerr << "ignore " << t  <<std::endl;
    };
    return(t == fd);
  }

};

int main(int argc, char* argv[])
{
#ifdef CHECK_RAY_QUERY
  if (argc != 2) {
    std::cerr << "Command ERROR! Usage: kdop_test_correctness mesh.off" << std::endl;
    return 0;
  }
#endif
#ifdef CHECK_DISTANCE_QUERY
  if (argc != 3) {
    std::cerr << "Command ERROR! Usage: kdop_test_correctness mesh.off points.xyz" << std::endl;
    return 0;
  }
#endif
  //---------------------------------------------------------------------------
  const char* filename = argv[1];
  std::ifstream input(filename);

  if (!input.is_open()) {
    std::cerr << "Mesh file could not be found!" << std::endl;
    return -1;
  }

  Mesh mesh;
  input >> mesh;

#ifdef CHECK_RAY_QUERY
  // create rays shooting from centroids of facets of the mesh, normal to the facets
  std::vector< Ray > rays;

  double d = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)?-1:1;

  for(face_descriptor fd : faces(mesh)){
    halfedge_descriptor hd = halfedge(fd,mesh);
    Point p = CGAL::centroid(mesh.point(source(hd,mesh)),
        mesh.point(target(hd,mesh)),
        mesh.point(target(next(hd,mesh),mesh)));

    Vector v = CGAL::Polygon_mesh_processing::compute_face_normal(fd,mesh);
    Ray ray(p, d*v);

    rays.push_back(ray);
  }
#endif

#ifdef CHECK_DISTANCE_QUERY
  const char* pointsFile = argv[2];
  std::ifstream pointsf(pointsFile);

  if (!pointsf.is_open()) {
    std::cerr << "Points file could not be found!" << std::endl;
    return -1;
  }

  std::cout << "read points from file" << std::endl;
  std::vector<Point> points;
  read_points(pointsf, points);
#endif

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

  tree_aabb.build();

  //===========================================================================
  // KDOP tree build
  //===========================================================================
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );

  // build the tree, including splitting primitives and computing k-dops
  tree_kdop.build();

  //===========================================================================
  // Ray intersection check using AABB tree and KDOP tree
  //===========================================================================
  int num_error = 0;

#ifdef CHECK_RAY_QUERY
  for (int i = 0; i < rays.size(); ++i) {
    const Ray& ray_query = rays[i];

#ifdef CHECK_DO_INTERSECT
    bool is_intersect_aabb = tree_aabb.do_intersect(ray_query);
    bool is_intersect_kdop = tree_kdop.do_intersect(ray_query);
    if (is_intersect_aabb != is_intersect_kdop) {
      std::cout << "ERROR: do_intersect!" << std::endl;
      num_error += 1;
    }
#endif
#ifdef CHECK_NUMBER_OF_INTERSECTED_PRIMITIVES
    int num_intersected_aabb = tree_aabb.number_of_intersected_primitives(ray_query);
    int num_intersected_kdop = tree_kdop.number_of_intersected_primitives(ray_query);
    if (num_intersected_aabb != num_intersected_kdop) {
      std::cout << "ERROR: number_of_intersected_primitives!" << std::endl;
      num_error += 1;
    }
#endif
#ifdef CHECK_ALL_INTERSECTED_PRIMITIVES
    std::list<Primitive_id_aabb> primitives_aabb;
    tree_aabb.all_intersected_primitives(ray_query, std::back_inserter(primitives_aabb));

    std::list<Primitive_id_kdop> primitives_kdop;
    tree_kdop.all_intersected_primitives(ray_query, std::back_inserter(primitives_kdop));

    if (primitives_aabb.size() > 0 && primitives_kdop.size() > 0) {
      if (primitives_aabb.size() == primitives_kdop.size()) {
        std::list<Primitive_id_aabb>::iterator it_aabb = primitives_aabb.begin();
        std::list<Primitive_id_kdop>::iterator it_kdop = primitives_kdop.begin();

        for (int ii = 0; ii < primitives_aabb.size(); ++ii) {
          if (*it_aabb != *it_kdop) {
            std::cout << "ERROR: all_intersected_primitives!" << std::endl;
            num_error += 1;
          }
          else {
            it_aabb++;
            it_kdop++;
          }
        }
      }
      else {
        std::cout << "ERROR: all_intersected_primitives!" << std::endl;
      }
    }
#endif
#ifdef CHECK_ANY_INTERSECTED_PRIMITIVE
    boost::optional<Primitive_id_aabb> primitive_id_aabb = tree_aabb.any_intersected_primitive(ray_query);
    boost::optional<Primitive_id_kdop> primitive_id_kdop = tree_kdop.any_intersected_primitive(ray_query);
    if (*primitive_id_aabb != *primitive_id_kdop) {
      std::cout << "ERROR: any_intersected_primitive!" << std::endl;
      num_error += 1;
    }
#endif
#ifdef CHECK_ALL_INTERSECTIONS
    std::list<Ray_intersection_aabb> intersections_aabb;
    tree_aabb.all_intersections(ray_query, std::back_inserter(intersections_aabb));

    std::list<Ray_intersection_kdop> intersections_kdop;
    tree_kdop.all_intersections(ray_query, std::back_inserter(intersections_kdop));

    if (intersections_aabb.size() > 0 && intersections_kdop.size() > 0) {
      if (intersections_aabb.size() == intersections_kdop.size()) {
        std::list<Ray_intersection_aabb>::iterator it_aabb = intersections_aabb.begin();
        std::list<Ray_intersection_kdop>::iterator it_kdop = intersections_kdop.begin();

        for (int ii = 0; ii < intersections_aabb.size(); ++ii) {
          if (*it_aabb != *it_kdop) {
            std::cout << "ERROR: all_intersections!" << std::endl;
            num_error += 1;
          }
          else {
            it_aabb++;
            it_kdop++;
          }
        }
      }
      else {
        std::cout << "ERROR: all_intersections!" << std::endl;
      }
    }
#endif
#ifdef CHECK_ANY_INTERSECTION
    Ray_intersection_aabb intersection_aabb = tree_aabb.any_intersection(ray_query);
    Ray_intersection_kdop intersection_kdop = tree_kdop.any_intersection(ray_query);

    if ( intersection_aabb || intersection_kdop ) {
      const Point* p_aabb = boost::get<Point>( &(intersection_aabb->first) );
      const Point* p_kdop = boost::get<Point>( &(intersection_kdop->first) );

      bool is_same = K().equal_3_object()(*p_aabb, *p_kdop);

      if (is_same == false) {
        std::cout << "ERROR: any_intersection!" << std::endl;
        num_error += 1;
      }
    }
#endif
#ifdef CHECK_FIRST_INTERSECTION
    Ray_intersection_aabb intersection_aabb = tree_aabb.first_intersection(ray_query);
    Ray_intersection_aabb intersection_kdop = tree_kdop.first_intersection(ray_query);

    if ( intersection_aabb || intersection_kdop ) {
      const Point* p_aabb = boost::get<Point>( &(intersection_aabb->first) );
      const Point* p_kdop = boost::get<Point>( &(intersection_kdop->first) );

      bool is_same = K().equal_3_object()(*p_aabb, *p_kdop);

      if (is_same == false) {
        std::cout << "ERROR: first_intersection!" << std::endl;
        num_error += 1;
      }
    }
#endif
#ifdef CHECK_FIRST_INTERSECTED_PRIMITIVE
    boost::optional<Primitive_id_aabb> primitive_id_aabb = tree_aabb.first_intersected_primitive(ray_query);
    boost::optional<Primitive_id_kdop> primitive_id_kdop = tree_kdop.first_intersected_primitive(ray_query);

    if ( primitive_id_aabb || primitive_id_kdop ) {
      if (*primitive_id_aabb != *primitive_id_kdop) {
        std::cout << "ERROR: first_intersected_primitive!" << std::endl;
        num_error += 1;
      }
    }
#endif
  }

  if (num_error == 0){
    std::cout << "The result of KDOP is the same as AABB." << std::endl;
  } else {
    std::cout << num_error << " differences for " << rays.size() << " queries" << std::endl;
    return -1;
  }
#endif

#ifdef CHECK_DISTANCE_QUERY
  for (int i = 0; i < points.size(); ++i) {
    Point point = points[i];

#ifdef CHECK_CLOSEST_POINT
    Point closest_point_aabb = tree_aabb.closest_point(point);
    Point closest_point_kdop = tree_kdop.closest_point(point);

    bool is_same = K().equal_3_object()(closest_point_aabb, closest_point_kdop);

    if (is_same == false) {
      std::cout << "ERROR: closest_point!" << std::endl;
      num_error += 1;
    }
#endif
#ifdef CHECK_SQUARED_DISTANCE
    double sq_dist_aabb = tree_aabb.squared_distance(point);
    double sq_dist_kdop = tree_kdop.squared_distance(point);

    if (std::abs(sq_dist_aabb - sq_dist_kdop) > 1.e-12) {
      std::cout << "ERROR: squared_distance!" << std::endl;
      num_error += 1;
    }
#endif
#ifdef CHECK_CLOSEST_POINT_AND_PRIMITIVE
    typename Tree_aabb::Point_and_primitive_id p_pr_aabb = tree_aabb.closest_point_and_primitive(point);
    typename Tree_kdop::Point_and_primitive_id p_pr_kdop = tree_kdop.closest_point_and_primitive(point);

    bool is_same = K().equal_3_object()(p_pr_aabb.first, p_pr_kdop.first);

    if ( (is_same == false) || (p_pr_aabb.second != p_pr_kdop.second) ) {
      std::cout << "ERROR: closest_point_and_primitive!" << std::endl;
      num_error += 1;
    }
#endif
  }

  if (num_error == 0){
    std::cout << "The result of KDOP is the same as AABB." << std::endl;
  } else {
    std::cout << num_error << " differences for " << rays.size() << " queries" << std::endl;
    return -1;
  }
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

