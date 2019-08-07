/*
 * An example of ray query to compute first intersections.
 *
 */

#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef K::Segment_3 Segment;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

// prescribe the number of directions in the k-dop
const unsigned int NUM_DIRECTIONS = 14;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive_kdop;
typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive_kdop> Traits_kdop;
typedef CGAL::KDOP_tree::KDOP_tree<Traits_kdop> Tree_kdop;

typedef boost::optional< Tree_kdop::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;

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

  const char* points_file = argv[2];
  std::ifstream pointsf(points_file);

  // read points from the file
  std::vector<Point> points;
  read_points(pointsf, points);

  // create a set of random rays, shooting from points read from the file to
  // centroids of primitives
  std::vector< Ray > rays;

  for(face_descriptor fd : faces(mesh)){
    halfedge_descriptor hd = halfedge(fd,mesh);
    Point p = CGAL::centroid(mesh.point(source(hd,mesh)),
        mesh.point(target(hd,mesh)),
        mesh.point(target(next(hd,mesh),mesh)));

    for (int i = 0; i < points.size(); ++i) {
      Ray ray(points[i], p);
      rays.push_back(ray);
    }

  }

  // read the mesh into the k-dop tree
  Tree_kdop tree_kdop( faces(mesh).first, faces(mesh).second, mesh );

  // build the tree, including splitting primitives and computing k-dops
  tree_kdop.build();

  // ray query to get first intersections
  for (int i = 0; i < rays.size(); ++i) {
    const Ray& ray_query = rays[i];

    Ray_intersection intersection = tree_kdop.first_intersection(ray_query);

    if (intersection) {
      const Point* p_kdop = boost::get<Point>( &(intersection->first) );
      std::cout << *p_kdop << std::endl;
    }
  }

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
