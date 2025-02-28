#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
 
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
 



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
 
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Point_set_3.h>
#include <iostream>
#include <string>
#include <vector>
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;
 
typedef CGAL::Point_set_3<Point>                              Point_set;

 
typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

using Traits = CGAL::Search_traits_3<K>;
using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Tree = Neighbor_search::Tree;
using Point_with_distance = Neighbor_search::Point_with_transformed_distance;
 
using Generator = CGAL::Random_points_in_sphere_3<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;
 
int main(int argc, char* argv[])
{
    
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("data/meshes/eight.off");
   
  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
     std::cerr << "Invalid input." << std::endl;
     return 1;
  }
   
  const int points_per_face = (argc > 2) ? std::stoi(argv[2]) : 30;
   
  std::vector<Point> pointz;
  PMP::sample_triangle_mesh(mesh,
                              std::back_inserter(pointz),
                              CGAL::parameters::number_of_points_per_area_unit(50000));
   
   
  Point_set point_set;
  PMP::sample_triangle_mesh(mesh,
                              point_set.point_back_inserter());
   
  std::cout << pointz.size() << std::endl;
    
  std::cout << mesh.number_of_faces() << std::endl;
  std::ofstream out("initial_sample.xyz");
  out << std::setprecision(17);
  std::copy(pointz.begin(), pointz.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();
    
    
  const unsigned int N = 1000;
  const unsigned int k = 6;
 
  // Generate N points in a sphere
  std::vector<Point> points;
  points.reserve (N);
  Generator generator;
  for (unsigned int i = 0; i < N; ++ i)
    points.push_back (*(generator++));
 
  // Build tree in parallel
  Tree tree(points.begin(), points.end());
  tree.build<CGAL::Parallel_tag>();
 
  // Query tree in parallel
  std::vector<std::vector<Point> > neighbors (points.size());
  tbb::parallel_for (tbb::blocked_range<std::size_t> (0, points.size()),
                     [&](const tbb::blocked_range<std::size_t>& r)
                     {
                       for (std::size_t s = r.begin(); s != r.end(); ++ s)
                       {
                         // Neighbor search can be instantiated from
                         // several threads at the same time
                         Neighbor_search search (tree, points[s], k);
                         neighbors[s].reserve(k);
 
                         // neighbor search returns a set of pair of
                         // point and distance <Point_3,FT>, here we
                         // keep the points only
                         for (const Point_with_distance pwd : search)
                           neighbors[s].push_back (pwd.first);
                       }
                     });
 
  return 0;
}
