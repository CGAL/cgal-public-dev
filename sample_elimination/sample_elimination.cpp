#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
//#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>

 
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
#include <cmath>
#include <queue>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;
 
typedef CGAL::Point_set_3<Point>                              Point_set;

 
typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Traits> Tree;

struct Weighted {
    Point point;
    double weight;

    Weighted(Point point, double weight) : point(point), weight(weight) {}
   
};

struct CompareWeighted {
    bool operator()(const Weighted& a, const Weighted& b) {
        return a.weight > b.weight; // Min-heap based on age
    }
};

const double minDistance = .02;
 
double weight (Point p, std::vector<Point> neighbors)
{
  double weight = 0;
  for(Point n : neighbors)
  {
      weight = weight + pow((1-sqrt(CGAL::squared_distance(p,n)))/(minDistance),8);
  }
    
  return weight;
}


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
   
  std::vector<Point> points;
  PMP::sample_triangle_mesh(mesh,
                              std::back_inserter(points),
                              CGAL::parameters::number_of_points_per_area_unit(50000));
   
   
  //Point_set point_set;
  //PMP::sample_triangle_mesh(mesh,
                           //   point_set.point_back_inserter());
   
  std::cout << points.size() << std::endl;
    
  std::cout << mesh.number_of_faces() << std::endl;
  std::ofstream out("initial_sample.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();
    
    

 
  // Build tree in parallel
  Tree tree(points.begin(), points.end());
  tree.build<CGAL::Parallel_tag>();
    
  Point query = points.front();
  Fuzzy_circle default_range(query, .02);
    
  std::optional<Point> any = tree.search_any_point(default_range);
    if(any)
       std::cout << *any << " is in the query circle\n";
     else
       std::cout << "Empty query circle\n";
    
     std::vector<Point> result;
     tree.search(std::back_inserter(result), default_range);
    
     std::cout << "\nPoints in circle with center " << query << " and radius 0.02" << std::endl;
    
 // std::list<Point>::iterator it;
  for (size_t i = 0; i < result.size(); ++i) {
            std::cout << result[i] << "\n ";
        }
  
    
  std::ofstream out1("ball_o_points.xyz");
  out1 << std::setprecision(17);
  std::copy(result.begin(), result.end(), std::ostream_iterator<Point>(out1, "\n"));
  out1.close();
    
  std::vector<Point> neighbors;
  std::priority_queue<Weighted, std::vector<Weighted>, CompareWeighted> weightedHeap;
  for (Point p : points)
  {
      Fuzzy_circle default_range(p,.02);
      tree.search(std::back_inserter(neighbors), default_range);
      double w = weight(p,neighbors);
      weightedHeap.push({p,w});
      neighbors.clear();
      
  }
    
/*
 
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
 */
  return 0;
}
