#include <CGAL/Simple_cartesian.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

#include "include/Saver.h"

using K = CGAL::Simple_cartesian<double>;
using Saver = CGAL::Levels_of_detail::Saver<K>;
using Color = CGAL::Color;

typedef K::FT               FT;
typedef K::Weighted_point_2 Weighted_point;
typedef K::Segment_2        Segment;
typedef K::Point_2          Point_2;
typedef K::Point_3          Point_3;

typedef CGAL::Regular_triangulation_vertex_base_2<K> Rvb;
typedef CGAL::Alpha_shape_vertex_base_2<K,Rvb>       Vb;
typedef CGAL::Regular_triangulation_face_base_2<K>   Rf;
typedef CGAL::Alpha_shape_face_base_2<K,Rf>          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>  Tds;

typedef CGAL::Regular_triangulation_2<K,Tds>      Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>      Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

using Pair = std::pair<Point_2, Color>;
using Point_map_2 = CGAL::First_of_pair_property_map<Pair>;
using K_neighbor_query =
CGAL::Levels_of_detail::internal::K_neighbor_query<K, std::vector<Pair>, Point_map_2>;

int main() {
  
  std::vector<Weighted_point> wpoints;
  const std::string testpath = "/Users/monet/Documents/lod/logs/input.xyz";
  std::ifstream file(testpath.c_str(), std::ifstream::in);
  file.precision(15);

  std::vector<Pair> pairs;
  Point_3 p; Color c;
  while (!file.eof()) {
    file >> p >> c;
    pairs.push_back(std::make_pair(Point_2(p.x(), p.y()), c));
    wpoints.push_back(Weighted_point(Point_2(p.x(), p.y()), FT(0)));
  }
  pairs.erase(pairs.begin() + pairs.size() - 1);
  wpoints.erase(wpoints.begin() + wpoints.size() - 1);

  std::vector<Point_2> tmp;
  Point_map_2 point_map_2;
  K_neighbor_query knq(pairs, FT(6), point_map_2);
  std::vector<std::size_t> neighbors;
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const auto& p = pairs[i];
    neighbors.clear();
    knq(p.first, neighbors);
    for (const std::size_t idx : neighbors) {
      if (int(p.second.red()) != int(pairs[idx].second.red())) {
        tmp.push_back(p.first);
        wpoints[i] = Weighted_point(p.first, FT(0));
        break;
      }
    }
  }

  Alpha_shape_2 alpha_shape(
    wpoints.begin(), wpoints.end(), 0.5, Alpha_shape_2::GENERAL);

  std::vector<Point_2> res;
  for (
    auto it = alpha_shape.alpha_shape_vertices_begin(); 
    it != alpha_shape.alpha_shape_vertices_end(); 
    ++it) {

    const auto& p = (*it)->point();
    res.push_back(p.point());
  }
  
  std::vector<Point_3> pts;
  for (const auto& p : res)
    pts.push_back(Point_3(p.x(), p.y(), FT(0)));
  
  Saver saver;
  saver.export_points(pts, Color(0, 0, 0), "/Users/monet/Documents/lod/logs/pts");

  return 0;
}
