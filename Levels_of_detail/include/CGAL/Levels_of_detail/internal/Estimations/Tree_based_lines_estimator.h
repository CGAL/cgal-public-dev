#ifndef CGAL_LEVELS_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H

// STL includes.
#include <map>

#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<class GeometricTraits, class Tree>
class Tree_based_lines_estimator {

public:
  using Traits = GeometricTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Line_2 = typename Traits::Line_2;

  Tree_based_lines_estimator(
    const std::vector<Point_2> &points, 
    const Tree &tree) : 
  m_tree(tree) {

    estimate_lines_2(points);
  }

  const std::map<Point_2, Line_2> &lines_2() const {
    return m_lines_2;
  }

  class Sorter {

  public:
    Sorter(
      const std::vector<Point_2> &points, 
      const std::map<Point_2, FT> &scores) :
    m_points(points), 
    m_scores(scores) 
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {
      return m_scores.at(m_points[i]) > m_scores.at(m_points[j]);
    }

  private:
    const std::map<Point_2, FT> &m_scores;
    const std::vector<Point_2> &m_points;
  };

  Sorter sorter(const std::vector<Point_2> &points) const {
    return Sorter(points, m_scores);
  }

private:
  const Tree &m_tree;

  std::map<Point_2, FT> m_scores;
  std::map<Point_2, Line_2> m_lines_2;

  void estimate_lines_2(const std::vector<Point_2> &points) {
                
    m_scores.clear();
    m_lines_2.clear();
                
    Line_2 line;
    std::vector<std::size_t> neighbors;

    for (std::size_t i = 0; i < points.size(); ++i) {
      
      const Point_2 &point = points[i];
      m_tree.search_2(point, neighbors);
      
      m_scores[point] = internal::fit_line_to_points_2(points, neighbors, line);
      m_lines_2[point] = line;
    }
  }

}; // Tree_based_lines_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
