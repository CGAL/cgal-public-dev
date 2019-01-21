#ifndef CGAL_LEVELS_OF_DETAIL_LINES_FROM_POINTS_2_H
#define CGAL_LEVELS_OF_DETAIL_LINES_FROM_POINTS_2_H

// STL includes.
#include <map>

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<class GeomTraits, class Tree>
class Lines_from_points_2 {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Line_2 = typename Traits::Line_2;

  Lines_from_points_2(
    const std::vector<Point_2> &points, 
    const Tree &tree) : 
  m_tree(tree) {

    estimate_lines(points);
  }

  const std::map<Point_2, Line_2> &lines() const {
    return m_lines;
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
  std::map<Point_2, Line_2> m_lines;

  void estimate_lines(const std::vector<Point_2> &points) {
                
    m_scores.clear();
    m_lines.clear();
                
    Line_2 line;
    std::vector<std::size_t> neighbors;

    for (std::size_t i = 0; i < points.size(); ++i) {
      
      const Point_2 &point = points[i];
      m_tree.search_2(point, neighbors);
      
      m_scores[point] = internal::fit_line_to_points_2(points, neighbors, line);
      m_lines[point] = line;
    }
  }

}; // Lines_from_points_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_LINES_FROM_POINTS_2_H
