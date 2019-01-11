#ifndef CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H

// STL includes.
#include <map>

#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Line_to_points_fitter.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class Tree>
class Tree_based_lines_estimator {

public:
  using Kernel    = GeomTraits;

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Line_2  = typename Kernel::Line_2;

  using Line_to_points_fitter   = LOD::Line_to_points_fitter<Kernel>;

  using Lines_2 = std::map<Point_2, Line_2>;
  using Scores  = std::map<Point_2, FT>;

  Tree_based_lines_estimator(const std::vector<Point_2>& points, const Tree& tree)
    : m_tree(tree)
  { 
    estimate_lines_2(points);
  }

  inline const Scores& scores() const {
    return m_scores;
  }

  inline const Lines_2& lines_2() const {
    return m_lines_2;
  }


  class Sorter {

    const std::vector<Point_2>& pts;
  public:
    Sorter (const std::vector<Point_2>& pts, const Scores& scores) :
      pts (pts), m_scores(scores) 
    { }

    bool operator() (const std::size_t& i, const std::size_t& j) const
    {
      return m_scores.at(pts[i]) > m_scores.at(pts[j]);
    }

  private:
    const Scores& m_scores;
  };

  Sorter sorter (const std::vector<Point_2>& pts) const
  {
    return Sorter (pts, m_scores);
  }


private:

  const Tree      &m_tree;

  Scores  m_scores;
  Lines_2 m_lines_2;

  void estimate_lines_2(const std::vector<Point_2>& points) {
                
    m_scores.clear();
    m_lines_2.clear();
                
    Line_2     line;
    std::vector<std::size_t> neighbors;

    const Line_to_points_fitter line_to_points_fitter;
    for (std::size_t i = 0; i < points.size(); ++ i)
    {
      const Point_2 &point = points[i];
      m_tree.search_2(point, neighbors);
      m_scores[point] = line_to_points_fitter.fit_line_2(neighbors, points, line);
      m_lines_2[point] = line;
    }
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_BASED_LINES_ESTIMATOR_H
