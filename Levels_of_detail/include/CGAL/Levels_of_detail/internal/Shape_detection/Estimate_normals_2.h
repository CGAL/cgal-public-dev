#ifndef CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_2_H
#define CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_2_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits, 
typename Connectivity>
class Estimate_normals_2 {

public:
  using Traits = GeomTraits;
  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Vector_2 = typename Traits::Vector_2;
  using Line_2 = typename Traits::Line_2;

  Estimate_normals_2(
    const std::vector<Point_2>& points, 
    const Connectivity& connectivity) : 
  m_points(points),
  m_connectivity(connectivity) {
    
    CGAL_precondition(m_points.size() > 0);
    estimate_normals();
  }

  const std::vector<Vector_2>& normals() const {
    
    CGAL_precondition(m_normals.size() == m_points.size());
    return m_normals;
  }

  class Sorter {

  public:
    Sorter(const std::vector<FT>& scores) :
    m_scores(scores) { 

      CGAL_precondition(m_scores.size() > 0);
    }

    bool operator()(const std::size_t i, const std::size_t j) const {
      
      CGAL_precondition(i >= 0 && i < m_scores.size());
      CGAL_precondition(j >= 0 && j < m_scores.size());

      return m_scores[i] > m_scores[j];
    }

  private:
    const std::vector<FT>& m_scores;
  };

  Sorter sorter() const {
    
    CGAL_precondition(m_scores.size() == m_points.size());
    return Sorter(m_scores);
  }

private:
  const std::vector<Point_2>& m_points;
  const Connectivity& m_connectivity;

  std::vector<FT> m_scores;
  std::vector<Vector_2> m_normals;

  void estimate_normals() {

    m_scores.clear();
    m_scores.resize(m_points.size());

    m_normals.clear();
    m_normals.resize(m_points.size());
                
    Line_2 line;
    std::vector<std::size_t> neighbors;
    for (std::size_t i = 0; i < m_points.size(); ++i) {

      m_connectivity.get_neighbors(i, neighbors);
      m_scores[i] = 
      internal::line_from_points_2(m_points, neighbors, line);

      m_normals[i] = line.to_vector();
      m_normals[i] = 
      m_normals[i].perpendicular(CGAL::COUNTERCLOCKWISE);

      const FT normal_length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(m_normals[i] * m_normals[i])));

      CGAL_precondition(normal_length > FT(0));
      m_normals[i] /= normal_length;
    }
  }

}; // Estimate_normals_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_2_H
