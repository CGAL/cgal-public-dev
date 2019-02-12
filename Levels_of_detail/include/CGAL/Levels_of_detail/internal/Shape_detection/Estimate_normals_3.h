#ifndef CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_3_H
#define CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_3_H

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
typename InputRange,
typename PointMap, 
typename Connectivity>
class Estimate_normals_3 {

public:
  using Traits = GeomTraits;
  using Input_range = InputRange;
  using Point_map = PointMap;

  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Plane_3 = typename Traits::Plane_3;

  Estimate_normals_3(
    const Input_range& input_range,
    const Point_map point_map, 
    const Connectivity& connectivity) : 
  m_input_range(input_range),
  m_point_map(point_map),
  m_connectivity(connectivity) {
    
    CGAL_precondition(m_input_range.size() > 0);
    estimate_normals();
  }

  const std::vector<Vector_3>& normals() const {
    
    CGAL_precondition(m_normals.size() == m_input_range.size());
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
    
    CGAL_precondition(m_scores.size() == m_input_range.size());
    return Sorter(m_scores);
  }

private:
  const Input_range& m_input_range;
  const Point_map m_point_map;
  const Connectivity& m_connectivity;

  std::vector<FT> m_scores;
  std::vector<Vector_3> m_normals;

  void estimate_normals() {

    m_scores.clear();
    m_scores.resize(m_input_range.size());

    m_normals.clear();
    m_normals.resize(m_input_range.size());
                
    Plane_3 plane;
    std::vector<std::size_t> neighbors;
    for (std::size_t i = 0; i < m_input_range.size(); ++i) {

      m_connectivity.get_neighbors(i, neighbors);
      m_scores[i] = 
      internal::plane_from_points_3(m_input_range, m_point_map, neighbors, plane);

      m_normals[i] = plane.orthogonal_vector();

      const FT normal_length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(m_normals[i] * m_normals[i])));

      CGAL_precondition(normal_length > FT(0));
      m_normals[i] /= normal_length;
    }
  }

}; // Estimate_normals_3

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ESTIMATE_NORMALS_3_H
