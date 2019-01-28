#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H

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

  template<typename GeomTraits>
  class Points_2_least_squares_line_fit_conditions {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    typename Traits::Compute_squared_length_2 squared_length_2;

    Points_2_least_squares_line_fit_conditions(
      const std::vector<Point_2>& points,
      const std::vector<Vector_2>& normals, 
      const FT noise_level, 
      const FT angle, 
      const FT min_length) :
    m_points(points),
    m_normals(normals),
    m_distance_threshold(noise_level),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_length(min_length) {

      CGAL_precondition(m_points.size() > 0);
      CGAL_precondition(m_points.size() == m_normals.size());

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));
      CGAL_precondition(m_min_length > FT(0));
    }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_points.size());

      const Point_2& query_point = m_points[query_index];
      const Vector_2& query_normal = m_normals[query_index];

      const FT distance_to_fitted_line = 
      internal::distance_2(query_point, m_line_of_best_fit);
      
      const FT cos_angle = 
      internal::cos_angle_2(query_normal, m_normal_of_best_fit);

      return (( distance_to_fitted_line <= m_distance_threshold ) && 
        ( cos_angle >= m_normal_threshold ));
    }

    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      
      const FT squared_min_length = m_min_length * m_min_length;
      return internal::points_squared_length_2(
        m_points, region, m_line_of_best_fit) >= squared_min_length;
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) {
                    
        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_points.size());

        const Point_2& point = m_points[region[0]];
        const Vector_2& normal = m_normals[region[0]];

        m_line_of_best_fit = Line_2(point, normal).perpendicular(point);
        m_normal_of_best_fit = normal;

      } else {

        internal::line_from_points_2(
          m_points, region, m_line_of_best_fit);
                    
        const Vector_2 normal = 
        m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
        const FT normal_length = static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              squared_length_2(normal))));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

  private:
        
    // Fields.
    const std::vector<Point_2>& m_points;
    const std::vector<Vector_2>& m_normals;
            
    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const FT m_min_length;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;

  }; // Points_2_least_squares_line_fit_conditions

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H
