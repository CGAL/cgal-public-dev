#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_3_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_3_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H

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
  typename PointMap>
  class Points_3_least_squares_plane_fit_conditions {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    typename Traits::Compute_squared_length_3 squared_length_3;

    Points_3_least_squares_plane_fit_conditions(
      const Input_range& input_range,
      const Point_map point_map,
      const std::vector<Vector_3>& normals, 
      const FT noise_level, 
      const FT angle, 
      const FT min_area) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_normals(normals),
    m_distance_threshold(noise_level),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_area(min_area) {

      CGAL_precondition(m_input_range.size() > 0);
      CGAL_precondition(m_input_range.size() == m_normals.size());

      CGAL_precondition(m_distance_threshold >= FT(0));
      CGAL_precondition(m_normal_threshold >= FT(0) && m_normal_threshold <= FT(1));
      CGAL_precondition(m_min_area > FT(0));
    }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      const Point_3& query_point = 
      get(m_point_map, *(m_input_range.begin() + query_index));
      
      const Vector_3& query_normal = 
      m_normals[query_index];

      const FT distance_to_fitted_plane = 
      internal::distance_3(query_point, m_plane_of_best_fit);
      
      const FT cos_angle = 
      internal::cos_angle_3(query_normal, m_normal_of_best_fit);

      return (( distance_to_fitted_plane <= m_distance_threshold ) && 
        ( cos_angle >= m_normal_threshold ));
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {

      return internal::points_area_3(
        m_input_range, m_point_map, region, m_plane_of_best_fit, m_min_area) 
        >= m_min_area;
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) {
                    
        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        const Point_3& point = 
        get(m_point_map, *(m_input_range.begin() + region[0]));

        const Vector_3& normal = 
        m_normals[region[0]];

        m_plane_of_best_fit = Plane_3(point, normal);
        m_normal_of_best_fit = normal;

      } else {

        internal::plane_from_points_3(
          m_input_range, m_point_map, region, m_plane_of_best_fit);
                    
        const Vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
        const FT normal_length = static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              squared_length_3(normal))));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

  private:
        
    // Fields.
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::vector<Vector_3>& m_normals;
            
    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const FT m_min_area;

    Plane_3 m_plane_of_best_fit;
    Vector_3 m_normal_of_best_fit;

  }; // Points_3_least_squares_plane_fit_conditions

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_3_LEAST_SQUARES_PLANE_FIT_CONDITIONS_H
