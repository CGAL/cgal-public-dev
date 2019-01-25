#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap, 
  typename NormalMap>
  class Points_2_least_squares_line_fit_conditions {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;
    using Vector_2 = typename GeomTraits::Vector_2;
    using Line_2 = typename GeomTraits::Line_2;

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_2 = typename Traits::Compute_squared_length_2;
    using Squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Scalar_product_2 = typename Traits::Compute_scalar_product_2;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;

    Points_2_least_squares_line_fit_conditions(
      const InputRange& input_range, 
      const FT distance_threshold = FT(1), 
      const FT normal_threshold = FT(9) / FT(10), 
      const std::size_t min_region_size = 2, 
      const PointMap point_map = PointMap(), 
      const NormalMap normal_map = NormalMap(), 
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(normal_threshold),
    m_min_region_size(min_region_size),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_2(traits.compute_squared_length_2_object()),
    m_squared_distance_2(traits.compute_squared_distance_2_object()),
    m_scalar_product_2(traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)) {

      CGAL_precondition(distance_threshold >= FT(0));
      CGAL_precondition(normal_threshold >= FT(0) && normal_threshold <= FT(1));
      CGAL_precondition(min_region_size > 0);
    }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      const auto& key = *(m_input_range.begin() + query_index);
      const Point_2& query_point = get(m_point_map, key);

      const Vector_2& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_2(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_2 query_normal = normal / normal_length;

      const FT distance_to_fitted_line = 
      m_sqrt(m_squared_distance_2(query_point, m_line_of_best_fit));
      
      const FT cos_angle = 
      CGAL::abs(m_scalar_product_2(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_line <= m_distance_threshold ) && 
        ( cos_angle >= m_normal_threshold ));
    }

    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and normal
                    
        CGAL_precondition(region[0] >= 0);
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit line will be a line through this point with its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_2& point = get(m_point_map, key);
        const Vector_2& normal = get(m_normal_map, key);
                    
        const FT normal_length = m_sqrt(m_squared_length_2(normal));
        CGAL_precondition(normal_length > FT(0));

        m_line_of_best_fit = Line_2(point, normal).perpendicular(point);
        m_normal_of_best_fit = normal / normal_length;

      } else { // update reference line and normal

        std::vector<Local_point_2> points(region.size());
        for (std::size_t i = 0; i < region.size(); ++i) {

          CGAL_precondition(region[i] >= 0);
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points[i] = m_to_local_converter(get(m_point_map, key));
        }
        CGAL_precondition(points.size() > 0);

        Local_line_2  fitted_line;
        Local_point_2 fitted_centroid;

        // The best fit line will be a line fitted to all region points with its normal being perpendicular to the line.
        #ifndef CGAL_EIGEN2_ENABLED
          linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Default_diagonalize_traits<Local_FT, 2>());
        #else 
          linear_least_squares_fitting_2(
            points.begin(), points.end(), 
            fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(), 
            Local_traits(), CGAL::Eigen_diagonalize_traits<Local_FT, 2>());
        #endif
                    
        m_line_of_best_fit = 
        Line_2(
          static_cast<FT>(fitted_line.a()), 
          static_cast<FT>(fitted_line.b()), 
          static_cast<FT>(fitted_line.c()));
                    
        const Vector_2 normal = 
        m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

  private:
        
    // Fields.
    const Input_range& m_input_range;
            
    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const std::size_t m_min_region_size;

    const Point_map m_point_map;
    const Normal_map m_normal_map;
            
    const Squared_length_2 m_squared_length_2;
    const Squared_distance_2 m_squared_distance_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;
  };

} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_2_LEAST_SQUARES_LINE_FIT_CONDITIONS_H
