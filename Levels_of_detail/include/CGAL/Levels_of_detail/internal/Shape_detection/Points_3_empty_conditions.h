#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_3_EMPTY_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_3_EMPTY_CONDITIONS_H

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Points_3_empty_conditions {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    Points_3_empty_conditions(
      const Input_range& input_range,
      const Point_map point_map) :
    m_input_range(input_range),
    m_point_map(point_map) {

      CGAL_precondition(m_input_range.size() > 0);
    }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>&) const {
      return true;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return region.size() >= 10;
    }

    void update(const std::vector<std::size_t>&) { }

  private:
        
    // Fields.
    const Input_range& m_input_range;
    const Point_map m_point_map;

  }; // Points_3_empty_conditions

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_3_EMPTY_CONDITIONS_H
