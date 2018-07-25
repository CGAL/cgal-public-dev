#ifndef CGAL_LEVEL_OF_DETAIL_AVERAGE_SPACING_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_AVERAGE_SPACING_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits>
class Average_spacing_estimator {
  
public:
  using Kernel    = GeomTraits;
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Local_kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Local_FT 	    = typename Local_kernel::FT;
  using Local_point_3 = typename Local_kernel::Point_3;

  Average_spacing_estimator(const size_t num_neighbours) :
    m_num_neighbours(num_neighbours)
  { }

  FT compute_on_segments_2(const std::vector<Segment_2>& segments) const {

    CGAL_precondition(segments.size() > 0);

    std::vector<Local_point_3> points;
    points.reserve(segments.size() * 2);

    for (std::size_t s = 0; s < segments.size(); ++ s) {
      const Segment_2 &segment = segments[s];

      const Point_2 &source = segment.source();
      const Point_2 &target = segment.target();

      const Local_FT sx = static_cast<Local_FT>(CGAL::to_double(source.x()));
      const Local_FT sy = static_cast<Local_FT>(CGAL::to_double(source.y()));

      const Local_FT tx = static_cast<Local_FT>(CGAL::to_double(target.x()));
      const Local_FT ty = static_cast<Local_FT>(CGAL::to_double(target.y()));

      points.push_back(Local_point_3(sx, sy, Local_FT(0)));
      points.push_back(Local_point_3(tx, ty, Local_FT(0)));
    }

    const Local_FT average_spacing
      = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, m_num_neighbours,
       CGAL::parameters::point_map(CGAL::Identity_property_map<Local_point_3>())
       .geom_traits (Local_kernel()));
                
    return static_cast<FT>(average_spacing);
  }

private:
  const size_t m_num_neighbours;
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_AVERAGE_SPACING_ESTIMATOR_H
