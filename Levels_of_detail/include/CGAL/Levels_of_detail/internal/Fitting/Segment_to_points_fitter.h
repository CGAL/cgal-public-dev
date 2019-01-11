#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Line_to_points_fitter.h>
#include <CGAL/Level_of_detail/internal/Estimations/End_points_estimator.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits>
class Segment_to_points_fitter {

public:
  using Kernel = GeomTraits;

  using FT        = typename Kernel::FT;
  using Line_2    = typename Kernel::Line_2;
  using Point_2   = typename Kernel::Point_2;
  using Vector_2  = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Line_to_points_fitter    = Line_to_points_fitter<Kernel>;
  using End_points_estimator     = End_points_estimator<Kernel>;

  using Projected_points   = std::vector<Point_2>;

  FT fit_segment_2(const std::vector<std::size_t>& elements,
                   const std::vector<Point_2>& points,
                   Segment_2 &segment) const {
    CGAL_precondition(elements.size() > 1);

    // Fit line to all points.
    Line_2 line;
    const Line_to_points_fitter line_to_points_fitter;
    const FT quality = line_to_points_fitter.fit_line_2(elements, points, line);

    // Project points onto the line.
    Projected_points projected_points;
    for (std::size_t i = 0; i < elements.size(); ++ i)
      projected_points.push_back(line.projection(points[elements[i]]));
                
    // Find end points of the segment;
    Point_2 a, b;

    const End_points_estimator end_points_estimator;
    end_points_estimator.find_end_points_wrt_barycentre_2(projected_points, a, b);

    segment = Segment_2(a, b);
    return quality;
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H
