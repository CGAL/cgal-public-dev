#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;
  // using Point_2 = typename Traits::Point_2;
  using FT = typename Traits::FT;

  template<typename Point_2>
  Point_2 compute_middle_point(const Point_2& source, const Point_2& target) {
    const FT half = FT(1) / FT(2);
    const FT x = half * (source.x() + target.x());
    const FT y = half * (source.y() + target.y());
    return Point_2(x, y);
  }

  template<typename Vector, typename Segment_2>
  Vector compute_direction(const Segment_2& segment) {
    Vector v = segment.to_vector(); 
    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
      v = -v;
    return v;
  }
  
  template<typename Vector>
  FT compute_orientation(Vector v) {
    const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
    FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
    if (orientation < FT(0)) 
      orientation += FT(180);
    return orientation;
  }



} // internal
} // Regularization
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H