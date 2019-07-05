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



} // internal
} // Regularization
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H