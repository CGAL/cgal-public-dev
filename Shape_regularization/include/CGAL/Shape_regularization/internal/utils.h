#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<typename Point_2>
  Point_2 compute_middle_point(const Point_2& source, const Point_2& target) {
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    const FT half = FT(1) / FT(2);
    const FT x = half * (source.x() + target.x());
    const FT y = half * (source.y() + target.y());
    return Point_2(x, y);
  }

  template<typename Vector, typename Segment_2>
  Vector compute_direction(const Segment_2& segment) {
    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;

    Vector v = segment.to_vector(); 
    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
      v = -v;
    return v;
  }
  
  template<typename Vector>
  typename Kernel_traits<Vector>::Kernel::FT
  compute_orientation(const Vector& v) {
    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;

    const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
    FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
    if (orientation < FT(0)) 
      orientation += FT(180);
    return orientation;
  }

  /* template<typename FT>
  FT compute_tij(const FT orientation_i, const FT orientation_j) {
      const FT mes_ij = orientation_i - orientation_j;
      const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

      const FT  t_ij = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;
      return t_ij;
  } */



} // internal
} // Regularization
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_UTILS_H
