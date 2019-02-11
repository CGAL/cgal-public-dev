#ifndef CGAL_LEVELS_OF_DETAIL_TREE_BOUNDARIES_2_H
#define CGAL_LEVELS_OF_DETAIL_TREE_BOUNDARIES_2_H

// STL includes.
#include <vector>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Tree_boundaries_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;

    void create_boundary_segments(
      const Point_2& center,
      const FT radius,
      const std::size_t min_faces_per_tree,
      std::vector<Segment_2>& segments) const {

      const std::size_t n = min_faces_per_tree;

      segments.clear();
      segments.resize(n);
      
      std::vector<Point_2> points(n);
      for (std::size_t i = 0; i < n; ++i) {

        const FT angle = 
        FT(2) * static_cast<FT>(CGAL_PI) * 
        (static_cast<FT>(i) / static_cast<FT>(n));

        points[i] = center + Vector_2(
          radius * static_cast<FT>(std::cos(CGAL::to_double(angle))),
          radius * static_cast<FT>(std::sin(CGAL::to_double(angle))));
      }

      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;

        const Point_2& p1 = points[i];
        const Point_2& p2 = points[ip];

        segments[i] = Segment_2(p1, p2);
      }
    }

  }; // Tree_boundaries_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_TREE_BOUNDARIES_2_H
