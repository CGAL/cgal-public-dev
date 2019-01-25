#ifndef CGAL_LEVELS_OF_DETAIL_THINNING_2_H
#define CGAL_LEVELS_OF_DETAIL_THINNING_2_H

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename KdTree>
  class Thinning_2 {

  public:
    using Traits = GeomTraits;
    using Tree = KdTree;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    Thinning_2(const Tree& tree, const FT scale) :
    m_tree(tree),
    m_scale(scale)
    { }

    void apply(std::vector<Point_2>& range) const {


    }

  private:
    const Tree& m_tree;
    const FT m_scale;

  }; // Thinning_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_THINNING_2_H
