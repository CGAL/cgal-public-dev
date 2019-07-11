#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2

// #include <CGAL/license/Shape_regularization.h>

// #include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits, 
    typename InputRange>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;

    Tree(
      InputRange& input_range) :
    m_input_range(input_range) {
      CGAL_precondition(input_range.size() > 0);
    }

  private:
    Input_range& m_input_range;
    

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2