#ifndef CGAL_SHAPE_REGULARIZATION_SEGMENT_DATA_2
#define CGAL_SHAPE_REGULARIZATION_SEGMENT_DATA_2

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits>
  struct Segment_data_2 {

  public:

    using Traits = GeomTraits;
    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Vector  = typename GeomTraits::Vector_2;
    using FT = typename GeomTraits::FT;

    Segment_data_2(
      GeomTraits& segment):
    m_segment(segment) {}

  private:
    size_t  m_index;
    const Segment& m_segment;
    FT      m_orientation;
    FT      m_difference;
    FT      m_length;
    Point   m_barycentre;
    Vector  m_direction;

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2