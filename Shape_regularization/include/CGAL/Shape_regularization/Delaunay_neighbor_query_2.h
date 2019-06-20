#ifndef CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2
#define CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange, 
    typename SegmentMap> //property map - read!
  class Delaunay_neighbor_query_2 {

  public:

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    Delaunay_neighbor_query_2(
      const InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(input_range.size() > 0);

    }

    // void operator() { 
      // returns std::vector indicies of neighbors
      // Use Delaunay triangulation to find neighbors
    // }

    //should be very similar to the K_neighbor_query class
    // in constructor:
    // 1) Check all segments in a "for" loop
    // 2) Call SegmentMap which returns a segment (center)
    // 3) Calculate a sentral point for each segment
    // 4) Send central points from each segment to Delaunay triangulation 
    // (from the old code); read Delaunay triangulation documentation
    // 5) Save Delaunay triangulation

  private:

    // Fields.
    const Input_range& m_input_range;
    // const std::size_t m_number_of_neighbors;
    const Segment_map m_segment_map;

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2