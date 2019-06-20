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

    //should be very similar to the K_neighbor_query class
    // in constructor:
    // 1) Check all segments in a "for" loop
    // 2) Call SegmentMap which returns a segment (center)
    // 3) Calculate a sentral point for each segment
    // 4) Send central points from each segment to Delaunay triangulation 
    // (from the old code); read Delaunay triangulation documentation
    // 5) Save Delaunay triangulation

  private:

    void operator() { // returns std::vector indicies of neighbors
      // Use Delaunay triangulation to find neighbors
    }
  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2