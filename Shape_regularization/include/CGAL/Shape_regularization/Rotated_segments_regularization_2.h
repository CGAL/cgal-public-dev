#ifndef CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange,
    typename SegmentMap>
  class Rotated_segments_regularization_2 {
  public:
    Rotated_segments_regularization_2 (); 
    // FT target_value(const int i, const int j) {return FT value} // takes indices of 2 segments and returns angle value; look up: regular segment in the old code
    // calculate t_ij and return it (like in Delaunay_neighbours_graph_builder)
    // we also need r_ij
    // void update(std::vector<FT> & result) {} // reorients (rotates) segments
    // class Tree from the old code
  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2