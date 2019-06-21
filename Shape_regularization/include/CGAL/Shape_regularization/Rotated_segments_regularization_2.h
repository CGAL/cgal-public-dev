#ifndef CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange,
    typename SegmentMap>
  class Rotated_segments_regularization_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Segment = typename GeomTraits::Segment_2;
    using Vector  = typename GeomTraits::Vector_2;

    Rotated_segments_regularization_2 (
      const InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(input_range.size() > 0);

    }

    FT target_value(const int i, const int j) {
      
      Vector vector_i = m_input_range[i].to_vector(); 
      Vector vector_j = m_input_range[j].to_vector(); 
      // compute_direction
      if (vector_i.y() < FT(0) || (vector_i.y() == FT(0) && vector_i.x() < FT(0))) 
        vector_i = -vector_i;
      if (vector_j.y() < FT(0) || (vector_j.y() == FT(0) && vector_j.x() < FT(0))) 
        vector_j = -vector_j;

      //compute_orientation
      const FT atan_i = static_cast<FT>(std::atan2(CGAL::to_double(vector_i.y()), CGAL::to_double(vector_i.x())));
      FT orientation_i = atan_i * FT(180) / static_cast<FT>(CGAL_PI);
      if (orientation_i < FT(0)) 
        orientation_i += FT(180);

      const FT atan_j = static_cast<FT>(std::atan2(CGAL::to_double(vector_j.y()), CGAL::to_double(vector_j.x())));
      FT orientation_j = atan_j * FT(180) / static_cast<FT>(CGAL_PI);
      if (orientation_j < FT(0)) 
        orientation_j += FT(180);

      const FT mes_ij = orientation_i - orientation_j;
      const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

      const FT  t_ij = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

  // not sure if we need r_ij
      int      r_ij;
      if (CGAL::abs(to_lower) < CGAL::abs(to_upper))
          r_ij = ((90 * static_cast<int>(mes90)) % 180 == 0 ? 0 : 1);
      else
          r_ij = ((90 * static_cast<int>(mes90 + 1.0)) % 180 == 0 ? 0 : 1);

      return t_ij;
    }

    // FT target_value(const int i, const int j) {return FT value} // takes indices of 2 segments and returns angle value; look up: regular segment in the old code
    // calculate t_ij and return it (like in Delaunay_neighbours_graph_builder)
    // we also need r_ij
    // void update(std::vector<FT> & result) {} // reorients (rotates) segments
    // class Tree from the old code

  private:
    // Fields.
    const Input_range& m_input_range;
    const Segment_map  m_segment_map;

    /*
    void compute_direction(int i, int j) {
    }

    void compute_orientation() {
    }
    */

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2