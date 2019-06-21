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
      // const Segment s_i = m_input_range[i];
      // const Segment s_j = m_input_range[j];

      // const FT mes_ij    = s_i.get_orientation() - s_j.get_orientation();

      // std::cout << "mes_ij = " << mes_ij << std::endl; 

      return FT(1);
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
    Vector  m_direction;
    FT      m_orientation;

/*
    void compute_direction(int i, int j) {
      m_direction.push_back(m_input_range[i]);
      m_direction.push_back(m_input_range[j]);
      if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) m_direction = -m_direction;
    }

    void compute_orientation() {
      compute_direction();

      const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(m_direction.y()), CGAL::to_double(m_direction.x())));
      m_orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);

      if (m_orientation < FT(0)) m_orientation += FT(180);
    }
    */

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2