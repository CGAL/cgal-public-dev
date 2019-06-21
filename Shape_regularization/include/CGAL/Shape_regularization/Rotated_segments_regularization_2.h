#ifndef CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <utility>

// use std::map where key-> pair and value t_ij
// the same is for r_ij and mu_ij
// make private functions for t_ij and r_ij calculations

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
      
      Vector v_i = compute_direction(i);
      Vector v_j = compute_direction(j);

      //compute_orientation
      const FT mes_ij = compute_orientation(v_i) - compute_orientation(v_j);
      const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

      const FT  t_ij = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

      // m_t_ijs.insert(std::pair<std::pair<int, int>, FT> (std::make_pair(i, j), t_ij));
      m_t_ijs[std::make_pair(i, j)] = t_ij;

      // we will need r_ij in update();  
      int      r_ij;
      if (CGAL::abs(to_lower) < CGAL::abs(to_upper))
          r_ij = ((90 * static_cast<int>(mes90)) % 180 == 0 ? 0 : 1);
      else
          r_ij = ((90 * static_cast<int>(mes90 + 1.0)) % 180 == 0 ? 0 : 1);
      
      //  m_r_ijs.insert(std::pair<std::pair<int, int>, FT> (std::make_pair(i, j), r_ij));
      m_r_ijs[std::make_pair(i, j)] = r_ij;
  
      return t_ij;
    }

    // FT target_value(const int i, const int j) {return FT value} // takes indices of 2 segments and returns angle value; look up: regular segment in the old code
    // calculate t_ij and return it (like in Delaunay_neighbours_graph_builder)
    // we also need r_ij
    void update(std::vector<FT> & result) {
      // reoirent segments from regularize angles (old code)
    } // reorients (rotates) segments
    // class Tree from the old code

    void debug_trmu_ijs() {
      std::cout << std::endl << "m_t_ijs: " << std::endl;
      for (typename std::map<std::pair<int, int>, FT>::iterator it = m_t_ijs.begin(); it!=m_t_ijs.end(); ++it)
        std::cout << "(" << it->first.first << ", " << it->first.second << ") => " << it->second << std::endl;
      std::cout << std::endl << "m_r_ijs: " << std::endl;
      for (typename std::map<std::pair<int, int>, FT>::iterator it = m_r_ijs.begin(); it!=m_r_ijs.end(); ++it)
        std::cout << "(" << it->first.first << ", " << it->first.second << ") => " << it->second << std::endl;
      std::cout << std::endl << "m_mu_ij = " << m_mu_ij << std::endl;
    }

  private:
    // Fields.
    const Input_range& m_input_range;
    const Segment_map  m_segment_map;
    std::map <std::pair<int, int>, FT> m_t_ijs;
    std::map <std::pair<int, int>, FT> m_r_ijs;
    const FT m_mu_ij = FT(4) / FT(5);

    
    Vector compute_direction(const int i) {
      Vector v = m_input_range[i].to_vector(); 
      if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
        v = -v;
      return v;
    }
    
    FT compute_orientation(Vector v) {
      const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
      FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
      if (orientation < FT(0)) 
        orientation += FT(180);
      return orientation;
    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2