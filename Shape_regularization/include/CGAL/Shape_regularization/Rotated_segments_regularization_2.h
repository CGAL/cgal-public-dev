#ifndef CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Shape_regularization/internal/Tree.h>
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

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
    using Tree = internal::Tree<Traits, Input_range>;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Rotated_segments_regularization_2 (
      InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_mu_ij(FT(4) / FT(5)) {

      CGAL_precondition(input_range.size() > 0);
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + i));
        const Segment_data seg_data(seg, i);
        m_segments.push_back(seg_data);
      }

    }

    // ~Rotated_segments_regularization_2() {
    //   delete m_tree_pointer;
    // }

    FT target_value(const std::size_t i, const std::size_t j) {

      //compute_orientation
      const FT mes_ij = m_segments[i].m_orientation - m_segments[j].m_orientation;
      const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

      const FT  t_ij = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

      m_t_ijs[std::make_pair(i, j)] = t_ij;

      // we will need r_ij in update();  
      int      r_ij;
      if (CGAL::abs(to_lower) < CGAL::abs(to_upper))
          r_ij = ((90 * static_cast<int>(mes90)) % 180 == 0 ? 0 : 1);
      else
          r_ij = ((90 * static_cast<int>(mes90 + 1.0)) % 180 == 0 ? 0 : 1);
      
      m_r_ijs[std::make_pair(i, j)] = r_ij;
  
      return t_ij;
    }

    // FT target_value(const int i, const int j) {return FT value} // takes indices of 2 segments and returns angle value; look up: regular segment in the old code
    // calculate t_ij and return it (like in Delaunay_neighbours_graph_builder)
    // we also need r_ij
    void update(std::vector<FT> & result) {
      // reoirent segments from regularize angles (old code)
      // reorients (rotates) segments
      // class Tree from the old code

      m_tree_pointer = new Tree(m_input_range, m_t_ijs, m_r_ijs, m_mu_ij, result/* m_final_orientations, m_qp_problem_data, m_parameters */);
      m_tree_pointer->apply_new_orientations();
      // delete m_tree_pointer;


      //segments = input_range
      //get_targets_matrix is my t_ijs vector => save them as a matrix
      //get_targets_matrix is r_ijs => save them as a matrix
      //orientations = me vector from QP solver 

    }

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
    Input_range& m_input_range;
    const Segment_map  m_segment_map;
    std::vector<Segment_data> m_segments;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_t_ijs;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_r_ijs;
    const FT m_mu_ij;
    Tree *m_tree_pointer;

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ROTATED_SEGMENTS_REGULARIZATION_2
