#ifndef CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2
#define CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <vector>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits,
    typename InputRange,
    typename SegmentMap>
  struct Parallel_groups_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Segment = typename GeomTraits::Segment_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;

    Parallel_groups_2 (
      const InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(m_input_range.size() > 0);

      build_segment_data();
      make_parallel_groups();

    }

    // void set_parallel_groups(const Groups_type & parallel_groups) {
    //   m_parallel_groups = parallel_groups;
    // }

    // Groups_type parallel_groups() {
    //   return m_parallel_groups;
    // }
    std::vector <std::vector <std::size_t>> get_parallel_groups() {
      CGAL_precondition(m_parallel_groups_angle_map.size() > 0);
      std::vector <std::vector <std::size_t>> parallel_groups;

      for(const auto & mi : m_parallel_groups_angle_map) {
        const std::vector <std::size_t> & group = mi.second;
        parallel_groups.push_back(group);
      }
      CGAL_postcondition(parallel_groups.size() > 0);
      return parallel_groups;
    }


  private:
    const Input_range& m_input_range;
    const Segment_map  m_segment_map;
    std::vector<Segment_data> m_segments;
    // std::vector<std::vector<std::size_t>> m_parallel_groups;
    std::map <FT, std::vector<std::size_t>> m_parallel_groups_angle_map;
    // Groups_type m_parallel_groups;

    void build_segment_data() {
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + i));
        const Segment_data seg_data(seg, i);
        m_segments.push_back(seg_data);
      }
      CGAL_postcondition(m_segments.size() > 0);
    }

    void make_parallel_groups() {
      for (const auto & seg : m_segments) {
        const FT angle = floor(seg.m_orientation);
        const std::size_t seg_index = seg.m_index;
        std::cout << seg_index << "). m_orientation = " << seg.m_orientation << "; angle = " << angle << std::endl;
        m_parallel_groups_angle_map[angle].push_back(seg_index);
      }
    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2