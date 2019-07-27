#ifndef CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>
#include <CGAL/Shape_regularization/internal/Conditions_angles_2.h>


namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange,
    typename SegmentMap>
  class Angle_regularization_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Conditions = typename internal::Conditions_angles_2<Traits>;
    using Grouping = internal::Grouping_segments_2<Traits, Conditions>;
    using Vector  = typename GeomTraits::Vector_2;
    using Targets_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>>;
    using Relations_map = std::map <std::pair<std::size_t, std::size_t>, std::pair<int, std::size_t>>;

    Angle_regularization_2 (
      InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(Grouping()) {

      CGAL_precondition(m_input_range.size() > 0);

      build_segment_data_map();

      CGAL_postcondition(m_segments.size() > 0);
    }

    FT target_value(const std::size_t i, const std::size_t j) {
 
      CGAL_precondition(m_segments.find(i) != m_segments.end());
      CGAL_precondition(m_segments.find(j) != m_segments.end());

      const Segment_data & s_i = m_segments.at(i);
      const Segment_data & s_j = m_segments.at(j);

      const FT mes_ij = s_i.m_orientation - s_j.m_orientation;
      const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

      const FT tar_val = CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;
      if (CGAL::abs(tar_val) < bound(i) + bound(j)) {
        m_targets[std::make_pair(i, j)] = tar_val;
 
        int rel_val;
        if (CGAL::abs(to_lower) < CGAL::abs(to_upper))
            rel_val = ((90 * static_cast<int>(mes90)) % 180 == 0 ? 0 : 1);
        else
            rel_val = ((90 * static_cast<int>(mes90 + 1.0)) % 180 == 0 ? 0 : 1);
        
        m_relations[std::make_pair(i, j)] = rel_val;
      } 
  
      return tar_val;
    }

    FT bound(const std::size_t i) {
      FT theta_max;
      m_input_range.size() > 3 ? theta_max = FT(25) : theta_max = FT(10);
      // theta_max = FT(25);
      return theta_max;
    }

    std::map<FT, std::vector<std::size_t>> parallel_groups_angle_map() {
      CGAL_precondition(m_parallel_groups_angle_map.size() > 0);
      return m_parallel_groups_angle_map;
    }


    void update(std::vector<FT> & result) {

      Targets_map targets;
      Relations_map relations;

      CGAL_precondition(m_targets.size() > 0);
      CGAL_precondition(m_targets.size() == m_relations.size());

      build_grouping_data(targets, relations);

      m_grouping.make_groups(m_input_range.size(), m_segments, result, m_parallel_groups_angle_map, targets, relations);
      rotate_parallel_segments();
    }


  private:
    Input_range& m_input_range;
    const Segment_map  m_segment_map;
    std::map <std::size_t, Segment_data> m_segments;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_targets;
    std::map <std::pair<std::size_t, std::size_t>, int> m_relations;
    Grouping m_grouping;
    std::map <FT, std::vector<std::size_t>> m_parallel_groups_angle_map;

    void build_segment_data_map() {

      m_segments.clear();
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + i));
        const Segment_data seg_data(seg, i);

        m_segments.emplace(i, seg_data);
      }

    }

    void build_grouping_data (Targets_map & targets,
                              Relations_map & relations) {

      auto ri = m_relations.begin();
      std::size_t tar_index = 0;

      for (const auto & ti : m_targets) {
        const std::size_t seg_index_tar_i = ti.first.first;
        const std::size_t seg_index_tar_j = ti.first.second;
        const FT tar_val = ti.second;

        const std::size_t seg_index_rel_i = ri->first.first;
        const std::size_t seg_index_rel_j = ri->first.second;
        const int rel_val = ri->second;

        CGAL_precondition(seg_index_tar_i == seg_index_rel_i);
        CGAL_precondition(seg_index_tar_j == seg_index_rel_j);

        targets[std::make_pair(seg_index_tar_i, seg_index_tar_j)] = std::make_pair(tar_val, tar_index);
        relations[std::make_pair(seg_index_rel_i, seg_index_rel_j)] = std::make_pair(rel_val, tar_index);

        ++ri;
        ++tar_index;
      }

      CGAL_postcondition(targets.size() == relations.size());
      CGAL_postcondition(targets.size() == tar_index);

    }

    void rotate_parallel_segments() {

      for (const auto & mi : m_parallel_groups_angle_map) {
        const FT theta = mi.first;
        const std::vector<std::size_t> & group = mi.second;

        // Each group of parallel segments has a normal vector that we compute with alpha.
        const FT x = static_cast<FT>(cos(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));
        const FT y = static_cast<FT>(sin(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));

        Vector v_dir = Vector(x, y);
        const Vector v_ort = Vector(-v_dir.y(), v_dir.x());
        
        const FT a = v_ort.x();
        const FT b = v_ort.y();

        // Rotate segments with precision.
        for (std::size_t i = 0; i < group.size(); ++i) {

          std::size_t seg_index = group[i];

          // Compute equation of the supporting line of the rotated segment.
          const Segment_data & seg_data = m_segments.at(seg_index);
          const Point & barycentre = seg_data.m_barycentre;
          const FT c = -a * barycentre.x() - b * barycentre.y();

          set_orientation(seg_index, a, b, c, v_dir);
        }
      }
    }

    void set_orientation(const std::size_t i, const FT a, const FT b, const FT c, const Vector &direction) {

      Vector l_direction = direction;
      
      if (l_direction.y() < FT(0) || (l_direction.y() == FT(0) && l_direction.x() < FT(0))) 
        l_direction = -l_direction;
      FT x1, y1, x2, y2;
      const Segment_data & seg_data = m_segments.at(i);
      const Point barycentre = seg_data.m_barycentre;
      const FT length = seg_data.m_length;
      if (CGAL::abs(l_direction.x()) > CGAL::abs(l_direction.y())) { 
        x1 = barycentre.x() - length * l_direction.x() / FT(2);
        x2 = barycentre.x() + length * l_direction.x() / FT(2);

        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;
      } 
      else {
        y1 = barycentre.y() - length * l_direction.y() / FT(2);
        y2 = barycentre.y() + length * l_direction.y() / FT(2);

        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }
      const Point source = Point(x1, y1);
      const Point target = Point(x2, y2);

      m_input_range[i] = Segment(source, target);
    } 

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2
