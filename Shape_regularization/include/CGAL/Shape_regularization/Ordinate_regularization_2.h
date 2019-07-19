#ifndef CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>

// use std::map where key-> pair and value t_ij
// the same is for r_ij and mu_ij
// make private functions for t_ij and r_ij calculations

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits, 
    typename InputRange,
    typename SegmentMap>
  class Ordinate_regularization_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using FT = typename GeomTraits::FT;
    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Grouping = internal::Grouping_segments_2<Traits>;
    using Vector  = typename GeomTraits::Vector_2;

    Ordinate_regularization_2 (
      InputRange& input_range,
      const std::map<FT, std::vector<std::size_t>> & parallel_groups_angle_map,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_parallel_groups_angle_map(parallel_groups_angle_map),
    m_segment_map(segment_map),
    m_mu_ij(FT(4) / FT(5)) {

      CGAL_precondition(input_range.size() > 0);
      for (const auto & m_i : m_parallel_groups_angle_map) {
        if (m_i.second.size() > 1) {
          Point frame_origin;
          for(std::size_t i = 0; i < m_i.second.size(); ++i) {
            const std::size_t seg_index = m_i.second[i];
            const Segment& seg = get(m_segment_map, *(m_input_range.begin() + seg_index));
            Segment_data seg_data(seg, seg_index);
            if (i == 0) {
              frame_origin = seg_data.m_barycentre;
            }
            const Point reference_coordinates = internal::transform_coordinates(seg_data.m_barycentre, frame_origin, m_i.first);
            seg_data.set_reference_coordinates(reference_coordinates);
            m_segments.emplace(seg_index, seg_data);
          }
        }
      }

      // m_grouping_ptr = std::make_shared<Grouping>(m_segments);

    }

    FT target_value(const std::size_t i, const std::size_t j) {

      //compute_orientation
      const Segment_data & s_i = m_segments.at(i);
      const Segment_data & s_j = m_segments.at(j);

      const FT y_ij = s_i.m_reference_coordinates.y() - s_j.m_reference_coordinates.y();
  
      return y_ij;

    }

    FT bound(const std::size_t i) {
      FT theta_max = FT(0.1);
      return theta_max;
    }

    // FT target_value(const int i, const int j) {return FT value} // takes indices of 2 segments and returns angle value; look up: regular segment in the old code
    // calculate t_ij and return it (like in Delaunay_neighbours_graph_builder)
    // we also need r_ij
    void update(std::vector<FT> & result) {
      // reoirent segments from regularize angles (old code)
      // reorients (rotates) segments
      // class Tree from the old code
      // std::vector<std::vector<std::size_t>> parallel_segments_groups;

    /*
      std::cout << "final orientations after qp: " << result.size() << std::endl;
      for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << std::endl;
      }
      */
    /*  m_parallel_groups_angle_map.clear();
      m_grouping_ptr->make_groups(m_t_ijs, m_r_ijs, m_mu_ij, result, m_parallel_groups_angle_map);

     /* 
      std::cout << "m_parallel_groups_angle_map: " << std::endl;
      std::size_t counter = 0;
      std::size_t iterator = 0;
      for (auto pgi = m_parallel_groups_angle_map.begin(); pgi != m_parallel_groups_angle_map.end(); ++pgi) {
        std::cout << iterator << ") Angle = " << pgi->first << "; ";
        for (std::size_t j = 0; j < pgi->second.size(); ++j) {
          std::cout << pgi->second[j] << " ";
          counter++;
        }
        std::cout << std::endl;
        ++iterator;
      }
      std::cout << "Counter = " << counter << std::endl; 
      // */
/*
      for (auto it_ps = m_parallel_groups_angle_map.begin(); it_ps != m_parallel_groups_angle_map.end(); ++it_ps) {
        const FT theta = it_ps->first;
        const std::vector<std::size_t> &group = it_ps->second;

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
          const Point &barycentre = m_segments[seg_index].m_barycentre;
          const FT c = -a * barycentre.x() - b * barycentre.y();

          set_orientation(seg_index, theta - m_segments[seg_index].m_orientation, a, b, c, v_dir);
        }

      }
      */
    }

   /* void debug_trmu_ijs() {
      std::cout << std::endl << "m_t_ijs: " << std::endl;
      for (typename std::map<std::pair<int, int>, FT>::iterator it = m_t_ijs.begin(); it!=m_t_ijs.end(); ++it)
        std::cout << "(" << it->first.first << ", " << it->first.second << ") => " << it->second << std::endl;
      std::cout << std::endl << "m_r_ijs: " << std::endl;
      for (typename std::map<std::pair<int, int>, FT>::iterator it = m_r_ijs.begin(); it!=m_r_ijs.end(); ++it)
        std::cout << "(" << it->first.first << ", " << it->first.second << ") => " << it->second << std::endl;
      std::cout << std::endl << "m_mu_ij = " << m_mu_ij << std::endl;
    } */


  private:
    // Fields.
    Input_range& m_input_range;
    const Segment_map  m_segment_map;
    std::map <std::size_t, Segment_data> m_segments;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_t_ijs;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_r_ijs;
    const FT m_mu_ij;
    std::shared_ptr<Grouping> m_grouping_ptr;
    // std::map<FT, std::vector<std::size_t>> m_parallel_groups_angle_map;
    const std::map<FT, std::vector<std::size_t>> & m_parallel_groups_angle_map;

    void set_orientation(std::size_t i, const FT new_orientation, const FT a, const FT b, const FT c, const Vector &direction) {

    /*  FT m_orientation = new_orientation;
      FT m_a = a;
      FT m_b = b;
      FT m_c = c;
      Vector m_direction = direction;
      
      if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) 
        m_direction = -m_direction;
      FT x1, y1, x2, y2;
      Point m_barycentre = m_segments[i].m_barycentre;
      FT m_length = m_segments[i].m_length;
      if (CGAL::abs(m_direction.x()) > CGAL::abs(m_direction.y())) { 
        x1 = m_barycentre.x() - m_length * m_direction.x() / FT(2);
        x2 = m_barycentre.x() + m_length * m_direction.x() / FT(2);

        y1 = (-m_c - m_a * x1) / m_b;
        y2 = (-m_c - m_a * x2) / m_b;
      } else {
        y1 = m_barycentre.y() - m_length * m_direction.y() / FT(2);
        y2 = m_barycentre.y() + m_length * m_direction.y() / FT(2);

        x1 = (-m_c - m_b * y1) / m_a;
        x2 = (-m_c - m_b * y2) / m_a;
      }
      const Point source = Point(x1, y1);
      const Point target = Point(x2, y2);

      m_input_range[i] = Segment(source, target); 
      */

    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2
