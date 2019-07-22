#ifndef CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2
#define CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_ordinates_2.h>

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
    using Grouping = internal::Grouping_segments_ordinates_2<Traits, Input_range>;
    using Vector  = typename GeomTraits::Vector_2;

    Ordinate_regularization_2 (
      InputRange& input_range,
      const std::map<FT, std::vector<std::size_t>> & parallel_groups_angle_map,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_parallel_groups_angle_map(parallel_groups_angle_map),
    m_segment_map(segment_map),
    m_grouping(Grouping(input_range)),
    m_mu_ij(FT(4) / FT(5)) {

      CGAL_precondition(input_range.size() > 0);
      // std::size_t counter = 0;
      for (const auto & m_i : m_parallel_groups_angle_map) {
        // ++counter;
        // std::cout << std::endl << counter << ") Angle = " << m_i.first << "; Segments: ";
        if (m_i.second.size() > 1) {
          Point frame_origin;
          for(std::size_t i = 0; i < m_i.second.size(); ++i) {
            const std::size_t seg_index = m_i.second[i];
            // std::cout << seg_index << " ";
            const Segment& seg = get(m_segment_map, *(m_input_range.begin() + seg_index));
            Segment_data seg_data(seg, seg_index, m_i.first);
            if (i == 0) {
              frame_origin = seg_data.m_barycentre;
            }
            const Point reference_coordinates = internal::transform_coordinates(seg_data.m_barycentre, frame_origin, m_i.first);
            seg_data.set_reference_coordinates(reference_coordinates);
            m_segments.emplace(seg_index, seg_data);
          }
        }
      }
      // std::cout << std::endl;

    }

    FT target_value(const std::size_t i, const std::size_t j) {

      //compute_orientation
      // add check whether m_segments[i] and m_segments[j] exist
      const Segment_data & s_i = m_segments.at(i);
      const Segment_data & s_j = m_segments.at(j);

      const FT y_ij = s_i.m_reference_coordinates.y() - s_j.m_reference_coordinates.y();

      // if (i == 8) {
      //   std::cout << "Target value: (" << i << ", " << j << ") = " << y_ij << std::endl;
      // }
      if (CGAL::abs(y_ij) < bound(i) + bound(j)) {
        m_t_ijs[std::make_pair(i, j)] = y_ij;
      }
  
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
      std::map<FT, std::vector<std::size_t>> collinear_groups_by_ordinates;
      std::map <std::size_t, Segment_data> temp_segments;
      std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> temp_t_ijs;
      std::size_t target_index;
      std::size_t counter = 0;
      for (const auto & mi : m_parallel_groups_angle_map) {
        if (mi.second.size() > 1) {
          collinear_groups_by_ordinates.clear();
          temp_segments.clear();
          temp_t_ijs.clear();
          for (const std::size_t it : mi.second) {
            const std::size_t seg_index = it;
            const Segment_data& seg_data = m_segments.at(seg_index);
            temp_segments.emplace(seg_index, seg_data);
            target_index = 0;
            for(const auto & ti : m_t_ijs) {
              if (ti.first.first == seg_index) {
                temp_t_ijs[std::make_pair(ti.first.first, ti.first.second)] = std::make_pair(ti.second, target_index);
              }
              ++target_index;
            }
          }
          if (temp_segments.size() > 0) {
            m_grouping.make_groups(temp_t_ijs, temp_segments, result, collinear_groups_by_ordinates);
            translate_collinear_segments(collinear_groups_by_ordinates);
            /*
            for (const auto & temp_it : collinear_groups_by_ordinates) {
              ++counter;
              std::cout << counter << ") Ordinate = " << temp_it.first << ". ";
              for (std::size_t index = 0; index < temp_it.second.size(); ++index) {
                std::cout << temp_it.second[index] << " ";
              }
              std::cout << std::endl;
            }
            std::cout << std::endl;
            // */
            //compute and set new data for the segments.
          }
        }
      }
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
    const FT m_mu_ij;
    Grouping m_grouping;
    const std::map<FT, std::vector<std::size_t>> & m_parallel_groups_angle_map;

    void translate_collinear_segments(const std::map<FT, std::vector<std::size_t>> & collinear_groups_by_ordinates) {

      for (const auto & mi : collinear_groups_by_ordinates) {
        const FT dt = mi.first;
        // Get the longest segment.
        FT l_max = -FT(1000000000000);
        int l_index = -1;
        for (const std::size_t it : mi.second) {
          const FT seg_length = m_segments.at(it).m_length;
          if (l_max < seg_length) {
            l_max = seg_length;
            l_index = it;
          }
        }
        CGAL_postcondition(l_index >= 0);
        FT new_difference = dt - m_segments.at(l_index).m_reference_coordinates.y();
        set_difference(l_index, new_difference);
        // Translate the longest segment and get the line equation.
        // compute_line_coefficients
        const Segment_data & l_data = m_segments.at(l_index);

        const FT l_a = l_data.m_a;
        const FT l_b = l_data.m_b;
        const FT l_c = l_data.m_c;
        const Vector & l_direction = l_data.m_direction;

        // Translate the other segments, so that they rest upon the line ax + by + c = 0.
        for (const std::size_t it : mi.second) {
          // std::cout << "mi.second = " <<  mi.second << std::endl;
          if (it != l_index) {
            // std::cout << "it = " << it << ". l_index = " << l_index << std::endl;
            new_difference = dt - m_segments.at(it).m_reference_coordinates.y();
            set_difference(it, new_difference, l_a, l_b, l_c, l_direction);
          }
        }
      }

    }

    void set_difference(const int i, const FT new_difference) {

      // std::cout << "new_difference = " << new_difference << std::endl;
        
      const FT difference = new_difference;

      // std::cout << "direction = " << new_difference << std::endl;
      // const Vector & m_direction = m_segments.at(i).m_direction;
      // FT theta = FT(88.0964689526762);

      Segment_data & seg_data = m_segments.at(i);
      // FT theta = seg_data.m_angle;
      // const FT x = static_cast<FT>(cos(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));
      // const FT y = static_cast<FT>(sin(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));

      // const Vector v_dir = Vector(x, y);
      // const Vector & direction = Vector(x, y);
      const Vector & direction = seg_data.m_direction;
      // std::cout << "direction = " << m_direction << std::endl;
      const Vector final_normal = Vector(-direction.y(), direction.x());
      // std::cout << "final_normal = " << final_normal << std::endl;

      const Point &source = seg_data.m_segment.source();
      const Point &target = seg_data.m_segment.target();

      Point new_source = Point(source.x() + difference * final_normal.x(), source.y() + difference * final_normal.y());
      Point new_target = Point(target.x() + difference * final_normal.x(), target.y() + difference * final_normal.y());
      
      // std::cout << "m_a = " << m_segments.at(i).m_a << std::endl;
      const FT bx = (new_source.x() + new_target.x()) / FT(2);
      // std::cout << "bx = " << bx << std::endl;

      // std::cout << "m_b = " << m_segments.at(i).m_b << std::endl;
      const FT by = (new_source.y() + new_target.y()) / FT(2);
      // std::cout << "by = " << by << std::endl;

      m_input_range[i] = Segment(new_source, new_target);
      seg_data.m_barycentre = Point(bx, by);

      // std::cout << "m_c before = " << m_segments.at(i).m_c << std::endl;
      seg_data.m_c = -seg_data.m_a * bx - seg_data.m_b * by;
      // std::cout << "m_c after = " << m_segments.at(i).m_c << std::endl;
      seg_data.m_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_input_range[i].squared_length())));
    }

    void set_difference(const int i, const FT new_difference, const FT a, const FT b, const FT c, const Vector &direction) {

      // std::cout << "new_difference = " << new_difference << " a = " << a << " b = " << b << " c = " << c << " direction = " << direction << std::endl;

      // We translate the segment by the distance new_difference in the direction of the normal vector.
      FT difference = new_difference;
      Segment_data & seg_data = m_segments.at(i);

      // We update the equation of the support line.
      // const FT seg_a = a;
      // const FT seg_b = b;
      // const FT seg_c = c;

      seg_data.m_direction = direction;
      if (seg_data.m_direction.y() < FT(0) || (seg_data.m_direction.y() == FT(0) && seg_data.m_direction.x() < FT(0))) 
        seg_data.m_direction = -seg_data.m_direction;

      Vector final_normal = Vector(-seg_data.m_direction.y(), seg_data.m_direction.x());
      FT x1, x2, y1, y2;

      const Point &source = seg_data.m_segment.source();
      const Point &target = seg_data.m_segment.target();

      if (CGAL::abs(seg_data.m_direction.x()) > CGAL::abs(seg_data.m_direction.y())) {

        x1 = source.x() + difference * final_normal.x();
        x2 = target.x() + difference * final_normal.x(); 


        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;

      } 
      else {
        
        y1 = source.y() + difference * final_normal.y();
        y2 = target.y() + difference * final_normal.y();

        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }

      const Point new_source = Point(x1, y1);
      const Point new_target = Point(x2, y2);

      m_input_range[i] = Segment(new_source, new_target);
    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2
